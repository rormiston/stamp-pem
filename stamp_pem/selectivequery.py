from gwpy.segments import DataQualityFlag, Segment, SegmentList
import numpy as np
from progressbar import (ProgressBar, Percentage, AnimatedMarker)
import os
import sys
import argparse
from stamp_pem import coh_io
from stamp_pem.coherence_segment import PEMSubsystem
from datetime import datetime
from gwpy.time import tconvert
import ConfigParser
from matplotlib import use
use('agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import subprocess as sp


def parse_command_line():
    """
    parse_command_line
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--ini-file", "-i", help="pipeline ini file",
        default=None, type=str, dest='ini')
    parser.add_argument(
        "-st", "--start-time", dest='st', help="query start time",
        type=int, default=None)
    parser.add_argument(
        "-et", "--end-time", dest='et', help="query end time",
        type=int, default=None)
    parser.add_argument(
        "-f", "--frequency_band", dest="frequency_band", nargs='+',
        help="Define frequency band to query", default=None)
    parser.add_argument(
        "--subsystems", dest='subsystems', type=str, help='subsystem',
        nargs='+', default=None)
    parser.add_argument(
        "--new-df", help="coarse graining frequency bin width",
        default=0.5, type=float, dest="new_df")
    parser.add_argument("--output", "-o", help="output directory",
        default=None, dest="output_dir", type=str)

    params = parser.parse_args()
    return params


def parse_cmd_and_configs(params):
    configs = {}
    try:
        ini = params.ini
    except:
        print('Must supply an ini file')

    try:
        configs['st'] = params.st
        configs['et'] = params.et
    except:
        print('Must specify start and stop times')

    if params.frequency_band is None:
        configs['flow'] = 10.0
        configs['fhigh'] = 500.0
    else:
        try:
            configs['flow'] = float(params.frequency_band[0])
            configs['fhigh'] = float(params.frequency_band[-1])
        except:
            configs['flow'] = 10.0
            configs['fhigh'] = 500.0

    # Unpack the config file
    pipeline_dict = coh_io.read_pipeline_ini(ini)
    env_params, run_params = coh_io.check_ini_params(pipeline_dict)
    configs['basedir'] = env_params['base_directory']
    configs['jobdur'] = int(env_params['job_duration'])
    configs['recipients_file'] = env_params['recipients']
    configs['inc'] = env_params['online_increment']
    configs['channel_list'] = env_params['list']
    configs['flag'] = run_params['flag']
    configs['stride'] = float(run_params['stride'])
    configs['darm_channel'] = run_params['darm_channel']
    configs['reference_st'] = int(run_params['reference_st'])
    configs['subsystems'] = run_params['subsystems']
    configs['segmentDuration'] = int(run_params['segmentduration'])
    configs['new_df'] = params.new_df
    configs['subsystem_list'] = get_subsystems(params,
                                              configs['basedir'])
    plotSNR = run_params['plotsnr']
    if plotSNR == "True":
        configs['plotSNR'] = True
    else:
        configs['plotSNR'] = False

    if params.output_dir is None:
        configs['output_dir'] = env_params['output_directory']
    elif "None" in params.output_dir:
        configs['output_dir'] = env_params['output_directory']
    else:
        configs['output_dir'] = env_params['base_directory']

    return configs


def get_subsystems(params, basedir):
    if params.subsystems is None:
        ini = params.ini
        pipeline_dict = coh_io.read_pipeline_ini(ini)
        _, run_params = coh_io.check_ini_params(pipeline_dict)
        subcheck = run_params['subsystems']
        if subcheck.lower() == 'all':
            all_dirs = os.listdir(basedir)
            misc_dirs = ['HTML', 'DAGS', 'SEGMENTS', 'FailedJobs', 'Calibrated ht 1']
            for i in range(len(misc_dirs)):
                try:
                    all_dirs.remove(misc_dirs[i])
                except:
                    # That directory does not exist
                    pass
            subsystem_list = all_dirs
        elif subcheck is not None:
            subsystem_list = run_params['subsystems'].split(', ')
        else:
            pass

    elif params.subsystems[0] == 'all':
        all_dirs = os.listdir(basedir)
        misc_dirs = ['HTML', 'DAGS', 'SEGMENTS', 'FailedJobs', 'Calibrated ht 1']
        for i in range(len(misc_dirs)):
            try:
                all_dirs.remove(misc_dirs[i])
            except:
                # That directory does not exist
                pass
        subsystem_list = all_dirs

    else:
        all_dirs = os.listdir(basedir)
        subsystem_list = [d for d in all_dirs for sub in params.subsystems if sub in d]
    subsystem_list = [s.replace('External', '(External)') for s in subsystem_list]
    return sorted(subsystem_list)


def get_all_subsystems(basedir):
    all_dirs = os.listdir(basedir)
    misc_dirs = ['HTML', 'DAGS', 'SEGMENTS', 'FailedJobs']
    for i in range(len(misc_dirs)):
        try:
            all_dirs.remove(misc_dirs[i])
        except:
            # That directory does not exist
            pass
    subsystem_list = all_dirs
    return sorted(subsystem_list)


def cmd_to_string(param):
    if param is not None:
        return str(param)[1:-1].replace("'", "").replace(",","")
    else:
        return None


def group_subsystems(chanlist, subsystem_list):
    f = open(chanlist)
    sub_names = [line[1:-2] for line in f.readlines() if line.startswith('[')]
    f.close()
    group_dict = {}
    for name in sub_names:
        group_dict[name] = []
        for sub in subsystem_list:
            if name in sub:
                group_dict[name].append(sub)
    for k, v in group_dict.items():
        if len(group_dict[k]) == 0:
            group_dict.pop(k)
    return group_dict


# Align start time to the nearest half hour
def adjusted_window(st, et, jobdur):
    """
    The given start time may not land on the start of
    a two hour segment. This function adjusts the start time
    and end time to the nearest two hour increments.

    Parameters
    ----------
    st : `int`
        gps start time
    et : `int`
        gps end time
    jobdur : `int`
        Job duration

    Returns
    -------
    adj_time_list : `list`
        Time list with adjusted start and end times
    """
    starttime = tconvert(st)
    daystart = int(tconvert(datetime(starttime.year, starttime.month,starttime.day)))
    times = np.arange(daystart, et + 2 * jobdur, jobdur)
    index_nearest_st = np.argmin(np.abs(times - st))
    index_nearest_et = np.argmin(np.abs(times - et))
    if times[index_nearest_st] - st > 0:
        index_nearest_st -= 1
    if times[index_nearest_et] - et < 0:
        index_nearest_et += 1

    nearest_st = times[index_nearest_st]
    nearest_et = times[index_nearest_et]
    adjusted_times = np.arange(nearest_st, nearest_et + jobdur, jobdur)
    return adjusted_times


# Pretty printing of active segments
def print_segs(flag, st, et, jobdur):
    """
    Print the active segment info. This should be identical
    to what you find on the detchar summary pages.

    Parameters
    ----------
    flag : `str`
        Either H1:DMT-ANALYSIS_READY:1
        or L1:DMT-ANALYSIS_READY:1
    st : `int`
        gps start time
    et : `int`
        gps end time
    jobdur : `int`
        Job duration
    """
    proxyCheck = sp.Popen(['grid-proxy-info'], stdout=sp.PIPE).communicate()[0]
    if proxyCheck:
        active_seg_list = DataQualityFlag.query(flag, st, et).active
    else:
        print('\nCannot query segment database without a valid proxy. Most likely,')
        print('you need to run: ligo-proxy-init albert.einstein')
        return 1

    adj_time = adjusted_window(st, et, jobdur)
    tot_active = 0
    print('\nActive segments for flag {0}'.format(flag))

    if (adj_time[0] != st) or (adj_time[-1] != et):
        print('in adjusted window ({0}, {1})'.format(adj_time[0],
                                                     adj_time[-1]))
    else:
        print('in window ({0}, {1})'.format(st, et))

    print('{0:15} {1:15} {2:15}'.format('Start', 'End', 'Duration'))

    for aseg in active_seg_list:
        st, et = int(aseg[0]), int(aseg[-1])
        dur = et - st
        tot_active += dur
        print('{0:<15} {1:<15} {2:<15}'.format(st, et, dur))

    print('{0}'.format('.' * 40))
    print('{0:<30}  {1:<15}\n'.format('Total', tot_active))


def logo():
    """
    Print the stamp_pem logo
    """
    print(r'''
          __
   _____ / /_ ____ _ ____ ___   ____          ____   ___   ____ ___
  / ___// __// __ `// __ `__ \ / __ \        / __ \ / _ \ / __ `__ \
 (__  )/ /_ / /_/ // / / / / // /_/ /       / /_/ //  __// / / / / /
/____/ \__/ \__,_//_/ /_/ /_// .___/______ / .___/ \___//_/ /_/ /_/
                            /_/    /_____//_/
    ''')


def usage(full_file_name):
    cmd = "python {0} --help".format(full_file_name)
    os.system(cmd)
    print()


def is_first_run(basedir, timefile):
    basedir += '/Calibrated ht 1'
    with open(timefile) as f:
        currentTime = f.read()[:5]
    PATH = basedir + '/' + currentTime
    all_dirs = os.listdir(PATH)
    if len(all_dirs) > 1:
        return False
    else:
        return True


def read_segs(subsystem, adj_time, basedir, output_dir,
               jobdur, darm_channel, subsys=None):

    First = 1
    for i in range(len(adj_time)):
        window_start = adj_time[i]
        window_end = window_start + jobdur
        cohdir = coh_io.get_directory_structure(subsystem, window_start,
                                                directory=basedir)
        cohfile = coh_io.create_coherence_data_filename(darm_channel,
                                                        subsystem,
                                                        window_start,
                                                        window_end,
                                                        directory=cohdir)
        if First:
            try:
                subsys = PEMSubsystem.read(subsystem, cohfile)
                First = 0
            except IOError:
                # print "Couldn't load %s" % cohfile
                continue
        else:
            try:
                temp = PEMSubsystem.read(subsystem, cohfile)
            except IOError:
                # print "Couldn't load %s" % cohfile
                continue
            subsys.update(temp)
    return subsys


# Make the necessary folders
def create_folder_structure(subsystems, st, directory='./'):
    """
    Function that creates the necessary folder hierarchy for a
    given subsystem(s) and start time.

    Parameters
    ----------
    subsystems : `str` or `list`
        Subsystems to run over
    st : `int`
        gsp start time
    """
    subsystems.append('HTML')

    for subsystem in subsystems:
        subsys = subsystem.replace(' ', '\ ').replace('(', '').replace(')', '')
        if subsystem == 'HTML':
            cmd = 'mkdir -p %s/%s' % (directory, subsys)
        else:
            cmd = 'mkdir -p %s/%s/%s' % (directory, subsys, str(st)[0:5])
        os.system(cmd)


def webpage_info(st, et, output_dir, jobdur):
    starttime = tconvert(st)
    daystart = int(tconvert(datetime(starttime.year, starttime.month,
                                     starttime.day)))
    times = np.arange(daystart, daystart + 8640000 + jobdur, jobdur)
    nearest_st = list(abs(times - st)).index(min(abs(times - st)))
    st = times[nearest_st]
    nearest_et = list(abs(times - et)).index(min(abs(times - et)))
    et = times[nearest_et]
    save_webpage_str = '%d%02d%02d' % (starttime.year, starttime.month,
                                       starttime.day)
    datestrmdy = '%02d-%02d-%d' % (starttime.month, starttime.day, starttime.year)
    datestrdmy = '%02d-%02d-%d' % (starttime.day, starttime.month, starttime.year)
    datestrymd = '%d%02d%02d' % (starttime.year, starttime.month, starttime.day)
    webpage_output_dir = '%s/%s/day/%s' % (output_dir, 'HTML', save_webpage_str)
    webpage_info = {}
    webpage_info['mdy'] = datestrmdy
    webpage_info['dmy'] = datestrdmy
    webpage_info['ymd'] = datestrymd
    webpage_info['output_dir'] = webpage_output_dir
    return webpage_info


def get_subsystem_channels(subsystem, channel_list):
    config = ConfigParser.ConfigParser()
    config.read(channel_list)
    try:
        chans = config.get(subsystem, 'channels').split('\n')
        chans = [chan for chan in chans if len(chan) > 0]
        return chans
    except:
        print('No channels found for subsystem: {0}'.format(subsystem))


def get_subsystem_info(subsystem_name, basedir, output_dir, darm_channel,
                       jobdur, adj_time, ref_adj_time=None, flow=10, fhigh=500,
                       residual=False, rst=None, ret=None, new_df=None):
    subsystem_info = {}
    pemsub = read_segs(subsystem_name, adj_time, basedir, output_dir,
                        jobdur, darm_channel, subsys=None)
    if pemsub is None:
        return None

    if new_df is not None:
        pemsub.coarse_grain(flowy=flow, deltaFy=new_df)

    frequencies = pemsub.values()[0].psd1.frequencies.value
    # Need to round to 0.1Hz (due to floats, computer gives
    # annoying rounding errors. i.e, 0.5 -> 0.4999999999)
    frequencies = np.array([round(x, 1) for x in frequencies])
    flow_index = np.where(frequencies == flow)[0][0]
    fhigh_index = np.where(frequencies == fhigh)[0][0]
    frequencies = frequencies[flow_index:fhigh_index]
    chans = pemsub.keys()
    cohs = [cohseg.get_coh().value[flow_index:fhigh_index]
            for cohseg in pemsub.values()]
    if residual:
        pemsubref = read_segs(subsystem_name, ref_adj_time, basedir, output_dir,
                               jobdur, darm_channel, subsys=None)
        if new_df is not None:
            pemsubref.coarse_grain(flowy=flow, deltaFy=new_df)

        cohsref = [cohseg.get_coh().value[flow_index:fhigh_index]
                   for cohseg in pemsubref.values()]
        cohs = [abs(cohs[i] - cohsref[i]) for i in range(len(cohs))]

    subsystem_info = {chans[i]:cohs[i] for i in range(len(chans))
                      if len(cohs[i]) != 0}
    return subsystem_info


def plot_residuals(sub1, sub2, webpage_basedir, st, et, ref_st, ref_et,
                   flow=None, fhigh=None, plotSNR=True):
    # initialize coherence matrix
    freqs = sub1[sub1.keys()[0]].psd1.frequencies.value
    max_freq = np.max(
                    np.asarray([sub1[key].psd2.frequencies.value[-1]
                                for key in sub1.keys()]))
    coh_matrix = np.zeros((len(sub1.keys()),
                           sub1[sub1.keys()[0]].psd1.size))
    labels = []
    for label in sub1.keys():
        label = label.replace('_OUT_DQ','').replace('H1:','').replace('H1:','')\
                .replace('L1:','').replace('_','\_')
        labels.append(label)

    # fill in coherence matrix
    for ii, key in enumerate(sub1.keys()):
        coh1 = sub1[key].get_coh()
        if sub2 is None:
            coh2 = np.zeros(len(coh1))
        else:
            coh2 = sub2[key].get_coh()
        coh = np.abs(coh2 - coh1)
        if plotSNR:
            coh_matrix[ii, :coh.size] = coh.value * sub1[key].N
        elif not plotSNR:
            coh_matrix[ii, :coh.size] = coh.value

    Navg = sub1[key].N
    fig = plt.figure(figsize=(12,6))
    ax = plt.gca()
    plt.yticks(np.arange(1, len(labels) + 1) - 0.5, labels, fontsize=6)

    if plotSNR:
        plt.pcolormesh(freqs, np.arange(0, len(sub1.keys())+1), coh_matrix,
                       norm=LogNorm(vmin=1, vmax=1e3), cmap='viridis')
        cbar = plt.colorbar(label='coherence SNR')
        ax.set_xlabel(r'Frequency [Hz]')
    elif not plotSNR:
        plt.pcolormesh(freqs, np.arange(0, len(sub1.keys())+1), coh_matrix,
                       norm=LogNorm(vmin=1e-3, vmax=1e0), cmap='viridis')
        cbar = plt.colorbar(label='coherence')
        plt.xscale('log')
        ax.set_xlabel(r'Log Frequency [Hz]')

    if fhigh == None:
        fhigh = freqs[-1]
    if fhigh > max_freq:
        fhigh = max_freq
    if flow == None:
        flow = freqs[0]

    ax.set_xlim(flow, fhigh)
    ax.set_title(r'%s, %d averages' % (sub1.subsystem, Navg))
    darm = str(sub1.darm_channel).replace(':', '-')
    temp_name = sub1.subsystem.replace('(','').replace(')','').split(' ')
    name = ''.join([word[0].upper() for word in temp_name])
    if "Microphone" in sub1.subsystem:
        name = name.replace("PEMM", "PEMMic")
    name = '{0}-{1}_{2}-{3}_{4}-{5}'.format(darm, name, st, et, ref_st, ref_et)
    PATH = webpage_basedir + '/Residuals_Plots'
    os.system('mkdir -p {}'.format(PATH))
    plt.savefig('{0}/{1}'.format(PATH, name))
