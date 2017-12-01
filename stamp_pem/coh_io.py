import ConfigParser
import os
from gwpy.segments import DataQualityFlag
import coherence_functions as cf


def config_section_map(c, section):
    dict1 = {}
    options = c.options(section)
    for option in options:
        dict1[option] = c.get(section, option)
    return dict1


def config_list(c):
    dict1 = {}
    for section in c.sections():
        dict1[section] = config_section_map(c, section)
    return dict1


def config_pipeline_ini(c):
    dict1 = {}
    for section in c.sections():
        dict1[section] = config_section_map(c, section)
    return dict1


def read_pipeline_ini(file):
    """
    Read ini file for analysis and put
    it into a dict.

    Parameters
    ----------
    file : `str`
        ini file for coherence pipeline

    Returns
    -------
    dict1 : `dict`
        dictionary with parameters for pipeline
    """
    c = ConfigParser.ConfigParser()
    c.read(file)
    dict1 = config_pipeline_ini(c)
    return dict1


def read_list(list):
    """
    read channel list from file

    Parameters
    ----------
    list : `str`
        file with list of channels in ini format

    Returns
    --------
    channels : `dict`
        channels in dict formaat with first keys
        as subsystems, second as numerals that
        reference channel names
    """
    c = ConfigParser.ConfigParser()
    c.read(list)
    channels = config_list(c)
    return channels


def extract_subsystem(channel_dict, subsystem):
    """
    extract channels into list form from channel_dict

    Parameters
    ----------
    channel_dict : `dict`
        channel dictionary
    subsystem : `str`
        subsystem tag

    Returns
    -------
    channels : `list` (`str`)
        list of channels from that subsystem
    """
    channels = []
    for chan in channel_dict[subsystem].keys():
        channels.append(channel_dict[subsystem][chan])
    return channels


def create_coherence_data_filename(darm_channel, subsystem, st, et,
                                   directory='./', tag=None):
    """
    Create coherence data filename. Name is:
    `directory/channel-subsystem-TAG-st-duration`

    Parameters
    ----------
    darm_channel : `str`
        darm channel name
    subsystem : `str`
        subsystem name
    st : `int`
        start time
    et : `int`
        end time
    directory : `str`, optional, default='./'
        top level directory for filename
    tag : `str`, optional
        tag to include in creating filename

    Returns
    -------
    filename : `str`
        name of file
    """
    filename = darm_channel.replace(
        ':', '-') + '-' + subsystem + '-' +\
        str(st) + '-' + str(et - st)
    chan = darm_channel.replace(':', '-')
    spl = subsystem.replace('(', ' ').replace(')', '').split()
    subsys_shortkey = ''
    for s in spl:
        subsys_shortkey += s[0].upper()

    if tag:
        filename = '%s%s-%s-%s-%d-%d' % (directory,
                                          chan, subsys_shortkey, tag, st, et - st)
    else:
        filename = '%s%s-%s-%d-%d' % (directory, chan, subsys_shortkey, st, et - st)
    return filename


def get_directory_structure(subsystem, st, directory='./', specgram=False):
    """
    Returns directory structure:
    `directory/subsystem/str(time)[:5]/plots`
    only includes `plots` directory if `specgram` is `True`.


    Parameters
    ----------
    subsystem : `str`
        subsystem
    st : `int`
        start time
    directory : `str`, optional, default='./'
        top level directory
    specgram : `bool`, optional, default=False
        if you're saving a spectrogram it goes in the plots directory
        whereas if you're saving something else it doesn't

    Returns
    -------
    directory : `str`
        return directory sructure for stamp_pem run
    """
    if specgram:
        return '%s/%s/%s/%s' % (directory, subsystem.replace('(','').replace(')',''), str(st)[0:5], 'plots')
    else:
        return '%s/%s/%s/' % (directory, subsystem.replace('(','').replace(')',''), str(st)[0:5])


def create_directory_structure(subsystems, st, directory='./'):
    """
    Creates the directory structure to save data for stamp-pem run.
    Includes extra directories for HTML (where the web page is saved),
    SEGMENTS (where the segment lists are saved), and DAGS, where the dag
    and sub files are saved.

    creates directory structure for:
    `directory/subsystem/str(time)[:5]/plots`

    Parameters
    ----------
    subsystems : `list`
        List of subsystems we're running over.
    st : `int`
        start time
    directory : `str`, optional, default='./'

    Returns
    -------
    None
    """
    subsystems.append('SEGMENTS')
    subsystems.append('DAGS')
    subsystems.append('HTML')
    subsystems.append('FailedJobs')

    for subsystem in subsystems:
        subsys = subsystem.replace(' ','\ ').replace('(','').replace(')','')
        cmd = 'mkdir -p %s/%s/%s/%s' % (directory, subsys,
                                        str(st)[0:5], 'plots')
        os.system(cmd)


def write_segs(flag, st, et, directory='./'):
    """
    Write segments to a file

    saved in:
    `directory/flag-st-duration`

    Parameters
    ----------
    flag : `str`
        DQ flag to query database for
    st : `int`
        start time
    et : `int`
        end time
    directory : `str`, optional, default='./'
        top level directory for saving segments
    """
    segments = DataQualityFlag.query_dqsegdb(flag, st, et,
                                             url='https://segments.ligo.org')
    directory = get_directory_structure('SEGMENTS', st, directory=directory)
    segments.write('%s/%s-%d-%d.xml.gz' % (directory, flag, st, (et - st)))


def increment_datetime(time, increment):
    """
    Increment a gps time

    Parameters
    ----------
    time : `int`
        start gps time
    increment : `int`
        time in (s) to increment by

    Returns
    -------
    incremented_time : `int`
        time + increment
    """
    time += increment
    return time


def increment_datetime_in_file(f, increment):
    """
    increment gps time from file

    Parameters
    ----------
    f : `str`
        file to increment
    increment : `int`
        amount of time by which to increment
    """
    l = read_time_from_file(f)
    new_time = increment_datetime(l, increment)
    write_time_to_file(f, new_time)


def read_time_from_file(f):
    f = open(f, 'r')
    l = int(f.read().rstrip())
    f.close()
    return l


def write_time_to_file(f, time):
    f = open(f, 'w')
    f.write(str(time))
    f.close()


def get_channel_dict_from_ascii(ascii_file, subsystems):
    channels = {}
    for sub in subsystems:
        f = open(ascii_file, 'r')
        counter = 0
        for line in f:
            if line[0] == '#' or len(line.split(':')) < 2:
                continue
            chan = line.split(':')[1]
            sub_temp = chan.split('_')[0]
            if sub in sub_temp:
                counter += 1
                channels[sub][str(counter)] = line
        f.close()
    return channels


def check_channel_and_flag(channel, flag):
    """
    make sure the ifos are the same for the two channels
    """
    ifo1 = channel.split(':')[0]
    ifo2 = channel.split(':')[0]
    if ifo1 == ifo2:
        return True
    else:
        return False


def check_ini_params(pipeline_dict):
    """
    Check that the necessary input parameters
    actually exist.
    """
    env_params = pipeline_dict['env']
    run_params = pipeline_dict['run']
    try:
        stride = run_params['stride']
        darm_channel = run_params['darm_channel']
        flag = run_params['flag']
    except KeyError:
        raise KeyError(
            'stride, darm_channel, and dq flag must be in the ini file')

    try:
        segmentDuration = run_params['segmentduration']
    except KeyError:
        run_params['segmentDuration'] = stride
        print 'segmentDuration not set. setting it to stride'

    try:
        fhigh = run_params['fhigh']
    except KeyError:
        print 'fhigh not seg. setting it to 16kHz'
        run_params['fhigh'] = 16384

    try:
        spec_fhigh = run_params['spec_fhigh']
    except KeyError:
        run_params['spec_fhigh'] = None

    try:
        spec_flow = run_params['spec_flow']
    except KeyError:
        run_params['spec_flow'] = None
    return env_params, run_params


def check_channels(channels, st):
    """
    Check whether the channel actually exists

    Parameters
    ----------
    channels : `list`
        list of channels
    st : `int`
        start time
    """
    counter = 0
    failed_chans = []
    for channel in channels:
        try:
            data = cf._read_data(channel, st, st + 1)
        except ValueError:
            failed_chans.append(channel)
            del channels[counter]

        counter += 1
    return channels, failed_chans
