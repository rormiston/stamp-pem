from stamp_pem.coherence_segment import (PEMSubsystem,
        ChannelDict)
from stamp_pem import coh_io
from gwpy import time
from gwpy.detector import Channel
import numpy as np
from gwpy.spectrogram import Spectrogram
"""
This module needs to be updated at each site based on its online setup.
"""

# Global variables customized for online setup
def get_channels_subsystem(channel_name, channel_dict):
    our_subsys = None
    for key in channel_dict.keys():
        for chan in channel_dict[key]:
            if isinstance(chan, Channel):
                chan = chan.name
            if chan == channel_name:
                our_subsys = key
    if our_subsys is not None:
        return our_subsys
    else:
        raise ValueError('%s was not found in the online channel list' % channel_name)

def get_channel_online_data(channel, st, et, format='spectrogram',
        remove_nonlocked_times=False, normalize_coherence=False,
        config_file='/home/stochastic/config_files/ini_files/H1.ini'):
    """
    Returns a list of PEMCoherenceSegment
    objects.

    Parameters
    ----------
    channel : str
        channel name you want to load
    st : str or int
        start time (in string format) or gps time
    et : str or int
        end time (in string format) or gps time
    format : str, optional, default='spectrogram'
        format to return. either spectrogram or seglist. Spectrogram returns a
        `gwpy.spectrogram.Spectrogram` and seglist returns a list of
        `stamp_pem.coherence_segment.PEMCoherenceSegment`.
    remove_nonlocked_times: bool, optional, default=False
        Removes non locked times from a spectrogram
    normalize_coherence : bool, optional, default=False
        Normalizes each column of spectrogram by the number of averages

    Returns
    -------
    out : `gwpy.spectrogram.Spectrogram` or list of
    `stamp_pem.coherence_segment.PEMCoherencSegment` objects
        representation of data between start and end times for a given channel
    """
    pipeline_dict = coh_io.read_pipeline_ini(config_file)
    env_params, run_params = coh_io.check_ini_params(pipeline_dict)
    channel_dict = ChannelDict.read(env_params['list'])
    jobdur = int(env_params['job_duration'])
    darm_channel = run_params['darm_channel']
    basedir = env_params['base_directory']

    if isinstance(st, str):
        st = int(time.to_gps(st))
    if isinstance(et, str):
        et = int(time.to_gps(et))
    starttimes = np.arange(st, et, jobdur)
    subsys = get_channels_subsystem(channel, channel_dict)
    seglist = []
    for starttime in starttimes:
        cohdir = coh_io.get_directory_structure(subsys, starttime, basedir)
        cohfile = coh_io.create_coherence_data_filename(darm_channel, subsys,
                starttime, starttime+jobdur, directory=cohdir)
        try:
            subsystem_data = PEMSubsystem.read(subsys, cohfile)
        except IOError:
            print "No data found between %d and %d for %s" % (starttime, starttime + jobdur, channel)
            continue

        if np.isnan(subsystem_data[channel].psd1.value[0]):
            continue
        seglist.append(subsystem_data[channel])

    N = 1

    if format=='spectrogram':
        if remove_nonlocked_times:
            foundtimes = np.asarray([seglist[ii].starttime for ii in
                range(len(seglist))])
            data = np.zeros((len(seglist), seglist[0].psd1.size))
            for ii in range(len(seglist)):
                if normalize_coherence:
                    N=seglist[ii].N
                if seglist[ii].get_coh()[0] == np.nan:
                    continue
                data[ii,:] = seglist[ii].get_coh() * N
                specgram = Spectrogram(data, epoch=foundtimes[0], dt=jobdur,
                        df=seglist[0].psd1.df)
        else:
            foundtimes = np.asarray([seglist[ii].starttime for ii in
                range(len(seglist))])
            count = 0
            data = np.nan * np.zeros((starttimes.size, seglist[0].psd1.size))
            for ii, starttime in enumerate(starttimes):
                if np.any(foundtimes==starttime):
                    if normalize_coherence:
                        N = seglist[count].N
                    data[ii,:] = seglist[count].get_coh() * N
                    count += 1
            specgram = Spectrogram(data, dt=jobdur,
                    epoch=starttimes[0],df=seglist[0].psd1.df)
        return specgram
    elif format=='seglist':
        return seglist
    else:
        raise ValueError('format needs to be "spectrogram" or "seglist"')
