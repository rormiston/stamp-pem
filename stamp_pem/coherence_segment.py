from gwpy.frequencyseries import FrequencySeries
from gwpy.timeseries import TimeSeries
import time
# import stamp_pem.coherence_functions as cf
import coherence_functions as cf
from glue import datafind
import numpy as np
from matplotlib import rc
rc('text', usetex=True)
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy import stats
from collections import OrderedDict
import h5py
from gwpy.detector import (Channel, ChannelList)
import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)
plt.rcParams.update({'figure.max_open_warning': 0})


class PEMCoherenceSegment(object):
    """
    docstring for PEMCoherenceSegment

    Methods
    -------
    .. py:method:: coarse_grain
        Coarse grain a :py:class:`gwpy.frequencyseries` object
    .. py:method:: coherence
        Calculates coherence between two
        channels and returns a :py:class:`PEMCoherenceSegment`
    .. py:method:: comb_finder
        Returns a :py:class:`CoherenceCombMatrix` object
    .. py:method:: get_coh
        Calculate coherence and cut down psd's to make sure
        they're the same size as the csd
    .. py:method:: plot
        Gives a plot of :math:`\\frac{\\text{CSD}}{\\text{PSD1} \\times \\text{PSD2}}`
    .. py:method:: sum_coh
        Sum up coherences (with a mask applied if you so desire).

        For use in :py:meth:`stamp_pem.coherence_segment.PEMCoherenceSegment.comb_finder`
    .. py:method:: update
        Update csd12, psd1, psd2, starttime, endtime, N
    """

    def __init__(self, darm_channel, aux_channel, csd12, psd1, psd2, N, starttime, endtime):
        super(PEMCoherenceSegment, self).__init__()
        self.darm_channel = darm_channel
        self.aux_channel = aux_channel
        self.csd12 = csd12
        self.psd1 = psd1
        self.psd2 = psd2
        self.N = N
        self.starttime = starttime
        self.endtime = endtime

    @classmethod
    def coherence(cls, channel1, channel2, st, et, overlap=None, pad=False, stride=1,
                  resamplerate1=None, resamplerate2=None):
        """
        Class methond that calculates coherence between two channels
        and return a coherence segment.

        Parameters
        ----------
        channel1 : `str`
            Name of first channel
        channel2 : `str`
            Name of second channel
        st : `int`
            Start time
        et : `int`
            End time
        stride : `int`, optional, default=1 second
            Length of ffts
        overlap : `int`, optional, default=0.5 seconds
            Amount of overlap between ffts
        pad : `bool`
            Determines whether or not to zeropad the data
            when taking ffts

        Returns
        -------
        segment : :class:`PEMCoherenceSegment`
            Coherence segment for these two channels
        """
        # read in data
        if isinstance(channel1, Channel):
            data1 = TimeSeries.find(channel1.name, st, et, frametype=channel1.frametype)
        else:
            data1 = TimeSeries.get(channel1, st, et)
        if resamplerate1 is not None and resamplerate1 < data1.sample_rate.value:
            data1 = data1.resample(resamplerate1)
        if isinstance(channel2, Channel):
            data2 = TimeSeries.find(channel2.name, st, et, frametype=channel2.frametype)
        else:
            data2 = TimeSeries.get(channel2, st, et)
        if resamplerate2 is not None and resamplerate2 < data2.sample_rate.value:
            data2 = data2.resample(resamplerate2)

        # get fft spectrograms
        fftgram1 = cf.fftgram(data1, stride, overlap=overlap, pad=pad)
        fftgram2 = cf.fftgram(data2, stride, overlap=overlap, pad=pad)
        # cut things down if frequency arrays are too long
        # TODO: eventually port this over to `gwpy.spectrogram.Spectrogram.crop()`
        # in some way using the `gwpy.detector.Channel.frequency_range()` specified
        # but want to keep backwards compatability in case strings are supplied
        maxlen = min(fftgram1.shape[1], fftgram2.shape[1])
        # take csd
        csd12 = cf.csdgram(fftgram1[:,:maxlen], fftgram2[:,:maxlen], stride,
                           overlap=overlap, pad=pad)

        # get number of segments analyzed
        N = fftgram1.shape[0]
        # take mean of csd
        csd12 = np.mean(csd12, 0)
        # take mean of fftgrams, take abs to get psds
        psd1 = np.mean(np.abs(fftgram1)**2, 0)
        psd2 = np.mean(np.abs(fftgram2)**2, 0)
        # return the segment
        return PEMCoherenceSegment(channel1, channel2, csd12, psd1, psd2, N, st, et)

    def update(self, cohseg):
        """
        Update csd12, psd1, psd2, starttime, endtime, N

        Parameters
        ----------
        cohseg : :py:class:`PEMCoherenceSegment` object.

        Returns
        -------
        Gives an updated :py:class:`PEMCoherenceSegment` object.
        """
        self.csd12 = (self.N * self.csd12 + cohseg.N * cohseg.csd12) / (self.N + cohseg.N)
        self.psd1 = (self.N * self.psd1 + cohseg.N * cohseg.psd1) / (self.N + cohseg.N)
        self.psd2 = (self.N * self.psd2 + cohseg.N * cohseg.psd2) / (self.N + cohseg.N)
        self.starttime = min(self.starttime, cohseg.starttime)
        self.endtime = max(self.endtime, cohseg.endtime)
        self.N = self.N + cohseg.N
        return

    def plot(self, **kwargs):
        """
        Gives a plot of :math:`\\frac{\\text{CSD}}{\\text{PSD1} \\times \\text{PSD2}}`
        """
        coh = np.abs(self.csd12) ** 2 / (self.psd1 * self.psd2)
        return coh.plot(**kwargs)

    def get_coh(self, cutoff=None):
        """
        Calculate coherence and cut down psd's to make sure
        they're the same size as the csd.
        """
        coh = np.abs(self.csd12) ** 2 / (self.psd1[:self.csd12.size] * self.psd2[:self.csd12.size])
        if cutoff is not None:
            coh[coh > np.float(cutoff)/self.N] = np.float(cutoff)/self.N
        return coh

    def sum_coh(self, mask, cutoff=None):
        """
        Sum up coherences (with a mask applied if you so desire). For use in
        :meth:`stamp_pem.coherence_segment.PEMCoherenceSegment.comb_finder`
        """
        coh = np.copy(self.get_coh())
        if cutoff is not None:
            coh[coh > np.float(cutoff)/self.N] = np.float(cutoff)/self.N
        return np.sum(coh * mask)

    def comb_finder(self, offsets, comb_widths, cutoff=None):
        """
        Parameters
        ----------
        offsets : `int`
        comb_widths : `int`
        cutoff : `bool`, optional

        Returns
        -------
        :py:class:`CoherenceCombMatrix` object
        """
        comb_mat = np.zeros((len(offsets), len(comb_widths)))
        likelihood_mat = np.zeros((len(offsets), len(comb_widths)))
        f0 = self.psd1.frequencies[0].value
        df = self.psd1.df.value
        coh = FrequencySeries(self.get_coh(cutoff=cutoff), frequencies=self.psd1.frequencies)
        for ii, offset in enumerate(offsets):
            start_idx = max(0, ((offset - f0) / df))
            for jj, width in enumerate(comb_widths):
                mask = np.zeros(self.psd1.size, dtype=bool)
                comb_idx_width = width / np.float(df)
                mask[start_idx::comb_idx_width] = True
                comb_mat[ii,jj] = self.sum_coh(mask, cutoff=cutoff)
                likelihood_mat[ii, jj] = stats.gamma.logpdf(self.N*self.sum_coh(mask, cutoff=cutoff),
                                                         np.sum(np.int_(mask))) * self.N
        return CoherenceCombMatrix(comb_mat, likelihood_mat, offsets, comb_widths, coh)

    def coarse_grain(self, flowy=None, deltaFy=None, Ny=None):
        """
        Coarse grain the :py:class:`gwpy.frequencyseries.FrequencySeries`
        objects in the :py:class:`PEMCoherenceSegment`, ie, csd12, psd1, psd2.
        Returns an updated :py:class:`PEMCoherenceSegment` object.

        Parameters
        ----------
        flowy : `float`
            Starting frequency value for the coarse grained
            spectrum. Called *f0* in :py:class:`gwpy.frequencyseries.FrequencySeries`.
            The default value is the minimum allowed,
            :math:`\\text{flowy} = \\text{flowx} +\\frac{1}{2}(\Delta F_y - \Delta F_x)`
        deltaFy : `float`
            Step size of the coarse grained spectrum. The
            default is the minimum allowed,
            :math:`\\text{flowy} = \\text{flowx}`
        Ny : `int`
            Number of steps. The default is the maximum amount
            with the given step size and starting frequency.
        """
        self.psd1 = cf.coarseGrain(self.psd1, flowy, deltaFy, Ny);
        self.psd2 = cf.coarseGrain(self.psd2, flowy, deltaFy, Ny);
        self.csd12 = cf.coarseGrain(self.csd12, flowy, deltaFy, Ny);
        self.flowy = flowy
        self.deltaFy = deltaFy
        self.Ny = Ny
        return

class PEMSubsystem(dict):
    """
    PEM subsystem object. Contains information
    on a subsystem and if coherence has been
    taken then it contains a dictionary of
    :py:class:`PEMCoherenceSegment` values with
    the non-darm channel (channel2) as the key.

    Methods
    -------
    .. py:method:: coherence
        Calculates coherence between two
        channels and return a coherence segment.
    .. py:method:: write
        Write class info into a file
    .. py:method:: read
        Read class info from a file
    .. py:method:: plot
        Return a plot of the coherence matrix
    .. py:method:: update
        Update csd12, psd1, psd2, starttime, endtime, N
    .. py:method:: update_from_file
        Invoke :py:func:`PEMSubsystem.read()` to update the class
    .. py:method:: coarse_grain
        Coarse grain a :py:class:`gwpy.frequencyseries` object

    Parameters
    ----------
    darm_channel : `str`
        Name of first channel
    subsystem : `str`
       Name of subsystem we're running on
    channel_dict : `str`
        `dict` of channels to read in from detchar channel list
        to cross correlate with DARM channel
    """

    def __init__(self, subsystem, darm_channel, channel_dict):
        super(PEMSubsystem, self).__init__()
        self.subsystem = subsystem
        self.channel_dict = channel_dict
        self.darm_channel = darm_channel
        self.failed_channels=[]

    @classmethod
    def coherence(cls, darm_channel, subsystem, channel_dict, st, et,
                  stride=1, overlap=None, pad=False, resamplerate1=None, resamplerate2=None):
        """
        Class method that calculates coherence between darm and subsystems.

        Parameters
        ----------
        darm_channel : `str`
            Name of first channel
        subsystem : `str`
           Name of subsystem we're running on
        channel_dict : `str`
            `dict` of channels to read in from detchar channel list
            to cross correlate with DARM channel
        st : `int`
            Start time
        et : `int`
            End time
        stride : `int`, optional, default=1 second
            Length of ffts
        overlap : `int`, optional, default=0.5 seconds
            Amount of overlap between ffts
        pad : `bool`
            Determines whether or not to zeropad the data
            when taking ffts

        Returns
        -------
        subsystem : :class:`PEMSubsystem`
            `dict` with channel_list elements as keys
            and PEMCoherenceSegment objects as values
        """
        start = time.time()
        subsys = PEMSubsystem(subsystem, darm_channel, channel_dict)
        for channel in channel_dict[subsystem]:
            subsys[channel.name] =\
                PEMCoherenceSegment.coherence(darm_channel, channel,
                                              st, et, stride=stride,
                                              overlap=overlap,
                                              pad=pad,
                                              resamplerate1=resamplerate1,
                                              resamplerate2=resamplerate2)
            #except ValueError as VE:
            #    print '%s failed with %s' % (channel, VE)
            #    subsys.failed_channels.append(channel)
        end = time.time()
        print (end - start)
        return subsys

    def write(self, filename):
        """
        Write to a file

        Parameters
        ----------
        filename : `str`
            File to write to
        """
        # open file
        if isinstance(filename, str):
            f = h5py.File(filename, 'w')
        else:
            f = filename
        # create groups and datasets
        f.create_group(self.subsystem)
        # for ease of coding...
        ss = self.subsystem
        for key in self.keys():
            if key in self.failed_channels:
                continue
            metadata = [('N', self[key].N),
                       ('st', self[key].starttime),
                       ('et', self[key].endtime)]
            f[ss].create_group(key)
            self[key].psd2.write(f[ss][key])
            self[key].csd12.write(f[ss][key])
            f[ss][key].create_dataset('metadata',data=metadata)
        f.create_group('psd1')
        self[key].psd1.write(f['psd1'])
        f.create_dataset('failed_channels', data=self.failed_channels)
        f.close()

    @classmethod
    def read(cls, subsys, filename):
        """
        read from a file

        Parameters
        ----------
        subsys : `str`
            Channel subsystem
        filename : `str`
            File to read from
        """
        # open file
        if isinstance(filename,str):
            f = h5py.File(filename,'r')
        else:
            f = filename
        chandict = {subsys: f[subsys].keys()}
        darm_channel = f['psd1'].keys()[0].split()[0].strip()
        ss = PEMSubsystem(subsys, darm_channel,chandict)
        psd1 = FrequencySeries.read(f['psd1'][f['psd1'].keys()[0]])
        for channel in chandict[subsys]:
            N = int(f[subsys][channel]['metadata'][0,1])
            st = int(f[subsys][channel]['metadata'][1,1])
            et = int(f[subsys][channel]['metadata'][2,1])
            csd12 = FrequencySeries.read(f[subsys][channel]['csd mean'])
            psd2 = FrequencySeries.read(f[subsys][channel]['%s mean' % channel])
            ss[channel] = PEMCoherenceSegment(darm_channel, channel,
                                              csd12, psd1, psd2, N, st, et)
        ss.failed_channels = f['failed_channels'][:]
        f.close()
        return ss


    def plot(self, plotSNR=True, **kwargs):
        """
        Return a plot of the coherence matrix
        """
        # initialize coherence matrix
        freqs = self[self.keys()[0]].psd1.frequencies.value
        max_freq = np.max(
                        np.asarray([self[key].psd2.frequencies.value[-1] for key in self.keys()]))
        coh_matrix = np.zeros((len(self.keys()),
                               self[self.keys()[0]].psd1.size))
        labels = []
        for label in self.keys():
            label = label.replace('_OUT_DQ','').replace('H1:','').replace('H1:','').replace('L1:','').replace('_','\_')
            labels.append(label)
        # fill in coherence matrix
        fig = plt.figure(figsize=(12,6))
        for ii, key in enumerate(self.keys()):
            coh = self[key].get_coh()
            if plotSNR:
                coh_matrix[ii, :coh.size] = coh.value * self[key].N
                plt.pcolormesh(freqs, np.arange(0, len(self.keys())+1), coh_matrix,
                               norm=LogNorm(vmin=1, vmax=1e3), cmap='viridis')
            else:
                coh_matrix[ii, :coh.size] = coh.value
                plt.pcolormesh(freqs, np.arange(0, len(self.keys())+1), coh_matrix,
                               norm=LogNorm(vmin=1e-3, vmax=1), cmap='viridis')

        Navg = self[key].N
        ax = plt.gca()
        plt.yticks(np.arange(1, len(labels) + 1) - 0.5, labels, fontsize=6)

        if plotSNR:
            cbar = plt.colorbar(label='coherence SNR')
            ax.set_xlabel(r'Frequency [Hz]')
        else:
            cbar = plt.colorbar(label='coherence')
            plt.xscale('log')
            ax.set_xlabel(r'Log Frequency [Hz]')

        ax.set_xlim(freqs[0], max_freq)
        ax.set_title(r'%s, %d averages' % (self.subsystem, Navg))
        return plt


    def bandplot(self, flow, fhigh, plotSNR=True, **kwargs):
        """
        Return a plot of the coherence matrix in a
        specified frequency band
        """
        # initialize coherence matrix
        freqs = self[self.keys()[0]].psd1.frequencies.value
        max_freq = np.max(
                        np.asarray([self[key].psd2.frequencies.value[-1] for key in self.keys()]))
        coh_matrix = np.zeros((len(self.keys()),
                               self[self.keys()[0]].psd1.size))
        labels = []
        for label in self.keys():
            label = label.replace('_OUT_DQ','').replace('H1:','').replace('H1:','').replace('L1:','').replace('_','\_')
            labels.append(label)
        # fill in coherence matrix
        fig = plt.figure(figsize=(12,6))
        for ii, key in enumerate(self.keys()):
            coh = self[key].get_coh()
            if plotSNR:
                coh_matrix[ii, :coh.size] = coh.value * self[key].N
                plt.pcolormesh(freqs, np.arange(0, len(self.keys())+1), coh_matrix,
                               norm=LogNorm(vmin=1, vmax=1e3), cmap='viridis')
            else:
                coh_matrix[ii, :coh.size] = coh.value
                plt.pcolormesh(freqs, np.arange(0, len(self.keys())+1), coh_matrix,
                               norm=LogNorm(vmin=1e-3, vmax=1), cmap='viridis')

        Navg = self[key].N
        ax = plt.gca()
        plt.yticks(np.arange(1, len(labels) + 1) - 0.5, labels, fontsize=6)

        if plotSNR:
            cbar = plt.colorbar(label='coherence SNR')
            ax.set_xlabel(r'Frequency [Hz]')
        else:
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
        ax.set_title(r'%s, %d averages' % (self.subsystem, Navg))
        return plt


    def update(self, subsystemseg):
        """
        Update the :py:class:`PEMSubsystem` object
        """
        for channel in self.keys():
            self[channel].update(subsystemseg[channel])
        return

    def update_from_file(self, subsystem, filename):
        """
        Parameters
        ----------
        subsystem : `str`
            Channel subsystem
        filename : `str`
            File to read update info from
        """
        s = PEMSubsystem.read(subsystem, filename)
        self.update(s)

    def coarse_grain(self, flowy=None, deltaFy=None, Ny=None):
        """
        Coarse grain the :py:class:`gwpy.frequencyseries.FrequencySeries`
        objects in the :py:class:`PEMCoherenceSegment`, ie, csd12, psd1, psd2.
        Returns an updated :py:class:`PEMCoherenceSegment` object.

        Parameters
        ----------
        flowy : `float`
            Starting frequency value for the coarse grained
            spectrum. Called *f0* in :py:class:`gwpy.frequencyseries.FrequencySeries`.
            The default value is the minimum allowed,
            :math:`\\text{flowy} = \\text{flowx} +\\frac{1}{2}(\Delta F_y - \Delta F_x)`
        deltaFy : `float`
            Step size of the coarse grained spectrum. The
            default is the minimum allowed,
            :math:`\\text{flowy} = \\text{flowx}`
        Ny : `int`
            Number of steps. The default is the maximum amount
            with the given step size and starting frequency.
        """
        for key in self.keys():
            self[key].coarse_grain(flowy=flowy, deltaFy=deltaFy, Ny=Ny)
        return

class CoherenceCombMatrix(object):
    """
    docstring for CoherenceCombMatrix

    Methods
    -------
    .. py:method:: plot
        Plots the :py:class:`CoherenceCombMatrix` object
    """
    def __init__(self, comb, likelihood, offsets, widths, coh):
        super(CoherenceCombMatrix, self).__init__()
        self.comb = comb
        self.offsets = offsets
        self.widths = widths
        self.coherence = coh
        self.likelihood = likelihood
        self.likelihood[self.likelihood < -1000] = -1000
        print likelihood

    def plot(self):
        """
        Returns a plot of the :py:class:`CoherenceCombMatrix` object
        """
        f, axarr = plt.subplots(3, 1, figsize=(8,12))
        axarr[0].pcolormesh(self.offsets, self.widths, np.abs(self.likelihood.T), cmap='viridis')
        axarr[0].set_xlabel('offsets')
        axarr[0].set_ylabel('widths')
        axarr[1].plot(self.coherence.frequencies.value, self.coherence.value)
        axarr[1].set_yscale('log')
        axarr[1].set_xscale('linear')
        axarr[1].set_ylim(1e-4, 1)
        cts, bins = np.histogram(self.comb.reshape((self.comb.size,1)), bins=500)
        pdf =  cts / np.float(np.sum(cts))
        axarr[2].step(bins[:-1], pdf)
        axarr[2].set_yscale('log')
        return plt

class ChannelDict(OrderedDict):
    """
    docstring for ChannelDict

    Methods
    -------
    .. py:method:: read
        Read from a file
    """
    def __init__(self, list_file):
        super(ChannelDict, self).__init__()
        self.list_file = list_file

    @classmethod
    def read(cls, list_file, maxchans=10):
        """
        Read from a file. Returns a channel `dict`

        Parameters
        ----------
        list_file : `str`
            File list to read dict from
        maxchans : `int`, optional
            Maximum number of channels to read.
            Default set to 10 channels.

        Returns
        -------
        Gives a `dict`
        """
        chandict = ChannelDict(list_file)
        channels = ChannelList.read(list_file)
        for channel in channels:
            # add channel to it's sub group
            try:
                chandict[channel.group].append(channel)
            # if key doesn't exist...add it
            except KeyError:
                chandict[channel.group] = []
                chandict[channel.group].append(channel)
        # Now remake the dict based on maxchans
        chandict2 = {}
        for key in chandict.keys():
            nchans = 0
            ngroups = (len(chandict[key]) / maxchans)
            if not ((len(chandict[key]) % maxchans) == 0):
                ngroups += 1
            for group in range(ngroups):
                newkey = '%s %d' % (key, group+1)
                chandict2[newkey] = []
                totmax = min((group+1)*maxchans, len(chandict[key]))
                for ii in range(group*maxchans,totmax):
                    chandict2[newkey].append(chandict[key][ii])
        return chandict2
