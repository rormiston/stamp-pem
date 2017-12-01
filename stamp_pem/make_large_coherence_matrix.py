from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
plt.rc('text', usetex=True)


def chandict_to_ordered_dict(channel_dict):
    """
    Convert channel dict to ordered dict
    """
    newkeys = np.copy(channel_dict.keys())
    newkeys.sort()
    ordered_dict = OrderedDict()
    for key in newkeys:
        ordered_dict[key] = np.copy(channel_dict[key])
    return ordered_dict


def plot_coh_matrix(coh_tab, freqs, nchans_per_subsystem, webpage_basedir, N=1,
                    datestr='', new_df=0.5, plotSNR=True):
    """
    plot a coherence matrix given certain pieces of data

    Parameters
    ----------
        coh_tab : TODO
        freqs : TODO
        nchans_per_subsystem : TODO
        webpage_basedir : TODO

    Returns
    -------
    TODO

    """

    # Now we have our full coherence table, we need to plot it
    # get location for y labels
    # number of channels in each subsystem
    ticknums = np.asarray([nchans_per_subsystem[key] for key in nchans_per_subsystem.keys()])
    subsystems = [key for key in nchans_per_subsystem.keys()]
    # half of number of channels in each subsystem
    halves = ticknums / 2.
    # cumulative sum over ticknumbers, but start it from 0
    sumticknums = np.zeros(ticknums.size)
    sumticknums[1:] = np.cumsum(ticknums[:-1])
    # now shift them by half of their original values
    # so if ticknums = [14, 12, 6] and halves = [7, 6, 3]
    # then sumticknums = [0, 14, 26]
    # then it becomes sumticknums = [7, 14+6, 14+12+3]
    sumticknums += halves
    N = 10 ** np.ceil(np.log10(N))
    plt.figure(figsize=(24,12))

    if plotSNR:
        plt.pcolormesh(coh_tab, cmap='viridis', norm=LogNorm(vmin=1, vmax=max(N, 1e3)))
        plt.xlabel(r'Frequency [Hz]')
        cbar = plt.colorbar(label='coherence SNR')
        # plt.xticks(np.arange(0,freqs[-1]+100,100))
    elif not plotSNR:
        plt.xscale('log')
        plt.pcolormesh(coh_tab, cmap='viridis', norm=LogNorm(vmin=1, vmax=1e3))
        plt.xlabel(r'Log Frequency [Hz]')
        cbar = plt.colorbar(label='coherence')

    plt.yticks(sumticknums, np.unique(subsystems), fontsize=6)
    plt.xlim(freqs[0], freqs[-1])
    plt.title(r'Coherence Matrix for %s' % datestr)
    plt.tight_layout()
    plt.savefig('%s/daily_full_coherence_matrix' % webpage_basedir)
