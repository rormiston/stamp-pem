from stamp_pem.coherence_segment import PEMSubsystem
from stamp_pem.coherence_segment import ChannelDict
from gwpy.detector import Channel

# define our darm channel, start time, end time again
ch1 = Channel('L1:GDS-CALIB_STRAIN', frametype='L1_HOFT_C00')
st = 1157850017
et = st + 10
stride = 1

# Read in our channel list into a dict with keys that are
# the subsystems specified by the headings in the ligo-channel-list
# file. You can open this ini file to see them
chandict = ChannelDict.read('L1-O2-reduced.ini', maxchans=10)

# The keys will tell us what subsystems we can choose from
print chandict.keys()

# The subsystems are broken down into subsets of 10 channels
# automatically. This can be changed with the "maxchans" keyword.
# We'll choose a subsystem with only a few channels for right now
subsystem = 'Thermal Compensation 1'


# Now we can run subsystem coherence on one of these subsystems
# it produces as dict with keys as channel names and PEMCoherenceSegments
# as the values for those keys
subcoh = PEMSubsystem.coherence(ch1, subsystem, chandict, st, et, stride=stride)

print subcoh.keys()

# make a line plot with each channel on it
plot = subcoh[subcoh.keys()[0]].get_coh().plot(label=subcoh.keys()[0].replace('_','\_'))
ax = plot.gca()
for key in subcoh.keys():
    ax.plot(subcoh[key].get_coh(), label=key.replace('_','\_'))
ax.set_title('Coherence for %s with %s' %(subsystem, ch1.name.replace('_','\_')), fontsize=10)
ax.set_ylabel(r'$C(f)$')
ax.set_xlabel('Frequency [Hz]')
plot.savefig('subsystem_coherence_line_plot.png')

# Make a coherence matrix for this
# subsystem
plt = subcoh.plot()
# only go up to 1kHz in plot
plt.xlim(0,1024)
plt.savefig('subsystem_coherence_matrix_plot.png')
