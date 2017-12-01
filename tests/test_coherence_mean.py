from coherence_segment import PEMCoherenceSegment
import numpy as np
import coherence_functions as cf
from gwpy.timeseries import TimeSeries
from gwpy.detector import Channel

# Supply the two channels to take the coherence between
ch1 = 'H1:GDS-CALIB_STRAIN'
ch2 = 'H1:LSC-PRCL_IN1_DQ'

# Supply a start and end time (10 minute interval shown here)
start_time = 1152007217
end_time = start_time + 60
data1 = TimeSeries(np.random.randn(1638400),sample_rate=16384)
data2 = TimeSeries(np.random.randn(1638400),sample_rate=16384)
coh_gwpy = data1.coherence(data2)
print coh_gwpy.median()

# make them channel objects
ch1 = Channel(ch1, frametype='H1_HOFT_C00')
ch2 = Channel(ch2, frametype='H1_R')

# Call the function
coherence_result = PEMCoherenceSegment.coherence(ch1, ch2, st=start_time, et=end_time, resamplerate1=512,  resamplerate2=512,stride=1)
print "number of averages",coherence_result.N
print "mean of gaussian",coherence_result.get_coh().mean()
coh,N = cf.coherence(TimeSeries(np.random.randn(1638400),sample_rate=16384),TimeSeries(np.random.randn(1638400),sample_rate=16384),4,overlap=0,pad=True)
print N
print coh.mean()
print coh



coherence_mean = np.median(coherence_result.get_coh())

print(''''median: {0}
N: {1}'''.format(coherence_mean, coherence_result.N))

# Create a plot
plot = coherence_result.plot()
ax = plot.gca()
ax.set_xlim(0, 256)
ax.set_ylabel(r'Coherence between two H1 channels')
ax.plot([0,256], [1./coherence_result.N, 1./coherence_result.N])

plot.savefig('coherence.png')

# Coarse grain the data
coherence_result.coarse_grain(deltaFy=2)
coherence_cg_mean = np.median(coherence_result.get_coh())
print(''''median: {0}
N: {1}'''.format(coherence_cg_mean, coherence_result.N))

# Create a plot
plot = coherence_result.plot()
ax = plot.gca()
ax.set_xlim(0, 256)
ax.set_ylabel(r'Coherence between two H1 channels')
ax.plot([0,256], [1./(coherence_result.N * 2), 1./(coherence_result.N * 2)])

plot.savefig('cg_coherence.png')

