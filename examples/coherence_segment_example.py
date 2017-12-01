from stamp_pem.coherence_segment import PEMCoherenceSegment
from gwpy.detector import Channel


# supplying frame types speeds things up considerably
ch1 = Channel('H1:GDS-CALIB_STRAIN', frametype='H1_HOFT_C00')
ch2 = Channel('H1:IMC-F_OUT_DQ', frametype="H1_R")
st = 1157850017
et = st + 10
stride = 1
pad = False

# run coherence
coh = PEMCoherenceSegment.coherence(ch1, ch2, st, et, stride=1)

# print out the coherence spectrum
# we use get_coh() so that we don't
# have to save an extra spectrum or hold it in memory
# unless we want to...the coh object initially just has
# PSDs and the CSD between our channels
print coh.get_coh() # gwpy.FrequencySeries

# plot it!
plot = coh.get_coh().plot() 
ax = plot.gca()
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel(r'$C(f)$')
ax.set_title('Coherence between %s and %s' %
                (ch1.name.replace('_','\_'),
                 ch2.name.replace('_','\_')),
                 fontsize=10)
plot.savefig('coherence_segment_example_plot.png')
