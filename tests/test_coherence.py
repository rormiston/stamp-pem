from coherence_segment import PEMCoherenceSegment

# Supply the two channels to take the coherence between
ch1 = 'H1:GDS-CALIB_STRAIN'
ch2 = 'H1:LSC-PRCL_IN1_DQ'

# Supply a start and end time (10 minute interval shown here)
start_time = 1152007217
end_time = start_time + 600

# Call the function
coherence_result = PEMCoherenceSegment.coherence(ch1, ch2, st=start_time, et=end_time)

# Create a plot
plot = coherence_result.plot()
ax = plot.gca()
ax.set_xlim(0, 600)
ax.set_ylabel(r'Coherence between two H1 channels')
ax.plot([0,600], [1./coherence_result.N, 1./coherence_result.N])

# To view the plot rather than saving it, uncomment the following
# plot.show()

# To save the figure, 
plot.savefig('H1_Coherence_Example.png')

# Coarse grain the data and plot it
coherence_result.coarse_grain(deltaFy=2)
plot = coherence_result.plot()
ax = plot.gca()
ax.set_xlim(0, 600)
ax.set_ylabel(r'Coherence between two H1 channels')
ax.plot([0,600], [1./(coherence_result.N * 2), 1./(coherence_result.N * 2)])

# To view the plot rather than saving it, uncomment the following
# plot.show()

# To save the figure, 
plot.savefig('H1_Coherence_Example_coarsegrained.png')

