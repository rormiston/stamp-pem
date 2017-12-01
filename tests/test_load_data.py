import stamp_pem.coherence_functions as cf

channel1 = 'H1:GDS-CALIB_STRAIN'
channel2 = 'H1:IMC-F_OUT_DQ'
st = 1135641617
dur = 10
chan1 = cf._read_data(channel1, st, st+dur)
chan2 = cf._read_data(channel2, st, st+dur)
coh, N = cf.coherence(chan1, chan2, 1)
print coh
print N


