import numpy as np
def get_excluded():
    excluded = [
                'GDS-CALIB_STRAIN',
                'OAF-CAL_DARM_DQ',
                'OAF-CAL_YARM_DQ',
                'OMC-DCPD_A_OUT_DQ',
                'OMC-DCPD_B_OUT_DQ',
                'OMC-DCPD_SUM_OUT_DQ',
                'OMC-DCPD_NORM_OUT_DQ',
                'LSC-DARM_*',
                'LSC-YARM_*',
                'LSC-ASAIR_A_RF45_Q_ERR_DQ',
                'SUS-*LOCK*',
                'SUS-*LOCK*',
                'SUS-ETMX*MASTER*',
                'SUS-ETMY*MASTER*',
                'SUS-ETMX*NOISEMON*',
                'SUS-ETMY*NOISEMON*',
                'SUS-ETMX*VOLTMON*',
                'SUS-ETMY*VOLTMON*',
                'SUS-ETMX*ISCINF*',
                'SUS-ETMY*ISCINF*',
                'SUS-*ESDAMON*',
                'ODC*',
                'TOP*',
                'CDS*',
                'TCS*',
                '*_EXC_*',
                'CAL-DELTAL*',
                'CAL-DARM*',
                '*ODC*',
                'HPI*MON*',
                '*SPARE*',
                '*TEMPERATURE*',
                'LSC-MICH_IN1_DQ',
                'LSC-PRCL_IN1_DQ',
                'LSC-SRCL_IN1_DQ',
                '*_OPLEV_SEG*',
                ]
    return  excluded

def cohe_color(c):
    """
    return a hex color code depending continuously on
    input: 1=red, 0=white. 

    Please note: all credit for this awesome function
    goes to Gabriele Vajente.
    """
    if np.isnan(c):
	c = 0
    if c == 0:
        return '#ffffff'
    else:
        if c == 1:
            return '#ff0000'
        else:
            s = hex(int((1.0 - c) * 256)).split('x')[1]
            if len(s) == 1:
                s = '0' + s
            return '#ff' + s + s



