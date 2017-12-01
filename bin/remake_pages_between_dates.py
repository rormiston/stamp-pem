from gwpy.time import tconvert, to_gps, from_gps
import optparse
from datetime import timedelta
from stamp_pem.coh_io import read_pipeline_ini
import os
import time

def parse_command_line():
    """
    parse command line
    WARNING: THIS MANUALLY CHANGES THE TIME FILE SPECIFIED IN THE CONFIG FILE
    YOU GIVE. IT"S NOT RECOMMENDED TO USE THIS FUNCTION WITH THE ONLINE CONFIG
    FILE. YOU SHOULD USE A SEPARATE 'CATCHUP' CONFIG FILE.
    """
    parser = optparse.OptionParser()
    parser.add_option("--start-date", "-s",
        help="start date", default=None,
        dest="sd", type=str)
    parser.add_option("--end-date", "-e",
        help="end date", default=None,
        dest="ed", type=str)
    parser.add_option("--ini", "-i",
        help="config file", default=None,
        dest="ini", type=str)
    params, args = parser.parse_args()
    return params

print """Reminder that this changes the time file specified in your config
file...be careful using this function...sleeping for 3 seconds just to let you
kill this if you want"""

time.sleep(3)

params = parse_command_line()
sd = tconvert(to_gps(params.sd))
ed = tconvert(to_gps(params.ed))
time_file = read_pipeline_ini(params.ini)['env']['online_time_file']
print 'time file is: %s' % time_file

while sd <= ed:
    f = open(time_file,'w+')
    f.write(str(to_gps(sd)))
    f.close()
    cmd = 'make_online_webpage -i %s' % params.ini
    os.system(cmd)
    sd += timedelta(days=1)


