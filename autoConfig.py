import sys
import os
from subprocess import Popen, PIPE
import argparse


"""
This file is used for automatically setting up the necessary folders and files
needed for stamp-pem. You may set things up manually and ignore this if you
wish. Be sure to use the correct ifo! By default, it is set to Hanford and you
can simply run

$ python autoConfigure.py

If you are installing this at LLO, run the following instead

$ python autoConfigure.py -ifo 'L1'
"""


def parse_command_line():
    """
    parse_command_line
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--interferometer", "-ifo", help="Interferometer (H1/L1)",
        default='H1', type=str, dest='ifo')

    params = parser.parse_args()
    return params


# Get account user (could be for a group account)
(stdout, stderr) = Popen(["whoami"], stdout=PIPE).communicate()
user = stdout.strip('\n')
basedir = '/home/' + user + '/'

# Get individual user
individual_user = raw_input('Enter your ligo username (e.g. albert.einstein): ')

# Make sure ligo-channel-lists exists
path_exists = os.path.isdir(basedir + 'ligo-channel-lists')
if path_exists:
    pass
elif not path_exists:
    print('''You do not have the ligo-channel-lists installed. If you already
have an SSH key on this machine, the lists can be installed now. Otherwise,
please generate your SSH key and clone the lists and then run this script again.
Github page: https://git.ligo.org/detchar/ligo-channel-lists''')
    response = raw_input('Install ligo-channel-lists now? (y/n): ')
    while response not in ['y', 'yes', 'n', 'no']:
        print('Invalid entry.')
        response = raw_input('Install ligo-channel-lists now? (y/n): ')

    if response.lower() in ['y', 'yes']:
        try:
            os.system('git clone git@git.ligo.org:detchar/ligo-channel-lists.git {0}'.format(basedir + 'ligo-channel-lists'))
        except:
            print('Unable to clone repository')
    elif response.lower() in ['n', 'no']:
        sys.exit()

# Get the interferometer
params = parse_command_line()
ifo = params.ifo
print('Interferometer set to {0}'.format(ifo))

# Set up the directory structure
print('Setting up directory structure for {0}'.format(user))
dir1 = basedir + 'config_files/ini_files/'
dir2 = basedir + 'config_files/time_files/'
dir3 = basedir + 'config_files/recipients/'
dir4 = basedir + 'stamppemtest/'
dir5 = basedir + 'query_database/'
os.system('mkdir -p {0} {1} {2} {3} {4}'.format(dir1, dir2, dir3, dir4, dir5))

# Write the necessary templates
# Set up the config_file
config_file = """[env]
base_directory = /home/{0}/stamppemtest
output_directory = /home/{0}/query_database
list = /home/{0}/ligo-channel-lists/O2/{1}-O2-standard.ini
accounting_user = {0}
accounting_tag = ligo.prod.o2.detchar.syswide_coh.stamp_pem
user = {2}
executable = /home/{0}/opt/stamp_pem_soft/bin/stamp-pem-pipeline
combine_executable = /home/{0}/opt/stamp_pem_soft/bin/combine_coherence
online_time_file = /home/{0}/config_files/time_files/times.txt
online_increment = 7200
job_duration = 1800
recipients = /home/{0}/config_files/recipients/recipients.txt

[run]
darm_channel = {1}:GDS-CALIB_STRAIN
flag = {1}:DMT-ANALYSIS_READY:1
stride = 10
segmentDuration = 10
fhigh = 16384
resamplerate1 = 2048
resamplerate2 = 2048
subsystems = all
reference_st = 1171936818
plotsnr = False
""".format(user, ifo, individual_user)
fname = dir1 + '{0}.ini'.format(ifo)
with open(fname, 'w') as f:
    f.write(config_file)

# Set up the time file (set to Feb 5, 2017)
time_file = "1170288018"
fname = dir2 + 'times.txt'
with open(fname, 'w') as f:
    f.write(time_file)

# Set up the recipients file
print('Automatic job reports can be emailed to a list of recipients')
print('Enter recipients now (as user@example.com). Type "q" when done...')
recipients = []
recipient = ''
num = 1
while recipient != 'q':
    recipient = raw_input('email address #{0}: '.format(num))
    recipients.append(recipient)
    num += 1
fname = dir3 + 'recipients.txt'
with open(fname, 'w') as f:
    for name in recipients[:-1]:
        f.write(name + '\n')

print('''If you made any errors, you can run "rm -r config_files stamppemtest query_database"
to remove everything just created and then run the AutoConfig script again, or, make
the edits manually''')
print('\nConfiguration complete')
