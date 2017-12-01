import sys
import os
import re
import subprocess
import smtplib
from email import encoders
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart


def SendEmail(sender, password, recipients, attachments):
    """
    Use this to send an email to a list of recipients
    with a text file attachment.

    Parameters
    ----------
    sender : `str`
        Sender's email address
    password : `str`
        Password for sender's email account
    recipients: `list`
        list of recipients (comma separated strings)

    Example
    -------
    SendEmail('sender@gmail.com', 'secret_pwd', ['user1@example.com', user2@example.com'])
    """

    # Create the enclosing (outer) message
    outer = MIMEMultipart()
    outer['Subject'] = 'Stamp-pem Pipeline: Failed Jobs Report'
    outer['To'] = ', '.join(recipients)
    outer['From'] = "Stamp-pem Team"

    # List of attachments
    # cwd = os.getcwd()
    # attachments = [cwd + '/FailedJobsReport.txt']

    # Add the attachments to the message
    for attachment in attachments:
        try:
            with open(attachment, 'rb') as fp:
                msg = MIMEBase('application', "octet-stream")
                msg.set_payload(fp.read())
            encoders.encode_base64(msg)
            msg.add_header('Content-Disposition', 'attachment', filename=os.path.basename(attachment))
            outer.attach(msg)
        except:
            print("Unable to open one of the attachments. Error: ", sys.exc_info()[0])
            raise

    composed = outer.as_string()

    # Send the email
    try:
        s = smtplib.SMTP('smtp.gmail.com')
        s.ehlo()
        s.starttls()
        s.ehlo()
        s.login(sender, password)
        s.sendmail(sender, recipients, composed)
        s.close()
        print("Email successfully sent to:")
        for recipient in recipients:
            print('  [+] {0}'.format(recipient))
        print('\n')
    except:
        print("Unable to send the email. Error: ", sys.exc_info()[0])
        raise


def CombineReports(file1, file2, st, et, basedir):
    """
    Append the data FROM file2 TO file1. Write this
    to a new file called FailedJobsReport_st.txt where
    st is the gps start time

    Parameters
    ----------
    file1 : `str`
        absolute path to job report file 1
    file2 : `str`
        absolute path to job report file 2
    st : `int`
        gps start time
    et : `int`
        gps end time
    basedir : `str`
        directory given in config file where coherence data
        is stored
    """
    AddTitle(file1, st, et)
    RemoveBlanks(file1)
    f1 = open(file1, 'a')
    AddTitle(file2, st, et)
    RemoveBlanks(file2)
    f2 = open(file2, 'r')
    lines = f2.readlines()
    for line in lines:
        f1.write(line)
    f1.close()
    f2.close()
    newfilename = 'FailedJobsReport_' + str(st) + '.txt'
    newfile = basedir + '/FailedJobs/' + newfilename
    cmd = 'cp {0} {1}'.format(file1, newfile)
    os.system(cmd)
    EmptyReports(file1, file2)


def CheckJobsRemaining():
    """
    This function checks condor to see how many jobs are left. For future
    scheduling purposes, we want to make sure that all of the jobs ran before
    we send the job reports
    """
    cwd = os.getcwd()
    usr = cwd.split('/')[2]
    output = subprocess.Popen(["condor_q", "-submitter", usr], stdout=subprocess.PIPE).communicate()[0]
    output = output.split('\n')[:-1]
    jobline = output[-1].split(' ')
    jobsremaining = int(jobline[0])
    if jobsremaining == 0:
        return 0


def RemoveBlanks(fname):
    """
    Get rid of excess blank lines in the files before
    they are sent. The cron job continues to call the coherence
    functions and will generate blank lines in the reports.

    Parameters
    ----------
    fname : `str`
        absolute path to file we wish to clean
    """
    f = open(fname, 'r')
    lines = f.readlines()
    lines.reverse()
    blank = 0
    for line in lines:
        if len(line) <= 1:
            blank += 1
        else:
            break
    f.close()
    g = open(fname, 'w')
    total = len(lines) - blank
    i = 0
    lines.reverse()
    for line in lines:
        if i <= total:
            g.write(line)
            i += 1


def EmptyReports(file1, file2):
    """
    Empty the job report files for combine_coherence and
    make_online_webpage. This data will have already been
    written to a file called FailedJobsReport_(\d+).txt so we
    can reuse these templates

    Paramters
    ---------
    file1 : `str`
        absolute path to job report file 1
    file2 : `str`
        absolute path to job report file 2

    """
    with open(file1, 'w') as f1:
        f1.write('')
    with open(file2, 'w') as f2:
        f2.write('')


def AddTitle(fname, st, et):
    """
    Write a title to a file

    Parameters
    ----------
    fname : `str`
        file name including absolute path
    title : `str`
        title to be written to file
    st : `int`
        gps start time
    et : `int`
        gps end time
    """
    read_file =  open(fname, 'r')
    lines = read_file.readlines()
    read_file.close()
    write_file = open(fname, 'w')
    title = fname.split('/')[-1].split('.')[0]
    heading = '{0} {1}-{2}\n'.format(title, st, et)
    write_file.write(heading)
    write_file.write('{0}\n'.format('-' * len(heading.strip('\n'))))
    for line in lines:
        write_file.write(line)
    write_file.close()


def GetEmailAddresses(recipientfile):
    """
    Specify email recipients in a text file and read them in from the config
    file using recipients = /path/to/recipient/file.txt

    Parameters
    ----------
    recipientfile : `str`
        Absolute path to recipient text file.
        The file should have one email per line

    Returns
    -------
    recipients : `list`
        list containing strings of recipients
    """
    recipients = []
    with open(recipientfile, 'r') as f:
        addresses = f.readlines()
        for address in addresses:
            recipients.append(address.strip('\n'))
        return recipients


def GetAttachments(basedir):
    """
    Return list of files in path if they match the form
    "(letters)_(numbers).txt." Most likely, dirname = basename + '/FailedJobs/'
    These will be the job reports. Use this to only email reports when the list
    is populated enough

    Parameters
    ----------
    dirname : `str`
        directory name where the jobs reports are located

    Returns
    -------
    filepaths : `list`
        list of all job report files (including the absolute path)
    """

    # get list of files
    filepaths = []
    dirname = basedir + '/FailedJobs/'
    for basename in os.listdir(dirname):
        regex = re.compile(r'^[^s].*_(\d+).txt')
        f = regex.search(basename)
        if f is not None:
            filename = os.path.join(dirname, basename)
            if os.path.isfile(filename):
                filepaths.append(filename)

    return filepaths


def SendReports(sender, password, recipients, attachment_list, frequency=12):
    """
    Send job reports emails where the frequency is the number of reports that
    must be created before they will be sent. Once there are enough
    jobs and the report is sent, the files are prefixed with "sent_".

    Parameters
    ----------
    attachment_list : `list`
        output of `GetAttachments`
    sender : `str`
        Sender's email address
    password : `str`
        Password for sender's email account
    recipients: `list`
        list of recipients (comma separated strings)
    frequency : `int`
        number of reports that must be generated before the email sends.
        By default, frequency=12
    """
    if len(attachment_list) >= frequency:
        print('Preparing to send the job report...')
        SendEmail(sender, password, recipients, attachment_list) 
        # Give emailed files the prefix "sent" so they aren't sent again
        for attachment in attachment_list:
            path = attachment.split('/')
            fname = path[-1]
            newname = '/'.join(path[:-1]) + '/sent_' + fname
            os.system('mv {0} {1}'.format(attachment, newname))

