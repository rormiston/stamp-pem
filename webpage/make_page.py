import markup
import optparse
import glob
import os
import numpy as np
import sys
import glob
from gwpy.segments import DataQualityFlag
from stamp_pem import coh_io


def parse_command_line():
    """
    parse command parse command line
    """
    parser = optparse.OptionParser()
    parser.add_option(
        "--post-proc-dir", "-d", help="Post processing files directory",
        default=None, type=str, dest='directory')
    parser.add_option(
        "--start-time", "-s", help="start time",
        type=int, dest="st")
    parser.add_option(
        "--end-time", "-e", help="end time",
        type=str, dest="et")
    parser.add_option(
        "--ini-file", "-i", help="ini file for stamp pem",
        type=str, dest="ini", default=None)

    params, args = parser.parse_args()
    return params


def readJobsFile(jobsFile):
    f = open(jobsFile, 'r')
    nums = []
    starts = []
    ends = []
    durs = []
    for line in f:
        s = line.split()
        nums.append(int(s[0].rstrip()))
        starts.append(int(s[1].rstrip()))
        ends.append(int(s[2].rstrip()))
        durs.append(int(s[3].rstrip()))
    f.close()
    return starts, ends, durs


def make_title(page, title, description):
    page.div(class_="panel panel-primary")
    page.div(class_="panel-heading")
    page.h1(title)
    page.div.close()
    page.div(class_="panel-body")
    page.p(description, size='3em')
    page.div.close()
    page.div.close()
    return page


def readNotesFile(notesFile):
    f = open(notesFile, 'r')
    notes = []
    for line in f:
        if line.rstrip():
            notes.append(line.rstrip())
    f.close()
    return notes


def make_index_page(page_name, params):
    starts, ends, durs = readJobsFile(params.jobsFile)
    parameters, freqsToRemove, nBinsToRemove = readParamsFile(params)
    page = page_init('Stochastic Analysis Page')
    title = 'Stochastic Analysis Overview!'
    description = 'Stochastic run from %d-%d. This page contains most information regarding the whole run.' % (
        starts[0], ends[-1])
    page = make_title(page, title, description)
    page.div(class_="panel panel-primary")
    page.div(class_="panel-heading")
    page.h3("Overview")
    page.div.close()
    page.div(class_="panel-body")
    if parameters['doShift1'].rstrip() == 'true':
        page.div(class_="well well-sm")
        page.p("We ran with a time shift of %d second(s)" %
               int(parameters['ShiftTime1'].rstrip()))
        page.div.close()
    else:
        page.div(class_="well well-sm")
        page.p("We ran a zero-lag analysis")
        page.div.close()
    njobs = len(starts)
    totTime = np.sum(np.asarray(durs)) / 86400.
    page.div(class_="well well-sm")
    page.p("There were %d jobs totaling approximately %4.4f days" %
           (njobs, totTime))
    page.div.close()
    page.div(class_="well well-sm")
    page.p("The final integrated energy density was found to be: ")
    page.div.close()
    page.div(class_="well well-sm")
    page.p("The SNR of this was found to be: ")
    page.div.close()
    page.div(class_="well well-sm")
    page.p("The rough upper limit would come out to: ")
    page.div.close()
    page.div.close()
    page.div.close()
    if params.notes is not None:
        notes = readNotesFile(params.notes)
        page.div(class_="panel panel-primary")
        page.div(class_="panel-heading")
        page.h3("Additional notes")
        page.div.close()
        page.div(class_="panel-body")
        for note in notes:
            page.div(class_="well well-sm")
            page.p(note)
            page.div.close()
        page.div.close()
        page.div.close()
    page_end(page, page_name)


def addNotes(page, params):
    starts, ends, durs = readJobsFile(params.jobsFile)
    parameters, freqsToRemove, nBinsToRemove = readParamsFile(params)
    page.div(class_="col-sm-4")
    page.div(class_="panel panel-primary")
    page.div(class_="panel-heading")
    page.h3("Overview %d-%d" % (starts[0], ends[-1]))
    page.div.close()  # heading
    page.div(class_="panel-body")
    if parameters['doShift1'].rstrip() == 'true':
        page.div(class_="well well-sm")
        page.p("We ran with a time shift of %d second(s)" %
               int(parameters['ShiftTime1'].rstrip()), id="sidenote")
        page.div.close()  # well
    else:
        page.div(class_="well well-sm")
        page.p("We ran a zero-lag analysis", id="sidenote")
        page.div.close()  # well
    njobs = len(starts)
    totTime = np.sum(np.asarray(durs)) / 86400.
    page.div(class_="well well-sm")
    page.p("There were %d jobs totaling approximately %4.4f days" %
           (njobs, totTime), id="sidenote")
    page.div.close()  # well
    page.div(class_="well well-sm")
    page.p(
        "The final integrated energy density was found to be: ", id="sidenote")
    page.div.close()  # well
    page.div(class_="well well-sm")
    page.p("The SNR of this was found to be: ", id="sidenote")
    page.div.close()  # well
    page.div(class_="well well-sm")
    page.p("The rough upper limit would come out to: ", id="sidenote")
    page.div.close()  # well
    page.div.close()  # body
    page.div.close()  # primary
    if params.notes is not None:
        notes = readNotesFile(params.notes)
        page.div(class_="panel panel-primary")
        page.div(class_="panel-heading")
        page.h3("Additional notes")
        page.div.close()  # heading
        page.div(class_="panel-body")
        for note in notes:
            page.div(class_="well well-sm")
            page.p(note, id="sidenote")
            page.div.close()  # well
        page.div.close()  # body
        page.div.close()  # primary
    page.div.close()  # col
    return page


def make_diagnostic_page(page_name, params):
    starts, ends, durs = readJobsFile(params.jobsFile)
    parameters, freqsToRemove, nBinsToRemove = readParamsFile(params)
    page = page_init('Stochastic Analysis Diagnostic Plots')
    # title
    page.div(classs_="row")
    page.div(class_="panel panel-primary")
    page.div(class_="panel-heading")
    page.h1("Diagnostic Plots")
    page.div.close()  # panel-heading
    page.div(class_="panel-body")
    page.p("Useful plots for stochastic analysis")
    page.div.close()  # panel-body
    page.div.close()  # panel-primary
    page.div.close()  # panel-row
    # plots
    page.div(class_="row")
    page.div(class_="col-sm-8")
    page.div(class_="panel panel-primary")
    page.div(class_="panel-heading")
    page.h2("Plots")
    page.div.close()  # panel-heading
    page.div(class_="panel-body")
    for plot in plots:
        name = plot.split('/')[-1]
        description = getPlotDescription(name)
        page.div(class_="panel panel-success")
        page.div(class_="panel-heading")
        page.p(description, size='2em')
        page.div.close()  # panel-heading
        page.div(class_="panel-body")
        os.system('cp %s webpage/' % plot)
        page.img(src=name, width='100%')
        page.div.close()  # panel-body
        page.div.close()  # panel-success
    page.div.close()  # panel-body
    page.div.close()  # panel-primary
    page.div.close()  # col-sm-9
    # make col of size 4 with notes and overview in it.
    page = addNotes(page, params)
    page.div.close()  # row!
    page_end(page, page_name)


def getPlotDescription(plot):
    name_temp = plot.split('.')[0]
    name_temp = name_temp.split('_')
    name = name_temp[1]
    if name == "ptEstIntegrand":
        name = name + '_' + name_temp[2]
    if name == "sensIntegrand" and len(name_temp) == 3:
        name = name + '_' + name_temp[2]
    return{
        'PanelPlot1': '$\Delta\sigma$ cut comparison',
        'PanelPlot2': '$\Delta\sigma$ cut comparison',
        'PanelPlot4': '$(1 - \\frac{\sigma}{\sigma_{naive}})$ vs. $\sigma$ scatter plot',
        'PanelPlot5': '$(1 - \\frac{\sigma}{\sigma_{naive}})$ and residuals scatter',
        'PanelPlot6': 'Residuals with and without cut',
        'PanelPlot7': 'KS test optimization statistic vs. Sigma scaling factor',
        'PanelPlot8': 'KS test information',
        'PanelPlot9': 'General residuals information',
        'PanelPlot11': '$ \\frac{\sigma^2}{<\sigma^2>} $ histogram',
        'PanelPlot12': '$\Omega_i$ for each segment.',
        'runningPointEstimate': 'Running point estimate',
        'runningSigma': 'Running sensitivity plot',
        'runningPtEstIntegrand': 'Segment by segment point estimate integrand',
        'ptEstIntegrand_real': 'Real part of point estimate spectrum',
        'ptEstIntegrand_imag': 'Imaginary part of point estimate spectrum',
        'ptEstIntegrand_cum': 'Cumulative point estimate spectrum',
        'ptEstIntegrand_abs': 'Absolute value of point estimate spectrum',
        'sensIntegrand': 'Sensitivity integrand spectrum',
        'sensIntegrand_cum': 'Cumulative integrand spectrum'

    }.get(name, '')


def readParamsFile(params):
    f2 = open(params.paramsFile, "r")
    parameters = {}
    freqsToRemove = []
    nBinsToRemove = []
    for line in f2:
        if line[0] == '%':
            continue
        if line.rstrip() == '':
            continue
        if line.split(' ')[0].rstrip() == "freqsToRemove":
            freqs = line.split(' ')[1].rstrip()
            freqsToRemove = freqs.split(',')
            continue
        if line.split(' ')[0].rstrip() == "nBinsToRemove":
            nbinsTot = line.split(' ')[1].rstrip()
            nBinsToRemove = nbinsTot.split(',')
            continue
        try:
            parameters[line.split(' ')[0]] = line.split(' ')[1]
        except IndexError:
            print line
            raise IndexError('problem')
    f2.close()
    return parameters, freqsToRemove, nBinsToRemove


def page_init(title):
    page = markup.page()
    page.init(css=("main.css",
                   "https://maxcdn.bootstrapcdn.com/bootstrap/3.3.2/css/bootstrap.min.css"),
              title=title
              )
    page.div(class_="navbar navbar-custom")
    page.div(class_="container")
    page.a("Home", id="pat",
           class_="navbar-brand", href="index.html")
    page.a("    ", id="pat", class_="navbar-brand", href="#")
    page.a("", id="pat",
           class_="navbar-brand", href="coherence_matrices.html")
    page.div.close()
    page.div.close()
    page.div(class_="container")
    return page


def page_end(page, filename):
    page.div(class_="row")
    page.div(class_="panel panel-primary")
    page.div(class_="panel-heading")
    page.p("How did I make this page?")
    page.div.close()
    page.div(class_="panel-body")
    page.pre(class_="prettyprint lang-bash")
    page.p("python make_page.py -d %s -p %s -j %s -n %s" %
           (params.directory, params.paramsFile, params.jobsFile, params.notes))
    page.pre.close()
    page.div.close()
    page.div.close()
    page.div.close()
    page.div.close()
    page.script(
        src="https://cdn.rawgit.com/google/code-prettify/master/loader/run_prettify.js")
    page.script.close()
    page.script(type="text/javascript",
                src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-MML-AM_CHTML")
    page.script.close()
    page.script(
        "MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}});", type="text/x-mathjax-config")
    page.script.close()
    page.script('MathJax.Hub.Config({\
                              config: ["MMLorHTML.js"],\
                              jax: ["input/TeX","input/MathML","output/HTML-CSS","output/NativeMML"],\
                              extensions: ["tex2jax.js","mml2jax.js","MathMenu.js","MathZoom.js"],\
                              TeX: {\
                              extensions: ["AMSmath.js","AMSsymbols.js","noErrors.js","noUndefined.js"]\
                                }});')
    page.script.close()
    f = open(filename, 'w')

    print >> f, page
    f.close()


def make_params_page(params_page, params):
    parameters, freqsToRemove, nBinsToRemove = readParamsFile(params)
    # NAVBAR
    page = page_init('Stochastic Analysis Params Page')
    # title
    page.div(class_="row")
    page.div(class_="panel panel-primary")
    page.div(class_="panel-heading")
    page.h1("Parameters")
    page.div.close()
    page.div(class_="panel-body")
    page.p("Parameters for analysis that was run!")
    page.p("In alphabetical order")
    page.div.close()
    page.div.close()
    page.div.close()
    # param table
    page.div(class_="row")
    page.div(class_="col-sm-8")
    page.div(class_="panel panel-primary")
    page.div(class_="panel-heading")
    page.h1("Parameters")
    page.div.close()
    page.div(class_="panel-body")
    page.table(class_="table")
    # params
    keys = sorted(parameters.keys())
    page.tr()
    page.td('Parameter')
    page.td('Value')
    page.tr.close()
    for key in keys:
        page.tr()
        page.td(key)
        page.td(parameters[key])
        page.tr.close()
    page.table.close()
    page.div.close()
    page.div.close()
    page.div.close()
    page = addNotes(page, params)
    page.div.close()
    page_end(page, params_page)


params = parse_command_line()

st_dir = int(str(params.st)[0:5])
et_dir = int(str(params.et)[0:5])
dirs = np.arange(st_dir, et_dir + 1)
pipeline_dict = coh_io.read_pipeline_ini(params.ini)
env_params, run_params = coh_io.check_ini_params(pipeline_dict)
segs = []
for directory in dirs:
    seg_files = sorted(glob.glob('%s/SEGMENTS/%d' %
                                 (env_params['base_directory'], directory)))
    for f in seg_files:
        temps = DataQualityFlag.read(f)
        segs.append(temps)



plots = sorted(glob.glob('%s/*.png' % params.directory))
notches_page = 'webpage/notches.html'
params_page = 'webpage/params.html'
diagnostic_plots_page = 'webpage/diagnostic_plots.html'
index_page = 'webpage/index.html'
os.system('mkdir -p webpage')
os.system('touch %s' % notches_page)
os.system('touch %s' % params_page)
os.system('touch %s' % diagnostic_plots_page)
os.system('touch %s' % index_page)
os.system('cp main.css webpage/')
make_index_page(index_page, params)
make_diagnostic_page(diagnostic_plots_page, params)
make_params_page(params_page, params)
make_notches_page(notches_page, params)
