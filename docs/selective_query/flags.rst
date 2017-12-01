Command Line Arguments
======================
    
The output results from running ``selective-query -h`` are as follows and a brief overview
of some these flags is given below:

.. literalinclude:: cmdhelp.txt
   :language: bash


Time Selection
--------------
Start and end times are supplied at GPS time to the command line with the flags
``-st`` and ``-et``. Since the pipeline saves data every half hour (1800s), the time
intervals that may be queried are :math:`n\times 1800` seconds long where
:math:`n = 1,\,2,\,3\ldots` The supplied start and end times will be adjusted
automatically to the nearest 1800s interval and therefore any times may be given.

Coarse Grain
------------
The coarse graining has a max resolution of 0.1Hz and defaults to 0.5Hz if no
input is supplied. To implement a new frequency bin width, use
``-new-df <deltaF>``. Keep in mind that wider binning results in a loss of the
more subtle features. If coarse graining to 0.1Hz, it is recommended that you
also use a smaller frequency band, as the interactive plots may become
extremely large and consequently quite slow. 

Subsystem Selection
-------------------
Specific subsystems may be supplied after the ``--subsystems`` flag
(e.g, ``--subsystems "Thermal Compensation 1"``) but the arguments accepted are
quite flexible and only partial names need to be supplied. For example, in
order to run the magnetometers, it is not necessary to specify the complete name
"Physical Environment Monitor: Magnetometers", but more simply just giving
"Magnetometer" or even "Mag" is enough. As another example, if one wanted to
run all of the "Physical Environment Monitors," then ``--subsystems  "Physical"``
would include all such subsystems. 

It is possible to specify several different subsystems as well by using the
format ``--subsystems sub1 sub2 sub3`` for as many systems as desired. If every
subsystem is to be run, then simply use the keyword "all"
(e.g, ``--subsystems all``).

Frequency Band
--------------
A particular frequency band of interest may be specified using
``-f <flow> <fhigh>``. If no arguments are supplied, then the default frequency
band from 10-500 Hz will be used.
