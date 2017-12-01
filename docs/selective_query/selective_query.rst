+++++++++++++++
Selective Query
+++++++++++++++

When the online pipeline has been run and coherence data has been written, we
may then parse this data and choose how to do so. Several options are availble,
among them:

* Segment duration to query
* Subsystem selection
* Coarse grain to new frequency bin width
* Frequency band specification

By default, results will be stored in the ``output_directory`` given in the
config file. To override this, use ``-o path/to/output``.

Output Results
--------------
The `output <https://ldas-jobs.ligo-wa.caltech.edu/~rich.ormiston/Stamp-PEM/SelectiveQueryExample/HTML/day/20170123/>`_
from running the ``selective-query`` pipeline is essentially the same as
running the online pipeline. Namely, we produce combined coherence data files
and a webpage containing the following:

* Full output coherence matrix of all queried subsystems
* Bruco-style table
* Coherence matrices of each subsystem
* Residual coherence matrices of each subsystem 
* Sortable coherence matrices (by frequency and channel)
* Sortable residual coherence matrices (by frequency and channel)

The difference is that we most likely queried a specific frequency range, 
subsystems, times and chose a different resolution for our search.
The interactive residual plots show the absolute value of the difference in
coherence between channels of a subsystem at different times per frequency.
Therefore, the largest coherences represent the *most changed* frequencies in
the subsystem between those times. The interactive plots allow one to sort by
channel and see the frequencies that have changed the most, or to sort by
frequency and see the channels that have changed the most.

.. warning::
   If javascript is not enabled in your browser, you will *not* be able to use
   the interactive D3 plots! Make sure that javascript *is* enabled.


Details:
--------

.. toctree::
   :maxdepth: 2

   flags
   configs
   example
   newtonian

