.. stamppem documentation master file, created by
   sphinx-quickstart on Tue Jun  7 17:32:16 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to stamppem's documentation!
====================================

``stamp-pem`` is an online coherence calculator used to try to identify sources of noise in the LIGO detector. On this page you'll find information on how to run coherences offline, how to query archived data created from the online pipeline, and how to set-up your own version of the pipeline.

The recommended use for ``stamp-pem`` is on the CIT, LHO or LLO clusters. At ``CIT`` and ``LHO`` right now you can source a working virtual environment here: ``/home/meyers/opt/stamp_pem_soft/bin/activate``. After this you should be able to easily run any of the examples provided in the tutorial without issue.

A sample output result from selective querying og the stamp-pem coherence data can be found `here <https://ldas-jobs.ligo-wa.caltech.edu/~rich.ormiston/Stamp-PEM/SelectiveQueryExample/HTML/day/20170123/>`_. The commands
to generate this page are covered in the "Selective Query" section.

For updates, releases and version info, be sure to check the `Announcements <https://ldas-jobs.ligo.caltech.edu/~rich.ormiston/stamp-pem/docs/_build/html/announcements/announcements.html>`_ page.

Contents:

.. toctree::
   :maxdepth: 2

   tutorial/tutorial
   ve/ve
   coherence/coherence
   online_pipeline/online_pipeline
   online_pipeline/query_online
   selective_query/selective_query
   io/io
   announcements/announcements


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`

