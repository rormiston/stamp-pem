Announcements
=============
Below lists the releases (with the most current on top) and relevant info. Feedback, questions and suggestions are warmly welcomed. Feel free to email Pat (patrick.meyers@ligo.org) or Rich (rich.ormiston@ligo.org)

Release 0.0 (Feb. 2017)
-----------------------

Features Added 
^^^^^^^^^^^^^^
* Selective querying of subsystems and gps times 
* Tutorial for selective querying
* Automated job summarization and optional reporting via email
* File and folder auto-configuration scripts (AutoConfig.py in base directory)
* Virtual environment installation instructions and scripts on document page

Bugs Fixed
^^^^^^^^^^
* Pipeline runs smoothly for first run without complaining about data not existing yet

Warnings
^^^^^^^^
matplotlib/colors.py and atropy/units/quantity.py are suddenly rasing a RuntimeWarning. This does NOT adversely affect the data, output or pipeline. 


Release 0.1 (June. 2017)
------------------------

Features Added 
^^^^^^^^^^^^^^
* Generation of residual coherence plots
* Interactive coherence plots sortable by frequency and/or channel
* Interactive residual coherence plots sortable by frequency and/or channel
* Simple querying of subsystems by keyword
* Adjustable frequency range to query
* Adjustable coarse grain frequency bin width
* Updated autoConfig script

Bugs Fixed
^^^^^^^^^^
* "TypeError: Object dtype dtype('O') has no native HDF5 equivalent" - Updated ``LIGOTimeGPS`` to ``gpsSeconds``
* Numpy IndexError in ``stamp_pem.coherence_functions.fftgram``
* Matplotlib RuntimeWarning in ``selective-query-combine-coherence``

