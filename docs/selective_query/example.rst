Example Usage
=============
In order to query all "physical" subsystems (microphones, accelerometers,
magnetometers, voltage monitors, radio receivers and seismometers) and
"suspension" systems (noise monitors, OSEMs) in the range from 20-120Hz,
coarse grained to 0.1Hz and in the time interval from 1169181018 to
1169184618 we would run:

.. code-block:: bash
   :linenos:

   selective-query -i /home/albert.einstein/config_files/ini_files/H1.ini -st 1169181018 -et 1169184618 --subsystems "Physical" "Suspension" --new-df 0.1 -f 20 120

The ``selective-query`` routine calls four different scripts:

1. ``selective-query-combine-coherence`` 
2. ``selective-query-webpage``
3. ``selective-query-d3plots``
4. ``selective-query-residuals``

Any of these may be run on their own using the same command line flags. In
this way for example, the interactive D3 plots can be generated without
generating the residual coherence matrix plots or the webpage command. A webpage
will still be generated for just these plots.
The output wepages of the above command can be found `here <https://ldas-jobs.ligo-wa.caltech.edu/~rich.ormiston/Stamp-PEM/SelectiveQueryExample/HTML/day/20170123/>`_
