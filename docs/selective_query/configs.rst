Config File
===========
It is necessary to include a params file in the command line arguments. This is
done using the ``-i`` flag. A sample params file is shown here.

.. code-block:: ini
    :linenos:

    [env]
    ;; base directory where input data is read from
    base_directory = /home/albert.einstein/stamppemtest/
    ;; output directory where results are stored
    ;; this can be the same as the base directory
    output_directory = /home/albert.einstein/output_results/
    ;; channel list
    list = /home/albert.einstein/ligo-channel-lists/O2/H1-O2-reduced.ini
    accounting_user = albert.einstein
    accounting_tag = ligo.prod.o1.detchar.syswide_coh.stamp_pem
    user = albert.einstein
    executable = /home/albert.einstein/opt/stamp_pem_ve/bin/stamp-pem-pipeline
    combine_executable = /home/albert.einstein/opt/stamp_pem_ve/bin/combine_coherence
    ;; end of most recent tiems run
    online_time_file = /home/albert.einstein/config_files/time_files/times.txt
    ;; time we increment with each run
    online_increment = 7200
    ;; time we break online increment down into
    ;; i.e. each online run is broken down into 4 different
    ;; job_duration number segments (for speed and memory 
    ;; purposes
    job_duration = 1800
    recipients = /home/albert.einstein/config_files/recipients/recipients.txt

    [run]
    ;; darm channel
    darm_channel = H1:GDS-CALIB_STRAIN
    ;; lock flag
    flag = H1:DMT-ANALYSIS_READY:1
    ;; stride length for coherence
    stride = 10
    ;; resample rate (for speed and memory)
    resamplerate1 = 2048
    resamplerate2 = 2048
    ;; subsystems that you want to run over. You can override
    ;; this by specifying subsystems from the command line
    subsystems = all
    ;; plot coherence instead of coherence SNR
    plotsnr = False
    ;; reference time for residual calculations
    reference_st = 1171936818

