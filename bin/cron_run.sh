#! /bin/bash

unset X509_USER_PROXY
export X509_USER_CERT="/home/meyers/robot_cert/cert.pem"
export X509_USER_KEY="/home/meyers/robot_cert/robot.key.pem"
export LIGO_DATAFIND_SERVER="10.13.5.22:80"
source /home/meyers/opt/stamp_pem_soft/bin/activate
generate-and-write-segs --ini /home/meyers/config_files/ini_files/L1.ini
run-stamp-pem-pipeline --ini /home/meyers/config_files/ini_files/L1.ini
deactivate
