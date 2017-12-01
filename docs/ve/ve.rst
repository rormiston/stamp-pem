+++++++++++++++++++++++++++
Install stamp-pem yourself
+++++++++++++++++++++++++++

The recommended installation for ``stamp-pem`` is on the LDAS clusters. There are working installations and a working virtualenv at CIT, LHO, and LLO. You can always access those virtualenvs by typing ``source /home/stochastic/opt/stamp_pem_soft/bin/actiavte``. If you'd like to create your own installation, I'd still recommend you do it at one of these places. If not, then please be aware that there are extra things you may have to install to get ``GWpy`` to work. You can find the ``GWpy`` documentation here_.

If you are installing this on an LDAS cluster, then the following script should work. However, you will need to create an ssh key to be able to access the stamp-pem gitlab repository and make sure that it is working on whatever computer you are installing to. You can find instructions for this by going to the ligo gitlab_ repository, clicking on the icon in the upper right, choosing ``settings``, ``ssh keys``, and following the instructions for creating a new ssh key.

``WARNING:`` If code changes are made to ``stamp-pem`` then you should pull from master and then uninstall and reinstall. That would mean sourcing your virtual environment, moving to the top level ``stamp-pem`` directory and running ``git pull origin master && pip uninstall stamp-pem && pip install .``

.. code-block:: bash
    :linenos:

    #! /bin/bash

    ### This is an installation script
    ### It's useful for installing stamp-pem
    ### in a new virtual environment
    ### on one of the LDAS clusters

    # make git repository directory
    mkdir ~/git_repos
    mkdir -p ~/opt/stamp_pem_soft
    cd ~/git_repos
    # clone gwpysoft environment
    # which has nice scripts for installing gwpy
    git clone https://github.com/gwpy/gwpysoft
    cd gwpysoft
    # install gwpyinto a new virtualenv called "stamp_pem_soft"
    # this installation will live in ~/opt/stamp_pem_soft/
    ./gwpysoft-init ~/opt/stamp_pem_soft packages.txt
    cd ../
    # clone the stamp-pem repository
    git clone git@git.ligo.org:patrick-meyers/stamp-pem.git
    cd stamp-pem
    # activate the new virtualenv you just created
    source ~/opt/stamp_pem_soft/bin/activate
    # install stamp-pem into that virtualenv
    pip install -r requirements.txt && pip install .


Once you have ``stamp-pem`` installed, you may set up the necessary folders by running ``AutoConfig.py`` located in the top level directory. If you're deploying this at LLO, be sure to include the flag ``-ifo 'L1'``. The flag will be set to "H1" by default.


.. _here: https://gwpy.github.io
.. _gitlab: https://git.ligo.org/patrick-meyers/stamp-pem/
