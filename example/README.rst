
===================================================================
km : a software for RNA-seq investigation using k-mer decomposition
===================================================================
Here you can found some script that you can use as an example,
that you can use and modify to fit your aims

---------------
all_install.sh:
---------------

Description:
------------
This script will be install jellyfish with python binding and km in a virtual
environement. All the code source will be downloaded in your
$HOME/software directory.

Need:
-----
This script need to have `pip`_ and `virtualenv`_ installed.

.. _pip: https://pip.pypa.io/en/stable/installing/
.. _virtualenv: https://virtualenv.pypa.io/en/stable/installation/

Usage:
------
Copy/past this script in a terminal.

The virtual environment need to loaded each time you open a new terminal, with this command:

.. code:: shell

  $ source $HOME/.virtualenvs/km/bin/activate

-----------------
run_leucegene.sh:
-----------------

Description:
------------

This script will:

* download all fastq files of sample 03H041 from `Leucegene`_
* Create the a k-mer count table of 31bp length from fastq downloaded with jellyfish
* made a km analysis on the 7 mutations from the `catalog`_

.. _Leucegene: https://leucegene.ca/
.. _catalog: https://github.com/iric-soft/km/tree/master/catalog/GRCh38

Need:
-----
This script is designed to run on a computer with 8 Go and 4 thread and need
46GB available space. Also km need to be installed and directly accessible
with the command "*km*".

To go futher:
-------------
To optimize this script on computer with less or more 8Go of RAM, you need to
change "**-s**" parameter of "*jellyfish count*" command with this formule:

 * -s = (0.50 * (8 * 1073741824 * [RAM]) / ([k_len] + [-c]))
