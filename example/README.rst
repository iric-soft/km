
===================================================================
km examples script
===================================================================
These scripts are designed to work on general case, but you could have to
modify some part to fit with your settings.

--------
Scripts:
--------
* `all_install.sh`_ : To install jellyfish and km
* `run_leucegene.sh`_ : To run your first km analysis

.. _all_install.sh: https://github.com/iric-soft/km/tree/master/example#all_installsh
.. _run_leucegene.sh: https://github.com/iric-soft/km/tree/master/example#run_leucegenesh

---------------
all_install.sh:
---------------

Description:
------------
This script will install jellyfish with python binding and km in a virtual
environement. Without modification, all the code source will be downloaded
in your $HOME/software directory and all executable will be available in
the virtual environement directory: $HOME/.virtualenvs/km.

Need:
-----
This script need to have `pip`_ and `virtualenv`_ installed.

.. _pip: https://pip.pypa.io/en/stable/installing/
.. _virtualenv: https://virtualenv.pypa.io/en/stable/installation/

Usage:
------

* Copy/past each line (and modify) in a terminal.
* The virtual environment need to loaded each time you open a new terminal, with this command:

.. code:: shell

  $ source $HOME/.virtualenvs/km/bin/activate

-----------------
run_leucegene.sh:
-----------------

Description:
------------

This script will:

* Download all fastq files of sample 03H041 of `Leucegene`_ from `GEO`_.
* Create a k-mer count table of 31bp length from fastq files downloaded, using Jellyfish.
* Run km to annotate 03H041 on the 7 mutations from the `catalog`_.

.. _Leucegene: https://leucegene.ca/
.. _catalog: https://github.com/iric-soft/km/tree/master/data/catalog/GRCh38
.. _GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1203307

Need:
-----
* `sratoolkit`_ installed to download fastq from GEO
* 8 Go of RAM, 4 thread and 46GB available space (for fastq).
* Km and Jellyfish need to be installed and directly accessible (into your PATH).
* To be execute in the km parent directory, with:

.. code:: shell

  $ ./example/run_leucegene.sh

.. _sratoolkit: https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=show&f=software&m=software&s=software

To go futher:
-------------
To optimize this script on computer with less or more 8Go of RAM, you need to
change **-s** parameter of *jellyfish count*. This command can help to found
the value who fit with your setting:

 .. code:: shell

   $ jellyfish mem -m 31 -s 799063683 -c 12
   7158303656 (6G)
