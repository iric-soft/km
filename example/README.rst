
===================================================================
km examples script
===================================================================
These scripts are designed to work on general case, but you could have to
modify some part to fit with your settings.

---------
Contents:
---------
* `run_leucegene`_ : To run your first km analysis

.. _run_leucegene: https://github.com/iric-soft/km/tree/master/example#run_leucegene

--------------
run_leucegene:
--------------

Description:
------------

This script will:

* Download all fastq files of sample 03H041 of `Leucegene`_ from `GEO`_.
* Create a k-mer count table of 31bp length from fastq files downloaded, using Jellyfish.
* Run km to annotate 03H041 on the 9 mutations from the `catalog`_.

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
