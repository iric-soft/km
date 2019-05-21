
===================================================================
km examples script
===================================================================
These scripts are designed to work on general case, but you could have to
modify some part to fit with your settings.

--------------
run_leucegene:
--------------

Description:
------------

This script run a km analysis from fastq of one `Leucegene`_ sample. 
It will help you to:

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
This script is optimize to run with 8Go of RAM. 
In other case, you need to change **-s** parameter of *jellyfish count*. 
To found the value who fit with your setting, use this command:

 .. code:: shell

   $ jellyfish mem -m 31 -s 799063683 -c 12
   7158303656 (6G)
