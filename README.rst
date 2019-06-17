
===================================================================
km : a software for RNA-seq investigation using k-mer decomposition
===================================================================

+-------------------------------------------------------------+-----------------------------------------------------------------+-----------------------------------------------------------------------------+
| .. image:: https://img.shields.io/badge/python-2.7-blue.svg | .. image:: https://travis-ci.org/iric-soft/km.svg?branch=master | .. image:: https://codecov.io/gh/iric-soft/km/branch/master/graph/badge.svg |
|    :target: https://www.python.org/download/releases/2.7.6/ |    :target: https://travis-ci.org/iric-soft/km                  |    :target: https://codecov.io/gh/iric-soft/km/                             |
+-------------------------------------------------------------+-----------------------------------------------------------------+-----------------------------------------------------------------------------+

-------------
Introduction:
-------------

This tool was developed to identify and quantify the occurence of single
nucleotide variants, insertions, deletions and duplications in RNA-seq data.  Contrary to most tools that try to report all variants in a complete genome, here we instead propose to focus the analysis on small regions of interest.

Given a reference sequence (typically a few hundred base pairs) around a
known or suspected mutation in a gene of interest, all possible sequences
that can be be created between the two end k-mers according to the
sequenced reads will be reported. A ratio of variant allele vs WT will be
computed for each possible sequence constructed.

-------
Citing:
-------
* Target variant detection in leukemia using unaligned RNA-Seq reads. bioRxiv 295808; doi: https://doi.org/10.1101/295808

-------------
Install:
-------------

Easy install:
-------------

`easy_install.sh`_ will install jellyfish with python binding, km in a virtual
environement, and test it. Without modification, all the code source will be
downloaded in your $HOME/software directory and all executable will be available
in the virtual environement directory: $HOME/.virtualenvs/km.

Requirements:
*************
* Python 2.7.6 or later with `pip`_ and `virtualenv`_ installed.
.. _pip: https://pip.pypa.io/en/stable/installing/
.. _virtualenv: https://virtualenv.pypa.io/en/stable/installation/

Usage:
******

* Copy/past each line in a terminal.
* The virtual environment need to loaded each time you open a new terminal, with this command:

.. code:: shell

  $ source $HOME/.virtualenvs/km/bin/activate

Test:
*****

.. code:: shell

  $ km find_mutation ./data/catalog/GRCh38/NPM1_4ins_exons_10-11utr.fa ./data/jf/02H025_NPM1.jf | km find_report -t ./data/catalog/GRCh38/NPM1_4ins_exons_10-11utr.fa
  Sample	Region	Location	Type	Removed	Added	Abnormal	Normal	Ratio	Min_coverage	Exclu_min_cov	Variant	Target	InfoVariant_sequence	Reference_sequence
  ./data/jf/02H025_NPM1.jf	chr5:171410540-171410543	chr5:171410544	ITD	0	4 | 4	2870.6	3055.2	0.484	2428		/TCTG	NPM1_4ins_exons_10-11utr	vs_ref	AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA	AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA
  ./data/jf/02H025_NPM1.jf		-	Reference	0	0	0.0	2379.0	1.000	2379		-	NPM1_4ins_exons_10-11utr	vs_ref	

.. _easy_install.sh: https://github.com/iric-soft/km/blob/master/easy_install.sh

Setup install:
--------------

If you have already installed Jellyfish with Python `bindings`_, you can install km using setup.py.

Requirements:
*************
* Python 2.7.6 or later
* Jellyfish 2.2 or later **with** Python `bindings`_.
* (Optional) Matplotlib

Usage:
******

.. code:: shell

  $ python setup.py install


Without install:
----------------
km can be executed directly from source code.

Requirements:
*************
* Python 2.7.6 or later
* Jellyfish 2.2 or later **with** Python `bindings`_.
* (Optional) Matplotlib

Usage:
******

.. code:: shell

  $ cd [your_km_folder]
  $ python -m km find_mutation ./data/catalog/GRCh38/NPM1_4ins_exons_10-11utr.fa ./data/jf/02H025_NPM1.jf | km find_report -t ./data/catalog/GRCh38/NPM1_4ins_exons_10-11utr.fa

.. _bindings: https://github.com/gmarcais/Jellyfish#binding-to-script-languages

----------------------------
Design your target sequence:
----------------------------
* km is design to made targeted analysis based on **target sequences**. These target sequences **need to be design** and given as km's input.
* A target sequence is a nucleotide sequence saved in a fasta file. Some target sequences are provide in `catalog <https://github.com/iric-soft/km/tree/master/km/data/catalog>`_.
* To feet your specific needs, you will have to create your own target sequences. 
* On generic cases, you can follow some good practices describe below:

.. image:: https://github.com/iric-soft/km/blob/master/data/figure/doc_target_sequence.pdf

* There are different methods to extract nucleotide sequences from genome, if needed two of them are discribe below:

   - samtools faidx chr2:25234341-25234405 GRCh38/genome.fa
   - `ucsc <https://genome.ucsc.edu/cgi-bin/hgc?hgsid=730614743_K2u5W9UIMXrPzrUlC5KaXmWjzf4R&o=25234340&g=getDna&i=mixed&c=chr2&l=25234340&r=25234405&db=hg38&hgsid=730614743_K2u5W9UIMXrPzrUlC5KaXmWjzf4R>`_.


-------------
Display help:
-------------

.. code:: shell

  $ km -h
    usage: PROG [-h] {find_mutation,find_report,linear_kmin,min_cov} ...
  
    positional arguments:
      {find_mutation,find_report,linear_kmin,min_cov}
                            sub-command help
        find_mutation       Identify and quantify mutations from a target sequence
                            and a k-mer database.
        find_report         Parse find_mutation output to reformat it in tabulated
                            file more user friendly.
        linear_kmin         Find min k length to decompose a target sequence in a
                            linear graph.
        min_cov             Compute coverage of target sequences.
   
    optional arguments:
      -h, --help            show this help message and exit


--------------------
km's tools overview:
--------------------

For more detailed documentation click `here <https://github.com/iric-soft/km/tree/master/km/tools>`_.

find_mutation:
--------------

This is the main tool of km, to identify and quantify mutations from
a target sequence and a k-mer jellyfish database.

.. code:: shell

  $ km find_mutation -h
  $ km find_mutation [your_fasta_targetSeq] [your_jellyfish_count_table]
  $ km find_mutation [your_catalog_directory] [your_jellyfish_count_table]

find_report:
------------
This tool parse find_mutation output to reformat it in more user friendly
tabulated file.

.. code:: shell

  $ km find_report -h
  $ km find_report -t [your_fasta_targetSeq] [find_mutation_output]
  $ km find_mutation [your_fasta_targetSeq] [your_jellyfish_count_table] | km find_report -t [your_fasta_targetSeq]

min_cov:
--------

This tools display some k-mer's coverage stats of a target sequence and a list of jellyfish database.

.. code:: shell

  $ km min_cov -h
  $ km min_cov [your_fasta_targetSeq] [[your_jellyfish_count_table]...]

linear_kmin:
------------

Length of k-mers is a central parameter:

* To produce a linear directed graph from the target sequence.
* To avoid false-positive. find_mutation shouldn't be use on jellyfish count table build with k<21 bp (we recommand k=31 bp, by default)

linear_kmin tool is design to give you the minimun k length to allow a
decomposition of a target sequence in a linear graph.

.. code:: shell

  $ km linear_kmin -h
  $ km linear_kmin [your_catalog_directory]

-------------------------------------------------
Runing km on a real sample from downloaded fastq:
-------------------------------------------------
In the `example`_ folder you can find a script to help you to
run a km analysis on one Leucegene sample.

  .. _example: https://github.com/iric-soft/km/tree/master/example
