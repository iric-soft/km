
===================================================================
km : a software for RNA-seq investigation using k-mer decomposition
===================================================================
+-------------------------------------------------------------+-----------------------------------------------------------------+-----------------------------------------------------------------------------+
| .. image:: https://img.shields.io/badge/python-2.7-blue.svg | .. image:: https://travis-ci.org/iric-soft/km.svg?branch=master | .. image:: https://codecov.io/gh/iric-soft/km/branch/master/graph/badge.svg |
|    :target: https://www.python.org/download/releases/2.7.6/ |    :target: https://travis-ci.org/iric-soft/km                  |    :target: https://codecov.io/gh/iric-soft/km/                             |
+-------------------------------------------------------------+-----------------------------------------------------------------+-----------------------------------------------------------------------------+

---------
Contents:
---------
* `Introduction`_
* `Requirements`_
* `Install`_
* `Usage`_

  - `General`_
  - `Runing km on a real sample`_

* `Tools`_

  - `find_mutation`_
  - `find_report`_
  - `linear_kmin`_

.. _Introduction: https://github.com/iric-soft/km#introduction
.. _Requirements: https://github.com/iric-soft/km#requirements
.. _Install: https://github.com/iric-soft/km#install
.. _Usage: https://github.com/iric-soft/km#usage
.. _Tools: https://github.com/iric-soft/km#tools

.. _General: https://github.com/iric-soft/km#general
.. _Runing km on a real sample: https://github.com/iric-soft/km#runing-km-on-a-real-sample
.. _find_mutation: https://github.com/iric-soft/km#find_mutation
.. _find_report: https://github.com/iric-soft/km#find_report
.. _linear_kmin: https://github.com/iric-soft/km#linear_kmin

-------------
Introduction:
-------------

This tool was developed to identify and quantify the occurence of single
nucleotide variants, insertions, deletions and duplications in RNA-seq data.  Contrary to most tools that try to report all variants in a complete genome, here we instead propose to focus the analysis on small regions of interest.

Given a reference sequence (typically a few hundred base pairs) around a
known or suspected mutation in a gene of interest, all possible sequences
that can be be created between the two end k-mers according to the
sequenced reads will be reported.  A ratio of variant allele vs WT will be
computed for each possible sequence constructed.

-------------
Requirements:
-------------
* Python 2.7.6 or later
* Jellyfish 2.2 or later **with** Python `bindings`_.
* (Optional) Matplotlib

--------
Install:
--------
Before installing or using km, Jellyfish needs to be installed with Python
`bindings`_. Should you need it, a script is available in the `example`_
folder, to help you to install Jellyfish with bindings. When jellyfish is
installed, km may be executed directly from source code or installed using
the setup script:

.. code:: shell

  $ python setup.py install
  $ km -h

.. _bindings: https://github.com/gmarcais/Jellyfish#binding-to-script-languages

------
Usage:
------

General:
--------

From source:
****************

.. code:: shell

  $ cd [your_km_folder]
  $ python -m km -h

After setup install:
********************

.. code:: shell

  $ km -h

Runing km on a real sample:
---------------------------

In the `example`_ folder you can find a bash script which can help you to
run your first km analysis on a Leucegene sample.

.. _example: https://github.com/iric-soft/km/tree/master/example

------
Tools:
------

Overview of km's tools, for more details see the full documentation `here`_.

.. _here: https://github.com/iric-soft/km/tree/master/km/tools

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
This tool parse find_mutation output to reformat it in tabulated file
more user friendly.

.. code:: shell

  $ km find_report -h
  $ km find_report -t [your_fasta_targetSeq] [find_mutation_output]
  $ km find_mutation [your_fasta_targetSeq] [your_jellyfish_count_table] | km find_report -t [your_fasta_targetSeq]

linear_kmin:
------------

Length of k-mers is a central parameter:

* To produce a linear directed graph from the target sequence.
* To avoid false-positive. `find_mutation`_ shouldn't be use on jellyfish count table build with k<21 bp (we recommand k=31 bp, by default)

linear_kmin tool is design to give you the minimun k length to allow a
decomposition of a target sequence in a linear graph.

.. code:: shell

  $ km linear_kmin -h
  $ km linear_kmin [your_catalog_directory]
