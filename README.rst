
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

.. _Introduction: https://github.com/iric-soft/km#introduction
.. _Requirements: https://github.com/iric-soft/km#requirements
.. _Install: https://github.com/iric-soft/km#install
.. _Usage: https://github.com/iric-soft/km#usage
.. _Tools: https://github.com/iric-soft/km#tools

.. _General: https://github.com/iric-soft/km#general
.. _Runing km on a real sample: https://github.com/iric-soft/km#runing-km-on-a-real-sample
.. _find_mutation: https://github.com/iric-soft/km#find_mutation
.. _find_report: https://github.com/iric-soft/km#find_report

-------------
Introduction:
-------------

This tool was developed to identify and quantify the occurence of single nucleotide variants, insertions, deletions and duplications in RNA-seq data.  Contrary to most tools that try to report all variants in a complete genome, here we instead propose to focus the analysis on small regions of interest.

Given a reference sequence (typically a few hundred base pairs) around a known or suspected mutation in a gene of interest, all possible sequences that can be be created between the two end k-mers according to the sequenced reads will be reported.  A ratio of variant allele vs WT will be computed for each possible sequence constructed.

-------------
Requirements:
-------------
* Python 2.7.6 or later
* Jellyfish 2.2 or later **with** Python `bindings`_.
* (Optional) Matplotlib

--------
Install:
--------
Before installing or using km, Jellyfish needs to be installed with Python `bindings`_.
Should you need it, a script is available in the `example`_ folder, to help you to install Jellyfish with bindings.
When jellyfish is installed, km may be executed directly from source code or installed using the setup script:

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

find_mutation:
--------------

Usage:
******

.. code:: shell

  $ km find_mutation -h
  $ km find_mutation [your_target_seq].fa [your_count_table].jf

Output:
*******

Here we are looking for a common 4-bp duplication that occurs in some leukemias, and that is especially troublesome to detect since it occurs a few base pairs from the start of the last exon.  Most standard mapping techniques will miss this variant.  Running the find_mutation command takes a few seconds and returns an output similar to this:

.. code:: shell

  Database	Query	Type	Variant name	Ratio	Expression	Min coverage	Start offset  Sequence	Reference ratio	Reference expression	Reference sequence	Info
  ./data/jf/02H025_NPM1.jf	NPM1_4ins_exons_10-11utr	Insertion	45:/TCTG:45	0.484	2870.6	2428	0	AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA	0.516	3055.2	AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA	vs_ref
  ./data/jf/02H025_NPM1.jf	NPM1_4ins_exons_10-11utr	Reference		1.000	2379.0	2379	0	AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA	1.000	2379.0	AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA	vs_ref
  ./data/jf/02H025_NPM1.jf	NPM1_4ins_exons_10-11utr	Insertion	45:/TCTG:45	0.484	2972.6	2428	9	CGGATGACTGACCAAGAGGCTATTCAAGATCTCTGTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAG	0.516	3172.9	CGGATGACTGACCAAGAGGCTATTCAAGATCTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAG	cluster 1 n=1

which shows that:

* a TCTG insertion was found at position 45 of the target sequence: NPM1_exons_10-11utr.
* the target sequence was found (without mutations).

The last line is the same as the first one with local calculation of Ratio, Expression and Min coverage.
It's a try to allowed long target sequence which can found several variants.

Output description:
*******************
Each line represents a path of the local assembly constructed from the target sequence.

* Database: Name of the Jellyfish kmer table queried
* Query: Name of the target sequence examined
* Type: Type of mutation found (Insertion, Deletion or Substitution).  A Reference type used to identify path without mutation
* Variant name: A description of the modification in the format start_position:deleted_bases/inserted_bases:end_position
* Ratio: Estimated ratio for the mutated allele represented by this path
* Expression: Estimated expression level for the mutated allele (coverage)
* Min coverage: Min k-mer count of all k-mers in the path
* Start offset: Starting position of sequences. Usefull for cluster quantification method (see Info column).
* Sequence: Sequence of the mutated path
* Reference ratio: Estimated ratio of the target allele
* Reference expression: Estimated expression level for the target
* Reference sequence: Target sequence used
* Info: Supplementary information regarding the quantification method.

  - vs_ref: means that each alternate path is compared in expression with the whole target sequence.
  - cluster: indicates that all alternate path in a subregion extending by k bases on each side of all overlapping mutations are considered at once to evaluate the expression of each

Using the -g argument, one can also obtain a coverage graph for the two variants, for example:

.. image:: https://github.com/iric-soft/km/blob/master/data/figure/figure_1.png

find_report:
-------------------

Usage:
******

.. code:: shell

  $ km find_report -h
  $ km find_report -t [your_target_seq].fa [find_mutation_output]
  $ km find_mutation [your_target_seq].fa [your_count_table].jf | km find_report -t [your_target_seq].fa

Output:
*******

.. code:: shell

  Sample	Region	Location	Type	Removed	Added	Abnormal	Normal	Ratio	Min coverage	Variant	Target	Info	Variant sequence	Reference sequence
  ./data/jf/02H025_NPM1.jf	chr5:171410540-171410543	chr5:171410544	ITD	0	4 | 4	2870.6	3055.2	0.484	2428	/TCTG	NPM1_4ins_exons_10-11utr	vs_ref	AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA	AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA
  ./data/jf/02H025_NPM1.jf		-	Reference	0	0	0.0	2379.0	1.000	2379	-	NPM1_4ins_exons_10-11utr	vs_ref

which shows that an ITD variant (TCTG insertion) was found at position chr5:171410544

Output description:
*******************
Each line represents a path that was constructed from the target sequence.

* Sample: name of the Jellyfish kmer table queried
* Region: the variant chromosome region
* Location: the variant chromosome position
* Type: the variant type
* Removed: number of nucleotides removed
* Added: number of nucleotides added spliced | unspliced
* Abnormal: estimated expression level for the mutated allele (coverage)
* Normal: estimated expression level for the target
* Ratio: estimated ratio for the mutated allele represented by this path
* Min coverage: Min k-mer count of all k-mers in the path
* Variant: A description of the variant in the format: deleted_bases/inserted_bases
* Target: name of the target sequence examined
* Info: supplementary information regarding the quantification method.
* Sequence: sequence of the mutated path
* Reference sequence: target sequence used
