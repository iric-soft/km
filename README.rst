
===================================================================
km : a software for RNA-seq investigation using k-mer decomposition
===================================================================

-------------
Introduction:
-------------

This tool was developed to identify and quantify the occurence of single nucleotide variants, insertions, deletions and duplications in RNA-seq data.  Contrary to most tools that try to report all variants in a complete genome, here we instead propose to focus the analysis on small regions of interest.

Given a reference sequence of interest (typically a few hundred base pairs) around a known or suspected mutation in a gene of interest, all possible sequences that can be be created between the two end k-mers according to the sequenced reads will be reported.  A ratio of variant allele vs WT will be computed for each possible sequence constructed.

-------------
Requirements:
-------------
* Python 2.7.6 or later
* Jellyfish 2.1 or later (http://www.genome.umd.edu/jellyfish.html)
* (Optional) Matplotlib

--------
Install:
--------
Before to install or use km, you need to install Jellyfish with `binding`_.
When it's done, you can use the setup install command or execute the
code from the source:

.. code:: shell

  $ python setup.py install
  $ km -h

.. _binding: https://github.com/gmarcais/Jellyfish#binding-to-script-languages

------
Usage:
------

From the source:
----------------

.. code:: shell

  $ cd [your_km_folder]
  $ python -m km find_mutation data/catalog/GRCh38/NPM1*.fa [your_kmer_count_table].jf

After setup install:
--------------------

.. code:: shell

  $ km find_mutation [your_target_seq].fa [your_kmer_count_table].jf


-------
Output:
-------

find_mutation output:
---------------------
Here we are looking for a common 4-bp duplication that occurs in some leukemias, and that is especially troublesome to detect since it occurs a few base pairs from the start of the last exon.  Most standard mapping techniques will miss this variant.  Running the find_mutation command takes a few seconds and returns an output similar to this:

.. code:: shell

  Database	Query	Type	Variant name	Ratio	Expression	Sequence	Reference ratio	Reference expression	Reference sequence	Info
  sample.jf	NPM1	Insertion	93:/TCTG:93	0.528	9020.0	CCAAGAGGCTATTCAAGATCTCTGTCTGGCAGTGGAGGAAGTCTCTT	0.472	8076.8	CCAAGAGGCTATTCAAGATCTCTGGCAGTGGAGGAAGTCTCTT	cluster 1 n=1

  real	0m6.712s
  user	0m0.203s
  sys	0m0.312s

which shows that a TCTG insertion was found at position 93 of the refence sequence.

Output description:
*******************
Each line represents a path that was constructed from the reference sequence.

* Database: name of the Jellyfish kmer table queried
* Query: name of the reference sequence examined
* Type: type of mutation found (Insertion, Deletion or Substitution).  A Reference type used to identify path without mutation
* Variant name: A description of the modification in the format start_position:deleted_bases/inserted_bases:end_position
* Ratio: estimated ratio for the mutated allele represented by this path
* Expression: estimated expression level for the mutated allele (coverage)
* Sequence: sequence of the mutated path
* Reference ratio: estimated ratio of the reference allele
* Reference expression: estimated expression level for the reference
* Reference sequence: reference sequence used
* Info: supplementary information regarding the quantification method.

  - vs_ref: means that each alternate path is compared in expression with the whole reference sequence.
  - cluster: indicates that all alternate path in a subregion extending by k bases on each side of all overlapping mutations are considered at once to evaluate the expression of each

Using the -g argument, one can also obtain a coverage graph for the two variants, for example:

figure_1.png


report_mutation output:
-----------------------
Soon...
