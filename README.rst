
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
* Jellyfish 2.2 or later (http://www.genome.umd.edu/jellyfish.html)
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

General:
--------

From the source:
****************

.. code:: shell

  $ cd [your_km_folder]
  $ python -m km -h

After setup install:
********************

.. code:: shell

  $ km -h

find_mutation tool:
-------------------

Usage:
******

.. code:: shell

  $ km find_mutation -h
  $ km find_mutation [your_target_seq].fa [your_count_table].jf

Output:
*******

Here we are looking for a common 4-bp duplication that occurs in some leukemias, and that is especially troublesome to detect since it occurs a few base pairs from the start of the last exon.  Most standard mapping techniques will miss this variant.  Running the find_mutation command takes a few seconds and returns an output similar to this:

.. code:: shell

  Database	Query	Type	Variant name	Ratio	Expression	Min coverage	Sequence	Reference ratio	Reference expression	Reference sequence	Info
  02H025/kmers-2.2.3_31.jf	NPM1_exons_10-11utr	Insertion	45:/TCTG:45	0.481	2865.2	2436	AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA	0.519	3097.0	AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA	vs_ref
  02H025/kmers-2.2.3_31.jf	NPM1_exons_10-11utr	Reference		1.000	2436.0	2449	AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA	1.000	2436.0	AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA	vs_ref
  02H025/kmers-2.2.3_31.jf	NPM1_exons_10-11utr	Insertion	45:/TCTG:45	0.480	2975.4	2436	CGGATGACTGACCAAGAGGCTATTCAAGATCTCTGTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAG	0.520	3224.1	CGGATGACTGACCAAGAGGCTATTCAAGATCTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAG	cluster 1 n=1

which shows that:

* a TCTG insertion was found at position 45 of the target sequence: NPM1_exons_10-11utr.
* the target sequence was found (without mutations).

The last line is the same as the first one with local calculation of Ratio, Expression and Min coverage.
It's a try to allowed long target sequence which can found several variants.

Output description:
*******************
Each line represents a path of the local assembly constructed from the target sequence.

* Database: name of the Jellyfish kmer table queried
* Query: name of the target sequence examined
* Type: type of mutation found (Insertion, Deletion or Substitution).  A Reference type used to identify path without mutation
* Variant name: A description of the modification in the format start_position:deleted_bases/inserted_bases:end_position
* Ratio: estimated ratio for the mutated allele represented by this path
* Expression: estimated expression level for the mutated allele (coverage)
* Min coverage: Min k-mer count of all k-mers in the path
* Sequence: sequence of the mutated path
* Reference ratio: estimated ratio of the target allele
* Reference expression: estimated expression level for the target
* Reference sequence: target sequence used
* Info: supplementary information regarding the quantification method.

  - vs_ref: means that each alternate path is compared in expression with the whole target sequence.
  - cluster: indicates that all alternate path in a subregion extending by k bases on each side of all overlapping mutations are considered at once to evaluate the expression of each

Using the -g argument, one can also obtain a coverage graph for the two variants, for example:

TODO add viz_1.png


find_report tool:
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
  02H025/kmers-2.2.3_31.jf	chr5:171410540-171410543	chr5:171410544	ITD	0	4 | 4	2865.2	3097.0	0.481	2436	/TCTG	NPM1_exons_10-11utr	vs_ref	AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA	AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA
  02H025/kmers-2.2.3_31.jf		-	Reference	0	0	0.0	2436.0	1.000	2449	-	NPM1_exons_10-11utr	vs_ref

which shows that an ITD variant (TCTG insertion) was found at position chr5:171410544

Output description:
*******************
Each line represents a path that was constructed from the target sequence.

* Sample: name of the Jellyfish kmer table queried
* Region: The variant chromosome region
* Location: The variant chromosome position
* Type: The variant type
* Removed: ...
* Added: ...
* Abnormal: estimated expression level for the mutated allele (coverage)
* Normal: estimated expression level for the target
* Ratio: estimated ratio for the mutated allele represented by this path
* Min coverage: Min k-mer count of all k-mers in the path
* Variant: A description of the variant in the format: deleted_bases/inserted_bases
* Target: name of the target sequence examined
* Info: supplementary information regarding the quantification method.
* Sequence: sequence of the mutated path
* Reference sequence: target sequence used
