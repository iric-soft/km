
===================================================================
km's tools documentation
===================================================================

---------
Contents:
---------
* `find_mutation`_

  - |fm-usage|_
  - `find_mutation output`_
  - `find_mutation output description`_

* `find_report`_

  - `fr-usage`_
  - `find_report output`_
  - `find_report output description`_

.. _find_mutation: https://github.com/iric-soft/km/tree/master/km/tools#find_mutation
.. _find_report: https://github.com/iric-soft/km/tree/master/km/tools#find_report

.. _fm-usage: https://github.com/iric-soft/km/tree/master/km/tools#usage
.. _fr-usage: https://github.com/iric-soft/km/tree/master/km/tools#usage-1

.. |fm-usage| replace:: Usage

--------------
find_mutation:
--------------
This is the main tool of km, to identify and quantify mutations from
a target sequence and a k-mer database.

Usage:
------

.. code:: shell

  $ km find_mutation -h
  $ km find_mutation [your_fasta_targetSeq] [your_jellyfish_count_table]
  $ km find_mutation [your_catalog_directory] [your_jellyfish_count_table]

Output:
-------

Here we are looking for a common 4-bp duplication that occurs in some
leukemias, and that is especially troublesome to detect since it occurs
a few base pairs from the start of the last exon. Running the find_mutation
command takes a few seconds and returns an output similar to this:

.. code:: shell

  Database	Query	Type	Variant name	Ratio	Expression	Min coverage	Start offset  Sequence	Reference ratio	Reference expression	Reference sequence	Info
  ./data/jf/02H025_NPM1.jf	NPM1_4ins_exons_10-11utr	Insertion	45:/TCTG:45	0.484	2870.6	2428	0	AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA	0.516	3055.2	AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA	vs_ref
  ./data/jf/02H025_NPM1.jf	NPM1_4ins_exons_10-11utr	Reference		1.000	2379.0	2379	0	AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA	1.000	2379.0	AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA	vs_ref
  ./data/jf/02H025_NPM1.jf	NPM1_4ins_exons_10-11utr	Insertion	45:/TCTG:45	0.484	2972.6	2428	9	CGGATGACTGACCAAGAGGCTATTCAAGATCTCTGTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAG	0.516	3172.9	CGGATGACTGACCAAGAGGCTATTCAAGATCTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAG	cluster 1 n=1

which shows that:

* a TCTG insertion was found at position 45 of the target sequence: NPM1_exons_10-11utr.
* the target sequence was found (without mutations).

The last line is the same as the first one with local calculation of Ratio,
Expression and Min coverage. It's a try to allowed long target sequence
which can found several variants.

Output description:
-------------------

Each line represents a path of the local assembly constructed from the
target sequence.

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

Using the -g argument, one can also obtain a coverage graph for the two
variants, for example:

.. image:: https://github.com/iric-soft/km/blob/master/data/figure/figure_1.png

------------
find_report:
------------

Usage:
------

.. code:: shell

  $ km find_report -h
  $ km find_report -t [your_fasta_targetSeq] [find_mutation_output]
  $ km find_mutation [your_fasta_targetSeq] [your_jellyfish_count_table] | km find_report -t [your_fasta_targetSeq]

Output:
-------

.. code:: shell

  Sample	Region	Location	Type	Removed	Added	Abnormal	Normal	Ratio	Min coverage	Variant	Target	Info	Variant sequence	Reference sequence
  ./data/jf/02H025_NPM1.jf	chr5:171410540-171410543	chr5:171410544	ITD	0	4 | 4	2870.6	3055.2	0.484	2428	/TCTG	NPM1_4ins_exons_10-11utr	vs_ref	AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA	AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA
  ./data/jf/02H025_NPM1.jf		-	Reference	0	0	0.0	2379.0	1.000	2379	-	NPM1_4ins_exons_10-11utr	vs_ref

which shows that an ITD variant (TCTG insertion) was found at position
chr5:171410544

Output description:
-------------------

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
