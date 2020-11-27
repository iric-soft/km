
===================================================================
km's tools documentation
===================================================================

---------
Contents:
---------
* `find_mutation`_

  - |fm-usage|_
  - |fm-output|_
  - |fm-output-desc|_

* `find_report`_

  - |fr-usage|_
  - |fr-output|_
  - |fr-output-desc|_

* `min_cov`_

  - |mc-usage|_
  - |mc-output|_
  - |mc-output-desc|_

* `linear_kmin`_

  - |lk-usage|_
  - |lk-output|_
  - |lk-output-desc|_

* `rVAF`_

.. _find_mutation: https://github.com/iric-soft/km/tree/master/km/tools#find_mutation
.. _find_report: https://github.com/iric-soft/km/tree/master/km/tools#find_report
.. _min_cov: https://github.com/iric-soft/km/tree/master/km/tools#min_cov
.. _linear_kmin: https://github.com/iric-soft/km/tree/master/km/tools#linear_kmin
.. _rVAF: https://github.com/iric-soft/km/tree/master/km/tools#rVAF

.. _fm-usage: https://github.com/iric-soft/km/tree/master/km/tools#usage
.. _fr-usage: https://github.com/iric-soft/km/tree/master/km/tools#usage-1
.. _mc-usage: https://github.com/iric-soft/km/tree/master/km/tools#usage-2
.. _lk-usage: https://github.com/iric-soft/km/tree/master/km/tools#usage-3

.. _fm-output: https://github.com/iric-soft/km/tree/master/km/tools#output
.. _fr-output: https://github.com/iric-soft/km/tree/master/km/tools#output-1
.. _mc-output: https://github.com/iric-soft/km/tree/master/km/tools#output-2
.. _lk-output: https://github.com/iric-soft/km/tree/master/km/tools#output-3

.. _fm-output-desc: https://github.com/iric-soft/km/tree/master/km/tools#output-description
.. _fr-output-desc: https://github.com/iric-soft/km/tree/master/km/tools#output-description-1
.. _mc-output-desc: https://github.com/iric-soft/km/tree/master/km/tools#output-description-2
.. _lk-output-desc: https://github.com/iric-soft/km/tree/master/km/tools#output-description-3

.. |fm-usage| replace:: Usage
.. |fr-usage| replace:: Usage
.. |mc-usage| replace:: Usage
.. |lk-usage| replace:: Usage

.. |fm-output| replace:: Output
.. |fr-output| replace:: Output
.. |mc-output| replace:: Output
.. |lk-output| replace:: Output

.. |fm-output-desc| replace:: Output description
.. |fr-output-desc| replace:: Output description
.. |mc-output-desc| replace:: Output description
.. |lk-output-desc| replace:: Output description

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

  Database	Query	Type	Variant_name	rVAF	Expression	Min_coverage	Start_offset  Sequence	Reference_expression	Reference_sequence	Info
  ./data/jf/02H025_NPM1.jf	NPM1_4ins_exons_10-11utr	Reference		nan	2379.0	2379	0	AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA	2379.0	AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA	vs_ref
  ./data/jf/02H025_NPM1.jf	NPM1_4ins_exons_10-11utr	Insertion	45:/TCTG:45	0.484	2870.6	2428	0	AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA	3055.2	AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA	vs_ref
  ./data/jf/02H025_NPM1.jf	NPM1_4ins_exons_10-11utr	Insertion	45:/TCTG:45	0.484	2972.6	2428	9	CGGATGACTGACCAAGAGGCTATTCAAGATCTCTGTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAG	3172.9	CGGATGACTGACCAAGAGGCTATTCAAGATCTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAG	cluster 1 n=1

which shows that:

* a TCTG insertion was found at position 45 (1-based) of the target sequence: NPM1_exons_10-11utr.
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
* `rVAF`_: Estimated reference Variant Alelle Frequencies for the mutated allele represented by this path.
* Expression: Estimated expression level for the mutated allele (coverage)
* Min_coverage: Min k-mer count of all k-mers in the path
* Start_offset: Starting position of sequences. Usefull for cluster quantification method (see Info column).
* Sequence: Sequence of the mutated path
* Reference_expression: Estimated expression level for the target
* Reference_sequence: Target sequence used
* Info: Supplementary information regarding the quantification method.

  - vs_ref: means that each alternate path is compared in expression with the whole target sequence.
  - cluster: indicates that all alternate path in a subregion extending by k bases on each side of all overlapping mutations are considered at once to evaluate the expression of each

Using the -g argument, one can also obtain a coverage graph for the two	variants, for example:

.. image:: https://github.com/iric-soft/km/blob/master/data/figure/figure_1.png

------------
find_report:
------------
This tool parse find_mutation output to reformat it in more user friendly
tabulated file.

Usage:
------

.. code:: shell

  $ km find_report -h
  $ km find_report -t [your_fasta_targetSeq] [find_mutation_output]
  $ km find_mutation [your_fasta_targetSeq] [your_jellyfish_count_table] | km find_report -t [your_fasta_targetSeq]

Output:
-------

.. code:: shell

  Sample	Region	Location	Type	Removed	Added	Abnormal	Normal	rVAF	Min_coverage	Exclu_min_cov  Variant	Target	Info	Variant_sequence	Reference_sequence
  ./data/jf/02H025_NPM1.jf		-	Reference	0	0	0.0	2379.0	nan	2379	 -	NPM1_4ins_exons_10-11utr	vs_ref
  ./data/jf/02H025_NPM1.jf	chr5:171410540-171410543	chr5:171410544	ITD	0	4 | 4	2870.6	3055.2	0.484	2428 0	/TCTG	NPM1_4ins_exons_10-11utr	vs_ref	AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA	AATTGCTTCCGGATGACTGACCAAGAGGCTATTCAAGATCTCTGGCAGTGGAGGAAGTCTCTTTAAGAAAATAGTTTAAA

which shows that an ITD variant (TCTG insertion) was found at position
chr5:171410544.

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
* `rVAF`_: Estimated reference Variant Alelle Frequencies for the mutated allele represented by this path.
* Min_coverage: Min k-mer count of all k-mers in the path
* Exclu_min_cov: Min k-mer count of all k-mers in the variant sequence from the jf database given with "-e".
* Variant: A description of the variant in the format: deleted_bases/inserted_bases
* Target: name of the target sequence examined
* Info: supplementary information regarding the quantification method.
* Variant_sequence: sequence of the mutated path
* Reference_sequence: target sequence used

--------
min_cov:
--------

This tools display some k-mer's coverage stats of a target sequence and a list of jellyfish database.

Usage:
------
.. code:: shell

  $ km min_cov -h
  $ km min_cov [your_fasta_targetSeq] [[your_jellyfish_count_table]...]

Output:
-------

.. code:: shell

  DB  count length  min max mean  kmer_nb kmer_nb_0
  /dev/shm/02H053.jf  455387  8371  0 318 54.60 8341  4
  /dev/shm/05H094.jf  58582 8371  0 36  7.02  8341  674
  /dev/shm/05H143.jf  1302959 8371  7 450 156.21  8341  0

Which shows that the sample 05H094 have at least one part of the target sequence not covered by k-mer count.

Output description:
-------------------

* DB: name of the Jellyfish kmer table queried
* cout: sum of k-mer count
* length: number of nucleotide of target sequence
* min: Minimum of k-mer count in the target sequence
* max: Maximum of k-mer count in the target sequence
* mean: Mean of k-mer count in the target sequence
* kmer_nb: Number of kmer in the target sequence
* kmer_nb_0: Number of kmer with 0 count in the target sequence

------------
linear_kmin:
------------
Length of k-mers is a central parameter:

* To produce a linear directed graph from the target sequence.
* To avoid false-positive.

**Warning**: `find_mutation`_ shouldn't be use on jellyfish count table build with k<21 bp (we recommand k=31 bp, by default)

linear_kmin tool is design to report the minimun k length to allow a linear decomposition of a target sequence.

Usage:
------

.. code:: shell

  $ km linear_kmin -h
  $ km linear_kmin [your_catalog_directory]

Output:
-------

.. code:: shell

  $ km linear_kmin -s 5 ./data/catalog/GRCh38/
  target_name linear_kmin
  FLT3-TKD_exon_20  8
  MYC_T58A_P59R_exon2 7
  NSD1_exon6-NUP98_exon13 9
  NUP98_exon11-NSD1_exon7 7
  DNMT3A_R882_exon_23 6
  FLT3-ITD_exons_13-15  10
  KMT2A-PTD_8-2 7
  NPM1_4ins_exons_10-11utr 7

For this catalog of target sequences, this output shows that
`find_mutation`_ need to be run on jellyfish count tables build
with at least k >= 10 bp.
Which is under the threshold to avoid the detection of false-positive
mutations. This is not always the case, especially on large target sequence
(like a transcript), where linear_kmin could be more longer than sequenced
read length (100 bp, Like ENST00000621744_NBPF19 need a k >= 3472 pb).

Output description:
-------------------

* target_name: name of target sequence.
* linear_kmin: minimum k length to decompose the target sequence in linear graph.

-----
rVAF:
-----
rVAF reported by km currently only considers mutated events independently even
if the signal at a given point may be spread out between more than 2 alternatives.

So, for example, if km returns multiple potential mutations:

* M1: mutated path 1 with a coverage of 80
* M2: mutated path 2 with a coverage of 250
* R: reference path with a coverage of 100

We currently compute the following rVAF by only comparing to the reference signal:

* rVAF(M1) = 80/(80+100) = 0.44
* rVAF(M2) = 250/(250+100) = 0.71
* rVAF(R) = nan

For cases where M1 and M2 do not overlap, our rVAF are in fact Variant Alelle
Frequencies. However when there is overlap, deconvoluting the total signal
between all paths is an aspect that we have not yet developed.
