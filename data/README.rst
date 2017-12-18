
===================================================================
DATA
===================================================================

Some data used for km project

---------
Contents:
---------
* `Catalog`_
* `Figure`_
* `JF`_

.. _Catalog: https://github.com/iric-soft/km/tree/master/data#catalog
.. _Figure: https://github.com/iric-soft/km/tree/master/data#figure
.. _JF: https://github.com/iric-soft/km/tree/master/data#jf
        https://github.com/iric-soft/km/tree/master/data#jf
--------
Catalog:
--------

Description:
------------

Catalog of target sequences used for AML samples, which includes:

* Single-base mutation in DNMT3A and MYC
* Insertions in NPM1
* Internal Tandem Duplication (ITD) and mutations in the Tyrosine Kinase Domain (TKD) of FLT3
* Partial Tandem Duplication (PTD) in KMT2A
* Fusion between NUP98 and NSD1

Create a target sequence:
----------------------------

A target sequence is a fasta file with a header, which specify chromosomic position of each sequences.

The header format is: [chrom]:[start]-[end] (see `FLT3-ITD`_ composed of 3 exons).

.. _FLT3-ITD: https://github.com/iric-soft/km/blob/master/data/catalog/GRCh38/FLT3-ITD_exons_13-15.fa

-------
Figure:
-------
Figures used to illustrate km documentation.

---
JF:
---
Small fraction of Jellyfish files used to `test`_ the km code.

Don't use these files to analyse a sample!

.. _test: https://github.com/iric-soft/km/tree/master/km/tests
