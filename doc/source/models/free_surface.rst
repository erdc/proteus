Free Surface
************


There are two implementations for dealing with the free surface.


VOF - NCLS
==========

In this implementation, the following 4 models must be solved in order:

* Volume of Fluid (VOF): :py:mod:`proteus.mprans.VOF`
* Non-conservative level set (NCLS): :py:mod:`proteus.mprans.NCLS`
* Redistancing: :py:mod:`proteus.mprans.RDLS`
* Mass correction: :py:mod:`proteus.mprans.MCorr`


CLSVOF
======

In this implementation, a conservative level set is used and only the following
model must be solved: :py:mod:`proteus.mprans.CLSVOF`.
