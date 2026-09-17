# BATH 1.0 release notes (Feb 2025)

BATH 1.0 is the first release of BATH. It is essentially a fork of the
HMMER 3.4 code base, extended to perform translated search (protein
query against DNA target) directly. bathsearch takes a protein query
and a DNA target and reports hits in the coordinates of the target
DNA. There is no need to compute ORFs in all six translation frames
first, and no bookkeeping to map hits on those ORFs back to the
encoding genome.

BATH also introduces novel frameshift-aware algorithms to detect
frameshift-inducing nucleotide insertions and deletions (indels).

Programs in this release: bathbuild, bathconvert, bathfetch, bathstat,
and bathsearch.
