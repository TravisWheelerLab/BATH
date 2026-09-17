# BATH 2.0 release notes (May 2026)

This directory holds BATH release notes. The notes for HMMER, from
which BATH is derived, are in the HMMER/ subdirectory.


## new features:

* Spliced alignment. With the new `--splice` option, bathsearch
  models intron splice sites (GT-AG, GC-AG, and AT-AC, with empirical
  signal probabilities) and joins ORF hits across introns into
  gene-level exon chains. The new `--exontblout <f>` option writes a
  tabular file with one row per exon, and `--min_intron <n>` and
  `--max_intron <n>` bound the intron lengths considered. A single
  search cannot yet combine `--splice` and `--fs`.

* SIMD acceleration of the frameshift-aware and spliced dynamic
  programming kernels: AVX, AVX2, and AVX-512 implementations on x86,
  and NEON on ARM, with automatic fallback to SSE where newer
  intrinsics are unavailable.

* A local bias filtering stage that recomputes composition bias over
  just the aligned region, to better reject short repetitive false
  positives.

* A new Forward-stage filter (F4): at least one ORF in a window must
  pass a standard Forward P-value threshold before bathsearch
  promotes the window to the frameshift or spliced Forward stage.

* Reworked ORF and DNA window construction, including handling of
  degenerate nucleotides: ORFs no longer span long runs of degenerate
  nucleotides, and split codons in spliced alignments may contain
  them.


## changes:

* Frameshift-aware alignment is now opt-in. By default bathsearch
  performs translated (non-frameshifted) search; `--fs` enables the
  frameshift-aware pipeline. The frameshift probability is now a fixed
  internal constant set at model build/convert time, and is no longer
  a tunable option.

* BATH 2.0 does not read model files built by BATH 1.x; run
  bathconvert to update them. bathsearch prints this instruction when
  it reads an old file.

* Removed: PowerPC VMX/AltiVec support; MPI support in bathbuild; the
  bathsearch bias-filter tuning options `--B1`, `--B2`, `--B3`; and
  the envelope-coordinate columns in search output.


## bug fixes:

* Fixed a bug in local bias filtering.

* Reduced memory use by reusing dynamic programming matrices (the
  Backward matrix for posterior probabilities, and the Forward matrix
  for Viterbi alignment when posterior decoding fails).


For more information, see the
[git log](https://github.com/TravisWheelerLab/BATH/commits/main).
