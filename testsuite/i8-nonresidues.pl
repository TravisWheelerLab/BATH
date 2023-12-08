#! /usr/bin/perl

# Regression test of handling a nonresidue '*' character. By design,
# '*' residues score 0 in insert states and N,C,J; and -inf in match
# states. Test case verifies that an inserted * is handled within one
# alignment, and a consensus * breaks an alignment in two.
# Implemented as a regression test against specific scores (in domtbl
# output) of a small manually checked example.
#
# Usage:    ./i8-nonresidues.pl <builddir> <srcdir> <tmpfile prefix>
# Example:  ./i8-nonresidues.pl ..         ..        tmpfoo
#
# SRE, Thu Oct 29 08:38:09 2009

BEGIN {
    $builddir = shift;
    $srcdir   = shift;
    $tmppfx   = shift;
}

use lib "$srcdir/testsuite";
use h3;

$bathsearch = "$builddir/src/bathsearch";
$hmm20aa   = "$srcdir/testsuite/20aa.bhmm";

# Two test sequences, to be aligned to 20aa.bhmm
# First one will get parsed into two hits (nucleotide replaced by '*')
# Second one will get parsed into one hit because '*' is an insertion;
# this is a "feature" of insertions being scored as 0 regardless of 
# residue identity.
# 
if (! open(FP1, ">$tmppfx.1")) { print "FAIL: couldn't open $tmppfx.1 for writing"; exit 1; }
if (! open(FP2, ">$tmppfx.2")) { print "FAIL: couldn't open $tmppfx.2 for writing"; exit 1; }

print FP1 <<"EOF";
>test1
GCATGTGACGAGTTTGGCCATATAAAA*TTATGAATCCACAGCGCTCAACTGTATGGTAT
EOF

print FP2 <<"EOF";
>test2
GCATGTGACGAGTTTGGCCATATAAAAC*TTATGAATCCACAGCGCTCAACTGTATGGTAT
EOF

system("$bathsearch --nofs --tblout $tmppfx.tbl $hmm20aa $tmppfx.1 > $tmppfx.out 2>&1");
if ($? != 0) { die "FAIL: bathsearch failed on first test sequence with --nofs\n"; }
&h3::ParseTbl("$tmppfx.tbl");

# Verify.
if ($h3::ntbl        != 1)      { printf("FAIL: expected one line in tbl; saw %d\n",              $h3::ntbl);          exit 1; }
if ($h3::fullsc[0]   != "67.3") { printf("FAIL: expected score of 67.3 for first hit; saw %s\n",  $h3::fullsc[0]);     exit 1; }
if ($h3::fullbias[0] != "0.1")  { printf("FAIL: expected bias of 0.1 for first hit; saw %s\n",    $h3::fullbias[0]);   exit 1; }

system("$bathsearch --fsonly --tblout $tmppfx.tbl $hmm20aa $tmppfx.1 > $tmppfx.out 2>&1");
if ($? != 0) { die "FAIL: bathsearch failed on first test sequence with --fsonly\n"; }
&h3::ParseTbl("$tmppfx.tbl");

# Verify.
if ($h3::ntbl        != 1)      { printf("FAIL: expected one line in tbl; saw %d\n",              $h3::ntbl);          exit 1; }
if ($h3::fullsc[0]   != "68.9") { printf("FAIL: expected score of 68.9 for first hit; saw %s\n",  $h3::fullsc[0]);     exit 1; }
if ($h3::fullbias[0] != "0.0")  { printf("FAIL: expected bias of 0.0 for first hit; saw %s\n",    $h3::fullbias[0]);   exit 1; }

system("$bathsearch -l 10 --nofs --tblout $tmppfx.tbl $hmm20aa $tmppfx.2 > $tmppfx.out 2>&1");
if ($? != 0) { print "FAIL: bathsearch failed on second test sequence with --nofs"; }
&h3::ParseTbl("$tmppfx.tbl");

if ($h3::ntbl    != 2)          { printf("FAIL: expected two lines in tbl; saw %d\n",     $h3::ntbl);    exit 1; }
if ($h3::fullsc[0]   != "31.1") { printf("FAIL: expected score of 31.1; saw %s\n",        $h3::fullsc[0]);   exit 1; }
if ($h3::fullbias[0] != "2.8")  { printf("FAIL: expected bias of 2.8; saw %s\n",          $h3::fullbias[0]); exit 1; }
if ($h3::fullsc[1]   != "25.2") { printf("FAIL: expected score of 25.2; saw %s\n",        $h3::fullsc[1]);   exit 1; }
if ($h3::fullbias[1] != "0.6")  { printf("FAIL: expected bias of 0.6; saw %s\n",          $h3::fullbias[1]); exit 1; }

system("$bathsearch -l 10 --fsonly --tblout $tmppfx.tbl $hmm20aa $tmppfx.2 > $tmppfx.out 2>&1");
if ($? != 0) { print "FAIL: bathsearch failed on second test sequence with --fsonly"; }
&h3::ParseTbl("$tmppfx.tbl");

if ($h3::ntbl        != 1)      { printf("FAIL: expected two lines in tbl; saw %d\n",             $h3::ntbl);          exit 1; }
if ($h3::fullsc[0]   != "63.5") { printf("FAIL: expected score of 63.5 for first hit; saw %s\n",  $h3::fullsc[0]);     exit 1; }
if ($h3::fullbias[0] != "0.0")  { printf("FAIL: expected bias of 0.0 for first hit; saw %s\n",    $h3::fullbias[0]);   exit 1; }

print "ok\n";
unlink "$tmppfx.1";
unlink "$tmppfx.2";
unlink "$tmppfx.tbl";
unlink "$tmppfx.out";
exit 0;

