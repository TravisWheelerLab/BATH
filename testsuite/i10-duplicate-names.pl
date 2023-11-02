#! /usr/bin/perl

# Check how we handle duplicate names, for both queries and targets.
# Indexing HMM or sequence files does require unique names/accessions.
# Aside from SSI indexes, other BATH operations do not require unique
# names/accessions.
#
# Usage:    ./i10-duplicate-names.pl <builddir> <srcdir> <tmpfile prefix>
# Example:  ./i10-duplicate-names.pl ..         ..       tmpfoo
#
# SRE, Sun Dec 13 14:41:31 2009 [Yokohama, Japan]


BEGIN {
    $builddir = shift;
    $srcdir   = shift;
    $tmppfx   = shift;
}
use lib "$srcdir/testsuite";  # The BEGIN is necessary to make this work: sets $srcdir at compile-time
use h3;


# Verify that we have all the executables we need for the test.
if (! -x "$builddir/src/bathbuild")   { die "FAIL: didn't find bathbuild binary in $builddir/src\n";  }
if (! -x "$builddir/src/bathsearch")  { die "FAIL: didn't find bathsearch binary in $builddir/src\n"; }

# Create our test files
if (! open(ALI1, ">$tmppfx.sto")) { print "FAIL: couldn't open $tmppfx.sto for write";  exit 1; }
if (! open(SEQ1, ">$tmppfx.fa"))  { print "FAIL: couldn't open $tmppfx.fa for write";   exit 1; }

print ALI1 <<"EOF";
# STOCKHOLM 1.0
#=GF ID profile
#=GF AC XX01234.5
#=GF DE A test description
seq1 ACDEFGHIKLMNPQRSTVWY
seq2 ACDEFGHIKLMNPQRSTVWY
seq3 ACDEFGHIKLMNPQRSTVWY
//
# STOCKHOLM 1.0
#=GF ID profile
#=GF AC XX01234.5
#=GF DE A test description
seq1 ACDEFGHIKLLMNPQRSTVWY
seq2 ACDEFGHIKLLMNPQRSTVWY
seq3 ACDEFGHIKLLMNPQRSTVWY
//
EOF

print SEQ1 << "EOF";
>seq
GCATGTGACGAGTTTGGCCATATAAAACTTATGAATCCACAGCGCTCAACTGTATGGTAT
>seq
GCATGTGACGAGTTTGGCCATATAAAACTTATGAATCCACAGCGCTCAACTGTATGGTAT
EOF

close ALI1;
close SEQ1;

# Build profiles from the test alignments
@output = `$builddir/src/bathbuild $tmppfx.bhmm $tmppfx.sto 2>&1`;
if ($? != 0) { die "FAIL: bathbuild failed\n"; }

# bathsearch should show four results
$output = `$builddir/src/bathsearch --tblout $tmppfx.tbl $tmppfx.bhmm $tmppfx.fa 2>&1`;
if ($? != 0) { die "FAIL: bathsearch failed\n"; }

&h3::ParseTbl("$tmppfx.tbl");
if ($h3::ntbl != 4) { die "FAIL: on expected number of hits, bathsearch\n"; } 


print "ok\n";
unlink "$tmppfx.sto";
unlink "$tmppfx.fa";
unlink "$tmppfx.tbl";
unlink <$tmppfx.bhmm*>;
exit 0;





