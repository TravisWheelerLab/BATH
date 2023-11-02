#! /usr/bin/perl

# Test the ability of bathbuild to deal with crappy alignments
# of lots of sequence fragments.
#
# Usage:    ./i7-bathbuild-fragments.pl <bathbuild binary> <tmpfile prefix>
# Example:  ./i7-bathbuild-fragments.pl ../src/bathbuild   tmpfoo
# 
# SRE, Tue Jun 16 13:37:05 2009


$bathbuild = shift;
$tmppfx   = shift;

if (! -x "$bathbuild")           { die "FAIL: didn't find bathbuild binary $bathbuild\n";  }
open (TMPFILE, ">$tmppfx.sto") || die "FAIL: couldn't open $tmppfx.sto for writing\n";   
print TMPFILE << "EOF";
# STOCKHOLM 1.0

#=GF ID test

seq1 ACDEFGHIKL------------------------------
seq2 ----------MNPQRSTVWY--------------------
seq3 --------------------ACDEFGHIKL----------
seq4 ------------------------------MNPQRSTVWY
//
EOF
close TMPFILE;


$output = `$bathbuild -O $tmppfx.sto2 $tmppfx.bhmm $tmppfx.sto 2>&1`;
if ($? != 0)   { die "FAIL: bathbuild failed unexpectedly\n"; }

$output =~ /1\s+test\s+4\s+40\s+(\d+)/;
if ($1 != 40)  { die "FAIL: should've built a M=40 model\n"; }


$output = `$bathbuild --fragthresh 0.0 $tmppfx.bhmm $tmppfx.sto 2>&1`;
if ($? == 0)   { die "FAIL: bathbuild should have failed but didn't\n"; }


print "ok\n"; 
unlink "$tmppfx.sto";
unlink "$tmppfx.sto2";
unlink "$tmppfx.bhmm";
exit 0;
