#! /usr/bin/perl

# Look for any problems in bathalign that corrupt the input sequences.
# Based on HMMER's i13-msa-integrity.pl, adapted for bathemit + bathalign.
#
# Usage:   ./i13-align-integrity.pl <builddir> <srcdir> <tmpfile prefix>
# Example: ./i13-align-integrity.pl ..         ..       tmpfoo

$builddir  = shift;
$srcdir    = shift;
$tmppfx    = shift;

# Verify that we have all the executables we need for the test.
if (! -x "$builddir/src/bathalign")                 { die "FAIL: didn't find bathalign binary in $builddir/src";  }
if (! -x "$builddir/src/bathemit")                  { die "FAIL: didn't find bathemit binary in $builddir/src";  }
if (! -x "$builddir/easel/miniapps/esl-reformat")   { die "FAIL: didn't find esl-reformat binary in $builddir/easel/miniapps";  }
if (! -x "$builddir/easel/miniapps/esl-shuffle")    { die "FAIL: didn't find esl-shuffle binary in $builddir/easel/miniapps";  }

# Verify that we have the datafile we need.
if (! -e "$srcdir/testsuite/RRM_1.bhmm")  { die "FAIL: didn't find RRM_1.bhmm in $srcdir/testsuite";  }
$profile = "$srcdir/testsuite/RRM_1.bhmm";

foreach $trial (1..5)
{
    foreach $n (1, 10, 100)
    {
	# homologous sequence fragments: generated from local profile
	`$builddir/src/bathemit -o $tmppfx.fa -N $n -L 0 -p --unilocal $profile`;
	if ($? != 0) { die "FAIL: bathemit"; }

	&align_integrity_check("$tmppfx.fa", $profile);

	# random sequences
	`$builddir/easel/miniapps/esl-shuffle -G -N $n -L 50 --amino -o $tmppfx.fa`;
	if ($? != 0) { die "FAIL: esl-shuffle"; }

	&align_integrity_check("$tmppfx.fa", $profile);
    }
}

print "ok\n";
unlink "$tmppfx.sto";
unlink <$tmppfx.fa*>;
exit 0;


sub align_integrity_check
{
    my ($fafile, $hmmfile) = @_;

    `$builddir/src/bathalign -o $tmppfx.sto $hmmfile $fafile > /dev/null 2>&1`;
    if ($? != 0) { die "FAIL: bathalign failed"; }

    `$builddir/easel/miniapps/esl-reformat -u fasta $tmppfx.sto > $tmppfx.fa1 2>/dev/null`;
    if ($? != 0) { die "FAIL: first esl-reformat failed"; }

    `$builddir/easel/miniapps/esl-reformat -u fasta $fafile    > $tmppfx.fa2 2>/dev/null`;
    if ($? != 0) { die "FAIL: second esl-reformat failed"; }

    `diff -b $tmppfx.fa1 $tmppfx.fa2 > /dev/null 2>&1`;
    if ($? != 0) { die "FAIL: alignment corrupted\n"; }
    0;
}
