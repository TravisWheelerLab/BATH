#! /usr/bin/perl

# Regression test of the per-exon table (--exontblout) from a spliced
# bathsearch. Every column of the table is checked against values from
# a manually checked run, for a two-exon hit and a four-exon hit.
# This catches errors in exon scoring that leave the whole-hit
# score and E-value untouched: exon scores and P-values are computed
# from the Forward matrix, and reusing that matrix for another purpose
# before the exons are scored gives wrong exon P-values.
#
# Usage:    ./i18-splice-exons.pl <builddir> <srcdir> <tmpfile prefix>
# Example:  ./i18-splice-exons.pl ..         ..        tmpfoo

BEGIN {
    $builddir = shift;
    $srcdir   = shift;
    $tmppfx   = shift;
}
use lib "$srcdir/testsuite";
use h3;

$bathsearch = "$builddir/src/bathsearch";
if (! -x $bathsearch) { die "FAIL: didn't find bathsearch binary in $builddir/src\n"; }

# Columns of each expected exon: hit ID, target name, target accession, query
# name, query accession, hmm len, seq len, full-hit E-value, score, bias,
# exon number, number of exons, hmm from, hmm to, ali from, ali to,
# exon P-value, exon PID, splice signal.
#
# Floating point columns (E-values, scores, biases, P-values) are compared with a
# tolerance, because they are printed to only two or three digits and the last
# digit can differ between platforms. Everything else must match exactly.
$tmem = [
  [ 1, "test_splice", "-", "tmem-258", "-",  81,  545, "4.9e-52", "162.4", "4.2", 1, 2,   1,  42,   1, 126, "1.7e-27", "100.00", "gtag" ],
  [ 1, "test_splice", "-", "tmem-258", "-",  81,  545, "4.9e-52", "162.4", "4.2", 2, 2,  43,  81, 245, 361,   "2e-30", "100.00", "----" ],
];

$pthr = [
  [ 1, "seq1", "-", "PTHR37536", "-", 279, 1300, "2.8e-28", "87.9", "5.2", 1, 4,  11, 135, 119,  491, "9.6e-19", "30.16", "gtag" ],
  [ 1, "seq1", "-", "PTHR37536", "-", 279, 1300, "2.8e-28", "87.9", "5.2", 2, 4, 136, 143, 577,  600, "7.3e-05", "12.50", "gtag" ],
  [ 1, "seq1", "-", "PTHR37536", "-", 279, 1300, "2.8e-28", "87.9", "5.2", 3, 4, 144, 180, 687,  798, "4.1e-10", "26.32", "gtag" ],
  [ 1, "seq1", "-", "PTHR37536", "-", 279, 1300, "2.8e-28", "87.9", "5.2", 4, 4, 181, 251, 952, 1159, "4.9e-12", "34.72", "----" ],
];

# spliced search of $hmm against $seqdb; check every column of the exon table
sub check_exons {
    my ($label, $hmm, $seqdb, $expected) = @_;

    system("$bathsearch --splice --exontblout $tmppfx.extbl -o $tmppfx.out $hmm $seqdb > $tmppfx.log 2>&1");
    if ($? != 0) { die "FAIL: bathsearch --splice failed on $label\n"; }
    &h3::ParseExonTbl("$tmppfx.extbl");

    if ($h3::nex != scalar(@$expected)) { printf("FAIL: %s: expected %d exons in table; saw %d\n", $label, scalar(@$expected), $h3::nex); exit 1; }

    for ($i = 0; $i < $h3::nex; $i++) {
	@got = ($h3::exhit[$i],   $h3::extname[$i],    $h3::extacc[$i],   $h3::exqname[$i],   $h3::exqacc[$i],
		$h3::exhmmlen[$i], $h3::exseqlen[$i],  $h3::exfullE[$i],  $h3::exfullsc[$i], $h3::exfullbias[$i],
		$h3::exnum[$i],   $h3::exof[$i],       $h3::exhmmfrom[$i], $h3::exhmmto[$i],  $h3::exalifrom[$i],
		$h3::exalito[$i], $h3::exP[$i],        $h3::expid[$i],    $h3::exsplice[$i]);
	@want = @{$expected->[$i]};

	# column index => [ column name, tolerance type ]
	%fuzzy = ( 7 => ["full-hit E-value", "rel"], 8 => ["full-hit score", "abs"], 9 => ["full-hit bias", "abs"], 16 => ["exon P-value", "rel"] );

	for ($c = 0; $c < scalar(@want); $c++) {
	    if (exists $fuzzy{$c}) {
		($name, $type) = @{$fuzzy{$c}};
		if ($type eq "rel") { $ok = (abs($got[$c] - $want[$c]) <= 0.1 * abs($want[$c])); }
		else                { $ok = (abs($got[$c] - $want[$c]) <= 0.2); }
	    }
	    else {
		$name = ("hit ID", "target name", "target accession", "query name", "query accession", "hmm len", "seq len",
			 "", "", "", "exon number", "number of exons", "exon hmm from", "exon hmm to", "exon ali from",
			 "exon ali to", "", "exon PID", "splice signal")[$c];
		$ok = ($got[$c] eq $want[$c]);
	    }
	    if (! $ok) { printf("FAIL: %s: exon %d: expected %s of %s; saw %s\n", $label, $i+1, $name, $want[$c], $got[$c]); exit 1; }
	}
    }
}

&check_exons("tmem-258",  "$srcdir/testsuite/tmem-258.bhmm", "$srcdir/testsuite/tmem-258.fa", $tmem);
&check_exons("PTHR37536", "$srcdir/testsuite/PTHR37536.bhmm", "$srcdir/testsuite/PTHR37536-nt.fa", $pthr);

print "ok\n";
unlink "$tmppfx.extbl";
unlink "$tmppfx.out";
unlink "$tmppfx.log";
exit 0;
