#! /usr/bin/perl

# Regression test of the CIGAR strings that bathsearch --cigar puts in
# its tabular output, in all three alignment modes: translated search
# (no flag), --fs, and --splice.
#
# Rather than compare against fixed strings, every hit's CIGAR is
# checked against the other columns of its own table row. Counts are in
# nucleotides:
#
#   all modes:  M + D (+ B with --fs)  =  3 * (hmm to - hmm from + 1)
#               M + I (+ F with --fs, + N with --splice)
#                                      =  |ali to - ali from| + 1
#
# --fs adds F (insertion frameshift) and B (deletion frameshift) ops, and
# the "shifts" column must equal the number of F and B ops. --splice
# adds N (intron) ops; there must be one per intron, and each must equal
# the intron length implied by the neighboring rows of the --exontblout
# table. These checks catch a count that is wrong for any reason, such as
# a counter that is never initialized.
#
# Usage:    ./i19-cigar.pl <builddir> <srcdir> <tmpfile prefix>
# Example:  ./i19-cigar.pl ..         ..        tmpfoo

BEGIN {
    $builddir = shift;
    $srcdir   = shift;
    $tmppfx   = shift;
}

$bathsearch = "$builddir/src/bathsearch";
if (! -x $bathsearch) { die "FAIL: didn't find bathsearch binary in $builddir/src\n"; }

$oxyhmm = "$srcdir/testsuite/2OG-FeII_Oxy_3.bhmm";
$tmemhmm = "$srcdir/testsuite/tmem-258.bhmm";
$pthrhmm = "$srcdir/testsuite/PTHR37536.bhmm";

# Test cases: label, mode, extra bathsearch options, query, target.
#   2OG-FeII_Oxy_3: 10 hits, on both strands, with insertions and deletions;
#                   the -fs target also has frameshifts and a stop codon
#   tmem-258:       one hit with one intron
#   PTHR37536:      one hit with three introns, on the forward strand and
#                   (the -rc target is its reverse complement) the reverse strand
@cases = (
  [ "2OG-FeII_Oxy_3",      "std",    "",         $oxyhmm,  "$srcdir/testsuite/2OG-FeII_Oxy_3-nt.fa"    ],
  [ "2OG-FeII_Oxy_3 --fs", "fs",     "--fs",     $oxyhmm,  "$srcdir/testsuite/2OG-FeII_Oxy_3-nt-fs.fa" ],
  [ "tmem-258 --splice",   "splice", "--splice", $tmemhmm, "$srcdir/testsuite/tmem-258.fa"             ],
  [ "PTHR37536 --splice",  "splice", "--splice", $pthrhmm, "$srcdir/testsuite/PTHR37536-nt.fa"         ],
  [ "PTHR37536-rc --splice", "splice", "--splice", $pthrhmm, "$srcdir/testsuite/PTHR37536-nt-rc.fa"    ],
);

%allowed = ( "std" => "MID", "fs" => "MIDFB", "splice" => "MIDN" );

sub fail { print "FAIL: @_\n"; exit 1; }

# read a bathsearch table: returns a list of array refs, one per hit
sub read_table {
    my ($file) = @_;
    my @rows;
    open(TBL, $file) || fail("couldn't open $file");
    while (<TBL>) {
	if (/^\#/) { next; }
	s/\s+$//;
	if ($_ eq "") { next; }
	push @rows, [ split(' ', $_) ];
    }
    close TBL;
    return @rows;
}

$total_hits = 0;
foreach $case (@cases) {
    ($label, $mode, $opts, $hmm, $db) = @$case;

    $cmd = "$bathsearch $opts --cigar --tblout $tmppfx.tbl";
    if ($mode eq "splice") { $cmd .= " --exontblout $tmppfx.extbl"; }
    system("$cmd -o $tmppfx.out $hmm $db > $tmppfx.log 2>&1");
    if ($? != 0) { fail("bathsearch failed on $label"); }

    @hits = read_table("$tmppfx.tbl");
    if (scalar(@hits) < 1) { fail("$label: found no hits"); }

    # exon table, by hit ID, exons in order
    %exons = ();
    if ($mode eq "splice") {
	foreach $r (read_table("$tmppfx.extbl")) {
	    push @{$exons{$r->[0]}}, [ $r->[10], $r->[14], $r->[15] ];   # exon number, ali from, ali to
	}
    }

    $shifted = 0;
    $introns = 0;
    foreach $r (@hits) {
	$hitid    = $r->[0];
	$hmmfrom  = $r->[6];
	$hmmto    = $r->[7];
	$alifrom  = $r->[9];
	$alito    = $r->[10];
	$cigar    = $r->[-1];
	$where    = sprintf("%s: hit %d (CIGAR %s)", $label, $hitid, $cigar);

	if ($cigar !~ /^(\d+[A-Z])+$/) { fail("$where: not a CIGAR string"); }

	%len = (); %nops = (); @nlens = ();
	foreach $op (split('', "MIDNFB")) { $len{$op} = 0; $nops{$op} = 0; }
	while ($cigar =~ /(\d+)([A-Z])/g) {
	    ($n, $op) = ($1, $2);
	    if (index($allowed{$mode}, $op) < 0) { fail("$where: unexpected operation $op in $mode mode"); }
	    if ($n < 1)                          { fail("$where: zero-length $op operation"); }
	    $len{$op} += $n;
	    $nops{$op}++;
	    if ($op eq "N") { push @nlens, $n; }
	}

	$modelnt  = $len{M} + $len{D} + $len{B};
	$targetnt = $len{M} + $len{I} + $len{F} + $len{N};
	$wantmodel  = 3 * ($hmmto - $hmmfrom + 1);
	$wanttarget = abs($alito - $alifrom) + 1;

	if ($modelnt  != $wantmodel)  { fail("$where: M+D+B is $modelnt nucleotides; the model span (hmm from/to) needs $wantmodel"); }
	if ($targetnt != $wanttarget) { fail("$where: M+I+F+N is $targetnt nucleotides; the target span (ali from/to) is $wanttarget"); }

	if ($mode eq "fs") {
	    $shifts = $r->[15];
	    if ($shifts != $nops{F} + $nops{B}) { fail("$where: shifts column is $shifts; the CIGAR has " . ($nops{F} + $nops{B}) . " frameshift operations"); }
	    $shifted += $nops{F} + $nops{B};
	}

	if ($mode eq "splice") {
	    @ex = @{$exons{$hitid}};
	    if (scalar(@ex) != $r->[11])           { fail("$where: exon count is $r->[11]; the exon table has " . scalar(@ex) . " rows"); }
	    if (scalar(@nlens) != scalar(@ex) - 1) { fail("$where: " . scalar(@nlens) . " intron operations for " . scalar(@ex) . " exons"); }
	    for ($e = 0; $e < scalar(@nlens); $e++) {
		$intron = abs($ex[$e+1][1] - $ex[$e][2]) - 1;
		if ($nlens[$e] != $intron) { fail("$where: intron " . ($e+1) . " is $nlens[$e]N; the exon table gives $intron"); }
	    }
	    $introns += scalar(@nlens);
	}
	$total_hits++;
    }

    # make sure the case exercised what it is here to exercise
    if ($mode eq "fs"     && $shifted < 1) { fail("$label: no frameshift operations in any CIGAR; test data no longer exercises --fs"); }
    if ($mode eq "splice" && $introns < 1) { fail("$label: no intron operations in any CIGAR; test data no longer exercises --splice"); }
}

print "ok\n";
unlink "$tmppfx.tbl";
unlink "$tmppfx.extbl";
unlink "$tmppfx.out";
unlink "$tmppfx.log";
exit 0;
