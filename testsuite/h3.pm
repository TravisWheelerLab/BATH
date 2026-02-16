#! /usr/bin/perl

package h3;


sub ParseTbl {
    my ($tblfile)    = @_;
    my (@fields);

    $ntbl     = 0;
    @tname    = ();
    @tacc     = ();
    @qname    = ();
    @qacc     = ();
    @hmmlen   = ();
    @hmmfrom  = ();
    @hmmto    = ();
    @seqlen   = ();
    @alifrom  = ();
    @alito    = ();
    @envfrom  = ();
    @envto    = ();
    @fullE    = ();
    @fullsc   = ();
    @fullbias = ();
    @pid      = ();
    @tdesc    = ();

    if (! open(TBLFILE, $tblfile)) { print "FAIL: couldn't open table file"; exit 1 ; }
    while (<TBLFILE>)
    {
	if (/^\#/) { next; }
	chop;
	@fields = split(' ', $_, 18);

	$tname[$ntbl]     = $fields[1];
	$tacc[$ntbl]      = $fields[2];
	$qname[$ntbl]     = $fields[3];
	$qacc[$ntbl]      = $fields[4];
    $hmmlen[$ntbl]    = $fields[5];
    $hmmfrom[$ntbl]   = $fields[6];
    $hmmto[$ntbl]     = $fields[7];
    $seqlen[$ntbl]    = $fields[8];
    $alifrom[$ntbl]   = $fields[9];
    $alito[$ntbl]     = $fields[10];
    $envfrom[$ntbl]   = $fields[11];
    $envto[$ntbl]     = $fields[12];
	$fullE[$ntbl]     = $fields[13];
	$fullsc[$ntbl]    = $fields[14];
	$fullbias[$ntbl]  = $fields[15];
    $pid[$ntbl]       = $fields[16];
	$tdesc[$ntbl]     = $fields[17];
	$ntbl++;
    }
    close TBLFILE;
    1;
}

sub ParseFSTbl {
    my ($tblfile)    = @_;
    my (@fields);

    $ntbl     = 0;
    @tname    = ();
    @tacc     = ();
    @qname    = ();
    @qacc     = ();
    @hmmlen   = ();
    @hmmfrom  = ();
    @hmmto    = ();
    @seqlen   = ();
    @alifrom  = ();
    @alito    = ();
    @envfrom  = ();
    @envto    = ();
    @fullE    = ();
    @fullsc   = ();
    @fullbias = ();
    @pid      = ();
    @shifts   = ();
    @stops    = ();
    @tdesc    = ();

    if (! open(TBLFILE, $tblfile)) { print "FAIL: couldn't open table file"; exit 1 ; }
    while (<TBLFILE>)
    {
	if (/^\#/) { next; }
	chop;
	@fields = split(' ', $_, 20);
	$tname[$ntbl]     = $fields[1];
	$tacc[$ntbl]      = $fields[2];
	$qname[$ntbl]     = $fields[3];
	$qacc[$ntbl]      = $fields[4];
    $hmmlen[$ntbl]    = $fields[5];
    $hmmfrom[$ntbl]   = $fields[6];
    $hmmto[$ntbl]     = $fields[7];
    $seqlen[$ntbl]    = $fields[8];
    $alifrom[$ntbl]   = $fields[9];
    $alito[$ntbl]     = $fields[10];
    $envfrom[$ntbl]   = $fields[11];
    $envto[$ntbl]     = $fields[12];
	$fullE[$ntbl]     = $fields[13];
	$fullsc[$ntbl]    = $fields[14];
	$fullbias[$ntbl]  = $fields[15];
    $pid[$ntbl]       = $fields[16];
    $shifts[$ntbl]    = $fields[17];
    $stops[$ntbl]     = $fields[18];
	$tdesc[$ntbl]     = $fields[19];
	$ntbl++;
    }
    close TBLFILE;
    1;
}

1;
