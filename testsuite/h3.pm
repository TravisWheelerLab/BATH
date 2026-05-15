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
	@fields = split(' ', $_, 16);

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
    $fullE[$ntbl]     = $fields[11];
    $fullsc[$ntbl]    = $fields[12];
    $fullbias[$ntbl]  = $fields[13];
    $pid[$ntbl]       = $fields[14];
    $tdesc[$ntbl]     = $fields[15];
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
    $fullE[$ntbl]     = $fields[11];
    $fullsc[$ntbl]    = $fields[12];
    $fullbias[$ntbl]  = $fields[13];
    $pid[$ntbl]       = $fields[14];
    $shifts[$ntbl]    = $fields[15];
    $stops[$ntbl]     = $fields[16];
    $tdesc[$ntbl]     = $fields[17];
    $ntbl++;
    }
    close TBLFILE;
    1;
}

1;
