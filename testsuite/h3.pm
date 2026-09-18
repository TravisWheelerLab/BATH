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

# ParseExonTbl(): parse the per-exon table saved by bathsearch --exontblout.
# Fills @exhit, @extname, @extacc, @exqname, @exqacc, @exhmmlen, @exseqlen,
# @exfullE, @exfullsc, @exfullbias, @exnum, @exof, @exhmmfrom, @exhmmto,
# @exalifrom, @exalito, @exP, @expid, @exsplice; $nex is the number of exons.
sub ParseExonTbl {
    my ($tblfile)    = @_;
    my (@fields);

    $nex        = 0;
    @exhit      = ();
    @extname    = ();
    @extacc     = ();
    @exqname    = ();
    @exqacc     = ();
    @exhmmlen   = ();
    @exseqlen   = ();
    @exfullE    = ();
    @exfullsc   = ();
    @exfullbias = ();
    @exnum      = ();
    @exof       = ();
    @exhmmfrom  = ();
    @exhmmto    = ();
    @exalifrom  = ();
    @exalito    = ();
    @exP        = ();
    @expid      = ();
    @exsplice   = ();

    if (! open(TBLFILE, $tblfile)) { print "FAIL: couldn't open exon table file"; exit 1 ; }
    while (<TBLFILE>)
    {
	if (/^\#/) { next; }
	s/\s+$//;
	@fields = split(' ', $_, 19);

	$exhit[$nex]      = $fields[0];
	$extname[$nex]    = $fields[1];
	$extacc[$nex]     = $fields[2];
	$exqname[$nex]    = $fields[3];
	$exqacc[$nex]     = $fields[4];
	$exhmmlen[$nex]   = $fields[5];
	$exseqlen[$nex]   = $fields[6];
	$exfullE[$nex]    = $fields[7];
	$exfullsc[$nex]   = $fields[8];
	$exfullbias[$nex] = $fields[9];
	$exnum[$nex]      = $fields[10];
	$exof[$nex]       = $fields[11];
	$exhmmfrom[$nex]  = $fields[12];
	$exhmmto[$nex]    = $fields[13];
	$exalifrom[$nex]  = $fields[14];
	$exalito[$nex]    = $fields[15];
	$exP[$nex]        = $fields[16];
	$expid[$nex]      = $fields[17];
	$exsplice[$nex]   = $fields[18];
	$nex++;
    }
    close TBLFILE;
    1;
}

1;
