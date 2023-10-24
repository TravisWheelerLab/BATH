#! /usr/bin/perl

# Check that bathsearch  can deal with HMMs with no optional annotation
# Bug #h69 segfaulted on this test.
#
# Usage:   ./i9-optional-annotation.pl <builddir> <srcdir> <tmpfile prefix>
# Example: ./i9-optional-annotation.pl ..         ..       tmpfoo
#
# SRE, Sun Nov 29 11:49:39 2009

BEGIN {
    $builddir  = shift;
    $srcdir    = shift;
    $tmppfx    = shift;
}
use lib "$srcdir/testsuite";
use h3;

# Verify that we have all the executables we need for the test.
if (! -x "$builddir/src/bathbuild")   { die "FAIL: didn't find bathbuild binary in $builddir/src\n";  }
if (! -x "$builddir/src/bathsearch")  { die "FAIL: didn't find bathsearch binary in $builddir/src\n"; }


# Create our test files.
if (! open(ALI1, ">$tmppfx.sto")) { die "FAIL: couldn't open $tmppfx.sto for write\n";  }
if (! open(SEQ1, ">$tmppfx.seq")) { die "FAIL: couldn't open $tmppfx.seq for write\n";  }

print ALI1 <<"EOF";
# STOCKHOLM 1.0
#=GF ID ali1 
#=GF AC XX01234.5
#=GF DE A test description
seq1 ACDEFGHIKLMNPQRSTVWY
seq2 ACDEFGHIKLMNPQRSTVWY
seq3 ACDEFGHIKLMNPQRSTVWY
//
# STOCKHOLM 1.0
#=GF ID ali2
seq1 ACDEFGHIKLMNPQRSTVWY
seq2 ACDEFGHIKLMNPQRSTVWY
seq3 ACDEFGHIKLMNPQRSTVWY
//
EOF

print SEQ1 <<"EOF";
ID   test1   STANDARD;  PRT;  20 AA.
AC   AC00001;
DE   Sequence description
SQ   SEQUENCE   20 AA; 99999 MW;  FFFFFFFFFFFFFFFF CRC64;
     GCATGTGACGAGTTTGGCCATATAAAACTTATGAATCCACAGCGCTCAACTGTATGGTAT
//
ID   test2   STANDARD;  PRT;  20 AA.
SQ   SEQUENCE   20 AA; 99999 MW;  FFFFFFFFFFFFFFFF CRC64;
     GCATGTGACGAGTTTGGCCATATAAAACTTATGAATCCACAGCGCTCAACTGTATGGTAT
//
EOF

close ALI1;
close SEQ1;

@output = `$builddir/src/bathbuild $tmppfx.bhmm $tmppfx.sto 2>&1`;
if ($? != 0) { die "FAIL: bathbuild failed\n"; }

@output = `$builddir/src/bathsearch --tblout $tmppfx.tbl $tmppfx.bhmm $tmppfx.seq 2>&1`;
if ($? != 0) { die "FAIL: bathsearch failed\n"; }

&h3::ParseTbl("$tmppfx.tbl");
if ($h3::ntbl       != 4)                      { die "FAIL: on expected number lines\n";   }
if ($h3::tname[0]   ne "test1")                { die "FAIL: on line 0 target name\n";      }
if ($h3::tacc[0]    ne "AC00001")              { die "FAIL: on line 0 accession\n";        }
if ($h3::tdesc[0]   ne "Sequence description") { die "FAIL: on line 0 desc\n";             }
if ($h3::qname[0]   ne "ali1")                 { die "FAIL: on line 0 query name\n";       }
if ($h3::qacc[0]    ne "XX01234.5")            { die "FAIL: on line 0 query accession\n";  }
if ($h3::tname[1]   ne "test2")                { die "FAIL: on line 1 target name\n";      }
if ($h3::tacc[1]    ne "-")                    { die "FAIL: on line 1 accession\n";        }
if ($h3::tdesc[1]   ne "-")                    { die "FAIL: on line 1 desc\n";             }
if ($h3::qname[2]   ne "ali2")                 { die "FAIL: on line 2 query name\n";       }
if ($h3::qacc[2]    ne "-")                    { die "FAIL: on line 2 query accession\n";  }


print "ok\n";
unlink "$tmppfx.sto";
unlink "$tmppfx.seq";
unlink "$tmppfx.tbl";
unlink <$tmppfx.bhmm*>;
exit 0;
