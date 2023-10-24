#! /usr/bin/perl

# Test that programs accept and reject argument of '-' (for reading
# data from stdin, rather than from files) as they're supposed to.
#
# Usage:   ./i17-stdin.pl <builddir> <srcdir> <tmpfile prefix>
# Example: ./i17-stdin.pl ..         ..       tmpfoo
#
# SRE, Wed Oct 27 13:05:10 2010 [Janelia]


BEGIN {
    $builddir  = shift;
    $srcdir    = shift;
    $tmppfx    = shift;
}

$verbose = 0;

# The test makes use of the following files:
#
# $model1.bhmm          <hmmfile>  Single RRM_1 model - BATH format
# $model2.bhmm          <hmmfile>  Single Caudal_act model - BATH format
# $model1.hmm           <hmmfile>  Single RRM_1 model - HMMER format
# $model2.hmm           <hmmfile>  Single Caudal_act model - HMMER format
# $model1.sto           <msafile>  Single RRM_1 alignment
# $model2.sto           <msafile>  Single Caudal_act alignment
# $model1-10.fa         <seqfile>  reverse translation of 10 seqs from RRM_1 model
# $model2-10.fa         <seqfile>  reverse translation 10 seqs from Caudal_act model

# It creates the following files:
# $tmppfx.bhmm          <hmmdb>    2 models, RRM_1 + Caudal_act - BATH format
# $tmppfx.hmm           <hmmdb>    2 models, RRM_1 + Caudal_act - HMMER format
# $tmppfx.sto           <msadb>    2 MSAs, RRM_1 + Caudal_act
# $tmppfx.db            <seqdb>   10 RRM_1 + 10 Caudal_act + 100 random seqs 
# $tmppfx.key           <keyfile>  "Caudal_act" in a file, to test bathfetch -f 
#

# All models assumed to be in testsuite subdirectory.
$model1   = "RRM_1";
$model2   = "Caudal_act";

@h3progs =  ("bathbuild", "bathconvert", "bathfetch", "bathsearch", "bathstat");
@eslprogs = ("esl-shuffle");

# Verify that we have all the executables and datafiles we need for the test.
foreach $h3prog  (@h3progs) { if (! -x "$builddir/src/$h3prog")             { die "FAIL: didn't find $h3prog executable in $builddir/src\n";              } }
foreach $eslprog (@eslrogs) { if (! -x "$builddir/easel/miniapps/$eslprog") { die "FAIL: didn't find $eslprog executable in $builddir/easel/miniapps\n";  } }

if (! -r "$srcdir/testsuite/$model1.bhmm")  { die "FAIL: can't read profile $model1.bhmm in $srcdir/testsuite\n"; }
if (! -r "$srcdir/testsuite/$model2.bhmm")  { die "FAIL: can't read profile $model2.bhmm in $srcdir/testsuite\n"; }

if (! -r "$srcdir/testsuite/$model1.sto")  { die "FAIL: can't read msa $model1.sto in $srcdir/testsuite\n"; }
if (! -r "$srcdir/testsuite/$model2.sto")  { die "FAIL: can't read msa $model2.sto in $srcdir/testsuite\n"; }


`cat $srcdir/testsuite/$model1.bhmm  $srcdir/testsuite/$model2.bhmm  > $tmppfx.bhmm`; if ($?) { die "FAIL: cat\n"; }
`cat $srcdir/testsuite/$model1.hmm  $srcdir/testsuite/$model2.hmm  > $tmppfx.hmm`;    if ($?) { die "FAIL: cat\n"; }

`cat $srcdir/testsuite/$model1.sto $srcdir/testsuite/$model2.sto > $tmppfx.sto`;      if ($?) { die "FAIL: cat\n"; }

`cat $srcdir/testsuite/$model1-10.fa $srcdir/testsuite/$model2-10.fa > $tmppfx.db`;   if ($?) { die "FAIL: cat\n"; }
`$builddir/easel/miniapps/esl-shuffle -G -N100 -L 1200 --dna >> $tmppfx.db`;          if ($?) { die "FAIL: esl-shuffle\n"; }

`echo $model1    > $tmppfx.key`;                                                      if ($?) { die "FAIL: cat\n"; }
`echo $model2   >> $tmppfx.key`;                                                      if ($?) { die "FAIL: cat\n"; }


################################################################
# bathbuild
#    don't diff HMM files, they may fail because of DATE field
#    reject - for <hmmfile>: can't send it to stdout.
################################################################

$tag  = "bathbuild";      $prog = "$builddir/src/$tag";      
$tag1 = "<msafile>";      $arg1 = "$tmppfx.sto";    
if ($verbose) { print "$tag...\n"; }

`$prog $tmppfx.bhmm.out1 $arg1                              | grep -v "^#" > $tmppfx.out1`;   if ($?) { die "FAIL: $tag <hmmfile> $tag1 \n"; }
`cat $arg1 | $prog --informat stockholm $tmppfx.bhmm.out2 - | grep -v "^#" > $tmppfx.out2`;   if ($?) { die "FAIL: $tag <hmmfile> -\n"; }
`diff -b $tmppfx.out1     $tmppfx.out2     2>&1 > /dev/null`; if ($?) { die "FAIL: $tag results differ if $tag1 comes through stdin\n"; }

$output = `$prog - $arg1`;
if (!$?) { die "FAIL: $tag should reject - for <hmmfile_out>\n"; }

################################################################
# bathconvert
#   reject - for <hmmfile>
################################################################

$tag  = "bathconvert";     $prog = "$builddir/src/$tag";    
$tag1 = "<bhmmfile>";      $arg1 = "$tmppfx.bhmm";    
$tag2 = "<hmmfile>";       $arg2 = "$tmppfx.hmm";
if ($verbose) { print "$tag...\n"; }

`$prog $arg1 $arg2   > $tmppfx.out1`;     if ($?)  { die "FAIL: $tag $tag1\n"; }
`cat $arg2 | $prog $arg1 > $tmppfx.out2`; if (!$?) { die "FAIL: $tag -\n"; }

################################################################
# bathfetch 
#    need to check all three use modes, including -f and --index
#    --index rejects -
#    w/ -f, only one of <hmmfile>, <keyfile> can be -
#    -f fetches in different orders depending on whether file is
#      indexed or not, so <keyfile> must be constructed to give
#      same fetch order either way.
################################################################

$tag  = "bathfetch";    $prog = "$builddir/src/$tag";   
$tag1 = "<hmmfile>";    $arg1 = "$tmppfx.bhmm";              
$tag2 = "<keyfile>";    $arg2 = "$tmppfx.key";              
if ($verbose) { print "$tag...\n"; }

`$prog $arg1 Caudal_act         > $tmppfx.out1`;          if ($?) { die "FAIL: $tag $tag1\n"; }
`cat $arg1 | $prog - Caudal_act > $tmppfx.out2`;          if ($?) { die "FAIL: $tag -\n"; }
`diff -b $tmppfx.out1 $tmppfx.out2 2>&1 > /dev/null`;  if ($?) { die "FAIL: $tag results differ if $tag1 comes through stdin\n"; }

`$prog -f $arg1 $arg2           > $tmppfx.out1`;          if ($?) { die "FAIL: $tag -f $tag1 $tag2\n"; }
`cat $arg1 | $prog -f - $arg2   > $tmppfx.out2`;          if ($?) { die "FAIL: $tag -f - $tag2\n"; }
`cat $arg2 | $prog -f $arg1 -   > $tmppfx.out3`;          if ($?) { die "FAIL: $tag -f $tag1 -\n"; }
`diff -b $tmppfx.out1 $tmppfx.out2 2>&1 > /dev/null`;  if ($?) { die "FAIL: $tag -f results differ if $tag1 comes through stdin\n"; }
`diff -b $tmppfx.out1 $tmppfx.out3 2>&1 > /dev/null`;  if ($?) { die "FAIL: $tag -f results differ if $tag2 comes through stdin\n"; }

$output = `cat $arg1 $arg2 | $prog -f - - 2>&1`;
if (! $?) { die "FAIL: $tag should have failed on double - -\n"; }
if ($output !~ /^Either <hmmfile> or <keyfile>/) { die "FAIL: $tag didn't give expected error message for the - - case.\n"; }

`$prog --index $arg1            > $tmppfx.out1`;   if ($?)   { die "FAIL: $tag --index $tag1\n"; }
$output = `cat $arg1 | $prog --index - 2>&1`;      if (! $?) { die "FAIL: $tag should reject - for <hmmfile> when using --index\n"; }
if ($output !~ /^Can't use - with --index/) { die "FAIL: $tag didn't give expected error message for the - - case.\n"; }

#################################################################
# bathsearch
#      reject - - case 
#      reject <seqdb> as - on multiquery
#################################################################
# note that the grep -v "^#" removes lines that would make diffs fail,
# like query name and cpu time.

$tag  = "bathsearch";        $prog = "$builddir/src/$tag";   
$tag1 = "<hmmfile>";         $arg1 = "$srcdir/testsuite/$model1.bhmm";   $arg1b = "$tmppfx.bhmm";   
$tag2 = "<seqdb>";           $arg2 = "$tmppfx.db";      
if ($verbose) { print "$tag...\n"; }

`$prog $arg1 $arg2          | grep -v "^#" > $tmppfx.out1`;  if ($?) { die "FAIL: $tag $tag1 $tag2\n"; }
`cat $arg1 | $prog - $arg2  | grep -v "^#" > $tmppfx.out2`;  if ($?) { die "FAIL: $tag - $tag2\n"; }
`cat $arg2 | $prog $arg1 -  | grep -v "^#" > $tmppfx.out3`;  if ($?) { die "FAIL: $tag $tag1 -\n"; }

`diff -b $tmppfx.out1 $tmppfx.out2 2>&1 > /dev/null`;  if ($?) { die "FAIL: $prog results differ if $tag1 comes through stdin\n"; }
`diff -b $tmppfx.out1 $tmppfx.out3 2>&1 > /dev/null`;  if ($?) { die "FAIL: $prog results differ if $tag2 comes through stdin\n"; }

$output = `cat $arg1 $arg2 | $prog - - 2>&1`;    if (! $?) { die "FAIL: $prog should have failed on double - -\n"; }
if ($output !~ /^Either <hmmfile> or <seqdb>/) { die "FAIL: $prog didn't give expected error message for the - - case.\n"; }

$output = `cat $arg2 | $prog $arg1b - 2>&1`;     if (! $?) { die "FAIL: $prog should fail on multiquery $tag1, stdin $tag2.\n"; }


################################################################
# bathstat
################################################################

$tag  = "bathstat";       $prog = "$builddir/src/$tag";    
$tag1 = "<hmmfile>";      $arg1 = "$tmppfx.bhmm";    
if ($verbose) { print "$tag...\n"; }

`$prog $arg1         > $tmppfx.out1`;                     if ($?) { die "FAIL: $tag $tag1\n"; }
`cat $arg1 | $prog - > $tmppfx.out2`;                     if ($?) { die "FAIL: $tag -\n"; }

`diff -b $tmppfx.out1 $tmppfx.out2 2>&1 > /dev/null`;  if ($?) { die "FAIL: $tag results differ if $tag1 comes through stdin\n"; }

unlink <$tmppfx.out*>;
unlink <$tmppfx.bhmm*>;
unlink "$tmppfx.sto";
unlink "$tmppfx.db";
unlink "$tmppfx.key";

print "ok\n";
exit 0;
