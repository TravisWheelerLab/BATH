#! /bin/sh

# Usage: ./i1-bathbuild.sh <bathbuild binary> <MSA file> <tmpfile prefix for filenames>
#
# Verifies that bathbuild builds HMMs reproducibly (no stochastic run-to-run variation by default).
#
if test ! $# -eq 3; then 
  echo "Usage: $0 <bathbuild binary> <MSA file> <tmpfile prefix>"
  exit 1
fi

bathbuild=$1;  if test ! -x $bathbuild;   then echo "FAIL: $bathbuild not executable";  exit 1; fi
alifile=$2;    if test ! -r $alifile;     then echo "FAIL: $alifile not readable";      exit 1; fi

hmmfile1=$3.1.bhmm
hmmfile2=$3.2.bhmm
outfile1=$3.1.out
outfile2=$3.2.out
diffile1=$3.1
diffile2=$3.2

# By default: no stochastic variation
$bathbuild $hmmfile1 $alifile > $outfile1;  if test $? -ne 0; then echo "FAIL: crash"; exit 1; fi
$bathbuild $hmmfile2 $alifile > $outfile2;  if test $? -ne 0; then echo "FAIL: crash"; exit 1; fi

cat $outfile1 | grep -v "^#" > $diffile1
cat $outfile2 | grep -v "^#" > $diffile2
diff $diffile1 $diffile2 > /dev/null
if test $? -ne 0 
then 
   echo "FAIL: output files differ"
   exit 1
fi

cat $hmmfile1| grep -v "^DATE" > $diffile1
cat $hmmfile2| grep -v "^DATE" > $diffile2
diff $diffile1 $diffile2 > /dev/null
if test $? -ne 0 
then 
   echo "FAIL: HMM files differ"
   exit 1
fi

#With --fs and not
$bathbuild --fs $hmmfile1 $alifile > $outfile1;  if test $? -ne 0; then echo "FAIL: crash"; exit 1; fi

cat $outfile1 | grep -v "^#" > $diffile1
diff $diffile1 $diffile2 > /dev/null
if test $? -eq 0
then
   echo "FAIL: output files identical, despite only one using --fs"
   exit 1
fi

cat $hmmfile1| grep -v "^DATE" > $diffile1
diff $diffile1 $diffile2 > /dev/null
if test $? -eq 0
then
   echo "FAIL: HMM files identical, despite only one using --fs"
   exit 1
fi

# both with --fs
$bathbuild --fs $hmmfile2 $alifile > $outfile2;  if test $? -ne 0; then echo "FAIL: crash"; exit 1; fi

cat $outfile1 | grep -v "^#" > $diffile1
cat $outfile2 | grep -v "^#" > $diffile2
diff $diffile1 $diffile2 > /dev/null
if test $? -ne 0
then
   echo "FAIL: output files differ"
   exit 1
fi

cat $hmmfile1| grep -v "^DATE" > $diffile1
cat $hmmfile2| grep -v "^DATE" > $diffile2
diff $diffile1 $diffile2 > /dev/null
if test $? -ne 0
then
   echo "FAIL: HMM files differ"
   exit 1
fi

# With different seeds: HMM files will differ because of statistical fits
$bathbuild --seed 1 $hmmfile1 $alifile > $outfile1;  if test $? -ne 0; then echo "FAIL: crash"; exit 1; fi
$bathbuild --seed 2 $hmmfile2 $alifile > $outfile2;  if test $? -ne 0; then echo "FAIL: crash"; exit 1; fi

cat $hmmfile1| grep -v "^DATE" > $diffile1
cat $hmmfile2| grep -v "^DATE" > $diffile2
diff $diffile1 $diffile2 > /dev/null
if test $? -eq 0 
then 
   echo "FAIL: HMM files identical, despite different rng seeds"
   exit 1
fi

# With different seeds: HMM files will differ because of statistical fits
$bathbuild --seed 1 --fs $hmmfile1 $alifile > $outfile1;  if test $? -ne 0; then echo "FAIL: crash"; exit 1; fi
$bathbuild --seed 2 --fs $hmmfile2 $alifile > $outfile2;  if test $? -ne 0; then echo "FAIL: crash"; exit 1; fi

cat $hmmfile1| grep -v "^DATE" > $diffile1
cat $hmmfile2| grep -v "^DATE" > $diffile2
diff $diffile1 $diffile2 > /dev/null
if test $? -eq 0
then
   echo "FAIL: HMM files identical, despite different rng seeds"
   exit 1
fi
echo "ok"

rm $hmmfile1 $hmmfile2
rm $outfile1 $outfile2
rm $diffile1 $diffile2
exit 0

