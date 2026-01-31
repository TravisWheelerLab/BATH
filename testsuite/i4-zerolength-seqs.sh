#! /bin/sh

# Usage:
#   ./i4-zerolength-seqs.sh <builddir> <srcdir> <HMM database> <tmpfile prefix>
# 
# Example:
#
# Verifies that the search programs can take zero length sequences as
# input without crashing.
# 
# The <HMM database> must be press'ed, for hmmscan to work on it.

if test ! $# -eq 4; then 
  echo "Usage: $0 <builddir> <srcdir> <HMM database> <tmpfile prefix>"
  exit 1
fi

builddir=$1;
srcdir=$2;
hmmfile=$3;
tmppfx=$4;

bathsearch=$builddir/src/bathsearch; if test ! -x $bathsearch; then echo "FAIL: $bathsearch not executable"; exit 1; fi
                                     if test ! -r $hmmfile;    then echo "FAIL: $hmmfile not readable";      exit 1; fi

fafile=$tmppfx.fa

cat > $fafile <<EOF
>foo
>bar
ACGTACGTACGTACGTACGTACGT
EOF

$bathsearch $hmmfile $fafile > /dev/null 2>&1; if test $? -ne 0; then echo "FAIL: crash"; exit 1; fi

echo "ok"

rm $fafile
exit 0

