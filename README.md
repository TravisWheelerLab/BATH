## FATHMM - Frameshift Aware Traslated Hidden Markov Models for the Annotation of Protien Coding DNA

FATHMM searches protein sequence databases for
homologous sequences, using either single DNA sequences or DNA 
multiple sequence alignments as queries. FATHMM implements 
innovations in profile hidden Markov model (pHMM) architecture 
and dynamic programing algotithms to align protein coding DNA 
directly to amino acid sequences even in the presecense of 
frameshifts.

### to download FATHMM:
Please select the correct binaries from the relseases section of 
the FATHMM github page.

### to clone a copy of FATHMM source from github:

```bash
   % git clone https://github.com/TravisWheelerLab/fathmm
   % cd fathmm
   % git clone https://github.com/TravisWheelerLab/easel
   % (cd easel; git checkout develop)
   % autoconf
```

and to build:

```bash
   % ./configure
   % make
   % make install               # optional: install FATHMM programs, man pages
```

### to report a problem:
[issues tracking page at github](https://github.com/TravisWheelerLab/fathmm/issues).

### User Guide

Install and Build FATHMM using one of the methods above.

FATHMM requires a protien query and a DNA target.  The target file must include one or more DNA sequences in a recognizable format (see --tformat). The query file may be either a profile model built using fathmmbuild, a sequence alignment, or a single sequence. Sequence based queries can be in a number
of formats (see --qformat), and can typically be autodetected. Note that only Stockholm format supports queries made up of more than one sequence alignment. If you have a profile model built for HMMER you will need to run fathmmconvert to reformat it for FATHMM. If you anticipate using the query more than once it is highly recomended that you create an hmm file either by prebuilding it with fathmmbuild or by building it as you run FATHMM by using the -hmmout flag. 

To run FATHMM: 
fathmm [options] queryfile targetfile

**Options** 

