## FATHMM - Frameshift Aware Traslated Hidden Markov Models for the Annotation of Protien Coding DNA

[FATHMM] searches protein sequence databases for
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
```

### to report a problem:
[issues tracking page at github](https://github.com/TravisWheelerLab/fathmm/issues).

