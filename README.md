## FATHMM - Frameshift Aware Traslated Hidden Markov Models for the Annotation of Protien Coding DNA

[FATHMM] searches protein sequence databases for
homologous sequences, using either single DNA sequences or DNA 
multiple sequence alignments as queries. FATHMM implements 
innovations in profile hidden Markov model (pHMM) architecture 
and dynamic programing algotithms to align protein coding DNA 
directly to amino acid sequences even in the presecense of 
framethifts.

Binaries are available for download 
To participate in HMMER development, visit us at
[github](https://github.com/EddyRivasLab/hmmer).  HMMER development
depends on the Easel library, also at
[github](https://github.com/EddyRivasLab/easel).


### to download and build the current source code release:

```
   % wget http://eddylab.org/software/hmmer/hmmer.tar.gz
   % tar zxf hmmer.tar.gz
   % cd hmmer-3.3.2
   % ./configure --prefix /your/install/path
   % make
   % make check                 # optional: run automated tests
   % make install               # optional: install HMMER programs, man pages
   % (cd easel; make install)   # optional: install Easel tools
``` 

Executable programs will be installed in `/your/install/path/bin`. If
you leave this optional `./configure` argument off, the default prefix
is `/usr/local`.

Files to read in the source directory:

   * INSTALL - brief installation instructions.
   * Userguide.pdf - the HMMER User's Guide.
 
To get started after installation, see the Tutorial section in the
HMMER User's Guide (Userguide.pdf).



### to clone a copy of HMMER3 source from github:

The tarball way, above, is a better way to install HMMER (it includes
a precompiled Userguide.pdf, for example), but you can also clone our
github repo. You need to clone both the HMMER and Easel repositories,
as follows:


```bash
   % git clone https://github.com/TravisWheelerLab/hmmer 
   % cd hmmer
   % git checkout develop
   % git clone https://github.com/TravisWheelerLab/easel
   % (cd easel; git checkout develop)
   % autoconf
```

and to build:

```bash
   % ./configure
   % make
```

Our [git workflow](https://github.com/TravisWheelerLab/hmmer/wiki/Git-workflow)
focuses on two main branches:

 * **develop** is the HMMER3 development branch
 * **h4-develop** is the HMMER4 development branch.

To contribute to our fork of HMMER3 development, you want to be on the **translatedsearch**
branch. **If you want to contribute to the main HMMER codebase, go to the
[EddyRivasLab](https://github.com/EddyRivasLab) repo**


### to report a problem:

Visit the main HMMER repo
[issues tracking page at github](https://github.com/EddyRivasLab/hmmer/issues).

