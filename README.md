## HMMER - biological sequence analysis using profile HMMs

[HMMER](http://hmmer.org) searches biological sequence databases for
homologous sequences, using either single sequences or multiple
sequence alignments as queries. HMMER implements a technology called
"profile hidden Markov models" (profile HMMs). HMMER is used by many
protein family domain databases and large-scale annotation pipelines,
including [Pfam](http://pfam.xfam.org) and other members of the
>>>>>>> 5e98c5c2cd5aadf51bffe6644e0b81c7028d512a
[InterPro Consortium](http://www.ebi.ac.uk/interpro/).

To obtain HMMER releases, please visit [hmmer.org](http://hmmer.org).

To participate in HMMER development, visit us at
[github](https://github.com/EddyRivasLab/hmmer).  HMMER development
depends on the Easel library, also at
[github](https://github.com/EddyRivasLab/easel).


### to download and build the current source code release:

```
   % wget http://eddylab.org/software/hmmer/hmmer.tar.gz
   % tar zxf hmmer.tar.gz
   % cd hmmer-3.2.1
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

<<<<<<< HEAD
```bash
   % git clone https://github.com/TravisWheelerLab/hmmer
   % cd hmmer
   % git clone https://github.com/TravisWheelerLab/easel
   % git checkout develop
   % (cd easel; git checkout develop)
=======
```
   % git clone https://github.com/EddyRivasLab/hmmer
   % cd hmmer
   % git clone https://github.com/EddyRivasLab/easel
>>>>>>> 5e98c5c2cd5aadf51bffe6644e0b81c7028d512a
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

<<<<<<< HEAD
To contribute to our fork of HMMER3 development, you want to be on the **develop**
branch. **If you want to contribute to the main HMMER codebase, go to the
[EddyRivasLab](https://github.com/EddyRivasLab) repo**
=======
To build the most recent official release, leave both HMMER and Easel
on their default **master** branch.  To contribute to HMMER3
development, you want to be on the **develop** branches.
>>>>>>> 5e98c5c2cd5aadf51bffe6644e0b81c7028d512a


### to report a problem:

Visit the main HMMER repo
[issues tracking page at github](https://github.com/EddyRivasLab/hmmer/issues).

