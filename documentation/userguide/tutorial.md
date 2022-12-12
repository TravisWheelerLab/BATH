# Tutorial

Once you have built FraHMMER you will be able to use the main search tool - frahmmer - as well as the several support tools which help create and format profile hidden Markov model (pHMM) files. This tutorial will focus on getting you familiar with the FraHMMER tools and files that will allow you to perform frameshift-aware translated homology searches. The following is a list of those tools.

**frahmmstat***   - show summary statistics for a FraHMMER formated pHMM file
**frahmmbuild**   - build and save FraHMMER formated pHMM from an input multiple sequence alignment (MSA) or unaligned sequences file
**frahmmconvert** - convert HMMER formated pHMM files to FraHMMER formated pHMM files
**frahmmfetch**   - copy selected pHMMs from a HMMER or FraHMMER formated file (converting if necissary) to a new FraHMMER formated pHMM file
**frahmmer**      - search one or more protein pHMMs against a DNA sequence database

All files necessary to complete the practices below are located in the directory FraHMMER/tutorial/. Run the following command to insure you are in the tutorial directotry before you proceed with the practices bellow.

```bash
   % cd /your_install_path/FraHMMER/tutorial/
```
If you have not added FraHMMER executables to your path (e.g. by running 'make install') you will need to add the full path to the FraHMMER/src/ directory to the start of any FraHMMER commands. 

## Step 1 - File types

Before you begin using FraHMMER, it will be helpful to become familar with the file types that are required for each frahmmer search. To cunduct a frahmmer search you will need a query file and a target file. The target file must include one or more DNA sequences in a recognizable single sequence or multiple sequence alignment format. Common single sequence formats include: fasta, embl, and genbank. Common alignment formats include: stockholm, a2m, afa, psiblast, clustal, and phylip. 

The Easel software suite developed by the Eddy/Rivas Lab (https://github.com/EddyRivasLab/easel)inculdes several miniapps deignsed to easily perform a number opertations on MSA and unaligned sequence files (see the HMMER user guide http://eddylab.org/software/hmmer/Userguide.pdf page 145-204). If you have already installed HMMER (https://github.com/EddyRivasLab/hmmer) you will also have installed the Easel miniapps. To avoid overwiting such a pervious install the miniapps are built but not installed with FraHMMER. If you do not have, nor desire to have, HMMER installed you can still use the miniapps with FraHMMER by including the full path <~/your_install_directory/FraHMMER/easel/miniapps/> to each command.

The query file must contian the protiens you wish to search the against the target DNA. The prefered format for query files is a FraHMMER formated pHMM file (altough you may also use a multiple sequence alignment (MSA), or an unalgined sequence file - see practice # and #). Since a pHMM file may coanin any number of individual models it is usefull to be able to quickly sumerize the contents.  The tool frahmmstat is designed to provide such a summary for FraHMMER formated pHMM files.
 
**Practice 1 : summerizing a pHMM file with frahmmstat**

The file ggg.hmm contains three FraHMMER foramted pHMMs. By running the following comand we will get a set of facts about each of these pHMMs:

```bash
   % frahmmstat GRK.hmm
```
This comand should produce the following output:

```bash
   #
   # idx  name                 accession        nseq eff_nseq      M relent   info p relE compKL
   # ---- -------------------- ------------ -------- -------- ------ ------ ------ ------ ------
     1    Glucosamine_iso      PF01182.15         30     1.18    193   0.59   0.62   0.54   0.02
     2    Ribosomal_S19e       PF01090.14         21     0.73    139   0.59   0.59   0.53   0.02
     3    K_oxygenase          PF13434.1          14     0.70    337   0.59   0.57   0.52   0.01
```

The ouput should look like this:

```bash
   % frahmmstat ggg.hmm
```
You will use frahmmstat several times in Step 2 as you get used to creating and manipulating pHMM files. 

## Step 2 - Preparing pHMM files

The sensativity of FraHMMER is powered, in large part, by the use of pHMMs. The pHMM files used by FraHMMER and almost identical to ones used by HMMER, but contain additional vital statistics needed to provide reliable evalues. Three of FraHMMERs five tools (frahmmbuild, frahmmconvert, and frahmmfetch) are used mainly to create or manipulate FraHMMER formated pHMM files. 

# frahmmbuild

**Practice 1 : building a pHMM from and MSA using frahmmbuild 

The file aaa.sto contains two stokholm formated protein MSAs (note that stokholm is the only format which allows multiple MSAs in a single file). By running the following command the each of those MSAs will be used to build a pHMM and ouput to the file aaa.hmm. 

```bash
   % frahmmbuild aaa.hmm aaa.sto
```

**Practice 2 : building a pHMM from single sequences using frahmmbuild 

If you wish to build a seperate pHMM for each sequence in your MSA or unaligned sequence file you will need to use the flag "--singlemx".  The file bbb.fa continans the consencus seqeunces from the two MSAs in aaa.sto.  To build two pHMMs from these single sequences use the following command:

```bash
   % frahmmbuild --singlemx bbb.hmm bbb.fa
```

# frahmmconvert

**Practice 3 : converting a HMMER formated pHMM file to FraHMMER format using frahmmconvert

If you have an existing HMMER formated pHMM file, either one you built yourself or one you downloaded from a site such as PFAM, and want to use it to run a frahmmer search you will first need to covert it to the FraHMMER format. The file ccc.hmm contians three pHMMs in HMMER format. The following comand will create the FraHMMER formated file ddd.hmm containing the same three pHMMs:


```bash
   % frahmmconvert ccc.hmm ddd.hmm
```

# frahmmfetch

**Practice 4 : copying and converting a single pHMM from larger pHMM file using frahmmfetch 

If you only need to search with a single pHMM but it is located in file with multiple pHMMs, you can save time by copying the desireed pHMM to a new file using frahmmfetch. If the originial pHMM file is in HMMER format, frahmmfeatch will aslo convert the selected pHMM to FraHMMER format. Use the following command to copy and convert the pHMM XXX from the HMMER formated pHMM file ccc.hmmm and write it it the file XXX.hmm:

```bash
   % frahmmfetch -o xxx.hmm ccc.hmm XXX
```

**Practice 5 : copying multiple pHMMs from larger pHMM file using frahmmfetch 

You can aslo use frahmmfetch to multiple pHMMs. To do so you will need to create a file contian the names of all the pHMMs you wish to copy.  



**Practice 4 : copying a multiple pHMMs from larger pHMM file and converting to FraHMMR format using frahmmfetch 


Every frahmmer search requires two input files - a query and a target. The target file must include one or more DNA sequences in a recognizable format. These sequences may be aligned or unaligned. Common unaligned sequence formats include fasta, embl, and genbank. Common alignment formats include stockholm, a2m, afa, psiblast, clustal, and phylip.

The query file must be either a FraHMMER formatted protein pHMM file, a protein multiple sequence alignment (MSA) or an unaligned protein sequence file. If you use an MSA or an unaligned sequence file the pHMMs will need to be built before the search can proceed and it is recommended that you use the --hmmout flag to save these pHMMs to file for future use (see Practice 1). You can also prebuild and save the pHMMs from your MSA or unaligned sequence file using frahmmbuild (see Practice 2). If you already have a HMMER formated pHMM file you can convert it to a FraHMMER formated pHMM file using frahmmconvert (see Practice 3).



