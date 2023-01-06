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

## Part 1 - Input files 

Before you begin using FraHMMER, it will be helpful to become familar with the file types that are required for each frahmmer search. To cunduct a frahmmer search you will need a query file and a target file. The target file must include one or more DNA sequences in a recognizable single sequence or multiple sequence alignment format. Common single sequence formats include: fasta, embl, and genbank. Common alignment formats include: stockholm, a2m, afa, psiblast, clustal, and phylip. 

The Easel software suite developed by the Eddy/Rivas Lab (https://github.com/EddyRivasLab/easel)inculdes several miniapps deignsed to easily perform a number opertations on MSA and unaligned sequence files (see the HMMER user guide http://eddylab.org/software/hmmer/Userguide.pdf page 145-204). If you have already installed HMMER (https://github.com/EddyRivasLab/hmmer) you will also have installed the Easel miniapps. To avoid overwiting such a pervious install the miniapps are built but not installed with FraHMMER. If you do not have, nor desire to have, HMMER installed you can still use the miniapps with FraHMMER by including the full path <~/your_install_directory/FraHMMER/easel/miniapps/> to each command.

The query file must contian the protiens you wish to search the against the target DNA. The prefered format for query files is a FraHMMER formated pHMM file (altough you may also use a multiple sequence alignment (MSA), or an unalgined sequence file - see practice # and #). Since a pHMM file may coanin any number of individual models it is usefull to be able to quickly sumerize the contents.  The tool frahmmstat is designed to provide such a summary for FraHMMER formated pHMM files.

<details><summary>Practice 1 : summerizing a pHMM file with frahmmstat</summary>
<p>


The file GRK.hmm contains three FraHMMER foramted pHMMs. By running the following comand we will get a set of facts about each of these pHMMs:

```bash
   % frahmmstat GRK.hmm
```
This comand should produce the following output to stdout:

```bash
  #
  # idx    name                 accession        nseq eff_nseq      M relent   info p relE compKL
  # ------ -------------------- ------------ -------- -------- ------ ------ ------ ------ ------
    1      Glucosamine_iso      PF01182.15         30     1.18    193   0.59   0.62   0.54   0.02
    2      Ribosomal_S19e       PF01090.14         21     0.73    139   0.59   0.59   0.53   0.02
    3      K_oxygenase          PF13434.1          14     0.70    337   0.59   0.57   0.52   0.01
```

Some of the fields above will be more meaningful to you than others. A brief description of each field is provided below.

```
idx            Number, in order in the database.

name           Name of the profile.

accession      Accession (if present; else ’-’).

nseq           Number of sequences in the alignment this profile was built from.

eff_nseq       Effective sequence number. This was the “effective” number of independent sequences that hmmbuild’s default “entropy weighting” step decided on, given the phylogenetic similarity of the nseq sequences in the input alignment. 

M              Length of the profile in consensus residues (match states).

relent         Mean relative entropy of the match state emission probabilities, relative to default null background frequencies, in bits. This is the average bit score per aligned consensus residue. This quantity is the target of frahmmbuild’s entropy weighting procedure for determining eff_nseq.

info           Mean information content per match state emission probability vector, in bits. Probably not useful to you. Information content is just a slightly different calculation from relent.

p relE         Mean positional relative entropy, in bits. Also probably not useful to you. This is an average relative entropy per position that takes into account the transition (insertion/deletion) probabilities. It should be a more accurate estimation of the average bit score contributed per aligned model consensus position.

compKL         Kullback-Leibler (KL) divergence from the average composition of the profile’s consensus match states to the default background frequency distribution, in bits. The higher this number, the more biased the residue composition of the profile is. Highly biased profiles may produce more false positives in searches, and can also slow the acceleration pipeline, by causing too many nonhomologous sequences to pass the filters.

```

You will use frahmmstat several times in Step 2 as you get used to creating and manipulating pHMM files. 

</p>
</details>

## Part 2 - Preparing pHMM files

The sensativity of FraHMMER is powered, in large part, by the use of pHMMs. The pHMM files used by FraHMMER and almost identical to ones used by HMMER, but contains additional infomration needed to perform accurate translation and provide reliable e-values. Three of FraHMMERs five tools (frahmmbuild, frahmmconvert, and frahmmfetch) are used mainly to create or manipulate FraHMMER formated pHMM files. Practices 1 thru ? will cover the use of these tools.

<details><summary>Practice 2 : building a pHMM from and MSA using frahmmbuild</summary>
<p>

The file JB.stk contains two stokholm formated protein MSAs (note that stokholm is the only format which allows multiple MSAs in a single file). In this pracitce you will use the frahmmbuild command to build pHMMs from those MSAs and dave them to the file JB.hmm. Run the following comand... 

```bash
   % frahmmbuild JB.hmm JB.stk
```
...and compare the sumary output that is printed to your stdout to the text bellow (the exact CPU and Elapsed time will vary):

```bash
  # input alignment file:             JB.stk
  # output HMM file:                  JB.hmm
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # idx    name                  nseq  alen  mlen eff_nseq re/pos description
# ------ -------------------- ----- ----- ----- -------- ------ -----------
  1      Jag_N                   30    60    52     5.01  1.066 Jag N-terminus
  2      BMFP                    30   102    79     1.48  0.717 Membrane fusogenic activity

# CPU time: 0.67u 0.00s 00:00:00.67 Elapsed: 00:00:00.41
```

Some of the fields above will be more meaningful to you than others.  A brief description of each field is provided below.

```
idx            Number, in order in the database.

name           Name of the profile.

nseq           Number of sequences in the alignment this profile was built from.

alen           Length of alignment - number of columns in the MSA.

mlen           Length of the profile in consensus residues (match states).

eff_nseq       Effective sequence number. This was the “effective” number of independent sequences that hmmbuild’s default “entropy weighting” step decided on, given the phylogenetic similarity of the nseq sequences in the input alignment. 

re/pos         Mean positional relative entropy, in bits. 

description    Despcriptoin on protien family - may be blank.
```

To check that the pHMMs were built and writen correctly, run frahmmstat on JB.hmm and compare your output to the text below:

```bash
  % frahmmstat JB.hmm
#
# idx    name                 accession        nseq eff_nseq      M relent   info p relE compKL
# ------ -------------------- ------------ -------- -------- ------ ------ ------ ------ ------
  1      Jag_N                PF14804.1          30     5.01     52   1.07   1.16   1.04   0.07
  2      BMFP                 PF04380.8          30     1.48     79   0.72   0.76   0.67   0.08
```
</p>
</details>

<details><summary>Practice 3 : building a pHMM from and MSA using frahmmbuild with a non-standard codon translation table</summary>
<p>

One of the fields that distingishes a FraHMMER fromated pHMM file from a HMMER formated pHMM file is a NCBI codon tranlsation table ID (for more information see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). The correct codon table depends on the origins of the target DNA you intend to search you pHMMs against. Matching the codon table of yout target sequence to the the one used to build your query will have an impact on the accuracy of the reported e-values when running a frahmmer search. By default frahmmbuild will use the standard code used by eukoriotic nuclear DNA.  To use an altenate codon translation table include the option -c followed by a table ID from the list bellow:

```bash
id  description
--- -----------------------------------
  1 Standard
  2 Vertebrate mitochondrial
  3 Yeast mitochondrial
  4 Mold, protozoan, coelenterate mitochondrial; Mycoplasma/Spiroplasma
  5 Invertebrate mitochondrial
  6 Ciliate, dasycladacean, Hexamita nuclear
  9 Echinoderm and flatworm mitochondrial
 10 Euplotid nuclear
 11 Bacterial, archaeal; and plant plastid
 12 Alternative yeast
 13 Ascidian mitochondrial
 14 Alternative flatworm mitochondrial
 16 Chlorophycean mitochondrial
 21 Trematode mitochondrial
 22 Scenedesmus obliquus mitochondrial
 23 Thraustochytrium mitochondrial
 24 Pterobranchia mitochondrial
 25 Candidate Division SR1 and Gracilibacteria
```

</p>
</details>

<details><summary>Practice 4 : converting a HMMER formated pHMM file to FraHMMER format using frahmmconvert</summary>
<p>
   

If you have an existing HMMER formated pHMM file, either one you built yourself or one you downloaded from a site such as PFAM, and want to use it to run a frahmmer search you will first need to covert it to the FraHMMER format. The file ccc.hmm contians three pHMMs in HMMER format. The following comand will create the FraHMMER formated file ddd.hmm containing the same three pHMMs:


```bash
   % frahmmconvert ccc.hmm ddd.hmm
```

</p>
</details>

<details><summary>Practice 5 : copying and converting a single pHMM from larger pHMM file using frahmmfetch </summary>
<p>


If you only need to search with a single pHMM but it is located in file with multiple pHMMs, you can save time by copying the desireed pHMM to a new file using frahmmfetch. If the originial pHMM file is in HMMER format, frahmmfeatch will aslo convert the selected pHMM to FraHMMER format. Use the following command to copy and convert the pHMM XXX from the HMMER formated pHMM file ccc.hmmm and write it it the file XXX.hmm:

```bash
   % frahmmfetch -o xxx.hmm ccc.hmm XXX
```
</p>
</details>

<details><summary>Practice 6 : copying multiple pHMMs from larger pHMM file using frahmmfetch </summary>
<p>


You can aslo use frahmmfetch to multiple pHMMs. To do so you will need to create a file contian the names of all the pHMMs you wish to copy.  

</p>
</details>


## Part 3 - Running frahmmer searches





