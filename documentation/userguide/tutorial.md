# Tutorial

Once you have built FraHMMER you will be able to use the main search tool - frahmmer - as well as several support tools. This tutorial will focus on a basic introduction to using the FraHMMER tools and files that will allow you to perform frameshift-aware translated homology searches. The following is a list of those tools.

**frahmmstat**   - show summary statistics for a FraHMMER formated pHMM file 

**frahmmbuild**   - build and save FraHMMER formated pHMMs from multiple sequence alignment (MSA) files

**frahmmconvert** - cconvert HMMER formated pHMM files to FraHMMER formated pHMM files

**frahmmfetch**   - copy selected pHMMs from one pHMM file to a new pHMM file (converting to FraHMMER format if necessary) 

**frahmmemit**    - generate sample sequences from a FraHMMER formated pHMM 

**frahmmer**      - search one or more protein pHMMs against a DNA sequence database using frameshift aware translation

If you desire more information on these tools than is provided in this tutoria, see the man pages. All files necessary to complete the practices below are located in the directory FraHMMER/tutorial/. Please cd into that directory before you proceed with the practices bellow.  If you have not added FraHMMER executables to your path (e.g. by running 'make install') you will need to include the full path to the FraHMMER/src/ to the start of any FraHMMER tool commands. 

<details><summary>Step 1 - File types</summary>
<p>
   
Before you begin using FraHMMER, it will be helpful to become familiar with the file types that are required for each frahmmer search. To conduct a frahmmer search you will need a query file and a target file. The target file must include one or more DNA sequences in a recognizable unaligned single sequence or MSA format. Common single sequence formats include fasta, embl, and genbank. Common alignment formats include stockholm, a2m, afa, psiblast, clustal, and phylip. 

The Easel software suite developed by the Eddy/Rivas Lab (https://github.com/EddyRivasLab/easel) includes several miniapps designed to easily perform a number of operations on MSA and unaligned single sequence files (see the HMMER user guide http://eddylab.org/software/hmmer/Userguide.pdf page 145-204). If you have already installed HMMER (https://github.com/EddyRivasLab/hmmer) you will also have installed the Easel miniapps. To avoid overwriting such a previous install, the miniapps are built but not installed with FraHMMER. If you do not have, nor desire to have, HMMER installed you can still use the miniapps with FraHMMER by including the full path to FraHMMER/easel/miniapps/ to each command.

The query file must contain the proteins you wish to search against the target DNA. The preferred format for query files is a FraHMMER formated pHMM file (although you may also use a multiple sequence alignment (MSA), or an unaligned sequence file - see practice #TBD). Since a pHMM file may contain any number of individual models it is useful to be able to quickly summarize the contents.  The tool frahmmstat is designed to provide such a summary for FraHMMER formated pHMM files.  To try using frahmmstat, and learn how to interpret its output, click on Practice 1 below and follow the instructions. 

<details><summary>Practice 1 : summarizing a pHMM file with frahmmstat</summary>
<p>
   
```bash
   Usage: frahmmstat [-options] <hmmfile>
```
   
The file GRK.hmm contains three FraHMMER formated pHMMs. By running the following command we will get a set of facts about each of these pHMMs:
   
```bash
   % frahmmstat GRK.hmm
```
This command should produce the following output to stdout:

```bash
  #
  # idx    name                 accession        nseq eff_nseq   mlen fs prob codon tbl relent   info p relE compKL
  # ------ -------------------- ------------ -------- -------- ------ ------- --------- ------ ------ ------ ------
    1      Glucosamine_iso      PF01182.15         30     1.18    193 0.01000         1   0.59   0.62   0.54   0.02
    2      Ribosomal_S19e       PF01090.14         21     0.73    139 0.01000         1   0.59   0.59   0.53   0.02
    3      K_oxygenase          PF13434.1          14     0.70    337 0.01000         1   0.59   0.57   0.52   0.01

Some of the fields above will be more meaningful to you than others. A brief description of each field is provided below.

```
idx            Number, in order in the database.

name           Name of the profile.

accession      Accession (if present; else ’-’).

nseq           Number of sequences in the alignment this profile was built from.

eff_nseq       Effective sequence number. This was the “effective” number of independent sequences that hmmbuild’s default “entropy weighting” step decided on, given the phylogenetic similarity of the nseq sequences in the input alignment. 

mlen           Length of the profile in consensus residues (match states).
   
fs prob        The probability of a single nucleotide indel - resulting in a frameshift - used to calculate important E-value parameters. This will need to match the frameshift probability used by any frahmmer search with this pHMM as the query.  
   
codon tbl      The NCBI codon translation table ID is used to calculate important E-value parameters. This will need to match the codon table used by any frahmmer search with this pHMM as the query.

relent         Mean relative entropy of the match state emission probabilities, relative to default null background frequencies, in bits. This is the average bit score per aligned consensus residue. This quantity is the target of frahmmbuild’s entropy weighting procedure for determining eff_nseq.

info           Mean information content per match state emission probability vector, in bits. Probably not useful to you. Information content is just a slightly different calculation from relent.

p relE         Mean positional relative entropy, in bits. Also probably not useful to you. This is an average relative entropy per position that takes into account the transition (insertion/deletion) probabilities. It should be a more accurate estimation of the average bit score contributed per aligned model consensus position.

compKL         Kullback-Leibler (KL) divergence from the average composition of the profile’s consensus match states to the default background frequency distribution, in bits. The higher this number, the more biased the residue composition of the profile is. Highly biased profiles may produce more false positives in searches, and can also slow the acceleration pipeline, by causing too many nonhomologous sequences to pass the filters. 

```
</p>
</details>
   
</p>
</details>

<details><summary>Step 2 - Preparing pHMM files</summary>
<p>

The sensitivity of FraHMMER is powered, in large part, by the use of pHMMs. The pHMM files used by FraHMMER and almost identical to the ones used by HMMER, but they contain additional information needed to perform accurate translations and provide reliable e-values. Three of FraHMMERs five tools (frahmmbuild, frahmmconvert, and frahmmfetch) are used mainly to create or manipulate FraHMMER formated pHMM files. Practices 2 thru #TBD will cover the use of these tools.

<details><summary>Practice 2 : building pHMMs from MSAs using frahmmbuild</summary>
<p>
   
```bash
   Usage: frahmmbuild [-options] <hmmfile_out> <msafile>
```   

The file met.stk contains two stokholm formated protein MSAs (note that stokholm is the only format which allows multiple MSAs in a single file). In this pracitce you will use the frahmmbuild command to build pHMMs from those MSAs and save them to the file JB.hmm. Run the following comand... 

```bash
   % frahmmbuild met.hmm met.stk
```
...and compare the summary output that is printed to your stdout to the text below (the exact CPU and Elapsed time will vary):
   
```bash
   # input alignment file:             met.stk
   # output HMM file:                  met.hmm
   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   # idx    name                  nseq  alen  mlen fs prob codon tbl eff_nseq re/pos description
   # ------ -------------------- ----- ----- ----- ------- --------- -------- ------ -----------
     1      metC                    11   487   409 0.01000         1     0.60  0.588
     2      metH                     8  1214  1204 0.01000         1     0.57  0.589

   # CPU time: 8.04u 0.01s 00:00:08.04 Elapsed: 00:00:06.01
```

Some of the fields above will be more meaningful to you than others.  A brief description of each field is provided below.

```
idx            Number, in order in the database.

name           Name of the profile.

nseq           Number of sequences in the alignment this profile was built from.

alen           Length of alignment - number of columns in the MSA.

mlen           Length of the profile in consensus residues (match states).
   
fs prob        The probability of a single nucleotide indel - resulting in a frameshift - used to calculate important E-value parameters. This will need to match the frameshift probability used by any frahmmer search with this pHMM as the query.
   
codon tbl      The NCBI codon translation table ID is used to calculate important E-value parameters. This will need to match the codon table used by any frahmmer search with this pHMM as the query.

eff_nseq       Effective sequence number. This was the “effective” number of independent sequences that hmmbuild’s default “entropy weighting” step decided on, given the phylogenetic similarity of the nseq sequences in the input alignment. 

re/pos         Mean positional relative entropy, in bits. 

description    Description of the protein family - may be blank.
```

To check that the pHMMs were built and writen correctly, run frahmmstat on met.hmm and compare your output to the text below:

```bash
  % frahmmstat met.hmm
```
   
```bash
   #
   # idx    name                 accession        nseq eff_nseq   mlen fs prob codon tbl relent   info p relE compKL
   # ------ -------------------- ------------ -------- -------- ------ ------- --------- ------ ------ ------ ------
     1      metC                 -                  11     0.60    409 0.01000         1   0.59   0.60   0.52   0.02
     2      metH                 -                   8     0.57   1204 0.01000         1   0.59   0.60   0.52   0.02
```
</p>
</details>
 
<details><summary>Practice 3 : building pHMMs from MSAs using frahmmbuild with a non-standard codon translation table</summary>
<p>

```bash
   Usage: frahmmbuild [-options] <hmmfile_out> <msafile>
```  
   
One of the fields that distinguishes a FraHMMER formatted pHMM file from an HMMER formated pHMM file is an NCBI codon translation table ID (for more information see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). The correct codon table depends on the origins of the target DNA you intend to search the pHMMs against. Matching the codon table of your target sequence to the one used to build your query will have an impact on the accuracy of the reported e-values when running a frahmmer search. By default, frahmmbuild will use the standard code used by eukaryotic nuclear DNA.  To use an alternate codon translation table include the option --ct followed by a table ID from the list below:
   
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

Since we did not use the --ct flag in Practice 2 the pHMMs in met.hmm were built with codon translation table 1, but these proteins actually come from endosymbiotic bacterial genomes which uses codon translation table 4.  In this practice, we will again build pHMMs from the MSAs in met.stk, but this time with the correct codon table via the --ct flag.  Run the following command and compare the output:
   
```bash
   % frahmmbuild --ct 4 met-C4.hmm met.stk
```
   
```bash
   # input alignment file:             met.stk
   # output HMM file:                  met-C4.hmm
   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   # idx    name                  nseq  alen  mlen fs prob codon tbl eff_nseq re/pos description
   # ------ -------------------- ----- ----- ----- ------- --------- -------- ------ -----------
     1      metC                    11   487   409 0.01000         4     0.60  0.588
     2      metH                     8  1214  1204 0.01000         4     0.57  0.589

   # CPU time: 8.03u 0.01s 00:00:08.03 Elapsed: 00:00:06.01
```
   
You can see that the codon tbl column now says 4. Using the correct codon translation table improves the accuracy of alignments and e-values.  We will see the effect of this difference on a frahmmer search in practice #TBD. 
   
</p>
</details>
   
</p>
</details>

