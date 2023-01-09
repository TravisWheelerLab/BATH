# Tutorial

This tutorial will focus on getting you familiar with FraHMMER. It is split into three sections: (1) Input files, (2) Running searches, and (3) Output files.  The included practices will provide you with specific instructions on using the FraHMMER tools below. 

**frahmmstat**   - show summary statistics for a FraHMMER formated pHMM file

**frahmmbuild**   - build and save FraHMMER formated pHMM from an input multiple sequence alignment (MSA) or unaligned sequences file

**frahmmconvert** - convert HMMER formated pHMM files to FraHMMER formated pHMM files

**frahmmfetch**   - copy selected pHMMs from a HMMER or FraHMMER formatted file (converting if necessary) to a new FraHMMER formated pHMM file

**frahmmer**      - search one or more protein pHMMs against a DNA sequence database

All the necessary files to complete the practices are located in the directory FraHMMER/tutorial/. You should cd into this directory before running the practice commands. If you have not run 'make install' you will need to add the path to the FraHMMER/src/ directory to the executables.

## Section 1 - Input files 

Before you begin using FraHMMER, it will be helpful to become familiar with the file types that are required as inputs. Each frahmmer search needs a protein query file and a DNA target file. The target file must include one or more DNA sequences in a recognizable unaligned sequence or multiple sequence alignment (MSA) format. Common unaligned sequence formats include fasta, embl, and genbank. Common MSA formats include stockholm, a2m, afa, psiblast, clustal, and phylip. 

FraHMMER includes the [Easel](https://github.com/EddyRivasLab/easel) software suite developed by the Eddy/Rivas Lab.  The Easel miniapps are a set of tools designed to perform a number of operations on MSA and unaligned sequence files.  To familiarize yourself with those tools see the [HMMER user guide](http://eddylab.org/software/hmmer/Userguide.pdf) (pages 145-204). 

The query file contains the proteins you wish to search for in the target DNA. The preferred format for query files is a FraHMMER formated pHMM file (although you may also use an MSA or unaligned sequence file - see practice # and #). The rest of this section will focus on practices to get you acquainted with the FraHMMER tools which are used to create and manipulate these pHMM files.

<details><summary>Practice 1: building a pHMM from an MSA using frahmmbuild</summary>
<p>

The sensitivity of FraHMMER is powered, in large part, by the use of pHMMs. The pHMM files used by FraHMMER and almost identical to the ones used by HMMER, but contains additional information needed to perform accurate frameshift-aware translations and provide reliable e-values.   FraHMMER formated pHMMs can be created from MSA files using the tool frahmmbuild. 

The file JB.stk contains two stockholm formatted protein MSAs (note that stockholm is the only MSA format that allows multiple MSAs in a single file). To build pHMMs from those MSAs and save them to the file JB.hmm. Run the following command:

```bash
   % frahmmbuild JB.hmm JB.stk
```
Now compare the summary output that is printed to your stdout to the text below (the exact CPU and elapsed time will vary):

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
The following is a brief description of each of the above fields. 

```
idx            Number, in order in the database.

name           Name of the profile.

nseq           Number of sequences in the alignment this profile was built from.

alen           Length of alignment - number of columns in the MSA.

mlen           Length of the profile in consensus residues (match states).

eff_nseq       Effective sequence number. This was the “effective” number of independent sequences that frahmmbuild’s default “entropy weighting” step decided on, given the phylogenetic similarity of the nseq sequences in the input alignment. 

re/pos         Mean positional relative entropy, in bits. 

description Description of the protein family - may be blank.
```
</p>
</details>

<details><summary>Practice 2: building a pHMM from an MSA using frahmmbuild with a non-standard codon translation table</summary>
<p>

One of the fields that distinguishes a FraHMMER formatted pHMM file from a HMMER formated pHMM file is an [NCBI codon translation table ID](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). The correct codon table depends on the origins of the target DNA you intend to search the pHMMs against. Matching the codon table of the target sequence to the one used to build your query will have an impact on the accuracy of the reported e-values when running a frahmmer search. By default, frahmmbuild will use the standard code employed by eukaryotic nuclear DNA.  To use an alternate codon translation table include the option -c followed by a table ID from the list below:

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

The proteins that we used to build pHMMs in Practice 1 come from the genomes of endosymbiotic bacteria.  In a later practice, you will search them against the genome of another endosymbiotic bacteria which uses codon translation table 4. To rebuild the pHMMs using the correct run the following command:

```bash
   % frahmmbuild -c 4 JB4.hmm JB.stk
```
Now compare the summary output that is printed to your stdout to the text below (the exact CPU and elapsed time will vary):

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
</p>
</details>

<details><summary>Practice 3: summerizing pHMM files with frahmmstat</summary>
<p>

Since a pHMM file may contain any number of individual models it is useful to be able to quickly summarize the contents.  The tool frahmmstat is designed to provide such a summary for FraHMMER formated pHMM files.  In this practice, you will compare the summaries of the two pHMM files made in practices 1 and 2.  First, run the following command to summarize the file built with the standard codon table. 
```bash
   % frahmmstat JB1.hmm
```
This command should produce the following output to stdout:

```bash
  #
  # idx    name                 accession        nseq eff_nseq      M relent   info p relE compKL
  # ------ -------------------- ------------ -------- -------- ------ ------ ------ ------ ------
    
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

</p>
</details>
   






