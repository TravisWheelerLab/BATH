# Tutorial

FraHMMER was built on top of of the existing HMMER3 code-base, and users who are familiar with HMMER will find that FraHMMER uses many of the same conventions. This tutorial will focus on getting you familiar with the five FraHMMER tools listed below and, to avoid redundency, will link to the HMMER user guide where applicable. There are three sections in this tutorial: (1) Input files, (2) Running searches, and (3) Output files.  The included practices will provide you with specific instructions on how to use the FraHMMER tools.  All the necessary files to complete the practices are located in the directory FraHMMER/tutorial/. **You should cd into this directory before running the practice commands**. If you have not run 'make install' you will need to add the path to the FraHMMER/src/ directory to the executables.

**Tools**
---

**frahmmbuild**   - build and save FraHMMER formated pHMM file from an input multiple sequence alignment (MSA) file
```
Usage: frahmmbuild [-options] <hmmfile_out> <msafile_in>
```
**frahmmstat**   - show summary statistics for a FraHMMER formated pHMM file
```
Usage: frahmmstat [-options] <hmmfile_in>
```
**frahmmconvert** - convert HMMER formated pHMM files to FraHMMER formated pHMM files
```
Usage: frahmmconvert [-options] <hmmfile_out> <hmmfile_in>
```
**frahmmfetch**   - copy selected pHMMs from a HMMER or FraHMMER formatted file (converting if necessary) to a new FraHMMER formated pHMM file
```
Usage: frahmmfetch [options] <hmmfile_in> <key>         (retrieves HMM named <key>)
Usage: frahmmfetch [options] -f <hmmfile_in> <keyfile>  (retrieves all HMMs in <keyfile>)
Usage: frahmmfetch [options] --index <hmmfile_in>       (indexes <hmmfile>)
```
**frahmmer**      - search one or more protein pHMMs against a DNA sequence database
```
Usage: frahmmer [options] <protein-queryfile> <DNA-targetfile>
```


## Section 1 - Input files 

Before you begin using FraHMMER, it will be helpful to become familiar with the file types that are required as inputs. Each frahmmer search needs a protein query file and a DNA target file. The target file must include one or more DNA sequences in a recognizable unaligned sequence or multiple sequence alignment (MSA) format. Common unaligned sequence formats include fasta, embl, and genbank. Common MSA formats include stockholm, a2m, afa, psiblast, clustal, and phylip. 

FraHMMER includes the [Easel](https://github.com/EddyRivasLab/easel) software suite developed by the Eddy/Rivas Lab.  The Easel miniapps are a set of tools designed to perform a number of operations on MSA and unaligned sequence files.  To familiarize yourself with those tools see the [HMMER user guide](http://eddylab.org/software/hmmer/Userguide.pdf) (pages 145-204). 

The query file contains the proteins you wish to search for in the target DNA. The preferred format for query files is a FraHMMER formated pHMM file (although you may also use an MSA or unaligned sequence file - see practice # and #). The rest of this section will focus on practices to get you acquainted with the FraHMMER tools which are used to create and manipulate these pHMM files.

<details><summary>Practice 1: building a pHMM from an MSA using frahmmbuild</summary>
<p>

The sensitivity of FraHMMER is powered, in large part, by the use of pHMMs. The pHMM files used by FraHMMER and almost identical to the ones used by HMMER, but contains additional information needed to perform accurate frameshift-aware translations and provide reliable e-values. If you would like more information on the pHMM files see the [HMMER user guide](http://eddylab.org/software/hmmer/Userguide.pdf) (page 208). FraHMMER formated pHMMs can be created from MSA files using the tool frahmmbuild. 

The file MET.msa contains two stockholm formatted protein MSAs (note that stockholm is the only MSA format that allows multiple MSAs in a single file). To build pHMMs from those MSAs and save them to the file JB.hmm. Run the following command:
```bash
   % frahmmbuild MET.hmm MET.msa
```
Now compare the summary output that is printed to your stdout to the text below (the exact CPU and elapsed time will vary):

```bash
# input alignment file:             MET.msa
# output HMM file:                  MET.hmm
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# idx    name                  nseq  alen  mlen fs_prob codon_tbl eff_nseq re/pos description
# ------ -------------------- ----- ----- ----- ------- --------- -------- ------ -----------
  1      metC                    11   487   409 0.01000         1     0.60  0.588 Cystathionine beta-lyase
  2      metH                     8  1214  1204 0.01000         1     0.57  0.589 Methionine synthase
 

# CPU time: 8.08u 0.00s 00:00:08.08 Elapsed: 00:00:06.04
```
The following is a brief description of each of the above fields. 

```
idx            Number, in order of the MSA file.

name           Name of the pHMM.

nseq           Number of sequences in the alignment this pHMM was built from.

alen           Length of alignment - number of columns in the MSA.

mlen           Length of the pHMM - number of match states.
   
fs_prob        The probability assigned to a nucleotide insertion that results in a frameshift

codon_tbl      The NCBI codon translation table ID assumed for the target DNA

eff_nseq       Effective sequence number. This was the “effective” number of independent sequences that frahmmbuild’s default “entropy weighting” step decided on, given the phylogenetic similarity of the nseq sequences in the input alignment. 

re/pos         Mean positional relative entropy, in bits. This can be ignored by most users. 
   
description    Description of the protein family - may be blank.
```
</p>
</details>

<details><summary>Practice 2: building a pHMM from an MSA using frahmmbuild with an alternate codon translation table</summary>
<p>

One of the fields that distinguishes a FraHMMER formatted pHMM file from a HMMER formated pHMM file is an [NCBI codon translation table ID](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). The correct codon table depends on the origins of the target DNA you intend to search the pHMMs against. When you run a frahmer search selecting the correct codon table will the highest quality alignments. Ensuring that the pHMMs were built with that same the codon table will produce the most accurate e-values for those alignments. By default, frahmmbuild will use the standard code employed by eukaryotic nuclear DNA. To use an alternate codon translation table include the option --ct followed by a table ID from the list below:

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

In a later practice, you will search the pHMMs in MET.msa against a target sequence from the genome of an endosymbiotic bacteria that uses codon table 4. Run the following command to rebuild the pHMMs using the correct codon table for that target:
```bash
   % frahmmbuild --ct 4 MET-ct4.hmm MET.msa
```
The summary output should be nearly identical to that in Pracitce 1, except for the output file name and the codon table field which should now say 4 for both pHMMs. 

</p>
</details>

<details><summary>Practice 3: summerizing pHMM files with frahmmstat</summary>
<p>

Since a pHMM file may contain any number of individual models, it is useful to be able to quickly summarize the contents. The tool frahmmstat is designed to provide such a summary for FraHMMER formated pHMM files.  In this practice, you will compare the summaries of the two pHMM files made in practices 1 and 2.  First, run the following command to summarize the file built with the standard codon table. 
```bash
   % frahmmstat MET.hmm
```
This command should produce the following output to stdout:

```bash
#
# idx    name                 accession        nseq eff_nseq   mlen fs prob codon tbl re/pos
# ------ -------------------- ------------ -------- -------- ------ ------- --------- ------
  1      metC                 -                  11     0.60    409 0.01000         1   0.52
  2      metH                 -                   8     0.57   1204 0.01000         1   0.52
    
```

The fields are mainly the same as those produced by frahmbuild, and detailed in practice 1, with the exception on the accession field which may contian an alphanumeric idetifier for the protein family. 

</p>
</details>
   






