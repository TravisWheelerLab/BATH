# Tutorial

FraHMMER was built on top of of the existing HMMER3 code-base, and users who are familiar with HMMER will find that FraHMMER uses many of the same conventions. This tutorial will focus on getting you familiar with the five FraHMMER tools listed below and, to avoid redundency, will link to the HMMER user guide where applicable. There are three sections in this tutorial: (1) Input files, (2) Running frahmmer searches, and (3) Output files.  The included practices will provide you with specific instructions on how to use the FraHMMER tools.  All the necessary files to complete the practices are located in the directory FraHMMER/tutorial/. **You should cd into this directory before running the practice commands**. If you have not run 'make install' you will need to add the path to the FraHMMER/src/ directory to the executables.

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

<details><summary>Practice 4: converting a HMMER formated pHMM file to FraHMMER format using frahmmconvert</summary>
<p>

If you have an existing HMMER formatted pHMM file and want to use it to run a frahmmer search, you will first need to convert it to the FraHMMER format using frahmmconvert. The file GRK-hmmer.hmm contains three pHMMs in HMMER3 format. The following command will create the FraHMMER formatted file GRK-frahmmer.hmm containing the same three pHMMs:

```bash
   % frahmmconvert GRK-hmmer.hmm GRK-frahmmer.hmm
```
Your summary output should match that shown below.
```
# input HMM file:                   GRK-hmmer.hmm
# output HMM file:                 GRK-frahmmer.hmm
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# idx    name                  nseq  mlen fs_prob codon_tbl eff_nseq re/pos description
# ------ -------------------- ----- ----- ------- --------- -------- ------ -----------
  1      Glucosamine_iso         30   193 0.01000         1     1.18  0.590 Glucosamine-6-phosphate isomerases/6-phosphogluconolactonase
  2      Ribosomal_S19e          21   139 0.01000         1     0.73  0.591 Ribosomal protein S19e
  3      K_oxygenase             14   337 0.01000         1     0.70  0.589 L-lysine 6-monooxygenase (NADPH-requiring)
# CPU time: 2.75u 0.00s 00:00:02.75 Elapsed: 00:00:02.76
```
You can also use frahmmconvert to change the codon table of an existing FraHMMER pHMM file using the --ct flag.  This will be faster than rebuilding from the original MSA. 

</p>
</details>

<details><summary>Practice 5: indexing a pHMM file and copying a single pHMM using frahmmfetch </summary>
<p>

If you only need to search with a single pHMM but it is located in a file with multiple pHMMs, you can save time by copying the desired pHMM to a new file using frahmmfetch. If the original file contains a large number of pHMMs, you may want to create an index file to speed up the fetch process.  The following command will index the create the index file GRK-frahmmer.ssi for the FraHMMER pHMM file we created in Practice 4. 
```bash
   % frahmmfetch --index GRK-frahmer.hmm 
```
The summary output should read as follows:
```
Working...    done.
Indexed 3 HMMs (3 names and 3 accessions).
SSI index written to file GRK-frahmmer.hmm.ssi
```
Whether or not you choose to create an index you will need the name of the pHMM you wish to copy to use as a key. The command below will copy the pHMM Ribosomal_S19e from the GRK-frahmmer.hmm.  The -o flag will direct the copied pHMM to the specified output file (RIB.hmm in this case). Otherwise, it will be printed to standard out. 
```bash
   % frahmmfetch -o RIB.hmm GRK-frahmer.hmm Ribosomal_S19e
```
The summary output should simply read as:
```
Retrieved HMM Ribosomal_S19e.
```
</p>
</details>

<details><summary>Practice 6: converting and copying pHMMs using frahmmfetch </summary>
<p>

You can also use frahmmfetch to copy multiple pHMMs. To do so you will need to create a key file that contains the names of all the pHMMs you wish to copy with one name per line and use the -f flag.  If the original pHMM file is in HMMER format frahmmfetch will automatically convert it to FraHMMER format. The following command will copy both of the pHMMs listed in the key file GK.txt (Glucosamine_iso and K_oxygenase) from a HMMER formated pHMM file, convert them to FraHMMER format, and print them to the output file GK.hmm
```bash
   % frahmmfetch -f -o GK.hmm GRK-hmmer.hmm GK.txt
```
The summary output should simply read as:
```
Retrieved 2 HMMs.
```
As with frahmmconvert, you can also use the --ct flag to change the codon table
</p>
</details>

## Section 2 - Running frahmmer searches

This section of the tutorial will focus on the tool frahmmer, which performs frameshift aware translated homology search using pHMMs and dynamic programming.  This tool allows the user to perform translated annotate of protein-coding DNA even when mutations or sequencing errors have introduced frameshifts.   Each of the practices in this section will involve running a frahmmer search with a different set of input formats,  options, and outputs. 

<details><summary>Practice 7: running a simple frahmmer search</summary>
<p>

Every frahmmer search requires two inputs - the query and the target.  In this practice, you will use the single pHMM you copied to its own file in Practice 5 (RIB.hmm) as the query.  For the target, you will use a single DNA sequence in the file seq1.fa, and the -o flag will be used to direct the hit data and alignment to the file RIB.out.  Run the following command:
```bash
   % frahmmer -o RIB.out RIB.hmm seq1.fa
```
The file RIB.out should contain a single hit between the Ribosomal_S19e protein and the DNA sequence. In Section 3 we will examine this output in detail.

</p>
</details>
 
<details><summary>Practice 8: running a frahmmer search on a target with an alternate codon translation table</summary>
<p>

As discussed in Practice 2, some DNA sequences use alternate codon translation tables. For frahmmer searches that use such DNA as the target, the best results are achieved by specifying the correct codon table both during the actual search and when building the pHMMs. In this  Practice you will first attempt to run a frahmmer search with a mismatch between the codon table specified for the target and the codon table used to build the pHMM, resulting in an error message.  

Run the following command to search the pHMMs in MET.hmm, which you built in Practice 1 using the standard codon translation, against the target DNA in the file seq2.fa that you will specify as using codon table 4 with the --ct flag

```bash
   % frahmmer --ct 4 -o MET.out MET.hmm seq2.fa
```
   
This will result in the following error message:

```bash
   Error: Requested codon tranlsation tabel ID 4 does not match the codon tranlsation tabel ID of the HMM file MET.hmm. Please run frahmmcovert with option '--ct 4'.
```

To avoid this error we need to use the pHMM file with the correct codon translation table by running the following command:

```bash
   % frahmmer --ct 4 -o MET.out MET-ct4.hmm seq2.fa
```
The file MET.out should contain a single hit between each of the pHMMS in MET-ct4.hmm and the DNA sequence. In Section 3 we will examine this output in detail.
   
</p>
</details>

