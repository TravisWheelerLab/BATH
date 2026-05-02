# Tutorial

BATH was built on top of the existing HMMER3 code base. Users who are familiar with HMMER will find that BATH uses many of the same conventions. This tutorial will focus on getting you familiar with the five BATH tools listed below and, to avoid redundancy, will link to the HMMER user guide where applicable. This tutorial is for BATH 2.0 and later releases. 

There are two sections in this tutorial. The first section - Input files - will cover the tools that are used to prepare your data before you begin a search. The second section - Running bathsearch - will focus on using BATH to perform translated homology search, and on interpreting the search results. All the necessary files to complete the practices are located in the directory BATH/tutorial/. **You should cd into this directory before running the practice commands**. Output files from the practice command are already present in the tutorial directory if you want to skip running the commands and jump to looking at the output. 

**Tools**
---

**bathsearch**      - search a DNA sequence database (or genome) for instances of one or more query proteins. The query can consist of a file of pHMMs (produced using bathbuild or bathconvert - see Practices 1 and 5 below) or a file containing sequences or sequence alignments (see Practice 11). The alignments bathsearch procedures are translated (codons to amino acids) and can be frameshift aware or spliced. 
```
Usage: bathsearch [options] <protein-queryfile> <DNA-targetfile>
```

**bathbuild**   - build BATH formatted profile hidden Markov models (pHMMs) from input multiple sequence alignments (MSAs) or unaligned sequences and save to file
```
Usage: bathbuild [-options] <hmmfile_out> <msa_or_seq_file_in>
```
**bathstat**   - show summary statistics for a BATH-formatted pHMM file
```
Usage: bathstat [-options] <hmmfile_in>
```
**bathconvert** - convert HMMER-formatted or older BATH-formatted pHMM files to the current BATH format
```
Usage: bathconvert [-options] <hmmfile_out> <hmmfile_in>
```
**bathfetch**   - copy selected pHMMs from an HMMER or BATH formatted file (converting if necessary) to a new BATH formatted pHMM file
```
Usage: bathfetch [options] <hmmfile_in> <key>         (retrieves HMM named <key>)
Usage: bathfetch [options] -f <hmmfile_in> <keyfile>  (retrieves all HMMs in <keyfile>)
Usage: bathfetch [options] --index <hmmfile_in>       (indexes <hmmfile>)
```


## Section 1 - Input files 

Before you begin using BATH, it will be helpful to become familiar with the file types it uses. Running bathsearch requires two input files: a protein query file and a DNA target file. The target file must contain one or more DNA sequences in a recognizable unaligned sequence or multiple sequence alignment (MSA) format. Accepted unaligned sequence formats include fasta, embl, and genbank. Accepted MSA formats include stockholm, a2m, afa, psiblast, clustal, and phylip. 

BATH's installation includes a branch of the [Easel](https://github.com/EddyRivasLab/easel) software suite developed by the Eddy/Rivas Lab.  The Easel miniapps are a set of tools designed to perform several operations on MSA and unaligned sequence files.  To familiarize yourself with those tools, see the [HMMER user guide](http://eddylab.org/software/hmmer/Userguide.pdf) (pages 145-204). 

A bathsearch query file contains the proteins you wish to search for in the target DNA. The preferred format for query files is a BATH-formatted pHMM file (although you may also use an MSA or unaligned sequence file - see practice 11). The current version of BATH performs translated (non-framshifted, non-spliced) alignment by default.  To enable the frameshift-aware algorithms, you can use the --fs flag.  To enable the splicing algorithm, you can use the --splice flag.  At this time, BATH does not support using --fs and --splice for the same search. The rest of this section will focus on practices to get you acquainted with the BATH tools used to create and manipulate BATH-formatted pHMM files.

<details><summary>Practice 1: building a pHMM from an MSA using bathbuild</summary>
<p>

The sensitivity of BATH is powered, in large part, by the use of pHMMs. The pHMM files used by BATH are nearly identical to the ones used by HMMER, but contain additional information needed to perform codon translations and provide e-values for frameshift-aware alignments. If you would like more information on pHMM files see the [HMMER user guide](http://eddylab.org/software/hmmer/Userguide.pdf) (page 208). 
   
BATH formated pHMMs can be created from MSA files or unaligned sequence files using the tool bathbuild. The file MET.msa contains two stockholm formatted protein MSAs (note that stockholm is the only MSA format that allows multiple MSAs in a single file). You can build pHMMs from those MSAs and save them to the file MET.bhmm by running the following command: 
   
```bash
% bathbuild MET.bhmm MET.msa
```
(note the file suffix '.bhmm' - this can help distinguish BATH-formatted pHMM files from HMMER-formatted ones, which often have the suffix '.hmm')

The summary output that is printed to your stdout should resemble the text below (the exact CPU and elapsed time will vary):

```bash
# bathbuild :: profile HMM construction from MSAs or unaligned sequences
# BATH 2.0 (May 2026); https://github.com/TravisWheelerLab/BATH
# Freely distributed under the BSD open source license.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# input file:                       MET.msa
# output HMM file:                  MET.bhmm
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# idx    name                  nseq   len  mlen ctbl eff_nseq re/pos description
# ------ -------------------- ----- ----- ----- ---- -------- ------ -----------
  1      metC                    11   483   409    1     0.60  0.591 Cystathionine beta-lyase
  2      metG                    24   494   458    1     0.62  0.589 Methionine--tRNA ligase

# CPU time: 4.67u 0.00s 00:00:04.67 Elapsed: 00:00:02.50
```
   
The following is a brief description of each of the above fields. 
   
```
idx            Number, in order of the MSA file.

name           Name of the pHMM.

nseq           Number of sequences in the alignment this pHMM was built from.

alen           Length of alignment - number of columns in the MSA.

mlen           Length of the pHMM - number of match states.

codon_tbl      The NCBI codon translation table ID assumed for the target DNA

eff_nseq       Effective sequence number. This is the “effective” number of independent sequences that bathbuild’s default “entropy weighting” step decided on, given the phylogenetic similarity of the nseq sequences in the input alignment. The higher the number the more diversity there is among the sequences in the MSA. 

re/pos         Mean positional relative entropy, in bits. This can be ignored by most users. 
   
description    Description of the protein family - may be blank.
```
</p>
</details>

<details><summary>Practice 2: building a pHMM from an MSA using bathbuild with an alternate codon translation table</summary>
<p>

One of the fields that distinguishes a BATH formatted pHMM file from an HMMER formated pHMM file is an [NCBI codon translation table ID](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). The correct codon table depends on the origins of the target DNA you intend to search the pHMMs against. When you run bathsearch, selecting the correct codon table will produce the highest quality alignments. Ensuring that the pHMMs were built with that same codon table will produce the most accurate e-values for those alignments. 
   
By default, bathbuild will use the standard code employed by eukaryotic nuclear DNA. To use an alternate codon translation table include the option --ct followed by a table ID from the list below:

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

In practice 10 you will search pHMMs built from MET.msa against a target sequence from the genome of an endosymbiotic bacterium that uses codon table 4. Running the following command will build the pHMMs using the correct codon table for that target:
   
```bash
% bathbuild --ct 4 MET-ct4.bhmm MET.msa
```
   
The summary output that is printed to your stdout should resemble the text below (the exact CPU and elapsed time will vary):

```bash
# bathbuild :: profile HMM construction from MSAs or unaligned sequences
# BATH 2.0 (May 2026); https://github.com/TravisWheelerLab/BATH
# Freely distributed under the BSD open source license.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# input file:                       MET.msa
# output HMM file:                  MET-ct4.bhmm
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# idx    name                  nseq   len  mlen ctbl eff_nseq re/pos description
# ------ -------------------- ----- ----- ----- ---- -------- ------ -----------
  1      metC                    11   483   409    4     0.60  0.591 Cystathionine beta-lyase
  2      metG                    24   494   458    4     0.62  0.589 Methionine--tRNA ligase

# CPU time: 4.66u 0.01s 00:00:04.67 Elapsed: 00:00:02.48
```

</p>
</details>

<details><summary>Practice 3: building a pHMM from an unaligned sequences file</summary>
<p>

If your queries are single unaligned sequences rather than MSAs you can still build HMMs using bathbuild. The file three_seqs.fa contains three unaligned protein sequences. To build three HMMs, one for each sequence, use the following command. 
   
```bash
% bathbuild three_seqs.bhmm three_seqs.fa
```

The summary output that is printed to your stdout should resemble the text below (the exact CPU and elapsed time will vary):

```bash
# bathbuild :: profile HMM construction from MSAs or unaligned sequences
# BATH 2.0 (May 2026); https://github.com/TravisWheelerLab/BATH
# Freely distributed under the BSD open source license.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# input file:                       three_seqs.fa
# output HMM file:                  three_seqs.bhmm
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# idx    name                  nseq   len  mlen ctbl eff_nseq re/pos description
# ------ -------------------- ----- ----- ----- ---- -------- ------ -----------
  1      AT1G01010.1              1   429   429    1     1.00  0.591
  2      AT1G01020.1              1   245   245    1     1.00  0.552
  3      AT1G01030.1              1   358   358    1     1.00  0.600

# CPU time: 4.58u 0.01s 00:00:04.59 Elapsed: 00:00:02.72
```

</p>
</details>

<details><summary>Practice 4: summerizing pHMM files with bathstat</summary>
<p>

Since a pHMM file may contain any number of individual models, it is useful to be able to quickly summarize the contents. The tool bathstat is designed to provide such a summary for BATH-formatted pHMM files.  The following command will summarize the pHMM file built in practice 1:
   
```bash
% bathstat MET.bhmm
```
   
This command should produce the following output to stdout:

```bash
# bathstat :: display summary statistics for a profile file
# BATH 2.0 (May 2026); https://github.com/TravisWheelerLab/BATH
# Freely distributed under the BSD open source license.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# idx    name                 accession        nseq eff_nseq   mlen codon_tbl re/pos
# ------ -------------------- ------------ -------- -------- ------ --------- ------
  1      metC                 -                  11     0.60    409         1   0.53
  2      metG                 -                  24     0.62    458         1   0.53
```

The fields are mainly the same as those produced by bathbuild, and detailed in practice 1, except for the accession field, which may contain an alphanumeric identifier for the protein family or be left blank if no accession is listed for the HMM. 

</p>
</details>

<details><summary>Practice 5: converting an HMMER-formatted or old BATH-formatted HMM file to the current BATH format using bathconvert</summary>
<p>

If you have an old BATH HMM file or a HMMER-formatted HMM file, you can use bathconvert to produce a new file with the current BATH format. The file tRNA-proteins.hmm contains 12 HMMs in HMMER3 format. The following command will create the BATH formatted file tRNA-proteins.bhmm containing the same twelve HMMs:

```bash
% bathconvert  tRNA-proteins.bhmm  tRNA-proteins.hmm
```
Your summary output should match that shown below.
   
```
# bathconvert :: convert HMMER or older BATH formatted HMM to current BATH format
# BATH 2.0 (May 2026); https://github.com/TravisWheelerLab/BATH
# Freely distributed under the BSD open source license.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# input HMM file:                   tRNA-proteins.hmm
# output HMM file:                  tRNA-proteins.bhmm
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# idx    name                  nseq  mlen codon_tbl eff_nseq re/pos description
# ------ -------------------- ----- ----- --------- -------- ------ -----------
  1      ATE_N                   30    78         1     1.11  0.726 Arginine-tRNA-protein transferase, N terminus
  2      GlutR_N                 12   152         1     0.87  0.590 Glutamyl-tRNAGlu reductase, N-terminal domain
  3      PTH2                    10   116         1     0.74  0.589 Peptidyl-tRNA hydrolase PTH2
  4      RtcB                    30   459         1     0.83  0.590 tRNA-splicing ligase RtcB
  5      TGT                     15   238         1     0.80  0.589 Queuine tRNA-ribosyltransferase
  6      Thg1                    30   131         1     0.69  0.589 tRNAHis guanylyltransferase
  7      Trm56                   11   121         1     0.64  0.590 tRNA ribose 2'-O-methyltransferase, aTrm56
  8      tRNA-synt_1_2           30   185         1     0.91  0.590 Leucyl-tRNA synthetase, Domain 2
  9      tRNA-synt_1c_C          14   192         1     0.81  0.591 tRNA synthetases class I (E and Q), anti-codon binding domain
  10     tRNA-synt_2d            19   247         1     0.73  0.592 tRNA synthetases class II core domain (F)
  11     tRNA-Thr_ED             12   136         1     0.63  0.590 Archaea-specific editing domain of threonyl-tRNA synthetase
  12     TruB_C                  11    56         1     1.64  0.994 tRNA Pseudouridine synthase II, C terminal
# CPU time: 7.86u 0.01s 00:00:07.87 Elapsed: 00:00:07.90
```

</p>
</details>

<details><summary>Practice 6: changing the codon table of a BATH-formatted pHMM using bathconvert</summary>
<p>
 
You can also use bathconvert to change the codon table of an existing BATH pHMM file using the --ct flag. This is faster than rebuilding from the original MSA/sequence. The following command will create the file tRNA-proteins-ct11.bhmm containing the same 12 pHMMs as tRNA-proteins.bhmm, but modified to use NCBI codon translation table 11:
   
```bash
% bathconvert --ct 11 tRNA-proteins-ct11.bhmm tRNA-proteins.bhmm
```

This should produce the following output:
 
```
# bathconvert :: convert HMMER or older BATH formatted HMM to current BATH format
# BATH 2.0 (May 2026); https://github.com/TravisWheelerLab/BATH
# Freely distributed under the BSD open source license.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# input HMM file:                   tRNA-proteins.bhmm
# output HMM file:                  tRNA-proteins-ct11.bhmm
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# idx    name                  nseq  mlen codon_tbl eff_nseq re/pos description
# ------ -------------------- ----- ----- --------- -------- ------ -----------
  1      ATE_N                   30    78        11     1.11  0.726 Arginine-tRNA-protein transferase, N terminus
  2      GlutR_N                 12   152        11     0.87  0.590 Glutamyl-tRNAGlu reductase, N-terminal domain
  3      PTH2                    10   116        11     0.74  0.589 Peptidyl-tRNA hydrolase PTH2
  4      RtcB                    30   459        11     0.83  0.590 tRNA-splicing ligase RtcB
  5      TGT                     15   238        11     0.80  0.589 Queuine tRNA-ribosyltransferase
  6      Thg1                    30   131        11     0.69  0.589 tRNAHis guanylyltransferase
  7      Trm56                   11   121        11     0.64  0.590 tRNA ribose 2'-O-methyltransferase, aTrm56
  8      tRNA-synt_1_2           30   185        11     0.91  0.590 Leucyl-tRNA synthetase, Domain 2
  9      tRNA-synt_1c_C          14   192        11     0.81  0.591 tRNA synthetases class I (E and Q), anti-codon binding domain
  10     tRNA-synt_2d            19   247        11     0.73  0.592 tRNA synthetases class II core domain (F)
  11     tRNA-Thr_ED             12   136        11     0.63  0.590 Archaea-specific editing domain of threonyl-tRNA synthetase
  12     TruB_C                  11    56        11     1.64  0.994 tRNA Pseudouridine synthase II, C terminal
# CPU time: 7.85u 0.01s 00:00:07.86 Elapsed: 00:00:07.93
```
 
</p>
</details>

<details><summary>Practice 7: indexing a pHMM file and copying a single pHMM using bathfetch </summary>
<p>

If you only need to search with a single pHMM, but it is located in a file with multiple pHMMs, you can save time by copying the desired pHMM to a new file using bathfetch. If the original file contains a large number of pHMMs, you may want to create an index file to speed up the fetch process.  The following command will create the index file tRNA-proteins.bhmm.ssi for the BATH pHMM file created in Practice 4. 
```bash
% bathfetch --index tRNA-proteins.bhmm 
```
The summary output should read as follows:
   
```
Working...    done.
Indexed 12 HMMs (12 names and 12 accessions).
SSI index written to file tRNA-proteins.bhmm.ssi
```
Whether or not you choose to create an index you will need the name of the pHMM you wish to copy to use as a key. The command below will copy the HMM PTH2 from the tRNA-proteins.bhmm.  The -o flag will direct the copied pHMM to the specified output file (PTH2.bhmm in this case). Otherwise, it will be printed to stdout. 
```bash
   %  bathfetch -o PTH2.bhmm tRNA-proteins.bhmm PTH2
```
The output should read as:
```
Retrieved HMM PTH2.
```
</p>
</details>

<details><summary>Practice 8: copying and converting multiple HMMs using bathfetch </summary>
<p>

You can also use bathfetch to copy multiple pHMMs. To do so you will need to create a key file containing the names of all the pHMMs you wish to copy, with one name per line, and use the -f flag. If the original pHMM file is in HMMER format, bathfetch will automatically convert it to BATH format. The following command will copy all 3 of the HMMs listed in the key file tRNA-synthetases-key.txt from an HMMER-formatted pHMM file, convert them to BATH format, and print them to the output file tRNA-synthetases.bhmm.
   
```bash
% bathfetch -f -o tRNA-synthetases.bhmm tRNA-proteins.hmm tRNA-synthetases-key.txt
```
   
The output should read as:
   
```
Retrieved 3 HMMs.
```
   
As with bathconvert, you can also use the --ct flag with bathfetch to change the codon table.
   
</p>
</details>

## Section 2 - Running bathsearch

This section of the tutorial will focus on the tool bathsearch. This tool allows the user to perform translated annotation of protein-coding DNA.  By using the --fs flag, bathsearch can align DNA to proteins even when mutations or sequencing errors have introduced frameshifts. By using the --splice flag, bathsearch can detect splice sites and join alignments across exon boundaries. Each of the practices in this section will involve running bathsearch with a different set of input formats, options, and outputs. 

<details><summary>Practice 9: running a simple bathsearch and reading the output</summary>
<p>

Bathsearch requires two inputs - the query and the target.  In this practice, you will use the single HMM in the file PTH2.bhmm as the query.  For the target, we will use a single DNA sequence in the file target-PTH2.fa. The -o flag directs the alignment output to the file PTH2.out (otherwise it will be printed to stdout).
   
```bash
% bathsearch -o PTH2.out PTH2.bhmm target-PTH2.fa
```
 
We can now open the file PTH2.out and see that the output is organized into the following sections:
     
   1) File Header - lines begin with '#' and contain basic information about the search parameters
```
# bathsearch :: search protein profile(s) against DNA sequence database
# BATH 2.0 (May 2026); https://github.com/TravisWheelerLab/BATH
# Freely distributed under the BSD open source license.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query HMM file:                                PTH2.bhmm
# target sequence database:                      target-PTH2.fa
# codon translation table:                       1
# output directed to file:                       PTH2.out
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
```
    
   2) Query Header - includes a summary of each query and a hits list sorted by E-value.  For each hit, the query header lists the E-value, bit score, and bias score adjustment (for more information on bias scores, see pages 60-61 of the [HMMER user guide](http://eddylab.org/software/hmmer/Userguide.pdf).  This is followed by the name of the target sequence where the hit was located, the target sequence position for the start and end of the alignment, and a target description (which may be blank).

```
Query:       PTH2  [M=116]
Accession:   PF01981.11
Description: Peptidyl-tRNA hydrolase PTH2
Scores for complete hits:
    E-value  score  bias  Sequence  start    end  Description
    ------- ------ -----  --------  -----  -----  -----------
      3e-35  110.6   0.3  seq1        672    325
    9.4e-28   86.4   0.0  seq1       1486   1731
    3.5e-12   36.2   0.0  seq1       2468   2343
    4.2e-12   36.0   0.3  seq1       1273   1359
```
   
   3) Annotation Lines- for each hit listed in the query header, bathsearch will produce an annotation line containing useful information about the hit. 
   
       As in the query header, the annotations line lists the score, bias, and E-value for each hit. It also lists three types of coordinates for the hit - the alignment start and end coordinates for both the query (hmm-from & hmm-to) and the target (ali-from & ali-to), as well as the envelope coordinates (env-from & env-to). The envelope is the region of the target that bathsearch has identified as containing the homology (the alignment coordinates are always contained within the envelope). It is the envelope coordinates that bound the target subsequence used to calculate the score, bias, and E-value. An explanation of the characters seen after the coordinates ('.','[', & ']') can be found on page 38 of the [HMMER user guide](http://eddylab.org/software/hmmer/Userguide.pdf). The annotation line also lists the full length of the target sequence (sq-len) and the alignment's accuracy score (acc), which is the average expected per-residue accuracy of the alignment. 
       
       Below is the annotation line for the first hit in the file PTH2.out
  
```
    score  bias    Evalue   hmm-from    hmm-to     ali-from    ali-to     env-from    env-to       sq-len    acc
   ------ ----- ---------   --------   -------    --------- ---------    --------- ---------    ---------   ----
 !  110.6   0.3     3e-35          2       116 .]       672       325 ..       675       325 ..      3000   0.92
```
   
   4) Alignment - Below each annotation line, bathsearch prints the alignment for that query-target hit. A typical bathsearch alignment will contain at least the following five rows (in order from top to bottom): (1) the query row, (2) the match row, (3) the translation row, (4) the target row, and (5) the posterior probability row. If the pHMM was built from an MSA containing consensus structure or reference annotations, those will be visible on separate CS and RF rows above the query row.  There are also three types of columns: (1) a match in which a query amino acid is aligned to a target codon or quasi-codon, (2) a deletion in which the query amino acid is aligned to target gap characters, or (3) an insertion in which the target codon is aligned to a query gap character. 
   
       The query row begins with the name of the query HMM followed by the coordinates of the first amino acid on that line of the alignment, and ends with the coordinates of the last amino acid on that line of the alignment. For each column, the query row shows either the query consensus letter, for matches and deletions, or a gap character ('.') for insertions. 
   
       The target row begins with the name of the target sequence, followed by the coordinates of the first nucleotide on that line of the alignment, and ends with the coordinates of the last nucleotide on that line of the alignment. For each column, the target row shows the target codon (or quasicodon if using frameshift search) that has been aligned to the query. In the case of a deletion, the target line will print three gap characters ('---') in place of the codon. 
   
       The translation row shows the amino acid translations of the codons in the target row.  The match row shows which columns in the alignment are positive scoring. Exact matches are shown as the matched amino acid residue (in lowercase), and positive-scoring mismatches are shown as a '+'. Finally, the posterior probability (PP) row gives the expected accuracy for each position of the alignment.
   
       Below is the first line of the first alignment in the file PTH2.out.  It contains a CS row in addition to the five basic rows detailed above.  All the columns in this line of the alignment are matches. 
```
  Alignment:
  score: 110.6 bits
                   E    E    E    E    E    E    E    E    S    C    C    S    S    -    H    H    H    H    H    H    H    H    H    H    H
      PTH2   2     l    k    q    v    i    v    v    r    t    d    l    k    m    g    k    G    k    l    a    a    q    v    a    h    a     26
                   +    k         v    +    v    v    r    t    d    l         m    +    k    G    k    +    a    a    q    +    +    h    a
                   V    K    L    V    L    V    V    R    T    D    L    G    M    T    K    G    K    I    A    A    Q    C    S    H    A
      seq1 672    GTG  AAG  CTT  GTG  CTG  GTT  GTG  AGG  ACA  GAT  CTG  GGC  ATG  ACC  AAA  GGC  AAA  ATC  GCC  GCC  CAG  TGC  TCG  CAT  GCA    598
                   8    9    9    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *     PP
```
      
   5) Query Footer - each query's output will conclude with a footer that provides information about the hit filtering process inside bathsearch.  The average user can ignore this data.  For those who are interested, more information on these data can be found on page 54 of the [HMMER user guide](http://eddylab.org/software/hmmer/Userguide.pdf).  There will also be a couple of lines listing run times and a line with '//', indicating the end of the output for the current query.
   
       Below is the query footer from PTH2.out.

```
Internal pipeline statistics summary:
-------------------------------------
Query model(s):                            1  (116 nodes)
Target sequence(s):                        1  (6000 residues searched)
Residues passing SSV filter:            1503  (0.251); expected (0.02)
Residues passing bias filter:           1503  (0.251); expected (0.02)
Residues passing Vit filter:            1401  (0.234); expected (0.001)
Residues passing Fwd filter:            1287  (0.214); expected (1e-05)
Total number of hits:                      4  (0.135)
# CPU time: 0.05u 0.01s 00:00:00.06 Elapsed: 00:00:00.04
# Mc/sec: 14.21
//
```  
    
</p>
</details>

<details><summary>Practice 10: running a bathsearch search on a target with an alternate codon translation table</summary>
<p>

As discussed in Practice 2, some DNA sequences use alternate codon translation tables, and the best results are achieved by specifying the correct codon table both when building the HMMs and when performing the search. For this reason, the codon table of your HMM file must match the codon table of your search. The following commands will convert the file MET.bhmm to codon table 4 and perform a codon table 4 search. 

```bash
% bathconvert --ct 4 MET-ct4.bhmm  MET.bhmm 
% bathsearch --ct 4 -o MET-ct4.out MET-ct4.bhmm  target-MET.fa
```   
</p>
</details>

<details><summary>Practice 11: running bathsearch with a sequence query</summary>
<p>

If you do not wish to build the query HMMs before searching, you can use a sequence file (MSA or unaligned) as the query, and bathsearch will build the HMMs on the fly. You can use the '--hmmout' flag to save the HMMs to a file, or the HMMs will be written in the /tmp/ directory. The following command uses the single unaligned sequence in the file AMP_N.fa as the query.
   
```bash
% bathsearch --hmmout AMP_N.bhmm -o AMP_N.out AMP_N.fa target-AMP_N.fa
```
   
The file AMP_N.bhmm will now contain a single HMM, and AMP_N.out will contain output for a single hit between that HMM and the DNA sequence in target-AMP_N.fa.
</p>
</details>

<details><summary>Practice 12: running bathsearch with frameshift-aware alignment</summary>
<p>

BATH no longer uses frameshift-aware alignment by default but instead requires the --fs flag to activate this alignment mode. The following command will run the same search as in Practice 11, but with frameshift-aware alignment and reusing the HMM file created with --hmmout. 
   
```bash
% bathsearch --fs -o AMP_N-fs.out AMP_N.bhmm target-AMP_N.fa
```
   
If you compare the alignments in AMP_N.out and AMP_N-fs.out you will see that bathsearch with --fs detected and aligned both frameshifts and a stop codon in the target sequence, resulting in a longer alignment with a better e-value. The frameshift and stop codon counts are reported in both the query header and the annotation lines. 

Below are some examples of what frameshift and stop codon columns will look like in the alignment. 
```
| one nucleotide insertion | two nucleotide deletion | stop codon as a match | stop codon as an insertion |
|            a             |           l             |          a            |              .             |
|            +             |           l             |          +            |                            |
|            S             |           L             |          S            |              *             |
|          TCaA            |          --G            |         TaA           |             TGA            |
|            9             |           9             |          9            |              3             |
```

You can use the --frameline flag to make it easier to spot frameshifts in alignment.  This flag will add a "FRAME" line to the alignment that will give the translation frame of each codon.

```bash
% bathsearch --fs --frameline -o AMP_N-frameline.out AMP_N.bhmm target-AMP_N.fa
```
On the following line from AMP_N-frameline.out, you can easily see the stop codon, with frame 0, as well as the frameshift where the frame switches from 3 to 1. 

```bash
     AMP_N  98     f    .    .    g    l    d    d    a    f    p    i    e    d    l    d    e    i    l    p    e    l    l    e    n    r     120
                   f              +         d    +    a    +         i    +         l    d    e         +         e         l         n    +
                   F    *    s    S    S    D    E    A    Y    A    I    D    L    L    D    E    E    I    L    E    K    L    L    N    K
      seq1 297    TTT  TGA  AGC  AGC  AGT  GAT  GAG  GCT  TAT  GCA  ATA  GAT  TTA  --G  GAT  GAA  GAA  ATC  CTT  GAA  AAA  TTG  CTC  AAT  AAA    369
                   3    0    3    3    3    3    3    3    3    3    3    3    3    1    1    1    1    1    1    1    1    1    1    1    1     FRAME
                   9    3    1    3    5    7    8    8    8    8    8    6    6    9    8    7    7    7    7    7    7    8    8    9    *     PP
```

</p>
</details>

<details><summary>Practice 13: running bathsearch with splice-aware alignment</summary>
<p>

BATH can perform spliced alignment for eukaryotic genes split across two or more exons by using the --splice flag. Spliced alignment requires one additional step - producing an index (ssi) file for the target file. The following commands will create that index and perform a splice-aware search:
   
```bash
% esl-sfetch --index target-PTHR37536.fa
% bathsearch --splice -o PTHR37536.out PTHR37536.bhmm target-PTHR37536.fa
```
The alignment in PTHR37536.out is split into four exons. Below are the alignments for the second and third exons. Each exon is labeled at the beginning of the translation row (e.g., "exon 2").  Acceptor and donor site signals are shown on the target line in lower case with a "||" in the PP line. 
```bash
  PTHR37536  136          n    g    r    v    y    n    d    a     143
                          +         +    +    y         +
     exon 2               D    Q    Q    L    Y    R    N    I
       seq1  575 ag  TA  GAT  CAA  CAA  CTC  TAC  CGG  AAT  A   gt 602
                 ||       9    9    *    *    *    *    *    *  || PP


  PTHR37536  144          w    f    e    w    f    p    d    y    a    y    d    f    .    n    l    a    i    n    t    g    d    v    i     165
                          w    +    +    w         p    +         a    y         +         n    +         +              g    d         +
     exon 3               W    W    Q    W    V    P    N    G    A    Y    T    I    t    N    I    P    V    F    A    G    D    W    F
       seq1  685 ag  TC  TGG  TGG  CAA  TGG  GTC  CCC  AAT  GGC  GCT  TAC  ACG  ATT  ACC  AAC  ATT  CCC  GTC  TTT  GCT  GGT  GAT  TGG  TTT    757
                 ||       *    *    *    *    *    *    *    *    *    *    8    6    5    9    *    *    *    *    *    *    *    *    9     PP

  PTHR37536  166     v    a    k    v    e    a    l    s    p    s    n    g    v    a    i     180
                                    +                   s         +         +              i
     exon 3          D    I    T    I    N    T    T    S    S    T    -    A    A    T    I
       seq1  758    GAC  ATC  ACC  ATC  AAT  ACA  ACC  TCC  TCC  ACA  ---  GCA  GCG  ACC  AT  gt 800
                     9    9    9    9    8    7    7    6    4    4    .    3    3    3    4  || PP
```

The number of exons for each alignment is listed in both the query header and the annotation lines.
</p>
</details>


<details><summary>Practice 14: producing and interpreting tabular output from bathsearch</summary>
<p>
   
In addition to the standard alignment output, bathsearch can also produce tabular summary files with the use of the following flags:
```bash
--tblout     : tabular hit list, works with all alignment modes
--exontblout : tabular exon list, works this spliced alignment only
```

<details><summary>Using --tblout</summary>
<p>
   
The following command will run the same search as in Practice 8, but with the addition of the '--tblout' flag directing the tabular output to the file PTH2.tbl.
   
```bash
% bathsearch -o PTH2.out --tblout PTH2.tbl PTH2.bhmm target-PTH2.fa
``` 

If you open the file PTH2.tbl you will see the following text:

```
# hit ID  target name         accession  query name           accession   hmm len  hmm from    hmm to   seq len  ali from    ali to  env from    env to    E-value  score  bias   PID  description of target
#------- ------------------- ---------- -------------------- ---------- --------- --------- --------- --------- --------- --------- --------- ---------  --------- ------ ----- ----- ---------------------
       1 seq1                 -          PTH2                 PF01981.11      116         2       116       3000       672       325       675       325     3e-35  110.6   0.3 42.74 -
       2 seq1                 -          PTH2                 PF01981.11      116        35       116       3000      1486      1731      1444      1731   9.4e-28   86.4   0.0 43.90 -
       3 seq1                 -          PTH2                 PF01981.11      116        71       113       3000      2468      2343      2483      2325   3.5e-12   36.2   0.0 60.47 -
       4 seq1                 -          PTH2                 PF01981.11      116         2        30       3000      1273      1359      1270      1383   4.2e-12   36.0   0.3 68.97 -
#
# Program:         bathsearch
# Query file:      PTH2.bhmm
# Target file:     target-PTH2.fa
# Option settings: /Users/genevievekrause/git/BATH/src/bathsearch -o PTH2.out --tblout PTH2.tbl PTH2.bhmm target-PTH2.fa
# Current dir:     /Users/genevievekrause/git/BATH/tutorial
# Date:            Fri May  1 18:49:03 2026
# [ok]
``` 

A brief description of each column header (from left to right) is provided below.

```
hit ID                  A numerical ID for each query's hits. Hit IDs will repeat with multiple queries.

target name             Name of the target sequence where the hit is located.

accession               Alphanumeric ID for the target sequence (if provided in the target file).

query name              Name of the query pHMM. 

accession               Alphanumeric ID for the query (if provided in the query file).

hmm len                 Length (number of match states) of the query pHMM.
   
hmm from                The start position of the alignment on the query pHMM. 
   
hmm to                  The end position of the alignment on the query pHMM. 
   
seq len                 The length of the target sequence (in nucleotides).
   
ali from                The start position of the alignment on the target sequence. If the hit is located on the reverse complement strand, ali from will be greater than ali to.
   
ali to                  The end position of the alignment on the target sequence.
   
env from                The start position of the hit envelope on the target sequence. 
   
env to                  The end position of the hit envelope on the target sequence. 

E-value                 The hit e-value. 
   
score                   The hit bit score.
   
bias                    The hit bias adjustment score.

PID                     The alignment percent identity - defined as the # of columns with exact matches over the total number of columns, including insertions and deletions. 
   
description of target   Description of the target sequence (if provided in the target file).
```

The --cigar flag can be used to replace the "description of target" column with an alignment CIGAR string. The numbers in the CIGAR strings are counted in nucleotides to make them compatible with frameshift and intron counts.  99M means 33 consecutive matches (3 nucleotides each), 6I means two consecutive insertions, and 3D means a single deletion. 

```bash
% bathsearch -o PTH2.out --cigar --tblout PTH2-cigar.tbl PTH2.bhmm target-PTH2.fa
```

```bash
# hit ID  target name         accession  query name           accession   hmm len  hmm from    hmm to   seq len  ali from    ali to  env from    env to    E-value  score  bias   PID CIGAR
#------- ------------------- ---------- -------------------- ---------- --------- --------- --------- --------- --------- --------- --------- ---------  --------- ------ ----- ----- ---------------------
       1 seq1                 -          PTH2                 PF01981.11      116         2       116       3000       672       325       675       325     3e-35  110.6   0.3 42.74 99M6I192M3D51M
```
<details><summary>Using --tblout with --fs</summary>
<p>

Using --tblout with --fs will produce two additional columns showing frameshift and stop codon counts. 

```bash
% bathsearch --fs --cigar --tblout AMP_N-fs.tbl -o AMP_N-fs.out AMP_N.bhmm target-AMP_N.fa
```

```bash
# hit ID  target name         accession  query name           accession   hmm len  hmm from    hmm to   seq len  ali from    ali to  env from    env to    E-value  score  bias   PID  shifts  stops CIGAR
#------- ------------------- ---------- -------------------- ---------- --------- --------- --------- --------- --------- --------- --------- ---------  --------- ------ ----- ----- ------- ------ ---------------------
       1 seq1                 -          AMP_N                -               134         1       131        411         1       402         1       411   1.9e-27   82.8   0.1 46.32       6      1 44M1F39M1B114M9I25M2B19M1B44M1B4M6I30M2B67M
```
The CIGAR string shows insertion frameshifts as either 1F or 2F and deletion frameshifts as either 1B or 2B.

</p>
</details>

<details><summary>Using --tblout with --splice</summary>
<p>

Using --tblout with --splice will produce an exon cnt column. The ali and env coordinates given are for the full gene, from the start of the first exon to the end of the last exon. To get start and end coordinates of individual exons, you need to use --exontblout. 

```bash
% bathsearch --splice --cigar --tblout PTHR37536.tbl -o PTHR37536.out PTHR37536.bhmm target-PTHR37536.fa
```

```bash

# hit ID  target name         accession  query name           accession   hmm len  hmm from    hmm to   seq len  ali from    ali to  exon cnt    E-value  score  bias   PID CIGAR
#------- ------------------- ---------- -------------------- ---------- --------- --------- --------- --------- --------- --------- ---------  --------- ------ ----- ----- ---------------------
       1 seq1                 -          PTHR37536            -               279        11       251       1300       119      1159        4    2.8e-28   87.9   5.2 30.33 24786M3D66M3I94M85N24M86N38M3I60M3D11M153N28M6D21M3D129M3I27M
```
The CIGAR string shows introns as N. 85N is an intron that is 85 nucleotides long. 

</p>
</details>

</p>
</details>

<details><summary>Using --exontblout</summary>
<p>

--exontblout give tabular per-exon data for splice-aware alignment.  Run the following command to produce tabular exon data:

```bash
% bathsearch --splice --tblout PTHR37536.tbl --exontblout PTHR37536.extbl  -o PTHR37536.out PTHR37536.bhmm target-PTHR37536.fa
```
The file PTHR37536.extbl should look like this:

```bash
#                                                                                            ------ full hit ------  ------------------------- this exon --------------------------
# hit ID  target name          accession  query name           accession   hmm len   seq len   E-value  score  bias   #  of  hmm from    hmm to  ali from    ali to   P-value   PID
#------- -------------------  ---------- -------------------- ---------- --------- --------- --------- ------ ----- --- --- --------- --------- --------- --------- --------- -----
       1 seq1                 -          PTHR37536            -                279      1300   2.8e-28   87.9   5.2   1   4        11       135       119       491   9.6e-19 30.16
       1 seq1                 -          PTHR37536            -                279      1300   2.8e-28   87.9   5.2   2   4       136       143       577       600   7.3e-05 12.50
       1 seq1                 -          PTHR37536            -                279      1300   2.8e-28   87.9   5.2   3   4       144       180       687       798   4.1e-10 26.32
       1 seq1                 -          PTHR37536            -                279      1300   2.8e-28   87.9   5.2   4   4       181       251       952      1159   4.9e-12 34.72
#
# Program:         bathsearch
# Query file:      PTHR37536.bhmm
# Target file:     target-PTHR37536.fa
# Option settings: /Users/genevievekrause/git/BATH/src/bathsearch --splice -o PTHR37536.out --tblout PTHR37536.tbl --exontblout PTHR37536.extbl PTHR37536.bhmm target-PTHR37536.fa
# Current dir:     /Users/genevievekrause/git/BATH/tutorial
# Date:            Fri May  1 19:45:18 2026
# [ok]
```

The layout for --exontblout output is similar to --domtblout from HMMER. There is the same basic hit information as the --tbltout, hit ID (matches the hit ID for the same hit in --tblout), target name and accession, query name and accession, and HMM and sequence length.  E-value, hit score, and bias score are then reported for the whole gene underneath "full hit". The rest of the columns are exon-specific and are explained below.

```
#                       The exon count for this exon

of                      the total number of exons in this alignment

hmm from                The HMM start position of the exon. 
   
hmm to                  The HMM end position of the exon.
   
ali from                The DNA sequence start position of the exon
   
ali to                  The DNA sequence end position of the exon

P-value                 The exon p-value. 

PID                     The exon percent identity 

```

</p>
</details>
   
</p>
</details>




