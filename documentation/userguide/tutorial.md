# Tutorial

FraHMMER\* was built on top of the existing HMMER3 code base. Users who are familiar with HMMER will find that FraHMMER uses many of the same conventions. This tutorial will focus on getting you familiar with the five FraHMMER tools listed below and, to avoid redundancy, will link to the HMMER user guide where applicable. 

There are two sections in this tutorial. The first section - Input files - will cover the tools that are used to prepare your data before you begin a frahmmer\*\* search. The second section - Running frahmmer searches - will focus on using frahmmer to perform frameshift aware translated homology search, and on interpreting the search results. All the necessary files to complete the practices are located in the directory FraHMMER/tutorial/. **You should cd into this directory before running the practice commands**. If you have not run 'make install' you will need to add the path to the FraHMMER/src/ directory to the commands.

<sup>\* 'FraHMMER' with capitalizations refers to the software package which includes a variety of tools used to support and perform frameshift aware translated homology search. </sup>

<sup>\*\* When written in all lowercase, 'frahmmer' refers to the frameshift aware translated homology search tool inside of FraHMMER.</sup>


**Tools**
---

**frahmmbuild**   - build FraHMMER formatted profile hidden Markov models (pHMMs) from input multiple sequence alignments (MSAs) and save to file
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
**frahmmfetch**   - copy selected pHMMs from an HMMER or FraHMMER formatted file (converting if necessary) to a new FraHMMER formated pHMM file
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

Before you begin using FraHMMER, it will be helpful to become familiar with the file types it uses. Each frahmmer search requires two input files; a protein query file and a DNA target file. The target file contains one or more DNA sequences in a recognizable unaligned sequence or multiple sequence alignment (MSA) format. Accepted unaligned sequence formats include fasta, embl, and genbank. Accepted MSA formats include stockholm, a2m, afa, psiblast, clustal, and phylip. 

FraHMMER's installation includes the [Easel](https://github.com/EddyRivasLab/easel) software suite developed by the Eddy/Rivas Lab.  The Easel miniapps are a set of tools designed to perform a number of operations on MSA and unaligned sequence files.  To familiarize yourself with those tools see the [HMMER user guide](http://eddylab.org/software/hmmer/Userguide.pdf) (pages 145-204). 

A frahmmer query file contains the proteins you wish to search for in the target DNA. The preferred format for query files is a FraHMMER formated pHMM file (although you may also use an MSA or unaligned sequence file - see practice 9). The rest of this section will focus on practices to get you acquainted with the FraHMMER tools which are used to create and manipulate FraHMMER pHMM files.

<details><summary>Practice 1: building a pHMM from an MSA using frahmmbuild</summary>
<p>

The sensitivity of FraHMMER is powered, in large part, by the use of pHMMs. The pHMM files used by FraHMMER and nearly identical to the ones used by HMMER, but contain additional information needed to perform accurate frameshift-aware translations and provide reliable e-values. This additional information includes the frameshift rate and codon translation table to be used in the frahmmer search as well as tau and lambda values that define the curve for the pHMMs score distribution from the frameshift-aware Forward algorithm. If you would like more information on the other information in pHMM files see the [HMMER user guide](http://eddylab.org/software/hmmer/Userguide.pdf) (page 208). 
   
FraHMMER formated pHMMs can be created from MSA files using the tool frahmmbuild. The file MET.msa contains two stockholm formatted protein MSAs (note that stockholm is the only MSA format that allows multiple MSAs in a single file). You can build pHMMs from those MSAs and save them to the file MET.fhmm by running the following command: (note the file suffix '.fhmm' - this can help distinguish FraHMMER formated pHMM files from HMMER formatted ones, which often have the suffix '.hmm')
   
```bash
   % frahmmbuild MET.fhmm MET.msa
```
   
The summary output that is printed to your stdout should reselmble the text below (the exact CPU and elapsed time will vary):

```bash
# input alignment file:             MET.msa
# output HMM file:                  MET.fhmm
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# idx    name                  nseq  alen  mlen fs_prob codon_tbl eff_nseq re/pos description
# ------ -------------------- ----- ----- ----- ------- --------- -------- ------ -----------
  1      metC                    11   483   409 0.01000         1     0.60  0.591 Cystathionine beta-lyase
  2      metG                    24   494   458 0.01000         1     0.62  0.589 Methionine--tRNA ligase

# CPU time: 4.38u 0.00s 00:00:04.38 Elapsed: 00:00:02.33
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

eff_nseq       Effective sequence number. This is the “effective” number of independent sequences that frahmmbuild’s default “entropy weighting” step decided on, given the phylogenetic similarity of the nseq sequences in the input alignment. The higher the number the more diversity there is amoung the sequences in the MSA. 

re/pos         Mean positional relative entropy, in bits. This can be ignored by most users. 
   
description    Description of the protein family - may be blank.
```
</p>
</details>

<details><summary>Practice 2: building a pHMM from an MSA using frahmmbuild with an alternate codon translation table</summary>
<p>

One of the fields that distinguishes a FraHMMER formatted pHMM file from an HMMER formated pHMM file is an [NCBI codon translation table ID](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). The correct codon table depends on the origins of the target DNA you intend to search the pHMMs against. When you run a frahmmer search, selecting the correct codon table will produce the highest quality alignments. Ensuring that the pHMMs were built with that same codon table will produce the most accurate e-values for those alignments. 
   
By default, frahmmbuild will use the standard code employed by eukaryotic nuclear DNA. To use an alternate codon translation table include the option --ct followed by a table ID from the list below:

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

In practice 8 you will search the pHMMs in MET.msa against a target sequence from the genome of an endosymbiotic bacteria that uses codon table 4. Running the following command will build the pHMMs using the correct codon table for that target:
   
```bash
   % frahmmbuild --ct 4 MET-ct4.fhmm MET.msa
```
   
The summary output should be nearly identical to that in Practice 1, except for the output file name and the codon table field which should now say 4 for both pHMMs. 

</p>
</details>

<details><summary>Practice 3: summerizing pHMM files with frahmmstat</summary>
<p>

Since a pHMM file may contain any number of individual models, it is useful to be able to quickly summarize the contents. The tool frahmmstat is designed to provide such a summary for FraHMMER formated pHMM files.  The following command will summarize the pHMM file built in practice 1:
   
```bash
   % frahmmstat MET.fhmm
```
   
This command should produce the following output to stdout:

```bash
#
# idx    name                 accession        nseq eff_nseq   mlen fs_prob codon_tbl re/pos
# ------ -------------------- ------------ -------- -------- ------ ------- --------- ------
  1      metC                 -                  11     0.60    409 0.01000         1   0.53
  2      metG                 -                  24     0.62    458 0.01000         1   0.53
```

The fields are mainly the same as those produced by frahmmbuild, and detailed in practice 1, except for the accession field which may contain an alphanumeric identifier for the protein family or be left blank if no accession is listed for the pHMM. 

</p>
</details>

<details><summary>Practice 4: converting an HMMER formated pHMM file to FraHMMER format using frahmmconvert</summary>
<p>

If you have an existing HMMER formatted pHMM file and want to use it to run a frahmmer search, you will first need to convert it to the FraHMMER format using frahmmconvert. The file tRNA-proteins.hmm contains 12 pHMMs in HMMER3 format. The following command will create the FraHMMER formatted file tRNA-proteins.fhmm containing the same three pHMMs:

```bash
   % frahmmconvert XXX.fhmm XXX.hmm
```
Your summary output should match that shown below.
   
```
# input HMM file:                   tRNA-proteins.hmm
# output HMM file:                  tRNA-proteins.fhmm
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# idx    name                  nseq  mlen fs_prob codon_tbl eff_nseq re/pos description
# ------ -------------------- ----- ----- ------- --------- -------- ------ -----------
  1      ATE_N                   30    78 0.01000         1     1.11  0.726 Arginine-tRNA-protein transferase, N terminus
  2      GlutR_N                 12   152 0.01000         1     0.87  0.590 Glutamyl-tRNAGlu reductase, N-terminal domain
  3      PTH2                    10   116 0.01000         1     0.74  0.589 Peptidyl-tRNA hydrolase PTH2
  4      RtcB                    30   459 0.01000         1     0.83  0.590 tRNA-splicing ligase RtcB
  5      TGT                     15   238 0.01000         1     0.80  0.589 Queuine tRNA-ribosyltransferase
  6      Thg1                    30   131 0.01000         1     0.69  0.589 tRNAHis guanylyltransferase
  7      Trm56                   11   121 0.01000         1     0.64  0.590 tRNA ribose 2'-O-methyltransferase, aTrm56
  8      tRNA-synt_1_2           30   185 0.01000         1     0.91  0.590 Leucyl-tRNA synthetase, Domain 2
  9      tRNA-synt_1c_C          14   192 0.01000         1     0.81  0.591 tRNA synthetases class I (E and Q), anti-codon binding domain
  10     tRNA-synt_2d            19   247 0.01000         1     0.73  0.592 tRNA synthetases class II core domain (F)
  11     tRNA-Thr_ED             12   136 0.01000         1     0.63  0.590 Archaea-specific editing domain of threonyl-tRNA synthetase
  12     TruB_C                  11    56 0.01000         1     1.64  0.994 tRNA Pseudouridine synthase II, C terminal
```

</p>
</details>

<details><summary>Practice 5: changing the codon table of a FraHMMER formated pHMM using frahmmconvert</summary>
<p>
 
You can also use frahmmconvert to change the codon table of an existing FraHMMER pHMM file using the --ct flag. This is faster than rebuilding from the original MSA.  The following command will create the file tRNA-proteins-ct11.fhmm containing the same 12 pHMMs as tRNA-proteins.fhmm but modified to use NCBI codon translation table 11:
   
```bash
   % frahmmconvert --ct 11 tRNA-proteins-ct11.fhmm tRNA-proteins.fhmm
```

This should produce the following output:
 
```
# input HMM file:                   tRNA-proteins.fhmm
# output HMM file:                  tRNA-proteins=ct11.fhmm
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# idx    name                  nseq  mlen fs_prob codon_tbl eff_nseq re/pos description
# ------ -------------------- ----- ----- ------- --------- -------- ------ -----------
  1      ATE_N                   30    78 0.01000        11     1.11  0.726 Arginine-tRNA-protein transferase, N terminus
  2      GlutR_N                 12   152 0.01000        11     0.87  0.590 Glutamyl-tRNAGlu reductase, N-terminal domain
  3      PTH2                    10   116 0.01000        11     0.74  0.589 Peptidyl-tRNA hydrolase PTH2
  4      RtcB                    30   459 0.01000        11     0.83  0.590 tRNA-splicing ligase RtcB
  5      TGT                     15   238 0.01000        11     0.80  0.589 Queuine tRNA-ribosyltransferase
  6      Thg1                    30   131 0.01000        11     0.69  0.589 tRNAHis guanylyltransferase
  7      Trm56                   11   121 0.01000        11     0.64  0.590 tRNA ribose 2'-O-methyltransferase, aTrm56
  8      tRNA-synt_1_2           30   185 0.01000        11     0.91  0.590 Leucyl-tRNA synthetase, Domain 2
  9      tRNA-synt_1c_C          14   192 0.01000        11     0.81  0.591 tRNA synthetases class I (E and Q), anti-codon binding domain
  10     tRNA-synt_2d            19   247 0.01000        11     0.73  0.592 tRNA synthetases class II core domain (F)
  11     tRNA-Thr_ED             12   136 0.01000        11     0.63  0.590 Archaea-specific editing domain of threonyl-tRNA synthetase
  12     TruB_C                  11    56 0.01000        11     1.64  0.994 tRNA Pseudouridine synthase II, C terminal
# CPU time: 8.65u 0.00s 00:00:08.65 Elapsed: 00:00:08.67
```
 
</p>
</details>

<details><summary>Practice 6: indexing a pHMM file and copying a single pHMM using frahmmfetch </summary>
<p>

If you only need to search with a single pHMM but it is located in a file with multiple pHMMs, you can save time by copying the desired pHMM to a new file using frahmmfetch. If the original file contains a large number of pHMMs, you may want to create an index file to speed up the fetch process.  The following command will index the create the index file tRNA-proteins.fhmm.ssi for the FraHMMER pHMM file created in Practice 4. 
```bash
   % frahmmfetch --index tRNA-proteins.fhmm 
```
The summary output should read as follows:
   
```
Working...    done.
Indexed 12 HMMs (12 names and 12 accessions).
SSI index written to file tRNA-proteins.fhmm.ssi
```
Whether or not you choose to create an index you will need the name of the pHMM you wish to copy to use as a key. The command below will copy the pHMM PTH2 from the tRNA-proteins.fhmm.  The -o flag will direct the copied pHMM to the specified output file (PTH2.hmm in this case). Otherwise, it will be printed to standard out. 
```bash
   % frahmmfetch -o PTH2.fhmm tRNA-proteins.fhmm PTH2
```
The summary output should simply read as:
```
Retrieved HMM PTH2.
```
</p>
</details>

<details><summary>Practice 7: copying and converting multiple pHMMs using frahmmfetch </summary>
<p>

You can also use frahmmfetch to copy multiple pHMMs. To do so you will need to create a key file that contains the names of all the pHMMs you wish to copy, with one name per line, and use the -f flag. If the original pHMM file is in HMMER format frahmmfetch will automatically convert it to FraHMMER format. The following command will copy all 3 of the pHMMs listed in the key file tRNA-synthetases-key.txt from an HMMER formated pHMM file, convert them to FraHMMER format, and print them to the output file tRNA-synthetases.fhmm.
   
```bash
   % frahmmfetch -f -o tRNA-synthetases.fhmm tRNA-proteins.hmm tRNA-synthetases-key.txt
```
   
The summary output should simply read as:
   
```
Retrieved 3 HMMs.
```
   
As with frahmmconvert, you can also use the --ct flag with frahmmfetch to change the codon table.
   
</p>
</details>

## Section 2 - Running frahmmer searches

This section of the tutorial will focus on the tool frahmmer. This tool allows the user to perform translated annotate of protein-coding DNA even when mutations or sequencing errors have introduced frameshifts. Each of the practices in this section will involve running a frahmmer search with a different set of input formats,  options, and outputs. 

<details><summary>Practice 8: running a simple frahmmer search and reading the output</summary>
<p>

Every frahmmer search requires two inputs - the query and the target.  In this practice, you will use the single pHMM in the file PTH2.fhmm as the query.  For the target, you will use a single DNA sequence in the file target-PTH2.fa. The -o flag is used to direct the standard output to the file PTH2.out. 
   
```bash
   % frahmmer -o PTH2.out PTH2.fhmm target-PTH2.fa
```
 
See Practice 9 for a breakdown of the frahmmer standard output in PTH2.out
  
    
</p>
</details>

<details><summary>Practice 9: interpreting frahmmer standard output</summary>
<p>

The file PTH2.out contains the standard frahmmer output for a search between the query file PTH2.fhmm and the target file target-PTH2.fa. If you open this file you will see that it is organized into the following sections:
     
   1) File Header - lines begin with '#' and contain basic information about the search parameters
```
# query HMM file:                  PTH2.fhmm
# target sequence database:        target-PTH2.fa
# frameshift probability:          0.010000
# codon translation table          1
# output directed to file:         PTH2.out
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
```
    
   2) Query Header - includes a summary of each query and a hits list sorted by E-value.  For each hit, the query header lists the E-value, bit score, and bias score adjustment (for more information on bias scores see pages 60-61 of the [HMMER user guide](http://eddylab.org/software/hmmer/Userguide.pdf).  This is followed by name of the target sequence where the hit was located, the target sequence position for the start and end of the alignment, the number of frameshifts and stop codons in that alignment, and finally a target description (which may be blank).

```
   Query:       PTH2  [M=116]
   Accession:   PF01981.11
   Description: Peptidyl-tRNA hydrolase PTH2
   Scores for complete hits:
    E-value  score  bias  Sequence     start    end  shifts  stops  Description
    ------- ------ -----  --------     -----  -----  ------  -----  -----------
    3.4e-34  110.1   0.3  PTH2-target    672    325       0      0
    4.2e-33  110.3   0.0  PTH2-target   1273   1731       0      0
    2.7e-27   91.6   0.2  PTH2-target   2659   2343       2      1
```
   
   3) Annotation Lines- for each hit listed in the query header, frahmmer will produce an annotation line containing useful information about the hit. After the line 'Annotation for each hit (and alignments):' these annotation lines (as well as the alignments) will appear, sorted first by target sequence and then by e-value.
   
       As in the query header, the annotations line lists the score, bias, and E-value for each hit. It also lists three types of coordinates for the hit - the alignment start and end coordinates for both the query (hmm-from & hmm-to) and the target (ali-from & ali-to), as well as the envelope coordinates (env-from & env-to). The envelope is the region of the target that frahmmer has identified as containing the homology (the hit alignment is always contained within the envelope). It is the envelope coordinates that bound the target subsequence used to calculate the score, bias, and E-value. An explanation of the characters seen after the coordinates ('.','[', & ']') can be found on page 38 of the [HMMER user guide](http://eddylab.org/software/hmmer/Userguide.pdf). The annotation line also lists the number of frameshifts and stop codons in the alignment (shift & stops), the full length of the target sequence (sq-len), and the alignment's accuracy score (acc) which is the average expected per residue accuracy of the alignment. 
       
       Below is the annotation line for the first hit in the file PTH2.out
  
```
Annotation for each hit (and alignments):
>> PTH2-target
    score  bias    Evalue   hmm-from    hmm-to     ali-from    ali-to     env-from    env-to    shifts  stops    sq-len    acc
   ------ ----- ---------   --------   -------    --------- ---------    --------- ---------    ------  ----- ---------   ----
 !  110.1   0.3   3.4e-34          2       116 .]       672       325 ..       675       325 ..      0      0      3000   0.92
```
   
   4) Alignment - Bellow each annotation line frahmmer prints the alignment for that query-target hit. A typical frahmmer alignment will contain at least the following five rows (in order from top to bottom): (1) the query row, (2) the match row, (3) the translation row, (4) the target row, and (5) the posterior probability row. If the pHMM was built from an MSA containing consensus structure or reference annotations those will be visible on separate CS and RF rows above the query row.  There are also three types of columns: (1) a match in which a query amino is aligned to a target codon or quasi-codon, (2) a deletion in which the query amino acid is aligned to target gap characters, or (3) an insertion in which the target codon is aligned to a query gap character. 
   
       The query row begins with the name of the query pHMM followed by the coordinates of the first amino acid on that line of the alignment and ends with the coordinates of the last amino acid on that line of the alignment. For each column, the query row shows either the query consensus letter, for matches and deletions, or a gap character ('.') for insertions. 
   
       The target row begins with the name of the target sequence followed by the coordinates of the first nucleotide on that line of the alignment and ends with the coordinates of the last nucleotide on that line of the alignment. For each column, the target row shows the target codon or pseudo-codon that has been aligned to the query. In the case of a deletion, the target line will print three gap characters ('---') in place of the codon. 
   
       The translation row shows the amino acid translations of the codons and pseudo-codons on the target row.  The match row shows which columns in the alignment are positive scoring. Exact matches are shown as the matched amino acid residue (in lowercase) and positive scoring mismatches are shown as a '+'. Finally, the posterior probability (PP) row gives the expected accuracy for each position of the alignment.
   
       Below is the first line of the first alignment in the file PTH2.out.  It contains a CS row in addition to the five basic rows detailed above.  All the columns in this line of the alignment are matches and, as there are no frameshifts, no quasi-codons are present. 

```
  Alignment:
  score: 110.1 bits
                    E    E    E    E    E    E    E    E    S    C    C    S    S    -    H    H    H    H    H    H    H    H    H    H    H  CS
         PTH2   2   l    k    q    v    i    v    v    r    t    d    l    k    m    g    k    G    k    l    a    a    q    v    a    h    a   26
                    +    k         v    +    v    v    r    t    d    l         m    +    k    G    k    +    a    a    q    +    +    h    a
                    V    K    L    V    L    V    V    R    T    D    L    G    M    T    K    G    K    I    A    A    Q    C    S    H    A
  PTH2-target 672  GTG  AAG  CTT  GTG  CTG  GTT  GTG  AGG  ACA  GAT  CTG  GGC  ATG  ACC  AAA  GGC  AAA  ATC  GCC  GCC  CAG  TGC  TCG  CAT  GCA  598
                    8    9    9    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *   PP
```
      
   5) Query Footer - each query's output will conclude with a footer that provides information about the hit filtering process inside frahmmer.  The average user can ignore this data.  For those who are interested, more information on these data can be found on page 54 of the [HMMER user guide](http://eddylab.org/software/hmmer/Userguide.pdf).  There will also be a couple of lines listing run times and a line with just '//', indicating the end of the output for the query.
   
       Bellow is the query footer from PTH2.out.

```
Internal pipeline statistics summary:
-------------------------------------
Query model(s):                            1  (116 nodes)
Target sequence(s):                        1  (6000 residues searched)
Residues passing SSV filter:            1503  (0.251); expected (0.02)
Residues passing bias filter:           1503  (0.251); expected (0.02)
Residues passing Vit filter:            1401  (0.234); expected (0.001)
Residues passing Fwd filter:            1961  (0.327); expected (1e-05)
Total number of hits:                      3  (0.187)
# CPU time: 0.04u 0.01s 00:00:00.05 Elapsed: 00:00:00.05
# Mc/sec: 12.20
//
```
   
   6) File Footer - If frahmmer did not encounter any errors the last line of the file will simply read '[ok]'
   
</p>
</details>


<details><summary>Practice 10: running a frahmmer search on a target with an alternate codon translation table</summary>
<p>

As discussed in Practice 2, some DNA sequences use alternate codon translation tables. For frahmmer searches that use such DNA as the target, the best results are achieved by specifying the correct codon table both during the actual search and when building the pHMMs. In this  Practice you will first attempt to run a frahmmer search with a mismatch between the codon table specified for the target and the codon table used to build the pHMM, resulting in an error message.  

Run the following command to search the pHMMs in MET.hmm, which you built in Practice 1 using the standard codon translation, against the target DNA in the file seq2.fa that you will specify as using codon table 4 with the --ct flag

```bash
   % frahmmer --ct 4 -o MET.out MET.hmm seq2.fa
```
   
This will result in the following error message:

```bash
   Error: Requested codon translation tabel ID 4 does not match the codon translation tabel ID of the HMM file MET.hmm. Please run frahmmcovert with option '--ct 4'.
```

To avoid this error we need to use the pHMM file with the correct codon translation table by running the following command:

```bash
   % frahmmer --ct 4 -o MET.out MET-ct4.hmm seq2.fa
```
The file MET.out should contain a single hit between each of the pHMMS in MET-ct4.hmm and the DNA sequence.
   
</p>
</details>

<details><summary>Practice 11: running a frahmmer search with a sequence query</summary>
<p>

If you do not wish to build the query pHMMs ahead of time, frahmmer can build them for you on the fly. However, depending on the number and length of the proteins, building pHMMs can be time-consuming.  If you chose to use a sequence query file (unaligned sequences or MSAs) it is recommended that you save the pHMMs to use in any subsequent searches.  The following command uses the unaligned sequences in the file Rib-Seqs.fa as the queries, building a pHMM for each one.   The '--hmmout' flag will direct frahmmer to print those pHMMs to the file Rib-Seqs.hmm.

```bash
   % frahmmer --hmmout Rib-Seqs.hmm -o Rib-Seqs.out Rib-Seqs.fa seq1.fa
```
</p>
</details>

<details><summary>Practice 12: producing and interpreting tabular output from a frahmmer search</summary>
<p>
   
</p>
</details>

<details><summary>Practice 13: finding frameshifts and stop codons in frahmmer alignments </summary>
<p>
   
</p>
</details>


