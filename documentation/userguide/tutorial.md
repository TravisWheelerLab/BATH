# Tutorial

BATH was built on top of the existing HMMER3 code base. Users who are familiar with HMMER will find that BATH uses many of the same conventions. This tutorial will focus on getting you familiar with the five BATH tools listed below and, to avoid redundancy, will link to the HMMER user guide where applicable. 

There are two sections in this tutorial. The first section - Input files - will cover the tools that are used to prepare your data before you begin a search. The second section - Running bathsearch - will focus on using BATH to perform translated homology search, and on interpreting the search results. All the necessary files to complete the practices are located in the directory BATH/tutorial/. **You should cd into this directory before running the practice commands**. If you have not run 'make install' you will need to add the path to the BATH/src/ directory to the commands.

**Tools**
---

**bathsearch**      - search a DNA sequence database (or genome) for instances of one or more query proteins. The query can consist of a file of pHMMs (produced using bathbuild or bathconvert - see Practices 1 and 5 below) or a file containing sequences or sequence alignments (see Practice 11).
```
Usage: bathsearch [options] <protein-queryfile> <DNA-targetfile>
```

**bathbuild**   - build BATH formatted profile hidden Markov models (pHMMs) from input multiple sequence alignments (MSAs) or unaligned sequences and save to file
```
Usage: bathbuild [-options] <hmmfile_out> <msa_or_seq_file_in>
```
**bathstat**   - show summary statistics for a BATH formated pHMM file
```
Usage: bathstat [-options] <hmmfile_in>
```
**bathconvert** - convert HMMER formated pHMM files to BATH formated pHMM files
```
Usage: bathconvert [-options] <hmmfile_out> <hmmfile_in>
```
**bathfetch**   - copy selected pHMMs from an HMMER or BATH formatted file (converting if necessary) to a new BATH formated pHMM file
```
Usage: bathfetch [options] <hmmfile_in> <key>         (retrieves HMM named <key>)
Usage: bathfetch [options] -f <hmmfile_in> <keyfile>  (retrieves all HMMs in <keyfile>)
Usage: bathfetch [options] --index <hmmfile_in>       (indexes <hmmfile>)
```


## Section 1 - Input files 

Before you begin using BATH, it will be helpful to become familiar with the file types it uses. Running bathsearch requires two input files; a protein query file and a DNA target file. The target file must contain one or more DNA sequences in a recognizable unaligned sequence or multiple sequence alignment (MSA) format. Accepted unaligned sequence formats include fasta, embl, and genbank. Accepted MSA formats include stockholm, a2m, afa, psiblast, clustal, and phylip. 

BATH's installation includes a branch of the [Easel](https://github.com/EddyRivasLab/easel) software suite developed by the Eddy/Rivas Lab.  The Easel miniapps are a set of tools designed to perform a number of operations on MSA and unaligned sequence files.  To familiarize yourself with those tools see the [HMMER user guide](http://eddylab.org/software/hmmer/Userguide.pdf) (pages 145-204). 

A bathsearch query file contains the proteins you wish to search for in the target DNA. The preferred format for query files is a BATH formated pHMM file (although you may also use an MSA or unaligned sequence file - see practice 11). If you are not interested in taking advantage of BATH's frameshift-aware algorithms, you can also use a HMMER formatted pHMM file along with the '--nofs' flag. The rest of this section will focus on practices to get you acquainted with the BATH tools that are used to create and manipulate BATH formated pHMM files.

<details><summary>Practice 1: building a pHMM from an MSA using bathbuild</summary>
<p>

The sensitivity of BATH is powered, in large part, by the use of pHMMs. The pHMM files used by BATH are nearly identical to the ones used by HMMER, but contain additional information needed to perform accurate frameshift-aware translations and provide reliable e-values for the alignments. This additional information includes the frameshift rate and codon translation table to be used by bathsearch as well as tau and lambda values that define the curve for the pHMMs score distribution from the frameshift-aware Forward algorithm. If you would like more information on pHMM files see the [HMMER user guide](http://eddylab.org/software/hmmer/Userguide.pdf) (page 208). 
   
BATH formated pHMMs can be created from MSA files or unaligned sequence files using the tool bathbuild. The file MET.msa contains two stockholm formatted protein MSAs (note that stockholm is the only MSA format that allows multiple MSAs in a single file). You can build pHMMs from those MSAs and save them to the file MET.bhmm by running the following command: 
   
```bash
   % bathbuild MET.bhmm MET.msa
```
(note the file suffix '.bhmm' - this can help distinguish BATH formated pHMM files from HMMER formatted ones, which often have the suffix '.hmm')

The summary output that is printed to your stdout should resemble the text below (the exact CPU and elapsed time will vary):

```bash
# input      file:                  MET.msa
# output HMM file:                  MET.bhmm
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# idx    name                  nseq  alen  mlen codon_tbl eff_nseq re/pos description
# ------ -------------------- ----- ----- ----- --------- -------- ------ -----------
  1      metC                    11   483   409         1     0.60  0.591 Cystathionine beta-lyase
  2      metG                    24   494   458         1     0.62  0.589 Methionine--tRNA ligase

# CPU time: 4.38u 0.00s 00:00:04.38 Elapsed: 00:00:02.33
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

In practice 10 you will search pHMMs built from MET.msa against a target sequence from the genome of an endosymbiotic bacteria that uses codon table 4. Running the following command will build the pHMMs using the correct codon table for that target:
   
```bash
   % bathbuild --ct 4 MET-ct4.bhmm MET.msa
```
   
The summary output that is printed to your stdout should resemble the text below (the exact CPU and elapsed time will vary):

```bash
# input file:                       MET.msa
# output HMM file:                  MET.bhmm
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# idx    name                  nseq  alen  mlen codon_tbl eff_nseq re/pos description
# ------ -------------------- ----- ----- ----- --------- -------- ------ -----------
  1      metC                    11   483   409         1     0.60  0.591 Cystathionine beta-lyase
  2      metG                    24   494   458         1     0.62  0.589 Methionine--tRNA ligase

# CPU time: 4.38u 0.00s 00:00:04.38 Elapsed: 00:00:02.33
```

</p>
</details>

<details><summary>Practice 3: building a pHMM from an unaligned sequences file</summary>
<p>

If your queries are single unaligned sequences rather than MSAs you can still build HMMs using bathbuild and the flag '--unali'. The file three_seqs.fa contains three unaligned protein sequences. To build three HMMs (one for each sequence) use the following command. 
   
```bash
   % bathbuild --unali three_seqs.bhmm three_seqs.fa
```

The summary output that is printed to your stdout should resemble the text below (the exact CPU and elapsed time will vary):

```bash
# input file:                       three_seqs.fa
# output HMM file:                  three_seqs.bhmm
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# idx    name                  nseq   len  mlen codon_tbl eff_nseq re/pos description
# ------ -------------------- ----- ----- ----- --------- -------- ------ -----------
  1      AT1G01010.1              1   429   429         1     1.00  0.591
  2      AT1G01020.1              1   245   245         1     1.00  0.552
  3      AT1G01030.1              1   358   358         1     1.00  0.600

# CPU time: 2.62u 0.00s 00:00:02.62 Elapsed: 00:00:02.14
```

</p>
</details>

<details><summary>Practice 4: summerizing pHMM files with bathstat</summary>
<p>

Since a pHMM file may contain any number of individual models, it is useful to be able to quickly summarize the contents. The tool bathstat is designed to provide such a summary for BATH formated pHMM files.  The following command will summarize the pHMM file built in practice 1:
   
```bash
   % bathstat MET.bhmm
```
   
This command should produce the following output to stdout:

```bash
#
# idx    name                 accession        nseq eff_nseq   mlen codon_tbl re/pos
# ------ -------------------- ------------ -------- -------- ------ --------- ------
  1      metC                 -                  11     0.60    409         1   0.53
  2      metG                 -                  24     0.62    458         1   0.53
```

The fields are mainly the same as those produced by bathbuild, and detailed in practice 1, except for the accession field which may contain an alphanumeric identifier for the protein family or be left blank if no accession is listed for the pHMM. 

</p>
</details>

<details><summary>Practice 5: converting an HMMER formated pHMM file to BATH format using bathconvert</summary>
<p>

If you have an existing HMMER formatted pHMM file and want to use it to run bathsearch with frameshift detectionG you will first need to convert it to the BATH format using bathconvert. The file tRNA-proteins.hmm contains 12 pHMMs in HMMER3 format. The following command will create the BATH formatted file tRNA-proteins.bhmm containing the same twelve pHMMs:

```bash
   % bathconvert  tRNA-proteins.bhmm  tRNA-proteins.hmm
```
Your summary output should match that shown below.
   
```
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
```

</p>
</details>

<details><summary>Practice 6: changing the codon table of a BATH formated pHMM using bathconvert</summary>
<p>
 
You can also use bathconvert to change the codon table of an existing BATH pHMM file using the --ct flag. This is faster than rebuilding from the original MSA.  The following command will create the file tRNA-proteins-ct11.bhmm containing the same 12 pHMMs as tRNA-proteins.bhmm but modified to use NCBI codon translation table 11:
   
```bash
   % bathconvert --ct 11 tRNA-proteins-ct11.bhmm tRNA-proteins.bhmm
```

This should produce the following output:
 
```
# input HMM file:                   tRNA-proteins.bhmm
# output HMM file:                  tRNA-proteins=ct11.bhmm
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
# CPU time: 8.65u 0.00s 00:00:08.65 Elapsed: 00:00:08.67
```
 
</p>
</details>

<details><summary>Practice 7: indexing a pHMM file and copying a single pHMM using bathfetch </summary>
<p>

If you only need to search with a single pHMM but it is located in a file with multiple pHMMs, you can save time by copying the desired pHMM to a new file using bathfetch. If the original file contains a large number of pHMMs, you may want to create an index file to speed up the fetch process.  The following command will create the index file tRNA-proteins.bhmm.ssi for the BATH pHMM file created in Practice 4. 
```bash
   % bathfetch --index tRNA-proteins.bhmm 
```
The summary output should read as follows:
   
```
Working...    done.
Indexed 12 HMMs (12 names and 12 accessions).
SSI index written to file tRNA-proteins.bhmm.ssi
```
Whether or not you choose to create an index you will need the name of the pHMM you wish to copy to use as a key. The command below will copy the pHMM PTH2 from the tRNA-proteins.bhmm.  The -o flag will direct the copied pHMM to the specified output file (PTH2.hmm in this case). Otherwise, it will be printed to standard out. 
```bash
   %  bathfetch -o PTH2.bhmm tRNA-proteins.bhmm PTH2
```
The summary output should simply read as:
```
Retrieved HMM PTH2.
```
</p>
</details>

<details><summary>Practice 8: copying and converting multiple pHMMs using bathfetch </summary>
<p>

You can also use bathfetch to copy multiple pHMMs. To do so you will need to create a key file that contains the names of all the pHMMs you wish to copy, with one name per line, and use the -f flag. If the original pHMM file is in HMMER format bathfetch will automatically convert it to BATH format. The following command will copy all 3 of the pHMMs listed in the key file tRNA-synthetases-key.txt from an HMMER formated pHMM file, convert them to BATH format, and print them to the output file tRNA-synthetases.bhmm.
   
```bash
   % bathfetch -f -o tRNA-synthetases.bhmm tRNA-proteins.hmm tRNA-synthetases-key.txt
```
   
The summary output should simply read as:
   
```
Retrieved 3 HMMs.
```
   
As with bathconvert, you can also use the --ct flag with bathfetch to change the codon table.
   
</p>
</details>

## Section 2 - Running bathsearch

This section of the tutorial will focus on the tool bathsearch. This tool allows the user to perform translated annotation of protein-coding DNA even when mutations or sequencing errors have introduced frameshifts. Each of the practices in this section will involve running bathsearch with a different set of input formats, options, and outputs. 

<details><summary>Practice 9: running a simple bathsearch and reading the output</summary>
<p>

Bathsearch requires two inputs - the query and the target.  In this practice, you will use the single pHMM in the file PTH2.bhmm as the query.  For the target, you will use a single DNA sequence in the file target-PTH2.fa. The -o flag is used to direct the standard output to the file PTH2.out (otherwise it will be printed to stdout).
   
```bash
   % bathsearch -o PTH2.out PTH2.bhmm target-PTH2.fa
```
 
We can now open the file PTH2.out and see that the output is organized into the following sections:
     
   1) File Header - lines begin with '#' and contain basic information about the search parameters
```
# query HMM file:                  PTH2.bhmm
# target sequence database:        target-PTH2.fa
# frameshift probability:          0.010000
# codon translation table          1
# output directed to file:         PTH2.out
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
```
    
   2) Query Header - includes a summary of each query and a hits list sorted by E-value.  For each hit, the query header lists the E-value, bit score, and bias score adjustment (for more information on bias scores see pages 60-61 of the [HMMER user guide](http://eddylab.org/software/hmmer/Userguide.pdf).  This is followed by the name of the target sequence where the hit was located, the target sequence position for the start and end of the alignment, the number of frameshifts and stop codons in that alignment, and finally a target description (which may be blank).

```
Query:       PTH2  [M=116]
Accession:   PF01981.11
Description: Peptidyl-tRNA hydrolase PTH2
Scores for complete hits:
    E-value  score  bias  Sequence  start    end  shifts  stops  Description
    ------- ------ -----  --------  -----  -----  ------  -----  -----------
    1.3e-34  110.1   0.3  seq1        672    325       0      0
    7.4e-27   90.7   0.0  seq1       2670   2343       3      1
    2.1e-25   86.0   0.0  seq1       1486   1731       0      0
    3.3e-11   40.2   0.0  seq1       1273   1387       1      0
```
   
   3) Annotation Lines- for each hit listed in the query header, bathsearch will produce an annotation line containing useful information about the hit. After the line 'Annotation for each hit (and alignments):' these annotation lines (as well as the alignments) will appear, sorted first by target sequence and then by e-value.
   
       As in the query header, the annotations line lists the score, bias, and E-value for each hit. It also lists three types of coordinates for the hit - the alignment start and end coordinates for both the query (hmm-from & hmm-to) and the target (ali-from & ali-to), as well as the envelope coordinates (env-from & env-to). The envelope is the region of the target that bathsearch has identified as containing the homology (the alignment coordinates are always contained within the envelope). It is the envelope coordinates that bound the target subsequence used to calculate the score, bias, and E-value. An explanation of the characters seen after the coordinates ('.','[', & ']') can be found on page 38 of the [HMMER user guide](http://eddylab.org/software/hmmer/Userguide.pdf). The annotation line also lists the number of frameshifts and stop codons in the alignment (shift & stops), the full length of the target sequence (sq-len), and the alignment's accuracy score (acc) which is the average expected per residue accuracy of the alignment. 
       
       Below is the annotation line for the first hit in the file PTH2.out
  
```
Annotation for each hit (and alignments):
>> seq1
    score  bias    Evalue   hmm-from    hmm-to     ali-from    ali-to     env-from    env-to    shifts  stops    sq-len    acc
   ------ ----- ---------   --------   -------    --------- ---------    --------- ---------    ------  ----- ---------   ----
 !  110.1   0.3   1.3e-34          2       116 .]       672       325 ..       675       325 ..      0      0      3000   0.92
```
   
   4) Alignment - Below each annotation line bathsearch prints the alignment for that query-target hit. A typical bathsearch alignment will contain at least the following five rows (in order from top to bottom): (1) the query row, (2) the match row, (3) the translation row, (4) the target row, and (5) the posterior probability row. If the pHMM was built from an MSA containing consensus structure or reference annotations those will be visible on separate CS and RF rows above the query row.  There are also three types of columns: (1) a match in which a query amino is aligned to a target codon or quasi-codon, (2) a deletion in which the query amino acid is aligned to target gap characters, or (3) an insertion in which the target codon is aligned to a query gap character. 
   
       The query row begins with the name of the query pHMM followed by the coordinates of the first amino acid on that line of the alignment and ends with the coordinates of the last amino acid on that line of the alignment. For each column, the query row shows either the query consensus letter, for matches and deletions, or a gap character ('.') for insertions. 
   
       The target row begins with the name of the target sequence followed by the coordinates of the first nucleotide on that line of the alignment and ends with the coordinates of the last nucleotide on that line of the alignment. For each column, the target row shows the target codon or pseudo-codon that has been aligned to the query. In the case of a deletion, the target line will print three gap characters ('---') in place of the codon. 
   
       The translation row shows the amino acid translations of the codons and pseudo-codons on the target row.  The match row shows which columns in the alignment are positive scoring. Exact matches are shown as the matched amino acid residue (in lowercase) and positive scoring mismatches are shown as a '+'. Finally, the posterior probability (PP) row gives the expected accuracy for each position of the alignment.
   
       Below is the first line of the first alignment in the file PTH2.out.  It contains a CS row in addition to the five basic rows detailed above.  All the columns in this line of the alignment are matches and, as there are no frameshifts, no quasi-codons are present. 

```
  Alignment:
  score: 110.1 bits
             E    E    E    E    E    E    E    E    S    C    C    S    S    -    H    H    H    H    H    H    H    H    H    H    H    H    H  CS
  PTH2   2   l    k    q    v    i    v    v    r    t    d    l    k    m    g    k    G    k    l    a    a    q    v    a    h    a    a    v   28
             +    k         v    +    v    v    r    t    d    l         m    +    k    G    k    +    a    a    q    +    +    h    a    +
             V    K    L    V    L    V    V    R    T    D    L    G    M    T    K    G    K    I    A    A    Q    C    S    H    A    T    L
  seq1 672  GTG  AAG  CTT  GTG  CTG  GTT  GTG  AGG  ACA  GAT  CTG  GGC  ATG  ACC  AAA  GGC  AAA  ATC  GCC  GCC  CAG  TGC  TCG  CAT  GCA  ACG  CTC  592
             8    9    9    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    8    7   PP
```
      
   5) Query Footer - each query's output will conclude with a footer that provides information about the hit filtering process inside bathsearch.  The average user can ignore this data.  For those who are interested, more information on these data can be found on page 54 of the [HMMER user guide](http://eddylab.org/software/hmmer/Userguide.pdf).  There will also be a couple of lines listing run times and a line with just '//', indicating the end of the output for the query.
   
       Below is the query footer from PTH2.out.

```
Internal pipeline statistics summary:
-------------------------------------
Query model(s):                            1  (116 nodes)
Target sequence(s):                        1  (6000 residues searched)
Residues passing SSV filter:            1503  (0.251); expected (0.02)
Residues passing bias filter:           1503  (0.251); expected (0.02)
Residues passing Vit filter:            1401  (0.234); expected (0.001)
Residues passing Fwd filter:            1961  (0.327); expected (1e-05)
Total number of hits:                      4  (0.173)
# CPU time: 0.04u 0.01s 00:00:00.05 Elapsed: 00:00:00.05
# Mc/sec: 12.20
//
```
   
   6) File Footer - If bathsearch did not encounter any errors the last line of the file will read '[ok]'
     
    
</p>
</details>

<details><summary>Practice 10: running a bathsearch search on a target with an alternate codon translation table</summary>
<p>

As discussed in Practice 2, some DNA sequences use alternate codon translation tables and the best results are achieved by specifying the correct codon table both when building the pHMMs and when performing the search. To prevent searches with mismatched codon tables, bathsearch responds to such searches with an error message. Running the following command will attempt a mismatches search by searching the pHMMs in MET.bhmm, built with the standard codon table, against the target DNA in the file target-MET.fa while specifying the use of the alternate codon table 4. 
   
```bash
   % bathsearch --ct 4 -o MET.out MET.bhmm target-MET.fa
```
   
This will result in the following error message:

```bash
   Error: Requested codon translation tabel ID 4 does not match the codon translation tabel ID of the HMM file MET.bhmm. Please either run bathsearch with option '--ct 1' or run bathconvert with option '--ct 4'.
```

In this case, we already have a pHMM file built with the correct codon table and can skip running bathconvert. The following command will use that pHMM file to run the same search without a codon table mismatch. 

```bash
   % bathsearch --ct 4 -o MET-ct4.out MET-ct4.bhmm  target-MET.fa
```
The file MET-ct4.out should contain a single hit between each of the pHMMS in MET-ct4.bhmm and the DNA sequence in target-MET.fa.
   
</p>
</details>

<details><summary>Practice 11: running bathsearch with a sequence query</summary>
<p>

If you do not wish to build the query pHMMs ahead of time you can use a sequence file (MSA or unaligned) as the query and bathsearch will build the pHMMs on the fly. However, depending on the number and length of the proteins, building pHMMs can be time-consuming. If you choose to use a sequence query file it is recommended that you use the '--hmmout' flag to save the pHMMs for use in any subsequent searches. The following command uses the single unaligned sequence in the file gidA.fa as the query, building a pHMM for that sequence and printing it the to file gidA.bhmm. The use of the '--ct' flag will determine the codon table used both to build the pHMM and to conduct the search. The standard output is directed to the file gidA.out using the '-o' flag. 
   
```bash
   % bathsearch --ct 4 --hmmout gidA.bhmm -o gidA.out gidA.fa target-gidA.fa
```
   
The file gidA.bhmm will now contain a single pHMM built with codon table 4 and gidA.out will contain output for a single git between that pHMM and the DNA sequence in target-gidA.fa.
</p>
</details>

<details><summary>Practice 12: producing and interpreting tabular output from bathsearch</summary>
<p>
   
In addition to the standard output, bathsearch can also produce a tabular summary file with a more easily parsable list of the hits found in a search. By using the '--tblout' flag you can direct bathsearch to create this tabular output and save it to the file of your choosing. The following command will run the same search as in Practice 8, but with the addition of the '--tblout' flag directing the tabular output to the file PTH2.tbl.
   
```bash
   % bathsearch -o PTH2.out --tblout PTH2.tbl PTH2.bhmm target-PTH2.fa
``` 

If you open the file PTH2.tbl you will see the following text (file directories and dates may vary):

```
# target name         accession  query name           accession   hmm len  hmm from    hmm to   seq len  ali from    ali to  env from    env to   E-value  score  bias  shifts  stops  pipe description of target
#------------------- ---------- -------------------- ---------- --------- --------- --------- --------- --------- --------- --------- --------- --------- ------ ----- ------- ------ ----- ---------------------
seq1                 -          PTH2                 PF01981.11       116         2       116      3000       672       325       675       325   1.3e-34  110.1   0.3       0      0   std -
seq1                 -          PTH2                 PF01981.11       116         1       113      3000      2670      2343      2677      2327   7.4e-27   90.7   0.0       3      1    fs -
seq1                 -          PTH2                 PF01981.11       116        35       116      3000      1486      1731      1426      1731   2.1e-25   86.0   0.0       0      0    fs -
seq1                 -          PTH2                 PF01981.11       116         2        40      3000      1273      1387      1269      1491   3.3e-11   40.2   0.0       1      0    fs -
#
# Program:         bathsearch
# Query file:      PTH2.bhmm
# Target file:     target-PTH2.fa
# Option settings: bathsearch -o PTH2.out --tblout PTH2.tbl PTH2.bhmm target-PTH2.fa
# Current dir:     BATH/tutorial
# Date:            Wed Mar  1 18:27:58 2023
# [ok]
``` 

A brief description of each column header (from left to right) is provided below.

```
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
   
env from                The start position of the hit envelope on the target sequence.  See practice 9 for an explanation of hit envelopes. 
   
env to                  The end position of the hit envelope on the target sequence. 

E-value                 The hit e-value. 
   
score                   The hit bit score.
   
bias                    The hit bias adjustment score. See pages 60-61 of the [HMMER user guide](http://eddylab.org/software/hmmer/Userguide.pdf) for an explanation of bias adjustment. 
   
shifts                  The number of frameshifts reported in the alignment.
   
stops                   The number of stop codons reported in the alignment.
   
pipe                    The translation pipeline used to produce the hit, either standard 'std' or frameshift aware 'fs'. To reduce runtimes and the potential for false frameshifts bathsearch uses standard (non frameshift aware) translation for hits it determines are unlikely to contain frameshifts.  Most users can ignore this column. 
   
description of target   Description of the target sequence (if provided in the target file).
```

</p>
</details>

<details><summary>Practice 13: locating frameshifts and stop codons in bathsearch alignments using '--frameline' </summary>
<p>
   
While both the standard and tabular outputs give the user the count of frameshifts and stop codons in an alignment, the user may also want to locate the quasi and stop codons.  Quasi-codons with deletions can be identified by looking for codons with one or two '-' characters in place of a nucleotide. Quasi-codons with insertions can be identified by looking for codons with more than 4 or 5 nucleotides (the nucleotides bathsearch determines to be the insertions will be shown in lowercase).  Stop codons can be identified by looking for codons with all two uppercase and one lowercase nucleotide (the nucleotide bathsearch determined to be a substitution will be in lowercase). Below are examples of quasi and stop codons taken from the alignment in gidA.out from Practice 11:
 
```
 | one nucleotide deletion | two nucleotide deletions | one nucleotide insertion | two nucleotide insertions | stop codon |
 |            w            |            g             |            k             |            v              |      a     |
 |            w            |            g             |            +             |            v              |      +     |
 |            W            |            G             |            Q             |            V              |      S     |
 |           -GA           |           --A            |          CtAA            |          GTtaT            |     TaA    |
 |            6            |            3             |            8             |            7              |      9     |
```

To make it easier to locate frameshifts and stop codons the '--frameline' flag can be used to add a row to the alignment that numbers the frame of each codon and quasi-codon. This line can be used to locate quasi-codons by looking for a change from one frame to another.  Stop codons can be identified on the frameline by a '0'. Note that, for hits on the reverse complement strand of a sequence, the frames will be negative (i.e. -1, -2, & -3). 
   
 Running the follwing comand will use the file gidA.bhmm, created in Practice 11, to search a single pHMM against the DNA sequence in the file target-gidA.fa using codon table 4. The '-o' flag will direct the standard output to the file gidA-frameline.out and the '--frameline' flag will add the frameline row to the alignment.
   
```bash
   % bathsearch --ct 4 --frameline -o gidA-frameline.out gidA.bhmm target-gidA.fa
```
   
The following is an excerpt of four lines from the alignment in gidA-frameline.out. This excerpt shows two frameshifts (one by deletion and one by insertion) as well as one stop codon (the frame is shown directly beneath each codon or quasi-codon). On the first line, the frame changes - from 2 to 1 -  due to the deletion of one nucleotide. There is a stop codon on the second line, with a 0 in the frameline.  On the third line, the frame changes again - from 1 to 2 - due to a single nucleotide insertion. 
 
```
  gidA   241   g    p    r    y    c    l    s    i    e    g    k    t    l    k    f    g    r    k    p    q    k    l    i    m    e    p   266
                    p         +    c                                                                               +         i    +    e    p
               H    P    P    F    C    N    Q    T    -    K    M    K    H    P    -    T    I    N    I    Y    R    P    I    I    E    P
  seq1 27209  CAT  CCA  CCC  TTT  T-C  AAC  CAA  ACT  ---  AAA  ATG  AAA  CAT  CCA  ---  ACC  ATA  AAC  ATC  TAC  AGA  CCC  ATC  ATA  GAA  CCA  27279
               2    2    2    2    1    1    1    1    .    1    1    1    1    1    .    1    1    1    1    1    1    1    1    1    1    1   FRAME
               7    7    7    7    1    1    1    1    .    0    0    0    0    0    .    4    4    5    6    7    8    8    9    *    *    *   PP

  gidA   267   e    a    t    g    s    s    s    v    y    v    n    g    l    s    t    s    m    .    p    i    e    l    q    l    q    l   291
               e    +              +         +    v    +    +    n    g         s         s                   +         +    q    l         +
               E    S    L    D    T    K    T    V    H    L    N    G    T    S    I    S    T    s    N    L    V    I    Q    L    N    I
  seq1 27280  GAA  TaA  CTT  GAT  ACT  AAA  ACT  GTA  CAC  CTA  AAT  GGT  ACC  TCT  ATC  TCA  ACC  TCC  AAT  CTG  GTA  ATC  CAA  CTT  AAC  ATA  27357
               1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1   FRAME
               *    9    9    9    *    *    *    *    *    *    *    *    9    9    7    7    5    1    5    8    8    9    *    *    *    *   PP

  gidA   292   l    k    f    t    k    a    f    r    g    a    k    i    i    k    a    g    y    a    i    e    y    d    c    v    c    s   317
               l    k         t    +                                  +    +    k    +         +         i    e    y    d              c    s
               L    K    P    T    Q    H    P    I    D    V    C    V    V    K    S    K    H    T    I    E    Y    D    V    T    C    S
  seq1 27358  CTA  AAA  CCC  ACA CtAA  CAC  CCA  ATT  GAT  GTC  TGT  GTT  GTT  AAG  TCC  AAA  CAC  ACC  ATA  GAA  TAC  GAT  GTA  ACT  TGT  TCA  27436
               1    1    1    1    2    2    2    2    2    2    2    2    2    2    2    2    2    2    2    2    2    2    2    2    2    2   FRAME
               *    *    *    *    8    7    7    8    9    9    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *   PP
```
 
</p>
</details>

<details><summary>Practice 14: locating frameshifts and stop codons with a tabular output </summary>
<p>
   
While the frameline makes it easier to find frameshifts and stop codons in individual alignments, some users may want to see the locations of frameshifts and stop codons across multiple alignments and in a  more parseable form. For this reason, the flag '--fstblout' allows the user to create a tabular output of frameshift and stop codon locations for all hits. The following command reruns the same search as in Practice 13, but rather than using '--frameline' it uses '--fstblout' to save frameshift and stop codon locations to the file gidA.fstbl.  
   
```bash
   % bathsearch --ct 4 -o gidA-fstbl.out --fstblout gidA.fstbl gidA.bhmm target-gidA.fa
``` 

If you open the file gidA.fstbl you will see the following text (file directories and dates may vary):

```
# target name         accession  query name           accession  E-value   ali from  ali to     I D S  length  seq start  ali start
#------------------- ----------- -------------------- ---------- --------- --------- ---------  -----  ------  ---------  ---------
 seq1                 -          gidA                 -              3e-13 26678     27759          I       2  26747             70
 seq1                 -          gidA                 -              3e-13 26678     27759          D       1  26815            138
 seq1                 -          gidA                 -              3e-13 26678     27759          I       1  26958            281
 seq1                 -          gidA                 -              3e-13 26678     27759          D       2  26995            318
 seq1                 -          gidA                 -              3e-13 26678     27759          D       1  27029            352
 seq1                 -          gidA                 -              3e-13 26678     27759          D       2  27076            399
 seq1                 -          gidA                 -              3e-13 26678     27759          D       2  27116            439
 seq1                 -          gidA                 -              3e-13 26678     27759          D       1  27141            464
 seq1                 -          gidA                 -              3e-13 26678     27759          D       1  27221            544
 seq1                 -          gidA                 -              3e-13 26678     27759          I       1  27370            693
 seq1                 -          gidA                 -              3e-13 26678     27759          I       1  27647            970
 seq1                 -          gidA                 -              3e-13 26678     27759          I       1  27696           1019
#
# Program:         bathsearch
# Query file:      gidA.bhmm
# Target file:     target-gidA.fa
# Option settings: /home/u3/gkrause/git/BATH/src/bathsearch -o gidA-fstbl.out --fstblout gidA.fstbl --ct 4 gidA.bhmm target-gidA.fa
# Current dir:     /home/u3/gkrause/git/BATH/tutorial
# Date:            Mon Mar 11 14:40:19 2024
# [ok]
```

A brief description of each column header (from left to right) is provided below.

```
target name             Name of the target sequence where the hit is located.

accession               Alphanumeric ID for the target sequence (if provided in the target file).

query name              Name of the query pHMM. 

accession               Alphanumeric ID for the query (if provided in the query file).

E-value                 The hit e-value. 
   
ali from                The start position of the alignment on the target sequence. If the hit is located on the reverse complement strand, ali from will be greater than ali to.
   
ali to                  The end position of the alignment on the target sequence.
   
I D S                   Type of codon/quasi-codons.  I - insertion.  D - deletion. S - stop codon.
   
length                  Number of nucleotides inserted or deleted (0 for stop codons).
   
seq start               Position (in the target sequence) of the start of the quasi-codon or stop codon.
   
ali start               Position (in the alignment) of the start of the quasi-codon or stop codon.  The first nucleotide in the alignment = 1.
```
</p>
</details>


