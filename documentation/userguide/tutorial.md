# Tutorial

FraHMMER<sup>\*</sup> was built on top of the existing HMMER3 code-base. Users who are familiar with HMMER will find that FraHMMER uses many of the same conventions. This tutorial will focus on getting you familiar with the five FraHMMER tools listed below and, to avoid redundency, will link to the HMMER user guide where applicable. 

There are two sections in this tutorial. The first section - Input files - will cover the tools that are used to prepare your data before you begin a frahmmer<sup>\*</sup> search. The second section - Running frahmmer searches - will focus on using frahmmer to perform frameshift aware tranlsated homology search, and on interpreting the search results. All the necessary files to complete the practices are located in the directory FraHMMER/tutorial/. **You should cd into this directory before running the practice commands**. If you have not run 'make install' you will need to add the path to the FraHMMER/src/ directory to the commands.

<sup>\* 'FraHMMER' with capitalizations refers to the software package which includes a varaity of tools used to support and preform frameshift aware tranalsted homology search. When writen in all lowercase, 'frahmmer' refers to the primiary search tool inside of FraHMMER.</sup>


**Tools**
---

**frahmmbuild**   - build FraHMMER formated profile hidden Markov models (pHMMs) from input multiple sequence alignments (MSAs) and save to file
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

Before you begin using FraHMMER, it will be helpful to become familiar with the file types it uses. Each frahmmer search reequires two input files; a protein query file and a DNA target file. The target file contains one or more DNA sequences in a recognizable unaligned sequence or multiple sequence alignment (MSA) format. Accepted unaligned sequence formats include fasta, embl, and genbank. Accepted MSA formats include stockholm, a2m, afa, psiblast, clustal, and phylip. 

FraHMMER's installation includes the [Easel](https://github.com/EddyRivasLab/easel) software suite developed by the Eddy/Rivas Lab.  The Easel miniapps are a set of tools designed to perform a number of operations on MSA and unaligned sequence files.  To familiarize yourself with those tools see the [HMMER user guide](http://eddylab.org/software/hmmer/Userguide.pdf) (pages 145-204). 

The query file contains the proteins you wish to search for in the target DNA. The preferred format for query files is a FraHMMER formated pHMM file (although you may also use an MSA or unaligned sequence file - see practice # and #). The rest of this section will focus on practices to get you acquainted with the FraHMMER tools which are used to create and manipulate these pHMM files.

<details><summary>Practice 1: building a pHMM from an MSA using frahmmbuild</summary>
<p>

The sensitivity of FraHMMER is powered, in large part, by the use of pHMMs. The pHMM files used by FraHMMER and nearly identical to the ones used by HMMER, but contain additional information needed to perform accurate frameshift-aware translations and provide reliable e-values. This additional information includes the frameshift rate and codon translation table to be used in the frahmmer search as well as tau and lambda values that define the curve for the pHMMs score distribution from the frameshift-aware Forward algorithm. If you would like more information on the other information in pHMM files see the [HMMER user guide](http://eddylab.org/software/hmmer/Userguide.pdf) (page 208). 
   
FraHMMER formated pHMMs can be created from MSA files using the tool frahmmbuild. The file MET.msa contains two stockholm formatted protein MSAs (note that stockholm is the only MSA format that allows multiple MSAs in a single file). You can build pHMMs from those MSAs and save them to the file MET.fhmm by running the following command: (note the file suffix '.fhmm' - this can help distinguish FraHMMER formated pHMM files from HMMER formatted ones, which often have the suffix '.hmm')
   
```bash
   % frahmmbuild MET.fhmm MET.msa
```
Now compare the summary output that is printed to your stdout to the text below (the exact CPU and elapsed time will vary):

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

eff_nseq       Effective sequence number. This was the “effective” number of independent sequences that frahmmbuild’s default “entropy weighting” step decided on, given the phylogenetic similarity of the nseq sequences in the input alignment. 

re/pos         Mean positional relative entropy, in bits. This can be ignored by most users. 
   
description    Description of the protein family - may be blank.
```
</p>
</details>

<details><summary>Practice 2: building a pHMM from an MSA using frahmmbuild with an alternate codon translation table</summary>
<p>

One of the fields that distinguishes a FraHMMER formatted pHMM file from a HMMER formated pHMM file is an [NCBI codon translation table ID](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). The correct codon table depends on the origins of the target DNA you intend to search the pHMMs against. When you run a frahmmer search, selecting the correct codon table will produce the highest quality alignments. Ensuring that the pHMMs were built with that same the codon table will produce the most accurate e-values for those alignments. By default, frahmmbuild will use the standard code employed by eukaryotic nuclear DNA. To use an alternate codon translation table include the option --ct followed by a table ID from the list below:

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

In practice 8 you will search the pHMMs in MET.msa against a target sequence from the genome of an endosymbiotic bacteria that uses codon table 4. Run the following command to build the pHMMs using the correct codon table for that target:
   
```bash
   % frahmmbuild --ct 4 MET-ct4.fhmm MET.msa
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
# idx    name                 accession        nseq eff_nseq   mlen fs_prob codon_tbl re/pos
# ------ -------------------- ------------ -------- -------- ------ ------- --------- ------
  1      metC                 -                  11     0.60    409 0.01000         1   0.53
  2      metG                 -                  24     0.62    458 0.01000         1   0.53
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

This section of the tutorial will focus on the tool frahmmer. This tool allows the user to perform translated annotate of protein-coding DNA even when mutations or sequencing errors have introduced frameshifts. Each of the practices in this section will involve running a frahmmer search with a different set of input formats,  options, and outputs. 

<details><summary>Practice 7: running a simple frahmmer search and reading the output</summary>
<p>

Every frahmmer search requires two inputs - the query and the target.  In this practice, you will use the single pHMM in the file Rib.hmm) as the query.  For the target, you will use a single DNA sequence in the file seq1.fa, and the -o flag will be used to direct the hit data and alignment to the file RIB.out. 
   
```bash
   % frahmmer -o RIB.out RIB.hmm seq1.fa
```
 
 The file RIB.out should now contain a single hit between the Ribosomal_S19e protein family and the DNA sequence. If you open this file you will see that it contains the following information:
     
   1) File Header - lines begin with '#' and contain basic information about the search parameters
```
   # query HMM file:                  Rib.hmm
   # target sequence database:        seq1.fa
   # frameshift probability:          0.010000
   # codon translation table          1
   # output directed to file:         Rib.out
   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
```
    
   2) Query Header - includes a summary of each query and a hits list sorted by E-value.  For each hit, the query header lists the E-value, bit score, and bias score adjustment (for more information on bias scores see pages 60-61 of the [HMMER user guide](http://eddylab.org/software/hmmer/Userguide.pdf).  This is followed by name of the target sequence where the hit was located, the target sequence position for the start and end of the alignment, the number of frameshifts and stop codons in that alignment, and finally a target description (which may be blank).

```
   Query:       Ribosomal_S19e  [M=139]
   Accession:   PF01090.14
   Description: Ribosomal protein S19e
   Scores for complete hits:
      E-value  score  bias  Sequence   start     end  shifts  stops  Description
      ------- ------ -----  --------   -----   -----  ------  -----  -----------
      1.6e-28  110.7   0.1  seq1     3197980 3197584       1      0
```
   
   3) Annotation - for each hit between the query and target, frahmmer will produce an annotation line containing useful information about the hit. As in the query header, the annotations line lists the score, bias, and E-value for each hit. It also lists three types of coordinates for the hit. These include the alignment start and end coordinates for the query (hmm-from & hmm-to) and for the target (ali-from & ali-to), as well as the envelope coordinates (env-from & env-to).  The envelope is the region of the target that frahmmer has identified as homologous and the hit alignment is always contained within the envelope. The reported score, bias, and E-value are all calculated for the target subsequence bound by the envelope coordinates. The annotation line also lists the frameshift and stop codon counts, the full length of the target sequence, and the alignment's accuracy score (acc) which is the average expected per residue accuracy of the alignment.
  
```
   Annotation for each hit (and alignments):
   >> seq1
       score  bias    Evalue   hmm-from    hmm-to     ali-from    ali-to     env-from    env-to    shifts  stops    sq-len    acc
      ------ ----- ---------   --------   -------    --------- ---------    --------- ---------    ------  ----- ---------   ----
    !  110.7   0.1   1.6e-28          8       137 ..   3197980   3197584 ..   3198010   3197575 ..      1      0  30000000   0.89
```
   
   4) Alignment - Bellow each annotation line is the alignment (here just the first line of the alignment is shown). This alignment contains 5 rows which are, from top to bottom, (1) the query row, (2) the match row, (3) the translation row, (4) the target row, and (5) the posterior probability row. The query row contains the query consensus letter for matched and deletions and a '.' for insertions.  The target row shows the target codons and pseudo-codons which have been aligned to the target.  In this example, only codons are present in the first row (no frameshifts).  In the case of an amino acid deletion target line will print three dashed '---' in place of the codon. Both the query and target row also show the name of the query or target and the start and end coordinates of the residues sown on that line of the alignment. The translation line shows the amino acid translations of the codons on the target line.  The match line shows which positions in the alignment are positive scoring.  Exact matches are shown as the matched amino acid residue (in lowercase) and positive scoring mismatches are shown as a '+'.   Finally, the posterior probability (PP) row gives the expected accuracy for each position of the alignment.

```
     Alignment:
     score: 110.7 bits
     Ribosomal_S19e       8   a    d    k    l    i    e    k    v    a    e    e    l    k    e    k    d    k    i    k    p    p    e    W   30
                              +    +    +    +    i    +              a    +         l    k    +    +         +    +    +         p    +    W
                              P    Q    E    F    I    A    T    Y    A    R    F    L    K    K    T    G    R    V    Q    I    P    K    W   
               seq1 3197980  CCA  CAA  GAA  TTC  ATT  GCT  ACC  TAC  GCA  AGA  TTC  TTA  AAG  AAA  ACT  GGT  CGT  GTT  CAA  ATC  CCA  AAA  TGG  3197912
                              5    7    8    9    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *    *   PP
```
      
   5) Query Footer - each query's output will conclude with a footer that provides information about the hit filtering process inside frahmmer.  The average user can ignore this data.  For those who are interested, more information on these data can be found on page 54 of the [HMMER user guide](http://eddylab.org/software/hmmer/Userguide.pdf).  There will also be a couple of lines listing run times and a line with just '//' indicating the end of the output for the query.

```
   Internal pipeline statistics summary:
   -------------------------------------
   Query model(s):                            1  (139 nodes)
   Target sequence(s):                        1  (60040812 residues searched)
   Residues passing SSV filter:         1327090  (0.0221); expected (0.02)
   Residues passing bias filter:        1222149  (0.0204); expected (0.02)
   Residues passing Vit filter:          129784  (0.00216); expected (0.001)
   Residues passing Fwd filter:             653  (1.09e-05); expected (1e-05)
   Total number of hits:                      1  (6.61e-06) 
   # CPU time: 3.86u 0.04s 00:00:03.90 Elapsed: 00:00:01.91
   # Mc/sec: 4360.71
   // 
```
   
   6) File Footer - If frahmmer did not encounter any errors the last line of the file will simply read '[ok]'
    
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
   Error: Requested codon translation tabel ID 4 does not match the codon translation tabel ID of the HMM file MET.hmm. Please run frahmmcovert with option '--ct 4'.
```

To avoid this error we need to use the pHMM file with the correct codon translation table by running the following command:

```bash
   % frahmmer --ct 4 -o MET.out MET-ct4.hmm seq2.fa
```
The file MET.out should contain a single hit between each of the pHMMS in MET-ct4.hmm and the DNA sequence.
   
</p>
</details>

<details><summary>Practice 9: running a frahmmer search with a sequence query</summary>
<p>

If you do not wish to build the query pHMMs ahead of time, frahmmer can build them for you on the fly. However, depending on the number and length of the proteins, building pHMMs can be time-consuming.  If you chose to use a sequence query file (unaligned sequences or MSAs) it is recommended that you save the pHMMs to use in any subsequent searches.  The following command uses the unaligned sequences in the file Rib-Seqs.fa as the queries, building a pHMM for each one.   The '--hmmout' flag will direct frahmmer to print those pHMMs to the file Rib-Seqs.hmm.

```bash
   % frahmmer --hmmout Rib-Seqs.hmm -o Rib-Seqs.out Rib-Seqs.fa seq1.fa
```
</p>
</details>

