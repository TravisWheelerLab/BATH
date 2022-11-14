# Tutorial

Once you have installed and built FraHMMER you will ba able to use the main search tool - frahmmer - as well as several support tools to help create, format and anayisis profile hidden Markov model (pHMM) files. You will also have access to the suite of easel miniapps developed by the Eddy/Rivas Lab (for details on easel tools see the HMMER user guide http://eddylab.org/software/hmmer/Userguide.pdf). This tutorial will focus on getting you familiar with the FraHMMER specific tools that will allow you to perform frameshift-aware translated homology searchs.

## Step 1 - Preparing HMM files

Every frahmmer search requires two input files - a query and a target.  The target file must include one or more DNA sequences in a recognizable single sequence or multiple sequence alignment format. Common single sequence formats include: fasta, embl, and genbank. Common alignment formats include: stockholm, a2m, afa, psiblast, clustal, and phylip.

The query file may be either an (1) HMM, (2) multiple sequence alignment (MSA), or (3) unalgined sequence file. In case (2) and (3), each MSA or unalgined sequence will be coverted to a pHMM. Building pHMMs can be computationally expensive so it is recomended to save them to file to avoid having to reporoduce the same calculations in subsequent searches. build them with frahmmbuild before searching. 

If you prefer, you can also use an MSA or unalgined sequence file as the query and save the resulting HMMs to a file using the frahmmer flag '--hmmout'.


### Practice 1 - Build HMM with frahmmbuild.

**1)** The tool frahmmbuild takes two arguments.  
'''bash
   % Usage: frahmmbuild [-options] <hmmfile_out> <msafile>
'''
   
first the name of the file you would like the pHMM to printed to, and the name of the file containing the MSA(s) or the single sequence you would like ot build the pHMM form. 

After cding into your FraHMMER directory run:
```bash
   % frahmmbuild tutorial/xxx.hmm tutorial/xxx.msa
```
**2)** This will create the file xxx.hmm and output a FraHMMER foramted HMM to that file. To ensure that the file was built properly you can compare it to the pre-built file yyy.hmm.  If the diff comand bellow produces no output then frahmmbuild is working correctly. 
```bash
   % diff tutorial/xxx.hmm tutorial/yyy.hmm
```

### Practice 2 - Outpur HMM with --hmmout flag.

**1)** Build an HMM file from a single sequence file and print to an HMM file using the --hmmout flag for frahmmer
   Endsure you are still in the FraHMMER directory and run:
```bash
   % frahmmer --hmmout tutorial/aaa.hmm -o tutorial/aaa.aliout tutorial/aaa.fa tutorial/seq1.fa
```
   This will create the file aaa.hmm and output a FraHMMER foramted HMM to that file. To ensure that the file was built properly you can compare it to the pre-built file bbb.hmm.  If the diff comand bellow produces no output then the --hmmout flag is working correctly. 
```bash
   % diff tutorial/aaa.hmm tutorial/bbb.hmm
```      
**Optional:**. If you would like to understand the contents on HMM files follow this link : LINK GOES HERE

## Step 2 - running frahmmer

Now that you know all three methods for creating FraHMMER formated HMM files, you are ready to run a frahmmer homology search. The usage for frahmmer is as follows:
```bash
Usage: frahmmer [options] <hmm, msa, or seq query file> <seq target file>
```
We will cover a few basic optios, or flags, here but a full list can be viewed by running:
```bash
frahmmer -h
```
For a more detailed explination of all frahmmer options see: LINK GOES HERE


 You already ran frahmmer when you used the --hmmout flag to build the aaa.hmm file from a single sequence file.  The alignments for this search were directed to the file tutorial/aaa.aliout.  To check that this output is correct run the following:
```bash
   % diff tutorial/aaa.aliout tutorial/bbb.aliout
```    
**2)** Now you will run frahmmer using the HMM file you created with frahmmbuild.  The -o flag will be used to redirect the alignments to tutorial/xxx.aliout and the --tblout flag will be used to gernerate a table listing each hit with all relevant data. 
```bash
   % frahmmer -o tutorial/xxx.aliout --tblout tutorial/xxx.tblout tutorial/xxx.hmm tutorial/seq2.fa 
```

To cheack that the table output is correct run the following:
```bash
   % diff tutorial/xxx.tblout tutorial/xxx.tblout
``` 

**Options** 

Basic Options:<br>
|Flag  | Description |
|:---|:---|
|**-h** &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Help; print a brief reminder of command line usage and all available options|                                               
            
Options directing output:
| Flag | Description |
| :--- | :--- |
| **-o \<f>** | Direct the main human-readable output to a file \<f> instead of the default stdout. |
| **-A \<f>** | Save a multiple alignment of all significant hits (those satisfying "inclusion thresholds") to the file \<f>. |
| **--tblout \<f>** &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Save a simple tabular (space-delimited) file summarizing the per-target output, with one data line per homologous target sequence found. |
| **--fstblout \<f>** | Save a simple tabular (space-delimited) file \<f> of the frameshift locations for each alignment. |
| **--hmmout \<f>** | If queryfile is sequence-based, write the internally-computed HMM(s) to file \<f>. |
| **--acc** | Use accessions instead of names in the main output, where available for profiles and/or sequences. |
| **--noali** | Omit the alignment section from the main output. This can greatly reduce the output volume. |
| **--notrans** | Omit the translated DNA sequence in the alignment output. |
| **--frameline** | Include the frame of each codon in the alignment output. |
| **--notextw** &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Unlimit the length of each line in the main output. The default is a limit of 120 characters per line, which helps in displaying the output cleanly on terminals and in editors, but can truncate target profile description lines.|
|**--textw \<n>**| Set the main output’s line length limit to \<n> characters per line. The default is 120. |
            
Options controlling translation:
| Flag | Description |
| :--- | :--- |   
| **-c \<id>** | Choose alternative genetic code <id> where <id> is the numerical code of one of the NCBI transl_tables. |
| **-l \<n>** | Set the minimum reported ORF length to \<n> aa. |
| **-m**      | Require ORFs to start with an initiator codon AUG (Met). |
| **-M** &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;| Require ORFs to start with an initiator codon, as specified by the allowed initiator codons in the NCBI transl_table. In the default Standard code, AUG, CUG, and UUG are allowed as initiators. An initiation codon is always translated as Met even if it does not normally encode Met as an elongator. |
| **--fs \<x>** &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Set the frameshift probabilty.  \<x> equals the probabilty of a emitting a psuedo-codon with a single nucleotide indel (defaults to 0.01).  For psuedo-codons with a two nucleotide indel the probability will be \<x>/2. The accepted range is 0.001<=\<x><=0.05. |
| **--rosalind** | Only align the top strand. |
| **--franklin** | Only align the bottom strand. |
            
Options controlling reporting thresholds:
| Flag | Description |
| :--- | :--- |   
| **-E \<x>** &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Report target sequences with an E-value of <= \<x>. The default is 10.0, meaning that on average, about 10 false positives will be reported per query, so you can see the top of the noise and decide for yourself if it’s really noise. |
| **-T \<x>**  | Instead of thresholding output on E-value, instead report target sequences with a bit score of >= \<x>. |
            
Options controlling inclusion (significance) thresholds:
| Flag | Description |
| :--- | :--- |   
| **--incE \<x>** &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Use an E-value of <= \<x> as the inclusion threshold. The defaul is 0.01, meaning that on average, about 1 false positive would be expected in every 100 searches with different query sequences. |
| **--incT \<x>** &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Instead of using E-values for setting the inclusion threshold, use a bit score of >= \<x> as the inclusion threshold. By default this option is unset. |
            
Options controlling acceleration heuristics:
| Flag | Description |
| :--- | :--- |   
| **--max** &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;| Maximum sensitivity. Turn off all filters, including the bias filter, and run full Forward/Backward postprocessing on every target. This increases sensitivity slightly, at a large cost in speed. |
| **--F1 \<x>** &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;| Set the P-value threshold for the MSV filter step. The default is 0.02, meaning that roughly 2% of the highest scoring nonhomologous targets are expected to pass the filter. |
| **--F2 \<x>** | Set the P-value threshold for the Viterbi filter step. The default is 0.001. |
| **--F3 \<x>** | Set the P-value threshold for the Forward filter step. The default is 1e-5. |
| **--nobias** &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Turn off the bias filter. This increases sensitivity somewhat, but can come at a high cost in speed, especially if the query has biased residue composition (such as a repetitive sequence region, or if it is a membrane protein with large regions of hydrophobicity). Without the bias filter, too many sequences may pass the filter with biased queries, leading to slower than expected performance as the computationally intensive Forward/Backward algorithms shoulder an abnormally heavy load. |
| **--nonull2** | Turn off the null2 score corrections for biased composition. |
| **--fsonly** | Turn off standard Forward filter, using only the frameshift aware Forward filter and pipeline. |
            
Other expert options:
| Flag | Description |
| :--- | :--- |   
| **-Z \<x>**  | Assert that the total number of targets in your searches is \<x>, for the purposes of per-sequence E-value calculations, rather than the actual number of targets seen. |
| **--seed \<n>** &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Set the random number seed to \<n>. Some steps in postprocessing require Monte Carlo simulation. The default is to use a fixed seed (42), so that results are exactly reproducible. Any other positive integer will give different (but also reproducible) results. A choice of 0 uses a randomly chosen seed. |
| **--w_beta \<x>** &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Window length tail mass. The upper bound, W, on the length at which FraHMMER expects to find an instance of the model is set such that the fraction of all sequences generated by the model with length >= W is less than <x>. The default is 1e-7. This flag may be used to override the value of W established for the model by frahmmbuild, or when the query is sequence based. |
|**--w_length \<n>** &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Override the model instance length upper bound, W, which is otherwise controlled by --w_beta. It should be larger than the model length. The value of W is used deep in the acceleration pipeline, and modest changes are not expected to impact results (though larger values of W do lead to longer run time). This flag may be used to override the value of W established for the model by hmmbuild, or when the query is sequence-based. |
| **--qformat \<s>** &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Assert that input queryfile is a sequence file (unaligned or aligned), in format \<s>, bypassing format autodetection. Common choices for \<s> include: fasta, embl, genbank. Alignment formats also work, and will serve as the basis for automatic creation of a profile HMM used for searching; common choices include: stockholm, a2m, afa, psiblast, clustal, phylip. |
|**--qsingle_seqs** &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Force queryfile to be read as individual sequences, even if it is in an msa format. For example, if the input is in aligned stockholm format, the --qsingle_seqs flag will cause each sequence in that alignment to be used as a seperate query sequence. |
| **--tformat \<s>** &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Assert that target sequence file seqfile is in format \<s>, bypassing format autodetection. Common choices for \<s> include: fasta, embl, genbank. Alignment formats also work; common choices include: stockholm, a2m, afa, psiblast, clustal, phylip. For more information, and for codes for some less common formats, see main documentation. The string \<s> is case-insensitive (fasta or FASTA both work). |
| **--block_length \<n>** | The input sequence is broken into blocks of size \<n> million letters. |
| **--cpu \<n>** &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Set the number of parallel worker threads to \<n>. On multicore machines, the default is 2. You can also control this number by setting an environment variable, HMMER_NCPU. There is also a master thread, so the actual number of threads that frahmmsearch spawns is \<n>+1. This option is not available if FraHMMER was compiled with POSIX threads support turned off. |
   
Available NCBI genetic code tables (for -c \<id>):
| id | Description |
| :--- | :--- |   
| 1 | Standard |
| 2 | Vertebrate mitochondrial |
| 3 | Yeast mitochondrial |
| 4 | Mold, protozoan, coelenterate mitochondrial; Mycoplasma/Spiroplasma |
| 5 | Invertebrate mitochondrial |
| 6 | Ciliate, dasycladacean, Hexamita nuclear |
| 9 | Echinoderm and flatworm mitochondrial |
| 10 | Euplotid nuclear |
| 11 | Bacterial, archaeal; and plant plastid |
| 12 | Alternative yeast |
| 13 | Ascidian mitochondrial |
| 14 | Alternative flatworm mitochondrial |
| 16 | Chlorophycean mitochondrial |
| 21 | Trematode mitochondrial |
| 22 | Scenedesmus obliquus mitochondrial |
| 23 | Thraustochytrium mitochondrial |
| 24 | Pterobranchia mitochondrial |
| 25 | Candidate Division SR1 and Gracilibacteria |

