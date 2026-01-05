# Basic ISCro4 target scoring (BITS) 
BITS is a simple tool for analyzing next-generation sequencing (NGS) data from CUT&RUN experiments on ISCro4 recombinase to nominate and score potential off-target sites.

# publication
BITS is associated with the following publication:

>**Programmable genome editing in human cells using RNA-guided bridge recombinases**
>Pelea et al., *Science* (2026)

Please see the [Schwank Lab GitHub](https://github.com/Schwank-Lab/ISCro4_Pelea_2026) for other code and analyses associated with this manuscript.

# usage
BITS is intended to work with genomic coverage data resulting from NGS alignment.

1. From raw NGS data (e.g., Read1/2 FASTQ datafiles), perform alignment, deduplication, and coverage calculations as described in the publication Methods.

2. Thereafter, take the resulting bamCoverage ".bw" bigWig files and deposit them in:
"./_step1_precompute/_input/". 

3. Run the Rnotebook "./step1_combineCasoffinderBigwigs.Rmd". This generates ".Rds.gz" datafiles that will be deposited in "./_step1_precompute/_output/".

4. Modify the analysis manifest file "./_step2_getHits/manifest.csv", if necessary. A pre-filled example is already included. The 3 columns in order are: 
	(i) Step1 output for the untransfected background
	(ii) Step1 output for the transfected sample
	(iii) gRNA sequence to use for analysis. This should be the same as used in CasOFFinder ("./_resources/casoffinder_14bp_sequences_4mm.out.gz")

5. Run the Rnotebook "./step2_findHits.Rmd". The results will be deposited in "./_step2_getHits/_output".

# extra considerations
BITS is intended only for initial nomination of ISCro4 target sites which should be further validated by orthogonal methods, such as targeted amplicon sequencing. BITS only evalautes targets with 0-4 mismatches to the target sequence; targets with >4 mismatches may exist. Both false positives and false negatives may be possible. 

BITS was created by Charles D Yeh and Lilly van de Venn. 

