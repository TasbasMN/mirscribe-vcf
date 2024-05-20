Download date: 14.05.2024

﻿# IntOGen RELEASE 20230531
# Headers of the driver files 
# For further information about methods please visit the documentation website: https://intogen.readthedocs.io/en/latest/index.html
# To access the source code of intogen please visit: https://bitbucket.org/intogen/intogen-plus/src/master/


Compendium_Cancer_Genes.tsv

1   SYMBOL: Hugo Symbol of the gene. 
2   TRANSCRIPT: ENSEMBL TRANSCRIPT identifier. 
3   COHORT: Name of the cohort where the gene has been detected as driver. 
4   CANCER_TYPE: Tumor type of the cohort. 
5   METHODS: Methods that detect this gene as driver (q-value <0.1). When the method is combination it means that the combination of p-values renders the gene as a driver. 
6   MUTATIONS: Number of somatic mutations (including short indels) of the gene in the cohort. 
7   SAMPLES: Number of mutated samples of the gene in the cohort. 
8   QVALUE_COMBINATION: Significance of the combined output. 
9   ROLE: Consensus mode of action of the gene. Derived from intOGen and from literature.  LoF (Loss of Function), Act (Activating) or amb (Ambiguous).  
10  CGC_GENE: Gene annotated in the Cancer Gene Census. 
11	CGC_CANCER_GENE: Gene annotated in the Cancer Gene Census as driver in the CANCER_TYPE of the COHORT. 	
12  DOMAINS: Comma-separated domains that are considered significant by smRegions (q-value <0.1).  The format is PFAM_ID:START_AA:END_AA. 
13  2D_CLUSTERS: Comma-separated linear clusters that are considered significant by OncodriveCLUSTL (p-value <0.05).  The format is START_AA:END_AA. 
14  3D_CLUSTERS: Comma-separated 3D clusters that are considered significant by HOTMAPS (q-value <0.05).  The format is AA_1, AA_2, AA_3, etc.
15  EXCESS_MIS: dNdScv missense excess rate.  
16  EXCESS_NON:  dNdScv nonsense excess rate. 
17  EXCESS_SPL:  dNdScv splicing excess rate.


