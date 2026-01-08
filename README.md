# Functional Annotation Guide

**Collection of command line protocols of functional annotation tools in plant genomics**

## Overview

Many tools are available for various tasks in plant genomics. Finding the right tool for a certain purpose can be challenging. This repository provides an overview of recommended tools specifically for functional annotation in plant genomics. This documentation is a guide elucidating step-by-step usage of some widely used functional annotation tools in plant genomics. The tools covered range from user-friendly databases to advanced command-line tools. This will be of use to biologists interested in functional annotation of genes of interest in their own research. Each of the commands and usage instructions here are based on an example file named 'sample' and need to adapted based on the user's own use-case. If the user is not familiar with linux environment and packages for software installation, it is first recommended to read these relevant guidelines specified in https://github.com/PuckerLab/PlantGenomicsGuide.


## Functional annotation database tools

### InterProScan5

Documentation: https://interproscan-docs.readthedocs.io/en/latest/

InterPro [https://doi.org/10.1093/bioinformatics/btu031] is a database integrating predictive information about protein function from a number of partner resources like CATH, CDD, PANTHER, and Pfam to name a few. It is hosted and maintained by EMBL-EBI. InterProScan5 is a software for functional annotation of proteins. It is integrated with the InterPro database and can be installed locally. Step by step instructions for using the tool can be found at https://github.com/PuckerLab/PlantGenomicsGuide 

### Mercator4

Mercator4 [https://doi.org/10.1007/978-1-0716-1609-3_9] [https://doi.org/10.1016/j.molp.2019.01.003] is an online tool for protein annotation of land plants and algae. It assigns protein functions in a hierarchical manner with each sub-node being more specific than the previous node and provides options to include Swiss-Prot annotations. The protein annotations can further be visualized as a tree structure or a heat map. After obtaining the functional annotation mapping file from Mercator4 it is also possible to do an enrichment analysis of genes of your interest using the database. Following is a step by step guide for using Mercator4:

```

- Prepare input FASTA file of protein or cDNA sequences in the correct FASTA format.

- Assess the correctness of FASTA file by uploading it to the Mercator4 FASTA validator.

- Upload the file in the protein annotation section of the database and submit the job. It is advisable to provide your email ID for long-running jobs.

- In case you have a list of genes of interest (GOI) say, differentially expressed genes for example, you could opt for the enrichment analysis as follows:

    - Upload the Mercator4 mapping file.
    - Here you can specify the list of your GOIs and background genes' list.
    - Make sure that the gene IDs you specify in this step and the gene IDs in the Mercator annotation mapping file are the same.
    - You can choose to perform either over-representation or under-representation analysis or both
    - Once you submit the job for the enrichmnet  analysis, you will get the results that can be downloaded as a TSV file
```

### KEGG

Kyoto Encyclopedia of Genes and Genomes (KEGG) [https://doi.org/10.1093/nar/gkae909] is a well-known database resource that helps in understanding biological functions of proteins from a molecular perspective. The annotation step using KEGG involves assigning KEGG Orthology (KO) identifiers to protein coding and RNA genes. The KEGG GENES dataset hosts a collection of genes and proteins from over 10000 complete genomes of cellular organisms and viruses and uses internal annotation tools like BlastKOALA, GhostKOALA, and KofamKOALA to assign the KO identifiers to obtain the functional orthologs. It integrates systemic, genomic, chemical and health information in the database and allows pathway mapping after the initial ortholog finding step using the assigned KO identifiers. This gives a better context-based information of the protein function in the organism and hence provides more cellular and organismal functional insights. KEGG can be used as follows for functional annotation of protein sequences:

```
- Upload your FASTA file of protein sequences to one of the internal annotation tools (BlastKOALA (https://www.kegg.jp/blastkoala/), GhostKOALA (https://www.kegg.jp/ghostkoala/), or KofamKOALA (https://www.genome.jp/tools/kofamkoala/)).

- Select your reference gene set from eukaryotes, prokaryotes or viruses.

- Provide your official email ID and request email confirmation.

- Click the link sent to your email  ID within 24 hours to submit your job.

- It is important to note that only one job canbe run at a time.

- After the annotation go to the KEGG mapper tool - 'Reconstruct' (https://www.kegg.jp/kegg/mapper/reconstruct.html) for pathway mapping.

- Upload the annotation results from the previous step into this tool and press exec. This file is a two column tab or space-separated file where the first column has the gene IDs of the query and the second column has the corresponding KO identifiers.

```


### BLAST2GO


## Ortholog finding tools

### OrthoFinder2

### SHOOT

### FASTOMA

## Sequence similarity based tools (aligners)

### BLAST

### DIAMOND

### MMSeqs

### MAFFT

### MUSCLE5

## Phylogeny and visualization tools

### IQ-TREE3

### FastTree

### iTOL

## Synteny analysis tools

### JCVI/ MCScan

### TBtools-II

### TOGA

### Orthorefine

### SOI

## Protein structure-based annotation tools

### Dali

### Foldseek

## Expression analysis-based annotation tools

### GENIE3

### WGCNA

## Annotation tools employing combined approaches

### CoExpPhylo

### KIPEs

### DupyliCate

## References



If you have questions about plant genomics that were not answered by any of these resources, please feel free to get in touch with the Plant Biotechnology and Biotechnology research group at the University of Bonn.
