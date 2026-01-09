# Functional Annotation Guide

**Collection of command line protocols of functional annotation tools in plant genomics**

## Overview

Many tools are available for various tasks in plant genomics. Finding the right tool for a certain purpose can be challenging. This repository provides an overview of recommended tools specifically for functional annotation in plant genomics. This documentation is a guide elucidating step-by-step usage of some widely used functional annotation tools in plant genomics. It is important to note that the installation method outlined for each of the tools is based on the ease-of-installation for users of all levels, and the possibility of errors and bugs that can be encountered while doing so. Hence, the installation method detailed can differ from the officially recommended installation methods in the respective tools' documentation pages. Further, the tools covered range from user-friendly databases to advanced command-line tools. This will be of use to biologists interested in functional annotation of genes of interest in their own research. Each of the commands and usage instructions here are based on an example file named 'sample' and need to adapted based on the user's own use-case. If the user is not familiar with linux environment and packages for software installation, it is first recommended to read these relevant guidelines specified in https://github.com/PuckerLab/PlantGenomicsGuide. 


## Functional annotation database tools

### 1. InterProScan5

Documentation: https://interproscan-docs.readthedocs.io/en/latest/

InterPro [https://doi.org/10.1093/bioinformatics/btu031] is a database integrating predictive information about protein function from a number of partner resources like CATH, CDD, PANTHER, and Pfam to name a few. It is hosted and maintained by EMBL-EBI. InterProScan5 is a software for functional annotation of proteins. It is integrated with the InterPro database and can be installed locally. Step by step instructions for using the tool can be found at https://github.com/PuckerLab/PlantGenomicsGuide 

### 2. Mercator4

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

### 3. KEGG

Kyoto Encyclopedia of Genes and Genomes (KEGG) [https://doi.org/10.1093/nar/gkae909] is a well-known database resource that helps in understanding biological functions of proteins from a molecular perspective. The annotation step using KEGG involves assigning KEGG Orthology (KO) identifiers to protein coding and RNA genes. The KEGG GENES dataset hosts a collection of genes and proteins from over 10000 complete genomes of cellular organisms and viruses and uses internal annotation tools like BlastKOALA, GhostKOALA, and KofamKOALA to assign the KO identifiers to obtain the functional orthologs. It integrates systemic, genomic, chemical and health information in the database and allows pathway mapping after the initial ortholog finding step using the assigned KO identifiers. This gives a better context-based information of the protein function in the organism and hence provides more cellular and organismal functional insights. KEGG can be used as follows for functional annotation of protein sequences:

```
- Upload your FASTA file of protein sequences to one of the internal annotation tools (BlastKOALA (https://www.kegg.jp/blastkoala/), GhostKOALA (https://www.kegg.jp/ghostkoala/), or KofamKOALA (https://www.genome.jp/tools/kofamkoala/)).

- Select your reference gene set from eukaryotes, prokaryotes or viruses.

- Provide your official email ID and request email confirmation.

- Click the link sent to your email  ID within 24 hours to submit your job.

- It is important to note that only one job can be run at a time.

- After the annotation go to the KEGG mapper tool - 'Reconstruct' (https://www.kegg.jp/kegg/mapper/reconstruct.html) for pathway mapping.

- Upload the annotation results from the previous step into this tool and press exec. This file is a two column tab or space-separated file where the first column has the gene IDs of the query and the second column has the corresponding KO identifiers.

```

### 4. Blast2GO

Blast2GO [https://doi.org/10.1155/2008/619832] is a user-friendly application for functional annotation. It is a suite of bioinformatics tools like BLAST and a number of databases like KEGG and InterProScan integrated with GO mapping and function assignment from these databases. It provides statistical summaries and nice features for visualization of the functional annotation results. However, only the basic version of the application is free for non-profit academic research purposes. The guideline below is a brief overview on using the basic version of Blast2GO for a simple functional annotation process.

```
- Register for the Blast2GO basic version at https://www.biobam.com/blast2go-basic/ You will receive the Blast2GO activation key via email.

- Download Blast2GO here https://www.biobam.com/blast2go-previous-versions/

- Select the executable depending on your operating system like Windows, Linux or Mac OS.

- After installation, you will be prompted to enter the activation key. Enter it and restart the application.

- Blast2GO provides options to run BLAST locally or online and interfaces with a number of databases like the nr, PIR, and Swiss-Prot to name a few.

- There is also an option to directly load the BLAST result files into the application.

- You can run a GO mapping where it interfaces with the GO database and maps your query sequences wwith relevant GO terms.

- Next, you can also opt to run InterProScan to obtain the functional annotation for your sequences or directly use the run annotation option in the 'annot' menu.

- After annotation, the analysis tab can be clicked to perform enrichment analysis and KEGG pathway mapping.

```

## Ortholog finding tools

### 5. OrthoFinder

OrthoFinder [https://doi.org/10.1101/2025.07.15.664860] is a widely used tool for ortholog identification. It identifies orthologous genes, and clusters them into orthogroups. These orthogroups are further analyzed to obtain gene and species trees. Gene duplication events are also reported using both these trees. The method is scalable to a large number of genomes. It helps in large scale comparative genomics across species and provides detailed statistics and summary of the results. 

**Installation:**

```
#clone the git repository in the destination folder of your choice

git clone https://github.com/OrthoFinder/OrthoFinder.git

# move to the OrthoFinder folder

cd OrthoFinder

# Create an virtual environment named of3_env

python3 -m venv of3_env

# Activate of3_env

. of3_env/bin/activate

# Insall the required dependencies with pip

pip install .

# Test your installation

# Print out help information

orthofinder --help

# Check the version

orthofinder --version

# Do a test run of OrthoFinder from the OrthoFinder folder

orthofinder -f ExampleData

If OrthoFinder and its dependencies were successfully installed, it will complete and display the time taken for the run without any errors.

```

**Sample OrthoFinder command:**

```
orthofinder -A mafft -T iqtree -f Sample_dir

```
In the above command the mandatory parameter is only the '-f' flag that takes in a directory of FASTA files for all the species to be analyzed. Each FASTA file must contain protein sequences. The flags -A and -T have famsa and fasttree respectively as defaults. They have been tweaked in the above sample command to show that the command needs to be tailored according to specific-use-cases. For example, iqtree takes longer to run than fasttree, but can provide more accurate results. So, such considerations must be done before constructing a command for a tool run.

**After completing the OrthoFinder run:**

```

# deactivate the virtual environment

deactivate

# Reactivate the environment when you want to use the tool again

#Move to the OrthoFinder folder

cd /path/to/OrthoFinder

#activate the environment from within the OrthoFinder folder

. of3_env/bin/activate

```

**Official documentation:** https://github.com/OrthoFinder/OrthoFinder


### 6. SHOOT

SHOOT [https://doi.org/10.1186/s13059-022-02652-8] is a phylogenetic search engine. It is an alternative to BLAST and an orthogonal approach to finding orthologs. Instead of searching for orthologs using sequence similarity, SHOOT directly relies on phylogenetic signal. It takes the query protein sequence, and searches it against a database of phylogenetic trees. It then places the sequence or more correctly grafts the sequence onto its suitable phylogenetic tree and returns the results along with the orthologs found for the specific sequence. 

**Using SHOOT:**

```
- Visit the SHOOT search engine page https://www.shoot.bio/

- Paste your query protein's amino acid sequence.

- Select a suitable phylogenetic database from all the domains of life like plants.

- You can also make the run more customized by tweaking the advanced options parameters like DIAMOND sensitivity.

- Press SHOOT to submit your job.

```

**Official documentation:** https://github.com/davidemms/SHOOT

### 7. FASTOMA
FASTOMA [https://doi.org/10.1038/s41592-024-02552-8] is another ortholog finding tool that combines ultrafast homology grouping with taxonomy-based sampling. It is an advanced update of the OMA algorithm and adopts a highly efficient parallel computing approach. It can take in multiple proteomes of species and also helps infer hierarchical orthologous groups (HOGs). Apart from the protein FASTA file for each of the species, it needs species tree in Newick format for all the species used in a particular run. 

**Installing pre-requisites:**

- FastOMA is a nextflow-based pipleine. nextflow needs Java 17 or above to be installed.

- Java 17 can be installed with
  ```
  sudo apt install openjdk-17-jdk
  ```
- Next nextflow needs to be installed for FastOMA to work and the steps to install it are as follows:
 
```
# move to the desired destination folder

cd /path/to/destination/folder

# download nextflow
    
curl -s https://get.nextflow.io | bash

# make nextflow executable

chmod +x executable

# If you are working on a UNIX-based system like Linux, update the bashrc file with the path of the destination folder containing the nextflow executable as follows:

nano ~/.bashrc

# add the following line at the end of the bashrc file

export PATH=$PATH:/path/to/destination/folder

# press ctrl+x and then when prompted y

# reload the bashrc file

source ~/.bashrc

```

**FastOMA installation**

```
# clone the FastOMA GitHub repository

git clone https://github.com/DessimozLab/FastOMA.git

# move into the FastOMA folder

cd FastOMA

# Make sure that docker is installed

For docker installation instructions please refer the official docker installation guide https://docs.docker.com/engine/install/

```
**Running FastOMA:**
```

# Test run of FastOMA using the nextflow docker profile

nextflow run FastOMA.nf -profile docker \
    --container_version "sha-$(git rev-list --max-count=1 --abbrev-commit HEAD)" \
    --input testdata/in_folder \
    --output_folder /path/to/output/folder/sample_fastoma_test

The in_folder in the above command has the protein FASTA files inside a folder called proteome. It also has the species tree file in the Newick format. Please note that the protein files must have the extension .fa

```

**Official documentation:** https://github.com/DessimozLab/FastOMA

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
