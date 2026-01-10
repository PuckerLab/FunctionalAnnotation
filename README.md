# Functional Annotation Guide

**Collection of command line protocols of functional annotation tools in plant genomics**

## Overview

Many tools are available for various tasks in plant genomics. Finding the right tool for a certain purpose can be challenging. This repository provides an overview of recommended tools specifically for functional annotation in plant genomics. This documentation is a guide elucidating step-by-step usage of some widely used functional annotation tools in plant genomics. It is important to note that the installation method outlined for each of the tools is based on the ease-of-installation for users of all levels, and the possibility of errors and bugs that can be encountered while doing so. Hence, the installation method detailed can differ from the officially recommended installation methods in the respective tools' documentation pages. Further, the tools covered range from user-friendly databases to advanced command-line tools. This will be of use to biologists interested in functional annotation of genes of interest in their own research. Each of the commands and usage instructions here are based on an example file named 'sample' and need to adapted based on the user's own use-case. If the user is not familiar with linux environment and packages for software installation, it is first recommended to read these relevant guidelines specified in https://github.com/PuckerLab/PlantGenomicsGuide. 


## Functional annotation database tools

### 1. InterProScan5

Documentation: https://interproscan-docs.readthedocs.io/en/latest/

InterPro [https://doi.org/10.1093/bioinformatics/btu031] is a database integrating predictive information about protein function from a number of partner resources like CATH, CDD, PANTHER, and Pfam to name a few. It is hosted and maintained by EMBL-EBI. InterProScan5 is a software for functional annotation of proteins. It is integrated with the InterPro database and can be installed locally. Step by step instructions for using the tool can be found at https://github.com/PuckerLab/PlantGenomicsGuide 

**Official documentation:** https://interproscan-docs.readthedocs.io/en/v5/

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

**Official documentation:** https://www.plabipd.de/mercator_main.html

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
**Official documentation:** https://www.genome.jp/kegg/

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

**Running OrthoFinder:**

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

### 8. BLAST

Basic Local Alignment Search Tool (BLAST) [https://doi.org/10.1016/S0022-2836(05)80360-2] is a widely used alignment tool. It compares the query sequences against a database of target sequences and can take in both nucleotide and protein sequences. It performs local alignment, meaning, it looks for highly similar regions between sequences, and does not perform an end-to-end alignment of the sequences. BLAST is available as a server as well as a standalone software. The standalone software supports analysis of a large number of sequences and is recommended to be used. 

**BLAST installation:**

```

# Move to the folder of your choice

cd /path/to/folder

# Fetch the latest BLAST tar ball

wget https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.17.0+-x64-linux.tar.gz

# Decompress the tar ball

tar -xzvf ncbi-blast-2.17.0+-x64-linux.tar.gz

# Move to the bin folder of the BLAST installation

cd ncbi-blast-2.17.0+/bin

```

**Making BLAST database:**

```
makeblastdb -in sampleA.pep.fasta -dbtype prot -parse_seqids -out sampleA_database

- In the above command, a BLAST database is created from the sampleA's peptide FASTA file.

```

**Running BLAST:**

```

./blastp -query sampleB.pep.fasta -db sampleA_database -evalue 1e-5 -out sampleBA.tsv -outfmt 6 -num_threads 10

- In the above command, the parameters for -evalue flag, the output format flag -outfmt and the -num_threads flag that takes in the number of cores for the analysis, can be tweaked based on the specific analysis at hand and the resources you have.

- There are also other available executables like blastx and blastn. More parameters of blastp and the other executables can be obtained by using the -h flag along with the specific executable's name.

```

**Official documentation:** https://blast.ncbi.nlm.nih.gov/doc/blast-help/

### 9. DIAMOND

DIAMOND [https://doi.org/10.1038/s41592-021-01101-x] is another aligner similar to BLAST. But it is 100x-10000x faster than BLAST. It facilitates both protein and translated DNA searches. It offers options to decide between the different levels of sensitivity, and thereby helps optimize the computational costs and accuracy. With its ability to cluster billions of proteins and low resource consumption, it is becoming a widely used tool both standalone, as well in integration with other bioinformatics tools and pipelines.

**DIAMOND installation:**

```
# Move to the folder of your choice

cd /path/to/folder

# Fetch the latest DIAMOND tar ball

wget http://github.com/bbuchfink/diamond/releases/download/v2.1.18/diamond-linux64.tar.gz

# Decompress the tar ball

tar xzf diamond-linux64.tar.gz
```

**Creating a DIAMOND database:**

```
./diamond makedb --in sampleA.fasta -d sampleA_database
```

**Running DIAMOND:**

```
./diamond blastp -d sampleA_database -q sampleB.fasta -o sampleBA.tsv

- To know about the other executables in diamond, the help flag can be used like this below

./diamond help
```

**Official documentation:** https://github.com/bbuchfink/diamond

### 11. MAFFT

MAFFT [https://doi.org/10.1093/nar/gkf436] [https://doi.org/10.1093/molbev/mst010] is a popular alignment tool for multiple sequence alignment (MSA). Being based on a fast Fourier transform method, it is a fast tool forMSA. It offers a number of modes and algorithms for MSA that can be tweaked according to the objective in hand, and the nature of the sequences. MAFFT is availbale as an online server as well as a standalone tool, but the latter is recommended for its ability to process a large number of sequences.

**MAFFT installation:**

```

# Fetch the mafft Debian package

wget https://mafft.cbrc.jp/alignment/software/mafft_7.526-1_amd64.deb

# Install the Debian package

sudo dpkg -i mafft_7.526-1_amd64.deb

```

**Running MAFFT:**

```
mafft --globalpair --maxiterate 1000 --amino sample.fasta > sample.aln

- The above sample command is based on the MAFFT G-INS-i method recommended when global alignment is needed using the Needleman-Wunsch algorithm.

- Details about other MAFFT algorithms can be found in https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html

- Tips for using the different MAFFT parameters and algorithms can be found in https://mafft.cbrc.jp/alignment/software/tips0.html

```
**MAFFT server and official documentation:** https://mafft.cbrc.jp/alignment/server/

### 12. Muscle5

Muscle5 [https://doi.org/10.1038/s41467-022-34630-w] is another orthogonal approach for multiple sequence alignment. Muscle5 is a major update of the original Muscle tool, and is more accurate, faster and much more scalable. It is reported to be approximately 20-30% more accurate than MAFFT. It is also available as an online server and a standalone software.

**Muscle5 installation:**

```
# Move to the folder of your choice

cd /path/to/folder

# Fetch the MUSCLE5 executable

wget https://github.com/rcedgar/muscle/releases/download/v5.3/muscle-linux-x86.v5.3

# Enable execute permission

chmod +x muscle-linux-x86.v5.3

# For ease of use, rename the executable as muscle

mv muscle-linux-x86.v5.3 muscle
```

**Running Muscle5:**
```
muscle -super5 sample.fa -output sample_aln.afa
```

**Muscle5 server:** https://www.ebi.ac.uk/jdispatcher/msa/muscle5

**Official documentation:** https://www.drive5.com/muscle5/

### 13. algntrim.py

After obtaining a multiple sequence alignment file, it is important to clean the file. This is because MSA columns can have a large number of gaps that turn into noise when the MSA file is used for a subsequnt phylogenetic analysis. This cleaning or trimming step helps decrease this noise, reduces computational costs of a phylogenetic tree building, and improves the signal from the sequences for phylogeny. Occupancy is the quantitative parameter used to determine the gaps present in a column. It is defined as the percentage of sequences in an aligned column that do not have gaps. Although tools like pxclsq from PHYX are used for this trimming step, here we recommend using an in-house developed Python script algntrim.py for its ease of use.

**Fetching and running algntrim.py:**
```
# Fetch the script from this GitHub repository in a folder of your choice

wget https://raw.githubusercontent.com/PuckerLab/FunctionalAnnotation/refs/heads/main/algntrim.py

# Run algntrim.py

python3 algntrim.py --in sample.aln --out sample_trimmed.aln
```

**Application note:** The trimmed MSA file obtained from this step can be used for phylogenetic tree building in the next subsequent steps.

## Phylogeny analysis tools

### 14. IQ-TREE3

IQ-TREE3 [https://ecoevorxiv.org/repository/view/8916/] is a well known phylogenetic tree building tool. It integrates ModelFinder, an algorithm for automatic selection of substitution model for the tree building, offering results with improved accuracy. Along with this, it offers a number of features like ultrafast bootstrap that performs bootstrapping comparatively faster than standard bootstrap without compromising the accuracy to a large extent. It also has options to perform topology testing of the phylogenetic tree using branch tests like SH-aLRT. Containing a wide variety of phylogenetic models, IQ-TREE3 is a highly recommended software for obtaining phylogenetic trees with good balance between accuracy and computtional times.

**IQ-TREE3 installation:**

```
# Move to the folder of your choice

cd /path/to/folder

# Fetch the IQ-TREE3 tar ball

wget https://github.com/iqtree/iqtree3/releases/download/v3.0.1/iqtree-3.0.1-Linux.tar.gz

# Decompress the IQ-TREE3 tar ball

tar -xvzf iqtree-3.0.1-Linux.tar.gz

# Move into the bin folder of the IQ-TREE3 installation

cd iqtree-3.0.1-Linux/bin
```

**Running IQTREE3:**

```
./iqtree3 -s sample_trimmed.aln -m MFP -wsr --alrt 1000 -B 1000 -pre /path/to/folder/sample -T 30

- The different parameters in the above command are explained below:

    m - MFP -> ModelFinder Plus (https://iqtree.github.io/doc/Substitution-Models)
    wsr -> write site rates
    alrt -> number of bootstrap replicates for SH-alrt test
    Shimodaira-Hasegawa approximate Likelihood Ratio Test -> how confident should you be about each branch in your tree
    B -> Ultrafast bootstrap
    pre -> full path to output folder with the name of your treefile
    T -> Number of cores for the phylogeny analysis
```

**Official documentation:** https://iqtree.github.io/doc/Home#documentation

### 15. FastTree

FastTree [https://doi.org/10.1371/journal.pone.0009490], as the name suggests is a tree building tool that can take in alignments with millions of sequences and perform tree building with a reasonable amount of memory and time. It is computationally fast and can take in both protein and nucleotide sequences. The accuracy of the tool is slightly traded for its computational efficiency. Hence, it is important to understand the goal of the tool, its merits and apply it in appropriate use cases.

**FastTree installation:**

```
# Move into a folder of your choice

cd /path/to/folder

# Fetch the FastTree executable

wget http://www.microbesonline.org/fasttree/FastTree

```

**Running FastTree:**

```
./FastTree -nopr -wag sample_trimmed.aln > sample.nwk

- The -nopr flag makes the program assume same evolutionary rate for all sites, reducing the computational time but compromising on the accuracy

- The -wag is a model suited for protein sequences. FastTree offers other models as well that can be found by using the ./FastTree -h command
```

**Official documentation:** https://morgannprice.github.io/fasttree/

## Synteny analysis tools

### JCVI/ MCScan

JCVI [https://doi.org/10.1002/imt2.211] is a versatile Python-based library. It offers a number of useful tools for analysing, wrangling genomic files, and for performing various aspects of genome annotation. It also facilitates comparative genomic studies across multiple genomes using tools like MCScan. MCScan is a specific utility of JCVI that is focussed on finding regions of synteny between genomes. Synteny is defined as conserved order of gene blocks between genomes. It helps obtain evolutionary insights about genomes and provides positional context in ortholog finding, making it more reliable. JCVI/ MCScan offers very good features to obtain the micro and macro-synteny plots and uses more sophisticated criteria and approach than BLAST for finding orthologs. 

**JCVI installation:**

```
# fetch JCVI toolkit from bioconda

conda create -p /sample/destination/folder/jcvi jcvi -c bioconda -c conda-forge

# activate the conda environment

conda activate /sample/destination/folder/jcvi

```

**Preparing input files for JCVI/ MCScan:**

- GFF files giving positional information and protein or coding sequence (CDS) FASTA files of the species to be analyzed for synteny are required for this analysis. We will take two example species - sampleA and sampleB as the species for analysing pairwise synteny in the example detailed below.

```
# convert the input GFF files to BED files

python -m jcvi.formats.gff bed --type=mRNA --key=Name --primary_only sampleA.gff3.gz -o sampleA.bed
python -m jcvi.formats.gff bed --type=mRNA --key=Name --primary_only sampleB.gff3.gz -o sampleB.bed

The parameters for the --type and --key flags in the above command must be adjusted according to the fields in your specific GFF file. If your species has many isoforms, it is recommended to use the --primary_only flag that retains only the primary transcript per gene.

# Reformat the input FASTA files

python -m jcvi.formats.fasta format sampleA.cds.fa.gz sampleA.cds
python -m jcvi.formats.fasta format sampleB.fa.gz sampleB.cds

```

**Synteny analysis with JCVI/ MCScan:**

```
# Perform the pairwise synteny search

python -m jcvi.compara.catalog ortholog sampleA sampleB --no_strip_names

- The --no_strip_names flag in the above command helps retain the original gene identifiers found in the input FASTA files in the MCScan output.

- This step produces the LAST file and anchors file. The LAST file is very similar to the BLAST output file and the anchors file is a list of high quality synteny blocks. This anchors file cna be used in subsequent steps for visualizing synteny at macro and micro scales.

```

**Synteny visualization with JCVI/ MCScan:**

```
# Visualizing synteny between two genomes using a dot plot

python -m jcvi.graphics.dotplot sampleA.sampleB.anchors

This step produces a dot plot in the PDF format and helps infer the genome-wide synteny pattern between the genomes.

# Macrosynteny visualization

- Karyotype plots help visualize synteny between specified chromosomes of the genomes being compared.

-First step for this is to produce a more concise form of the anchors file called the anchors.simple file

# Make a succint version of the anchors file

python -m jcvi.compara.synteny screen --minspan=30 --simple sampleA.sampleB.anchors sampleA.sampleB.anchors.new

The parameter of the --minspan flag in the above command specifies the minimum number of genes spanning a synteny block for the block to be retained. This step produces a .anchors.simple file.

- Next, two more input files need to be prepared - first the seqids file and then the layout file.

- seqids file is a simple TXT file named seqids. It should have the chromosome identifiers of sampleA in first row and sampleB in second row, each separated by commas, as shown below

chrA1,chrA2,chrA3,chrA4,chrA5
chrB1,chrB2,chrB3,chrB4

- layout file is a simple TXT design file named layout that specifies positional and style information of the plot. It is divided into two sections - the upper section and the lower section. The upper section has eight columns - y axis position, x axis start, x axis end, degree of rotation, colour (colour needs to be specified in hex code like #fc8d62), labels to be displayed, vertical alignment of the chromosomes per species, and the name of the BED file to be used. The lower section specifies the boundaries to draw the edges based on the information from the .anchors.simple file. A sample layout file (adopted from the official MCScan documentation) is shown below:

# y, xstart, xend, rotation, color, label, va,  bed
 .6,     .1,    .8,       0,      , sampleA, top, sampleA.bed
 .4,     .1,    .8,       0,      , sampleB, top, sampleB.bed
# edges
e, 0, 1, sampleA.sampleB.anchors.simple

# Plot the karyotype plot

python -m jcvi.graphics.karyotype seqids layout

# Microsynteny visualization

- JCVI/ MCScan offers great flexibility and wide range of synteny visualization options. Microsynteny visualization helps view the local synteny at the gene level of your desired region in the genome.

- After the pairwise synteny step, for a microsyteny visualization, a blocks file needs to be obtained.

# Obtain blocks file for miscrosynteny

python -m jcvi.compara.synteny mcscan sampleA.sampleB sampleA.sampleB.lifted.anchors --iter=1 -o sampleA.sampleB.i1.blocks

In the above command, the parameter of the --iter flag specifies the number of best hits to be extracted per gene from the LAST results to make the blocks file and can be tweaked according to your use case.

# Extract specific region of interest from the full blocks file

Let's say you want to visualize the first 50 genes in the blocks file, then you need to extract this region into another blocks file that will be used for the visualization.

head -50 sampleA.sampleB.i1.blocks > blocks

# Prepare the blocks layout file

The microsynteny visualization does not need the seqids file but it needs the TXT layout file similar to the one required for macrosynteny visualization. Let us name it blocks.layout, since this is the design file describing the synteny plot for the blocks file we extracted in the previous step. Everthing remains the same in this file like the macrosynteny lyout file except that, since no .anchors.simple file is not produced in microsynteny visualization, the lower section of the layout file does not have the name of the .anchors.simple file.

# Get a microsynteny plot

python -m jcvi.formats.bed merge sampleA.bed sampleB.bed -o sampleA_sampleB.bed

python -m jcvi.graphics.synteny blocks sampleA_sampleB.bed blocks.layout

# deactivate the conda environment after use

conda deactivate

- More styling and customization options for the synteny plots can be obtained in the official documentation page of JCVI/ MCScan.
```
**Official documentation:** https://github.com/tanghaibao/jcvi/wiki/                                                                                          

### SOI

SOI is another Python-based toolkit that helps in finding syntenic orthologs. It employs a method called Orthology Index (OI) using which it infers the proportion of pre-inferred orthologs within a syntenic block. It needs synteny results and ortholog results as mandatory input files and filters out orthologs in the syntenic region with a default OI value of 0.6. It accepts synteny outputs from JCVI/ MCScan, MCScanX and WGDI and ortholog outputs from OrthoFinder, and OrthoMCL. Apart from providing the filtered syntenic orthologs as results, it can also be used for other evolutionary analyses like obtaining a dotplot and clustering the syntenic orthologs into syntenic orthogroups to name a few.

**SOI Installation:**

```

# Create the conda environment in a folder of your choice and name it OrthoIndex

conda create -p /path/to/folder/OrthoIndex

# Install the SOI tool using conda

conda install -p /path/to/folder/OrthoIndex -c conda-forge -c bioconda soi

# Activate the conda environment

conda activate /path/to/folder/OrthoIndex

# Test SOI installation

soi -h

- The above command must display the usage instructions for SOI

```

**Running SOI:**

```
soi filter -s sample.collinearity.gz -o OrthoFinder/Sample_results -c 0.8 > sample.collinearity.ortho.results

The parameter of the -c flag is 0.6 by default and can be adjusted according to the strictness level of orthologs needed.

# Deactivate the conda environment

conda deactivate
```
**Official documentation:** https://github.com/zhangrengang/SOI

## Protein structure-based annotation tools

### Dali

#### Dali server

Dali is a popular webserver [https://doi.org/10.1093/nar/gkac387] used for finding protein homologs based on protein structure search against databases like the AlphaFoldDB and PDB. It offers three main options - (i) Comparing a query protein structure against the PDB, PDB25 or the AlphaFoldDB databases, (ii) Pairwise structure comparison amongst the user given list of protein structures that allows a maximum of 10 structures per job, (iii) All against all protein structure comparison amongst the user given list of prorein structures that allows a maximum of 64 structures per job. For all these options, the input needs to be the PDB identifier of the query protein(s) along with the chain identifier. PDB identifiers based on key word search can be accessed at http://www.rcsb.org/. The chain identifier denotes the chains in a protein structure and must be given along with the PDB identifier while submitting your job on the Dali server. While Dali is versatile, the web server is limited by the number of protein structures that can be analyzed concurrently. 

**Official documentation:** http://ekhidna2.biocenter.helsinki.fi/dali/DaliTutorial.pdf

#### DaliLite.v5 standalone software

To overcome the input limitations of the Dali server, DaliLite.v5 [https://doi.org/10.1093/bioinformatics/btz536] a standalone software is available.

**DaliLite.v5 installation:**

```

# move to your desired destination folder

cd /path/to/folder

# fetch the DaliLite tar file

wget http://ekhidna2.biocenter.helsinki.fi/dali/DaliLite.v5.tar.gz

# Decompress the tar ball

tar -xvzf DaliLite.v5.tar.gz

# move to the DaliLite bin folder

cd /path/to/folder/DaliLite.v5/bin

# compile the binaries

make clean

make

```

**Running Dali:**

```

- The DaliLite software needs the PDB file format to be converted into an internal format that DaliLite accepts. This is achieved with the import.pl script. The user's private protein file can be used when going for pairwise or all against all searches. But in case of database search, the public database structure fields must be mirrored and converted to the required format using the import.pl script.

# Import private protein structure files

cd /path/to/folder/DaliLite.v5/bin

./import.pl --pdbfile sample.pdb --pdbid samp --dat /path/to/DATA --clean

In the above command the PDB ID must alway be of 4 letters as it is hard coded. --dat flag specifies the output folder for the reformatted file. All structures for a comparison must be provided in a single directory.

# Import public database structure files 

./import.pl --rsync --pdbmirrordir /path/to/folder/pdb --dat /path/to/folder/DATA --clean

The PDB structures will be stored in the location specified to the --pdbmirrordir flag and the Dali formatted files will be stored in the path specified to the --dat flag.

# Make a BLAST database for the structure based database search

- BLAST is used to group together sequences that show sequence level similarity before the structural alignment step

 The following commands can be used to extract sequences from the imported structures into a FASTA file:

# Obtain a list of protein structures in the Dali format from the previous step(s)

ls /path/to/folder/DATA | perl -pe 's/\.dat//' > pdb.list

# Obtain FASTA sequences from the protein structure files

./dat2fasta.pl /path/to/folder/DATA < pdb.list | awk -v RS=">" -v FS="\n" -v ORS="" ' { if ($2) print ">"$0 } ' > pdb.fasta

# Create a BLAST database

makeblastdb -in pdb.fasta -out /path/to/folder/pdb.blast -dbtype prot

# Pairwise Dali alignment

./dali.pl --query sample.list --db target.list

# All-against-all Dali alignment

./dali.pl --matrix --query sample.list

# Database search

./dali.pl --hierarchical --repset pdb25.list --query sample.list --db pdb.list

```
**Official documentation:** http://ekhidna2.biocenter.helsinki.fi/dali/README.v5.html

### Foldseek

Foldseek [https://doi.org/10.1038/s41587-023-01773-0] is a fast structural alignment tool that is available as a server as well as a standalone tool. It is capable of ultra-sensitive searches using protein sequences without the need for structure by leveraging language models. It relies on something called a 'structural alphabet' which is a descriptor of protein tertiary interactions and can be thought of as 3D interaction alphabet that can be used as a proxy for protein structural information while drastically reducing the computation times compared to the other protein structural alignment softwares.

#### Foldseek server:

Foldseek is available as a web server. It can be used to search a protein structure or FASTA sequence against the AlphaFoldDB and PDB in the order of seconds. Under the databases and search settings tab it provides options to select the databases to be searched against along with the mode that needs to be used for the structural alignment process. There is also an option to choose an iterative search process. Further, it offers more customization by providing option to choose the taxonomy class to filter the search to only those groups, helping optimize the computational time and the number of false positive results.

**Official documentation:** https://search.foldseek.com/

#### Foldseek standalone software

Foldseek software can be used to analyze a large number of protein sequences or structures. 

**Foldseek installation:**

```

# Create a conda environment in the folder of your choice

conda create -p /path/to/folder/foldseek

# Install foldseek with conda

conda install -p /vol/data/tools/foldseek -c conda-forge -c bioconda foldseek

# Move to the bin folder of foldseek

cd /path/to/folder/foldseek/bin

# Test the installation

./foldseek

This should display the usage instructions for running foldseek

```

**Running foldseek:**

```
- There are a number of options available in foldseek starting with structural search of simple and complex proteins, clustering, and database creation to name a few. Here we will look at sample commands for structural search of proteins in the light of their relevance to functional annotation.

foldseek easy-search /path/to/query/structure/sample /path/to/target/structures/database aln /path/to/tmp

The easy-search command above is a part of the easy workflows of foldseek. It enables structural search of simple single chain protein structure files or FASTA files against a target database. sample_aln represents the resulting srtuctural alignment file and tmp is the temporary folder for storing intermediate files during the run.

foldseek search sampleDB targetDB resultDB tmp

The above search command is a part of the main workflows of foldseek. It enables searching a database of protein structures against another target database and also provides options for sensitive search and clustering with additional parameters. In the above command, sampleDB represents the sample database, targetDB represents the target database, and resultDB represents the result database of the run with temporary files being stored in the tmp folder.

```
Official documentation: https://github.com/steineggerlab/foldseek

Official video tutorial: https://www.youtube.com/watch?v=k5Rbi22TtOA

## Expression analysis-based annotation tools

### GENIE3

GENIE3 is an R package that performs gene regulatory network analysis from expression data. It is based on machine learning and uses tree-based ensemble methods like Random Forests for the coexpression analysis. It is available as Python, Matlab and R-based (R/C) implementations. However, the R/C implementation is stated as the fastest GENIE3 implementation by the developers. Hence the installation and usage instructions given below are for this R/C implementation of GENIE3. This requires R to be installed priorly. https://cran.r-project.org/ provides detailed steps on R installation for the different operating systems.

**GENIE3 installation:**
```
# The installation must be done in R

#Install BiocManager for Bioconductor packages

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install GENIE3 from Bioconductor
BiocManager::install("GENIE3")

```

**Running GENIE3:**

```
GENIE3 [https://doi.org/10.1371/journal.pone.0012776] requires expression data file in the following format - the rows should correspond to genes and the columns should correspond to samples. It is important to note that normalized expression data needs to be loaded for the GENEIE3 analysis.

# Load expression data file into R

expr_data <- read.table("sample_counts.tsv", header =TRUE, row.names = 1)

# Load the GENIE3 library into R

library(GENIE3)

# Set a fixed seed value to ensure reproducibility of results

set.seed(42)

# Obtain the weight matrix from GENIE3

weightMat <- GENIE3(expr_data)

- The weight matrix has the genes in rows and columns and specifies the correlation weight between a gene pair. The greater the weight between a gene pair, stronger is the regulatory link between them.

- By default all genes in your analysis are taken as regulatory candidates. In case you know apriori about which genes are transcription factors, that subset of genes can be specified as regulators.

# Specify a subset of regulator genes

regulators <- c("geneA", "geneB")

# GENIE3 analysis based on the new set of regulators

weightMat <- GENIE3(expr_data, regulators=regulators)

```
**Official documentation:** https://bioconductor.org/packages/release/bioc/html/GENIE3.html

## Annotation tools employing combined approaches

### CoExpPhylo

CoExpPhylo [https://doi.org/10.1186/s12864-025-12061-3] is a Python script that combines coexpression analysis and phylogeny for robust ortholog identification. It can be used to identify genes involved in biosynthetic pathways across a number of species. It requires bait sequences of genes known to be involved in the specific pathway of interest for facilitating the phylogeny analysis. It also provides options to upload the phylogenetic trees obtained in the analysis and view them directly on iTOL phylogenetic tree viewer. 

**CoExpPhylo installation:**
```
# Clone the GitHub repository in a folder of your choice

git clone https://github.com/bpucker/CoExpPhylo.git

# Move into the CoexpPhylo folder

cd CoExpPhylo

# Create the conda environment

conda env create -f environment.yml

# Activate the conda environment

conda activate env_CoExpPhylo
```

**Running CoExpPhylo:**

```
python3 coexp_phylo.py --config sample_config.txt --out sample_results

- The config file is a simple TXT file that needs four mandatory columns separated by commas and an optional fifth column.

- The four mandatory columns to be specified in order are name of the species to be analyzed, full path to the expression file of the species, full path to the coding sequence FASTA file of the species, full path to the baits file that has the gene identifiers of the bait sequences.

- The full path to be specified must also have the name of the respective file.

- The fifth optional column is full path to the peptide FASTA file of the species.

```

**Official documentation:** https://github.com/bpucker/CoExpPhylo

### KIPEs

Knowledge-based Identification of Pathway Enzymes (KIPEs) [https://doi.org/10.1371/journal.pone.0294342] is a Python tool that helps in automatic annotation of genes involved in flavonoid biosynthesis. Given enough knowledge of another biosynthetic pathway like well-known gene players, the tool can be adopted for such pathways as well. For instance, the tool was recently extended for carotenoid biosynthesis. The tool has a very comprehensive list of bait sequences of flavonoid biosynthesis and carotenoid biosynthesis genes in a number of plants that is used for the ortholog search step. The script starts with a local alignment using BLAST to identify orthologs and then looks for sequence level properties of the protein sequences like conserved residues. Since an enzyme's catalytic functions are heavily dependent on these conserved residues, the tool in a way combines sequence similarity with functional cues of the enzyme to identify the correct ortholog. Along with this, it also offers an option to perform a global alignment of top candidates from local alignment and infer orthologs from a phylogenetic tree, combining multiple levels of evidence to determine orthologs. 

**KIPEs installation:**

```
# Installing dependencies

- The tool needs dendropy, BLAST, MAFFT as mandatory depedencies and FastTree as optional dependencies. The installation instructions for BLAST and MAFFT can be found in the aligners section above and those for FastTree can be found in the phylogenetic analysis section above. dendropy can be installed as follows:

# Install the pip package manager

sudo apt install python3-pip

# Install dendropy using pip

python3 -m pip install git+https://github.com/jeetsukumaran/DendroPy.git

# Clone the KIPEs GitHub repository in a folder of your choice

git clone https://github.com/bpucker/KIPEs

# Move into the KIPEs folder

cd /path/to/folder/KIPEs

- In the above folder you will find the KIPEs3.py script, the flavonoid_baits.tar.gz which is the baits file for flavonoid biosynthesis genes and the carotenoid_baits.tar.gz which is the baits file for carotenoid biosynthesis genes.

```

**Running KIPEs:**

```

python3 KIPEs3.py --baits /path/to/folder/KIPEs/flavonoid_baits.tar.gz --out /path/to/output/folder --subject /path/to/sample.pep.fasta

- In the above command, sample.pep.fasta is the FASTA file of your query proteins that are to be annotated using KIPEs. It is also possible to specify a folder of peptide FASTA files to be annotated instead of a single file using the --subjectdir flag.

```  
**Official documentation:** https://github.com/bpucker/KIPEs

### DupyliCate

DupyliCate [https://doi.org/10.1101/2025.10.10.681656 ] is a Python tool for mining and analyzing gene duplications. It is able to identify gene duplications in a large number of species and is scalable. It combines gene expression data and helps study the expression divergence of the identified gene duplicates. Apart from this, it also offers an option to perform comparative analyses with respect to a reference species. In case, this reference-based analysis is opted for, apart from intra-species gene duplications, the tool also outputs orthologs across the species with respecto the reference species using a combination of local alignment, global alignment, synteny and phylogeny approaches. 

**DupyliCate installation:**

```

# docker installation

- For docker installation instructions please refer the official docker installation guide https://docs.docker.com/engine/install/

# Pull the latest DupyliCate docker image from docker hub

docker pull shakunthalan/dupylicate:latest

# Test the DupyliCate image

 docker run --rm -u $(id -u) -v /path/to/data:/data shakunthalan/dupylicate:latest

- In the above command /path/to/data:/data must be replaced by your specific host and container paths.

- The above command should display the DupyliCate usage indicating a successful pull of the docker image.

```

**Running DupyliCate:**

```

docker run --rm -u $(id -u) -v /path/to/data:/data shakunthalan/dupylicate:latest --gff path/to/folder/sample.gff --pep path/to/folder/sample.pep.fasta out /path/to/output/folder

- The above command is a basic sample command for running DupyliCate. Options to integrate the expression analysis and the ortholog search can be found in the official documentation page. In the above command it is important to note that the names of the gff and pep  file without the extension should be the same. The --gff and --pep flags can also take in folders in which case the folder path must include the folder name and include  a / at the end of the full path.

```

**Official documentation:** https://github.com/ShakNat/DupyliCate

## References



If you have questions about plant genomics that were not answered by any of these resources, please feel free to get in touch with the Plant Biotechnology and Biotechnology research group at the University of Bonn.
