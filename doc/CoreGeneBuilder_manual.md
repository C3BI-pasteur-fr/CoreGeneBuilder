Manual of CoreGeneBuilder
=========================

  I-CITATION <br>
  II-INTRODUCTION <br>
  III-ARCHITECTURE OF THE PIPELINE <br>
  IV-USAGE <br>
  V-QUICK START <br>

  ANNEXES <br>
    ANNEXE 1 - DETAILS OF EXPECTED INPUT AND OUTPUT FILES <br>
    ANNEXE 2 - DETAILS OF MODULES <br>


___

I-CITATION
----------
Please cite CoreGeneBuilder using the following doi:
<TODO: INSERT_DOI_OF_GITHUB_REPOSITORY>


___


II-INTRODUCTION
---------------
CoreGeneBuilder can be used to extract a core genome (or a persistent genome) from a given set of bacterial genomes. <br>

It is separated in three ordered modules: **DIVERSITY, ANNOTATION, COREGENOME**. <br>
It requires as input fasta files of the genomes, and optionally an annotation file of one of those genomes (genbank full format).
Warning: the genomes must have similar length (see section **DIVERSITY MODULE**). <br>
It is better to have genomes as contigs and not as scaffolds. <br>
Because `andi` requires closely related genomes, it is better to provide closely related genomes as input. The program is typically intended to build core genomes of strains within a bacterial species or closely related species. <br>
The optionally provided reference annotation must be linked to the provided reference genome, i.e. they must have the same accession ids (eg accession/locus in genbank files must be the same as sequence id in reference genome fasta file). <br>

CoreGeneBuilder outputs: <br>
- list of core genes: `core_genome/CoreGenome*.lst`
- fasta files of these genes and their variants:  `core_genome/core_genes_by_gene`
- fasta files of core genes for each genome:  `core_genome/core_genes_by_genome`

See an overview of CoreGeneBuilder in `doc/coregenebuilder_overview_fig.pdf`. <br>
<br>
To simplify, we will call the analysis directory `$DIR` <br>

<br>
#### DIVERSITY MODULE
This module starts by **filtering out** genome sequences based on **several criteria**: <br>
- **contig length** (`default: 500 bp`)
- **N50 length** (`default: 0` (no filter))
- **number of contigs** (`default: -1` (no filter))
- **genome size**: for next steps of the DIVERSITY and ANNOTATION modules, the external tools (`andi`, `ecamber`) require to have genomes with limited variation in size (number of cumulated bases). Therefore genomes with a genome size that is lower than [ 0.5 * median(genome_sizes) ] or higher than [ 1.5 * median(genome_sizes) ] are discarded from the analysis. Indeed, `ecamber`, by default, does not annotate variable size genomes ; and `andi` outputs pairwise distance values equals to 'NA' for genome pairs with too different sizes.
- **presence of scaffolding stretches of Ns in sequences**: by default, sequences having stretch of 10 ’N’s or more are split (default: 10 ; -1: no filter) ; minimum length of scaffolding stretch of Ns (`default: 10`).

These filters define a new set of clean genomes. <br>
Note that the program will exit if the optionally provided reference genome happens to be filtered out. <br>

The next step of this module analyzes the **phylogenetic diversity** of the set of remaining genomes. <br>
This step allows to optionally select a subset of genomes based on the phylogenetic diversity. This is meant to remove multiple genomes from e.g. a single outbreak or single clonal group. <br>
The output set of genomes of this module will be used for the subsequent steps: annotation and core genome building.
The option `-s` defines this parameter (default: all genomes are selected for the next steps). <br>
If your genome number is higher than 300 and depending on the size of your genomes, we recommend selecting a subset of those genomes. <br>
The program `andi` is used to compute pairwise genome distances and outputs the corresponding matrix. This matrix is then used to generate a hierarchical clustering (average linkage clustering, i.e. UPGMA). If option `-s` is supplied, the genome having the biggest N50 contig size in each cluster is selected. <br>

###### INPUT FILES
Requires **fasta** files as genome input files in a directory called `assemblies`. <br>
Optionally you can provide an annotation file (**genbank** format only) in a directory called `ref_gbk_annotation` <br>
###### OUTPUT FILES
The output files will be placed in the directory `$DIR/diversity`. <br>

<br>
#### ANNOTATION MODULE
The next module consists in the _**de novo**_ **syntactic annotation** of the remaining genomes (using `prodigal`), except for the reference genome if a reference genbank annotation is provided. <br>
The positions of gene starts of all genomes are then corrected (harmonizationn of start codons). <br>
These two steps are done by `ecamber` (`prodigal` is a dependency of `ecamber`). <br>
If a reference annotation is provided, based on orthologous gene clustering by `ecamber` (cluster ids), **the functional annotation** of the reference genome is transferred to other annotations. <br>

###### INPUT FILES
Input files came from `$DIR/diversity` <br>
###### OUTPUT FILES
The output files will be placed in the directories: <br>
~~~~
$DIR/genes
$DIR/proteins
$DIR/ecamber_output
~~~~

<br>
#### CORE_GENOME MODULE
The last module generates **a set of core (or persistent) genes** from input proteomes.

- The program `opscan` (written by **Joel Pothier**) is used to compute this set:  it identifies homologs as unique pairwise reciprocal best hits (eg bidirectional best hits) between the reference proteome as a pivot and each of the remaining proteomes, using end-gap free global alignment. It compares reference proteome against all other proteomes, but does not compare all proteomes versus all proteomes.

- Hits with **less than x% similarity** (`default:80%` ) **in amino acid sequence** and **more than x% difference in protein length** (`default:20%`) are discarded.

- For every pairwise comparison, the resulting list of homologs is then refined taking into account the conservation of gene neighborhood. <br>
To do this, synteny criteria (see **SYNTENY** options `-S`,`-R` of coregenebuilder) are used, offering the possibility to retain only orthologous genes that are located in blocks of synteny. Indeed, genomes from the same species typically show low levels of genome rearrangements and this information can be used to identify orthologs more accurately. This synteny constraint can be ignored using `-S 0`.
- The filtered set of genes (CDS only) is then post-processed to keep the set of homologous genes shared by at least p % of the genomes (see option `-p`, `default: 95%`). <br>
If p is equal to `100%`, the final gene set is **the true core genome**. <br>
Else if p is `lower than 100%` (and ideally greater than 70%), the final gene set can be considered as representing **the persistent genome**, rather than a core strictly speaking. <br>

###### INPUT FILES
Input files came from these directories: <br>
~~~~
$DIR/genes
$DIR/proteins
~~~~
###### OUTPUT FILES
The output files will be placed in the directory
~~~~
$DIR/core_genome
~~~~


___

III-ARCHITECTURE OF THE PIPELINE
--------------------------------

The pipeline is organized in several directories:
~~~~
config  data  doc  ext-tools  src
~~~~
~~~~
config    : contains config files used to define static parameters of external tools (ecamber, opscan) and to define environment variables (path of dependencies)
src       : contains main script (coregenebuilder) to launch the pipeline and secondary other scripts (bash, python, awk) called by this main script.
ext-tools : contains installed external tools
data      : directory dedicated to input and output files
doc       : contains this README.txt file and a figure showing the main steps of the pipeline and their associated parameters.
~~~~


___

IV-USAGE
--------

Use this command to display the usage: <br>
`$ coregenebuilder`

Or if you run CoreGeneBuilder using docker, prefix this command by: <br>

`$ docker run -v
 <your_local_directory_where_are_your_data>:/root/mydisk coregenebuilder`

Example: <br>
`$ docker run -v /home/dupont/cganalysis:/root/mydisk coregenebuilder`

~~~~
USAGE:
     coregenebuilder [options] -d <input_directory> -n <name_of_4_chars>
     -d <inDirectory>  directory where are stored input and output files, it must contain at least a directory called 'assemblies' where are genomic sequence fasta files
     -n <string>  four letter name (ex: esco (es:escherichia; co:coli))
     where 'options' are:
     ### PRE-PROCESSING OF GENOMIC SEQUENCES ###
     -f <int>     filtering out genomic sequences below this length value (bp) (default: 500); if 0, no filter
     -N <int>     cut sequences with mers of 'N' of size '-N' (scaffolds) into contigs (default : 10); if -1, no filter
     -y <int>     filtering out genomes having N50 below this N50 value (bp) (default: 0); if 0, no filter
     -c <int>     filtering out genomes having number of sequences above this value (default: -1); if -1, no filter
     ### DIVERSITY ###
     -s <int>     number of genomes to select for core genome construction (default : -1); if -1, no selection (all genomes are kept)
     ### CORE-GENOME ###
     -i <int>     similarity percent (default: 80) -- core genome construction step
     -l <float>   protein length ratio (default: 1.2) -- core genome construction step
     -S <int>     synteny - synteny threshold: minimal number of syntenic genes found in the neighborhood of a given homologous gene,
                         the boundaries of explored neighborood are defined by window size (option -R); (default: 4); if 0, no synteny criteria is applied
     -R <int>     synteny - radius size around each analyzed homologous gene (a number of genes) (default: 5) -- window_size=(radius_size * 2 + 1)
     -p <int>     core genes are present at least in p% of the genomes (default: 95) -- if '-p 100' is supplied, core gene set is output; else if 'p' is lower than 100, persistent gene set is output
     ### GENERAL OPTIONS ###
     -z <STEPS>   steps to run ; ex: '-z DAC' or '-z DA' or '-z D' ; D: diversity - A: annotation - C: coregenome (default: DAC)
     -g <inFile>  reference genome fasta file -- annotation and core genome construction steps
     -a <inFile>  reference genome genbank file -- annotation step
     -e <inFile>  prefix string of contig ids found in the reference genome fasta and genbank files (eg id in fields LOCUS and ACCESSION of genbank file) (ex: 'NC_'; 'NZ_' ; 'AKAC' ; etc)
     -t <int>     number of threads used during diversity and annotation steps (default: 1)
     -r <int>     if 1, remove files from precedent run of this pipeline inside directory <input_directory> (default: 0)
     -v           print version number
     -h           print help

DESCRIPTION:
     CoreGeneBuilder version 1.0
     CoreGeneBuilder extracts a core genome or a persistent genome from a set of bacterial genomes.
     The core genome is the set of homologous genes (protein families) shared by all genomes.
     The persistent genome is the set of homologous genes (protein families) shared by less than 100% of the genomes.

EXAMPLES:
     #provided reference genome and related genbank annotation, the functional annotation of reference genome will be transferred to the other genomes:
     coregenebuilder -d klpn5refannot -n klpn -g MGH78578_NC.fasta -a MGH78578_NC.gb -e NC_ -p 95 -t 4

     #provided reference genome but not provided related genbank annotation, the reference genome will be de novo annotated:
     coregenebuilder -d klpn5refannot -n klpn -g MGH78578_NC.fasta -p 95 -t 4

     #not provided reference genome and related genbank annotation, the reference genome will be the first on genome list sorted alphabetically:
     coregenebuilder -d klpn5refannot -n klpn -p 100 -t 4
~~~~

___

V-QUICK START
-------------

#### Run CoreGeneBuilder on EXAMPLE data
One dataset is provided (inputs only).
It can be founded into directory `data/klpn5refannot`. <br>
More precisely:
~~~~
# installation from ansible (on IFB cloud) or docker repository
/usr/local/share/coregenebuilder/data/klpn5refannot
# IFB cloud only:
cp -pr /usr/local/share/coregenebuilder/data/klpn5refannot /root/mydisk/

# installation from coregenebuilder git repository
CGPIPELINE=<where_coregenebuilder_distribution_is_installed_on_local_machine>
$CGPIPELINE/data/klpn5refannot
~~~~

We call `$DIR` the analysis directory. Here `DIR=klpn5refannot`. <br>
Run this command to test your pipeline installation: <br>
~~~~
#provided genbank annotation
$ coregenebuilder -d klpn5refannot -n klpn -g MGH78578_NC.fasta -a MGH78578_NC.gb -e "NC_" -p 95 -t 6  
~~~~



#### Run CoreGeneBuilder on YOUR data
To run new analyses from our dataset, you must create this directory/file architecture. <br>
We call `$DIR` the new analysis directory. Here `DIR=cg_analysis_ex`.

1.Move to data directory:
~~~~
# if you run coregenebuilder on IFB cloud:
$ cd /root/mydisk
$ mkdir cg_analysis_ex
~~~~
~~~~
# if you run the docker image:
$ cd <your_local_directory_where_are_your_data>
$ cd /home/dupont/cganalysis
$ mkdir cg_analysis_ex
~~~~
~~~~
## if coregenebuilder is installed on local machine:
$ cd $CGPIPELINE/data
$ mkdir cg_analysis_ex
~~~~

2.Create input directory to store genome fasta files, they must be already stored in a directory named `assemblies` <br>
~~~~
$ mkdir cg_analysis_ex/assemblies
~~~~
And import genomes:
~~~~
$ cp <PATH_OF_GENOME_FASTA>/* cg_analysis_ex/assemblies/.

# or create a link to the target input directory that contain fasta files
$ cd cg_analysis_ex/assemblies
$ ln -s <PATH_OF_GENOME_FASTA> assemblies
~~~~
If you provide a reference genome (option '-g'), it will be in the directory `assemblies`. <br>
Here we call it `ref.fasta`:
~~~~
$ ls assemblies/ref.fasta
~~~~

3.If you want to provide a genbank annotation, create the directory named `ref_gbk_annotation`.
Here we call it `ref.gb`:
~~~~
$ mkdir cg_analysis_ex/ref_gbk_annotation
$ cp ref.gb cg_analysis_ex/ref_gbk_annotation/.
~~~~

4.Then launch the pipeline with these parameters for example:
~~~~
$ coregenebuilder -d cg_analysis_ex -n klpn -g ref.fasta -a ref.gb -e NC_ -p 95 -t 6
~~~~
 Note:<br>
 `-n klpn` => value `klpn` to designate genomes of KLebsiella PNeumoniae <br>
 `-e NC_` => prefix of contig ids of files `ref.fasta` and `ref.gb` <br>




#### Architecture of `$DIR` when all steps of pipeline have been done
~~~~
assemblies        : contains genomes in fasta format (extensions .fasta .fas .fa .fna are only accepted)
ref_gbk_annotation: contains a reference genbank annotation and
logs              : contains log files of the 3 modules of the pipeline (DIVERSITY, ANNOTATION, COREGENOME)
progress_file.txt : finished steps and their status (status 'OK' if no errror)
diversity         : contains input and output files of module DIVERSITY
ecamber_output    : contains some output files of module ANNOTATION (ecamber output)
genes             : contains nucleic sequences of CDS
proteins          : contains amino-acid sequences of CDS
core_genome       : contains input and output files of module COREGENOME, contains core genes as nucleic and amino-acid sequences (fasta format)
~~~~


___

ANNEXES
---------

ANNEXE 1 - DETAILS OF EXPECTED INPUT AND OUTPUT FILES
-----------------------------------------------------

###### DIVERSITY FILES
Given these `CoreGeneBuilder` options:
~~~~
-n klpn: 2 first lowercase letters of Genus followed by two first letters of species ; example: klpn: kl for Klebsiella ; pn for pneumoniae
-s 5   : 5 genomes will be selected among the initial genomes
-N 1   : scaffolds are cut into contigs (polyN removing)
~~~~
The ouput files and directory will be:
~~~~
$ cd $DIR/diversity
$ ls

genomes                                 #directory 'genomes' contains all genomes:
## - remaining genomes after first filtering based on assembly statistics; they are renamed using prefix specified by option '-n' (here 'klpn') (see renamed_genomes_list.txt for details)
## discarded genomes after first filtering based on assembly statistics; they are not renamed like the others but by adding these extensions:
.n50_nbctg_problem  #if options '-y' and '-c' are used, and the genome is filtered out based on these two thresholds (N50, number of contigs)
.sizeproblem        #if genome have a size that is lower than (0.5 * median genome size) or higher than (1.5 * median genome size).

genomes_and_genomesize_list             #assembly statistics of all genomes of directory 'assemblies'
renamed_genomes_list.txt                #new and old names of genomes (old names are those found in directory 'assemblies', new names are those found in directory 'diversity/genomes')
klpn.5.mat                              #genome distance matrix output by andi, computed from remained genomes after first filtering of assemblies (filters based on assembly statistics: contig length (option '-f'), N5O length (option '-y'), number of contigs (option '-c'), after cutting scaffolds into contigs if '-N 1' or not (-N 0))
renamed_genomes_list.links.txt          #link names associated to new names of genomes
klpn.5.mat.renamed                      #genome distance matrix output by andi, but with old names of genomes
klpn.4.woref.mat                        #the same file as klpn.5.mat but without reference
renamed_genomes_list.woref.links.txt    #the same file as renamed_genomes_list.links.txt but without reference
klpn.4.woref.mat.renamed                #the same file as klpn.5.mat.renamed but without reference
diversity_klpn_plots.pdf                #hierarchical clustering (average linkage) computed from andi matrix 'klpn.5.mat'
diversity_klpn_plots_renamed.pdf        #the same file as 'diversity_klpn_plots.pdf' but with old names of genomes
~~~~

If `-s` option is supplied, there are also these additional files:
~~~~
#reference genome is removed from the clustering to keep it, so it is not in these files:
diversity_klpn_genome_cluster_list_wo_ref.txt
diversity_klpn_genome_cluster_list_wo_ref_renamed.txt
diversity_klpn_genome_cluster_list_wo_ref.txt.infos
diversity_klpn_genome_cluster_list_wo_ref.txt.infos.sort  #clusters of genomes, sorted by cluster then N50 length (fields: 1:genome_name, 2:cluster_id, 3:N50_length, 4:genome_size, 5:contig_number)
diversity_klpn_plots_wo_ref.pdf         #the same file as 'diversity_klpn_plots.pdf' but without reference genome, and with displayed clusters ; is computed from andi matrix 'klpn.4.woref.mat'
diversity_klpn_plots_wo_ref_renamed.pdf #the same file as 'diversity_klpn_plots_renamed.pdf' but without reference genome, and with displayed clusters
~~~~

The last output file generated by the module DIVERSITY is `selected_genomes_list.txt`. <br>
It contains selected genomes after steps of: <br>
1. (optional) cuts of scaffolds into contigs, <br>
2. filtering based on assembly statistics <br>
3. selection of N (= VALUE OF `-s` PARAMETER) genomes among filtering genomes: the tree output from hierarchical clustering is cut into N clusters, when cluster contain multiple genomes, the chosen genome is those that has the higher N50 length among genomes of this cluster. This method allows to select representative genomes that cover the genomic diversity of analyzed genomes. By default all filtered genomes are kept. <br>
This set of genomes is then annotated by the next module, the module ANNOTATION. <br>




###### ANNOTATION FILES:
The annotation files are separated in multiple directories
~~~~
$ cd $DIR
$ ls

#input and output files of program ecamber are in this directory:
$DIR/../../ext-tools/ecamber/datasets/$DIR

#ecamber output post-processed (new formats) are in these directories
$DIR/ecamber_ouput  # ortholog clusters of ecamber
$DIR/genes          # fasta files of CDS containing nucleic sequences (.gen), sequence headers in (.lst) files ; (.fasta) correspond to ecamber outputs
$DIR/proteins       # fasta files of genes containing proteic sequences (.prt), sequence headers in (.lst) files ; (.fasta) correspond to ecamber outputs
~~~~

Directories `$DIR/genes` and `$DIR/proteins`:
~~~~
#FIELDS OF FASTA HEADERS (files .gen .prt .lst of directories `$DIR/genes` and `$DIR/proteins`):
## $1  gene_id (gembase format)
## $2  strand (D/C)
## $3  start_codon
## $4  stop_codon
## $5  start
## $6  stop
## $7  feature_type (=CDS)
## $8  gene_name(s) or locus_tag(s) if no gene_name (several (transferred) gene_name(s) are separated by a ';' ; 'x' if not available)
## $9  gene_size
## $10 formatted contig_id (=fastaname if genome have only ONE contig)
## $11 locus_tag(s) (several (transferred) locus_tag(s) are separated by a ';' ; 'x' if not available)
## $12 gene_product(s) ('x' if not available)
## $13 strain (input fasta name)
~~~~

######  I - IF A GENBANK ANNOTATION IS PROVIDED
The _de novo_ syntactic annotation of the genomes is done (`prodigal`), except for the reference genome. <br>
The positions of gene starts of all genomes are then corrected (harmonization of start codons) (`ecamber`). <br>
Then, the functional annotation of the reference genome (genbank file) is transferred to other annotations.

~~~~
DIR=klpn5refannot

#reference genome - and genbank annotation provided (only one gene_name, locus_tag, product by gene)
$ less klpn5refannot/genes/klpn.0000001.c001.lst
KLPN0000001gAAA_000010 D ATG TAA 21 137 CDS KPN_RS00005 117 klpn_0000001_c001 KPN_RS00005 hypothetical_protein MGH78578_NC
KLPN0000001iAAA_000020 D ATG TAA 340 2802 CDS thrA 2463 klpn_0000001_c001 KPN_RS00010 bifunctional_aspartokinase_I/homoserine_dehydrogenase_I MGH78578_NC
KLPN0000001iAAA_031180 C ATG TAA 3116660 3117163 CDS x 504 klpn_0000001_c001 x insertion_element_IS1_protein_InsB;transposase MGH78578_NC
KLPN0000001iAAC_055650 C TTG TAA 28300 28803 CDS x 504 klpn_0000001_c003 x insertion_element_IS1_protein_InsB;transposase MGH78578_NC

#other genome - annotations transferred from reference genbank annotation
$ less klpn5refannot/genes/klpn.0000002.c001.lst
KLPN0000002iAAA_000020 D ATG TAA 531 1754 CDS KPN_RS08450 1224 klpn_0000002_c001 KPN_RS08450 transcriptional_regulator 03-9138_380
KLPN0000002gAAB_000070 C ATG TGA 255 1919 CDS mhpA 1665 klpn_0000002_c002 KPN_RS11435 3-(3-hydroxy-phenyl)propionate/3-hydroxycinnamic_acid_hydroxylase 03-9138_380
KLPN0000002gAAC_000080 C ATG TGA 9 1244 CDS x 1236 klpn_0000002_c003 x x 03-9138_380
KLPN0000002iACD_009510 D ATG TAA 6914 8431 CDS lysS 1518 klpn_0000002_c064 KPN_RS02675;KPN_RS17790 lysine--tRNA_ligase 03-9138_380
KLPN0000002iALX_051450 C ATG TAA 4502 5005 CDS KPN_RS04085;KPN_RS26420;KPN_RS26635;KPN_RS26735;KPN_RS26750 504 klpn_0000002_c372 KPN_RS04085;KPN_RS26420;KPN_RS26635;KPN_RS26735;KPN_RS26750 insertion_element_IS1_protein_InsB;transposase 03-9138_380
~~~~

######  II - IF NO GENBANK ANNOTATION IS PROVIDED
A _de novo_ syntactic annotation of the genomes is done (`prodigal`) (also for reference genome). <br>
The positions of gene starts of all genomes are then corrected (harmonization of start codons) (`ecamber`).

~~~~
DIR=klpn5

#reference genome - without genbank annotation provided (potentially several transferred gene_names, locus_tags, products by gene)
$ less klpn5/genes/klpn.0000001.c001.lst
KLPN0000001gAAA_000010 D ATG TAA 340 2802 CDS x 2463 klpn_0000001_c001 x x MGH_78578_25
KLPN0000001iAAA_000020 D ATG TAA 2804 3733 CDS x 930 klpn_0000001_c001 x x MGH_78578_25
KLPN0000001iAAA_000030 D ATG TAA 3737 5017 CDS x 1281 klpn_0000001_c001 x x MGH_78578_25

#other genome - without transferred annotations
$ less klpn5/genes/klpn.0000002.c001.lst
KLPN0000002gAAA_000010 D ATG TGA 18 410 CDS x 393 klpn_0000002_c001 x x 03-9138_380
KLPN0000002iAAA_000020 D ATG TAA 531 1754 CDS x 1224 klpn_0000002_c001 x x 03-9138_380
KLPN0000002iAAA_000030 C ATG TAA 2530 3825 CDS x 1296 klpn_0000002_c001 x x 03-9138_380
~~~~
The proteome `$DIR/proteins/*.prt` is then used as input for the last MODULE, module COREGENOME.


######  ECAMBER OUTPUTS:
~~~~
$DIR/ecamber_ouput     #contains ortholog clusters output by ecamber; come from directory $DIR/../../ext-tools/ecamber/datasets/$DIR/output/
$DIR/genes/*.fasta     #come from $DIR/../../ext-tools/ecamber/datasets/$DIR/output/genes_dna/
$DIR/proteins/*.fasta  #come from $DIR/../../ext-tools/ecamber/datasets/$DIR/output/genes_aa/
~~~~

Header of ecamber dna fasta and aa fasta, eg header of files `$DIR/genes/*.fasta` and files `$DIR/proteins/*.fasta`:
~~~~
# if reference functional annotation supplied
$ grep '>' klpn.0000001.c001.fasta |head
>KPN_RS00005|7944|KPN_RS00005|+|21|137|klpn_0000001_c001|KPN_RS00005|hypothetical_protein|MGH78578_NC
>KPN_RS00010|2848|thrA|+|340|2802|klpn_0000001_c001|KPN_RS00010|bifunctional_aspartokinase_I/homoserine_dehydrogenase_I|MGH78578_NC
>KPN_RS00015|6009|KPN_RS00015|+|2804|3733|klpn_0000001_c001|KPN_RS00015|homoserine_kinase|MGH78578_NC
>x|3844|x|+|19212|19421|klpn_0000001_c001|x|x|MGH78578_NC
#fields:
#gene_id|cluster_id|gene_name|strand|start|stop|contig_id
gene_id =locus_tag, if exist in genbank file
gene_id='x', if created by ecamber (transfer)
~~~~
~~~~
# if de novo annotation by ecamber
$ grep '>' klpn.0000005.c001.fasta |head
>1_1|4331|KPN_RS25580|+|481|1602|klpn_0000005_c001|KPN_RS25580|50S_ribosomal_protein_L16_arginine_hydroxylase|WGLW1_1561_contigs
>1_2|1934|KPN_RS25575|-|2183|1650|klpn_0000005_c001|KPN_RS25575|acetyltransferase|WGLW1_1561_contigs
>1_3|2038|KPN_RS25570|-|2460|2194|klpn_0000005_c001|KPN_RS25570|toxin-antitoxin_system_antitoxin_component|WGLW1_1561_contigs
#fields:
#gene_id|cluster_id|gene_name|strand|start|stop|contig_id
gene_id=prodigal_id, else:
gene_id='x', if created by ecamber (transfer)
~~~~

ecamber files of ortholog cluster: <br>
~~~~
$ ls $DIR/ecamber_ouput:
cluster_gene_names.txt
clusters_table.txt
~~~~

`cluster_gene_names.txt`: <br>
~~~~
$ less cluster_gene_names.txt    #contains only clusters having gene_name (if functional annotation is supplied (option '-a'))
cluster_id	gene_name
5989	cobT
3515	narZ
$ cat cluster_gene_names.txt     #only header line if no functional annotation is supplied
cluster_id	gene_name
~~~~

`clusters_table.txt`: <br>
~~~~
$ less clusters_table.txt       # all ortholog clusters output by ecamber (having gene_name or not (value 'x'))
cluster_id	cluster_gene_name	cluster_type	cluster_mg_count	cluster_mg_ann	cluster_strain_count	cluster_strain_ann	klpn.0000001.c001	klpn.0000002.c001	klpn.0000003.c001	klpn.0000004.c001	klpn.0000005.c001
5988	x	ANCHOR	1	1	1	1				1_276.272167.271832.-.klpn_0000004_c001
5989	cobT	ANCHOR	5	5	5	5	KPN_RS13285.2686757.2685690.-.klpn_0000001_c001	163_6.5088.4021.-.klpn_0000002_c163	277_10.9668.10735.+.klpn_0000003_c277	114_5.3622.4689.+.klpn_0000004_c114	10_52.62260.63327.+.klpn_0000005_c010
5980	x	ANCHOR	5	5	5	5	KPN_RS22920.4658435.4658794.+.klpn_0000001_c001	205_4.4855.4496.-.klpn_0000002_c205	391_7.7699.8058.+.klpn_0000003_c391	33_4.5138.4779.-.klpn_0000004_c033	4_4.5095.4736.-.klpn_0000005_c004
2257    narH    NON_ANCHOR      10      10      5       5       KPN_RS11930.2422546.2421011.-.klpn_0000001_c001;KPN_RS10125.2074657.2076201.+.klpn_0000001_c001 301_4.6516.8060.+.klpn_0000002_c301;227_53.55534.53999.-.klpn_0000002_c227      657_3.3205.1661.-.klpn_0000003_c657;677_5.3917.3042.-.klpn_0000003_c677 138_87.88468.86924.-.klpn_0000004_c138;35_109.1

=> the fields 8 until the last field correspond with each genome of the analysis (we will call them strain fields)
=> for each cluster of genes, each member gene is indicated in corresponding strain field,
=> for each <strain_name>, each gene is named by its gene_id/start/stop/strand/contig_id, several genes of the same strain are separated by ';'
details for each member gene:
example '10_52.62260.63327.+.klpn_0000005_c010' gene of the strain 'klpn.0000005.c001'
10_52               #gene_id (prodigal_id or locus_tag)
62260               #start (if strand +) / stop (if strand -)
63327               #stop (if strand +) / start (if strand -)
+                   #strand
klpn_0000005_c010   #contig_id

other fields:
=> cluster_mg_count (cluster_multigene_count) =  number of genes in the cluster
=> cluster_strain_count = number of strains found in the cluster
~~~~


If reference functional annotation is supplied, these three additional files are present in the directory `$DIR/ecamber_ouput`:
~~~~
- cluster_gene_names.ecamber_dna.txt   #only clusters having gene_name and/or locus_tag ; file 'cluster_gene_names.txt' with additional fields: locus_tag (cluster_gene_id) and gene function (cluster_gene_product)
- cluster_gene_names.ecamber_aa.txt    #exactly the same file as cluster_gene_names.ecamber_dna.txt
- clusters_table.ecamber.txt           # all ortholog clusters output by ecamber, it corresponds with file 'clusters_table.txt' post-processed to add these informations: when cluster has a functional annotation but no gene_name, a locus_tag stands for cluster_gene_name, if no functional annotation cluster_genename='x'
~~~~

~~~~
$ less cluster_gene_names.ecamber_dna.txt
cluster_id	cluster_gene_name	cluster_gene_id	cluster_gene_product
5989	cobT	KPN_RS13285	nicotinate-nucleotide--dimethylbenzimidazole_phosphoribosyltransferase
3515	narZ	KPN_RS10120	nitrate_reductase_A_subunit_alpha
~~~~

~~~~
$ less clusters_table.ecamber.txt
cluster_id	cluster_gene_name	cluster_type	cluster_mg_count	cluster_mg_ann	cluster_strain_count	cluster_strain_ann	klpn.0000001.c001	klpn.0000002.c001	klpn.0000003.c001	klpn.0000004.c001	klpn.0000005.c001
5988	x	ANCHOR	1	1	1	1				1_276.272167.271832.-.klpn_0000004_c001
5989	cobT	ANCHOR	5	5	5	5	KPN_RS13285.2686757.2685690.-.klpn_0000001_c001	163_6.5088.4021.-.klpn_0000002_c163	277_10.9668.10735.+.klpn_0000003_c277	114_5.3622.4689.+.klpn_0000004_c114	10_52.62260.63327.+.klpn_0000005_c010
5980	KPN_RS22920	ANCHOR	5	5	5	5	KPN_RS22920.4658435.4658794.+.klpn_0000001_c001	205_4.4855.4496.-.klpn_0000002_c205	391_7.7699.8058.+.klpn_0000003_c391	33_4.5138.4779.-.klpn_0000004_c033	4_4.5095.4736.-.klpn_0000005_c004
2257	narH	NON_ANCHOR	10	10	5	5	KPN_RS11930.2422546.2421011.-.klpn_0000001_c001;KPN_RS10125.2074657.2076201.+.klpn_0000001_c001	301_4.6516.8060.+.klpn_0000002_c301;227_53.55534.53999.-.klpn_0000002_c227	657_3.3205.1661.-.klpn_0000003_c657;677_5.3917.3042.-.klpn_0000003_c677	138_87.88468.86924.-.klpn_0000004_c138;35_109.115810.114275.-.klpn_0000004_c035	7_278.319231.320766.+.klpn_0000005_c007;14_335.338144.336600.-.klpn_0000005_c014
~~~~



###### CORE_GENOME FILES:
~~~~
Genomes-klpn5refannot.lst          : list of input proteomes used to build core genomes
CoreGenome-klpn5refannot.lst       : list of core genes present in all proteomes (equivalent to parameter -p 100) ; for each core gene, it contains as first field gene ids of the reference proteome and then gene ids of other proteomes
klpn.0000001.c001.histo            : for each core gene, the number of genomes having this gene.
klpn.0000001.c001.095.listprot     : list of core genes ; it contains gene ids of the reference proteome (genome and proteome named by number 1)
CoreGenome-klpn5refannot-095.lst   : list of core genes defined using parameter -p (here -p 95), eg core genes shared by 95% of the proteomes ; for each core gene, it contains as first field gene ids of the reference proteome and then gene ids of other proteomes ; if '-p' parameter value is lower than 100, this is not a strict core genome, but the persistent core genome (eg genes conserved in a majority of proteomes)
core_genes_by_genome               : directory containing output fasta files of core genes, they are grouped by genome, (see section 'ANNOTATION FILES' for .gen and .prt file format specifications)
core_genes_by_gene                 : directory containing output fasta files of core genes, they are grouped by gene, (see section 'ANNOTATION FILES' for .gen and .prt file format specifications)
temp                               : directory containing output files of program 'opscan' and its secondary scripts
~~~~


ANNEXE 2 - DETAILS OF MODULES
-----------------------------

###### MODULE COREGENOME - SYNTENY OPTIONS
~~~~
$ ./coregenebuilder
-S <int>     synteny - synteny threshold: minimal number of syntenic genes found in the neighborhood of a given homologous gene,
the boundaries of explored neighborood are defined by window size (option -R); (default: 4); if 0, no synteny criteria is applied
-R <int>     synteny - radius size around each analyzed homologous gene (a number of genes) (default: 5) -- window_size=(radius_size * 2 + 1)
~~~~

A threshold `-S` greater than zero triggers the search of synteny conservation in a set of homologous genes, to define a final subset of orthologous genes from this set of homologous genes. <br>
If this threshold is equal to zero, no search of synteny is done.
