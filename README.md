CoreGeneBuilder
================

CoreGeneBuilder can be used to extract a core genome (or a persistent genome) from a given set of bacterial genomes.


CITATION
--------
Please cite CoreGeneBuilder using the following DOI:
[![DOI](https://zenodo.org/badge/62309201.svg)](https://zenodo.org/badge/latestdoi/62309201)

CONTACT
-------
For help please contact:

Julien Guglielmini <br>
julien.guglielmini@pasteur.fr <br>
Institut Pasteur <br>
Bioinformatics and Biostatistics Hub <br>
C3BI, USR 3756 IP CNRS <br>
25-28 rue du docteur Roux <br>
75015 Paris, France

Sylvain Brisse <br>
sylvain.brisse@pasteur.fr <br>
Institut Pasteur <br>
Microbial Evolutionary Genomics <br>
CNRS, UMR 3525 <br>
25-28 rue du docteur Roux <br>
75015 Paris, France <br>


INSTALLATION
------------
An appliance 'CoreGeneBuilder' is available on the IFB cloud: https://www.france-bioinformatique.fr/fr/cloud. <br>
A docker container 'CoreGeneBuilder' is also hosted on the docker registry BioShadock: https://docker-ui.genouest.org. <br>
See `INSTALL.md` file for instructions.


USAGE
-----
Use this command to display the usage: <br>
`$ coregenebuilder`

Or if you run CoreGeneBuilder using docker, prefix this command by: <br>
`$ docker run -v <your_local_directory_where_are_your_data>:/root/mydisk coregenebuilder`

Example: <br>
~~~~
$ DOCKER_DATA=/home/dupont/cganalysis
$ docker run -v $DOCKER_DATA:/root/mydisk coregenebuilder
~~~~

QUICK START
-----------
#### Run CoreGeneBuilder on EXAMPLE data
One dataset is provided (inputs only).
It can be founded into directory `data/klpn5refannot`. <br>
More precisely:
~~~~
# installation on the IFB cloud
$ cp -pr /usr/local/share/coregenebuilder/data/klpn5refannot /root/mydisk/
$ cd /root/mydisk/klpn5refannot

# installation from the docker repository
### we suppose that your data will be into the following local directory:
$ DOCKER_DATA=/home/dupont/cganalysis
### then download the data:
$ cd $DOCKER_DATA
$ wget https://github.com/C3BI-pasteur-fr/CoreGeneBuilder/archive/v1.0.tar.gz
$ tar -zxf v1.0.tar.gz && mv ./CoreGeneBuilder-1.0/data/klpn5refannot . && rm -r ./CoreGeneBuilder-1.0 ./v1.0.tar.gz
### the quick start data are now here:
$ cd $DOCKER_DATA/klpn5refannot

# installation from the coregenebuilder git repository
$ CGPIPELINE=<where_coregenebuilder_distribution_is_installed_on_your_local_machine>
### for example:
$ CGPIPELINE=/home/dupont/CoreGeneBuilder
$ cd $CGPIPELINE/data/klpn5refannot
~~~~

We call `$DIR` the analysis directory. Here `DIR=klpn5refannot`. <br>
Run this command to test your pipeline installation: <br>
~~~~
#provided genbank annotation
$ coregenebuilder -d klpn5refannot -n klpn -g MGH78578_NC.fasta -a MGH78578_NC.gb -e "NC_" -p 95 -t 1 -s 3  
~~~~



#### Run CoreGeneBuilder on YOUR data
To run new analyses from our dataset, you must create this directory/file architecture. <br>
We call `$DIR` the new analysis directory. Here `DIR=cg_analysis_ex`.

1.Move to data directory:
~~~~
# if you run coregenebuilder on the IFB cloud :
$ cd /root/mydisk
$ mkdir cg_analysis_ex
~~~~
~~~~
# if you run the docker image:
$ cd $DOCKER_DATA
$ mkdir cg_analysis_ex
~~~~
~~~~
## if coregenebuilder is installed on a local machine from the git repository:
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

# or create a link to the target input directory that contains fasta files
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

#### For more information, please refer to the manual of CoreGeneBuilder :
 `doc/CoreGeneBuilder_manual.md`


ACKNOWLEDGMENTS
---------------
  We thank Bertrand Néron and Amandine Perrin of Institut Pasteur for their contribution to
the deployment of CoreGeneBuilder on the IFB cloud and as a docker image on the registry BioShadock. <br>
  This work was financially supported
by the French Institute of Bioinformatics (Grant ANR-11-INBS-0013)
and by the Pasteur International Bioresources Network (PIBnet) programme.


AUTHORS
-------
 Elise Larsonneur, Marie Touchon, Damien Mornico, Alexis Criscuolo, Sylvain Brisse, Eduardo P. C. Rocha

LICENSING
---------
CoreGeneBuilder is distributed under the terms of the
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0).
For further details see COPYING file.
