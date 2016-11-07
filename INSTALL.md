Installation of CoreGeneBuilder
===============================

I   - Installation of a docker image <br>
II  - Installation on the IFB cloud <br>
III - Installation from source files on Unix machine <br>



I   - Installation of a docker image
------------------------------------
A docker container 'CoreGeneBuilder' is hosted on the docker registry BioShadock, you can download it here: https://docker-ui.genouest.org

You can then display help of CoreGeneBuilder using this command : <br>
~~~~
$ docker run -v <your_local_directory_where_are_your_data>:/tmp coregenebuilder

#example :
$ docker run -v /home/dupont/cganalysis:/tmp coregenebuilder
~~~~




II  - Installation on the IFB cloud
-----------------------------------
Nothing to be done. An appliance 'CoreGeneBuilder' is already available on the cloud of the French Institute of Bioinformatics (IFB). <br>
You need to have a user account to use CoreGeneBuilder on the cloud:
https://cloud.france-bioinformatique.fr/cloud

Display help using this command :
~~~~
$ coregenebuilder
~~~~


III - Installation from source files (git repository) on Unix machine
---------------------------------------------------------------------
###### DOWNLOAD
Download the distribution of CoreGeneBuilder here :
https://github.com/C3BI-pasteur-fr/CoreGeneBuilder

###### REQUIREMENTS
`python` (version 2.x) <br>
biopython (http://biopython.org) <br>
`R` <br>

###### DEPENDENCIES
`ecamber` - v1.08 <br>
`andi` - v0.9.4 <br>
`opscan` <br>

___
Before starting to install the 3 dependencies, define this environment variable
~~~~
export CGPIPELINE=<ABSOLUTE_PATH_TO_CORE_GENOME_PIPELINE>

#example :
export CGPIPELINE=/home/username/programs/CoreGeneBuilder
~~~~

___

###### INSTALLATION OF ANDI

Installation of a stable version of andi under linux system : version v0.9.4 <br>

Create directory `andi` in the  directory `ext-tools` of CoreGeneBuilder :
~~~~
mkdir $CGPIPELINE/ext-tools/andi
cd $CGPIPELINE/ext-tools/andi
wget https://github.com/EvolBioInf/andi/releases/download/v0.9.4/andi-0.9.4.tar.gz
tar -xf andi-0.9.4.tar.gz
~~~~

###### Option 1 - Install andi with the C/C++ library `libdivsufsort`
Install the C/C++ dependent library
~~~~
cd $CGPIPELINE/ext-tools/andi
wget https://github.com/EvolBioInf/andi/releases/download/v0.9.4/libdivsufsort-2.0.2-1.tar.gz
tar -xf libdivsufsort-2.0.2-1.tar.gz
cd libdivsufsort-2.0.2-1
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE="Release" -DCMAKE_INSTALL_PREFIX="$CGPIPELINE/ext-tools/andi/libdivsufsort-2.0.2-1" ..
make install
#files .c and .so are into $CGPIPELINE/ext-tools/andi/libdivsufsort-2.0.2-1/lib and files .h into $CGPIPELINE/ext-tools/andi/libdivsufsort-2.0.2-1/include
~~~~
Install andi
~~~~
cd $CGPIPELINE/ext-tools/andi/andi-0.9.4
## prepare installation of andi defining path of the precedent library (variable $CGPIPELINE had to be defined)
./configure CPPFLAGS="-I${CGPIPELINE}/ext-tools/andi/libdivsufsort-2.0.2-1/include" CFLAGS="-I${CGPIPELINE}/ext-tools/andi/libdivsufsort-2.0.2-1/include" LDFLAGS="-L${CGPIPELINE}/ext-tools/andi/libdivsufsort-2.0.2-1/lib -Wl,-rpath,${CGPIPELINE}/ext-tools/andi/libdivsufsort-2.0.2-1/lib" --prefix=$CGPIPELINE/ext-tools/andi/andi-0.9.4
make
make install
~~~~

###### Option 2 - Install andi without the C/C++ library `libdivsufsort`
~~~~
cd $CGPIPELINE/ext-tools/andi/andi-0.9.4
./configure --prefix=$CGPIPELINE/ext-tools/andi/andi-0.9.4 --without-libdivsufsort
make
make install
make clean
~~~~

###### Test andi
~~~~
$CGPIPELINE/ext-tools/andi/andi-0.9.4/bin/andi -h
$CGPIPELINE/ext-tools/andi/andi-0.9.4/bin/andi <fasta> 1> result.matrix
~~~~

___

###### INSTALLATION OF ECAMBER DEPENDENCIES

###### Install prodigal
Requirements: <br>
`gcc` <br>
`make` <br>
`install` <br>

~~~~
cd $CGPIPELINE/ext-tools
mkdir prodigal
cd prodigal
wget https://github.com/hyattpd/Prodigal/archive/v2.6.3.tar.gz
tar -xf v2.6.3.tar.gz
cd Prodigal-2.6.3
mkdir bin
make install INSTALLDIR=$CGPIPELINE/ext-tools/prodigal/Prodigal-2.6.3/bin
#quick test :
./bin/prodigal -h
~~~~

###### install blast+
~~~~
mkdir $CGPIPELINE/ext-tools/blast+
cd $CGPIPELINE/ext-tools/blast+

#Download the latest release (binary)
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.4.0+-x64-linux.tar.gz
tar -xf ncbi-blast-2.4.0+-x64-linux.tar.gz
ncbi-blast-2.4.0+/bin/blastn -h
~~~~

###### install muscle
~~~~
mkdir $CGPIPELINE/ext-tools/muscle
mkdir $CGPIPELINE/ext-tools/muscle/3.8.31
cd $CGPIPELINE/ext-tools/muscle/3.8.31

#1a-download binary file only
wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz
tar -xf muscle3.8.31_i86linux64.tar.gz
mv muscle3.8.31_i86linux64 muscle
#display help to quickly test it
./muscle

#OR 1b-install from src files
wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_src.tar.gz
tar -xf muscle3.8.31_src.tar.gz
cd muscle3.8.31/src
make
./muscle
~~~~

###### INSTALLATION OF ECAMBER

Visit this website : <br>
http://bioputer.mimuw.edu.pl/ecamber/software.html

And follow this procedure of installation :
~~~~
mkdir $CGPIPELINE/ext-tools/ecamber
# so download the last stable version of ecamber (v1.08):
wget http://bioputer.mimuw.edu.pl/ecamber/software/ecamber_v1.08/ecamber_src.zip
unzip ecamber_src.zip
~~~~

The next steps are common to an installation starting from archive `ecamber_src.zip` or starting from provided src/bin files.

Ecamber files to configure the path of external tools :
~~~~
$CGPIPELINE/ext-tools/ecamber/config/config_resources.txt # ecamber output directory
$CGPIPELINE/ext-tools/ecamber/config/config_aln_paths.txt # MUSCLE, BLAST+, (PHYLIP)
$CGPIPELINE/ext-tools/ecamber/config/config_extpaths.txt #BLAST+ (blastn, blastp, makeblastdb), PRODIGAL
~~~~

Modify these config files :
~~~~
$ vim ext-tools/ecamber/config/config_resources.txt
#replace
ECAMBER_DATASETS:ECAMBER_PATH:datasets/
#by
ECAMBER_DATASETS:ECAMBER_PATH:../../data/
~~~~
~~~~
$ vim CGPIPELINE/ext-tools/ecamber/config/config_aln_paths.txt

#replace
TBLASTN_PATH:ECAMBER_PATH:ext-tools/tblastn
#by the absolute path of blast binary files, for example :
TBLASTN_PATH:/usr/local/bin/tblastn

#replace
MUSCLE_PATH:ext-tools/muscle#muscle#
#by the absolute path of muscle binary files, for example :
MUSCLE_PATH:/usr/local/bin/muscle#muscle#
~~~~
~~~~
$ vim $CGPIPELINE/ext-tools/ecamber/config/config_extpaths.txt
#replace :
BLASTN_PATH:ECAMBER_PATH:/ext-tools/blastn              # BLAST+ blastn      --- corresponding to blastall blastn
BLASTP_PATH:ECAMBER_PATH:/ext-tools/blastp              # BLAST+ blastp      --- corresponding to blastall blastn
MAKEBLASTDB_PATH:ECAMBER_PATH:/ext-tools/makeblastdb    # BLAST+ makeblastdb --- corresponding to formatdb
PRODIGAL_PATH:ECAMBER_PATH:/ext-tools/prodigal

#by the absolute path of blast and prodigal binary files, for example :
BLASTN_PATH:/usr/local/bin/blastn              # BLAST+ blastn      --- corresponding to blastall blastn
BLASTP_PATH:/usr/local/bin/blastp              # BLAST+ blastp      --- corresponding to blastall blastn
MAKEBLASTDB_PATH:/usr/local/bin/makeblastdb    # BLAST+ makeblastdb --- corresponding to formatdb
PRODIGAL_PATH:/usr/local/bin/prodigal
#note : blastall.exe and formatdb.exe are not used by ecamber
~~~~

Before using ecamber, correct the code of the file `camber_format_utils.py` using these commands:
~~~~
file=ecamber/src/soft/utils/camber_format_utils.py
sed -e 's/left_bound, right_bound = tokens\[1\]\.strip("complement(")\.split("\.\.")/left_bound, right_bound = tokens\[1\]\.strip("complement(")\.strip(")")\.split("..")/' $file > $file.tmp
mv $file.tmp $file
~~~~

Test ecamber
~~~~
cd $CGPIPELINE/ext-tools/ecamber/

#please visit the documentation and download data available here
http://bioputer.mimuw.edu.pl/ecamber/software.html

#then run this kind of command :
python ecamber.py -a pr -d klpn5refannot -w 1
~~~~

___

###### INSTALLATION OF OPSCAN

You need the linux binary file `opscan`.
Test it :
~~~~
cd $CGPIPELINE/src
./opscan -h
~~~~

###### Only for limited users:
If this binary file 'opscan' is not compatible with your system (Linux, MacOS), you need to recompile it.
To do this, you need to have this archive : `opscan_src_files.tar.gz`

~~~~
#compile opscan
cd $CGPIPELINE/src
tar -xvf opscan_src_files.tar.gz
cd opscan_src_files
make clean
make
~~~~

~~~~
#test opscan binary file :
cd opscan_src_files/tests
../bin/opscan -h  #display help if no error

# additional test
cd opscan_src_files/tests
./COMMANDE
diff log.out testout/log.out |wc -l   #if this last command returns '0', installation is successful
~~~~

Copy these files to `src`
~~~~
cp -p ext-tools/CoreGenome/src/opscan_src_files/bin/opscan src/.
cp -p ext-tools/CoreGenome/src/opscan_src_files/data/BLOSUM60 ext-tools/CoreGenome/data/.
~~~~

___

###### CONFIGURE COREGENEBUILDER PIPELINE

Once all dependencies are installed, please configure path of those programs in the configuration file `config_env.txt` :
~~~~
$ vim $CGPIPELINE/config/config_env.txt
## variables to edit :
CGPIPELINE=<path_where_pipeline_is_installed>
ANDI=$CGPIPELINE/ext-tools/andi/andi-0.9.4/bin/andi     #andi binary file
ECAMBER=$CGPIPELINE/ext-tools/ecamber                   #directory where is 'ecamber.py'
COREGENOME=$CGPIPELINE/ext-tools/CoreGenome             #directory where is 'opscan' binary file
OPSCAN_MATRIX=$(readlink -f $COREGENOME/data/BLOSUM60)  #opscan matrix path
~~~~
