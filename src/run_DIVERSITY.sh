#!/bin/bash

#####################################################################################################################
# CoreGeneBuilder - extracts a core genome or a persistent genome from a set of bacterial genomes.                  #
# Authors: Elise Larsonneur, Marie Touchon, Damien Mornico, Alexis Criscuolo, Sylvain Brisse, Eduardo P. C. Rocha   #
# Copyright © 2016 IFB, CNRS, Institut Pasteur                                                                      #
#                                                                                                                   #
# Please read README file for contact information.                                                                  #
#                                                                                                                   #
# This file is part of CoreGeneBuilder.                                                                             #
#                                                                                                                   #
# CoreGeneBuilder is free software: you can redistribute it and/or modify                                           #
# it under the terms of the GNU General Public License as published by                                              #
# the Free Software Foundation, either version 3 of the License, or                                                 #
# (at your option) any later version.                                                                               #
#                                                                                                                   #
# CoreGeneBuilder is distributed in the hope that it will be useful,                                                # 
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                                    # 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                     #
# GNU General Public License for more details.                                                                      #
#                                                                                                                   #
# You should have received a copy of the GNU General Public License                                                 #
# along with CoreGeneBuilder.  If not, see <http://www.gnu.org/licenses/>.                                          #
#                                                                                                                   #
# run_DIVERSITY.sh: computes genomic sequence classification based on estimated genomic evolutionary distances      #
#                   outputed by andi.                                                                               #
# Author: Elise Larsonneur                                                                                          #
#####################################################################################################################


set -o pipefail


########################################################################################
# GLOBAL VARIABLES & SOURCE LIBRARY AND CONFIGURATION FILES                            #
########################################################################################

CGB_CONFIG='../config/'
source "${CGB_BIN}utils.sh"

## parameters initialization
DIRECTORY='N.O.D.I.R'      ## -d
NAME='N.O.N.A.M.E'         ## -n
REFGENOME='N.O.R.E.F'      ## -g
SELECTEDGENOMENB='-1'      ## -s
CONTIGLENGTH=500           ## -f
N50THRESHOLD=0             ## -y
CUTN=10                    ## -N
NBCTGTHRESHOLD='-1'        ## -c
THREADS=1                  ## -t
REFANNOTATION='N.O.R.E.F'  ## -a
REFIDPATTERN='N.O.P.A.T.T.E.R.N'   ## -e

Rexe=R



########################################################################################
# FUNCTIONS                                                                            #
########################################################################################

# display_usage 
# This function displays the usage of this program.
# No parameters
function display_usage { 
  echo '';
  echo 'USAGE :';
  echo "     $0 [options] -d <input_directory> -n <name_of_4_chars>";
  echo '   -d <inDirectory>  directory where are stored input and output files, it must contain at least a directory called 'assemblies' where are genomic sequence fasta files';
  echo '   -n <string>  four letter name (ex : esco (es:escherichia; co:coli))';
  echo "  where 'options' are :";
  echo '   ### PRE-PROCESSING OF GENOMIC SEQUENCES ###';
  echo '   -f <int>     filtering out genomic sequences below this length value (bp) (default : 500); if 0, no filter';
  echo "   -N <int>     cut sequences with mers of 'N' of size '-N' (scaffolds) into contigs (default : 10); if -1, no filter";
  echo '   -y <int>     filtering out genomes having N50 below this N50 value (bp) (default : 0); if 0, no filter';
  echo '   -c <int>     filtering out genomes having number of sequences above this value (default : -1); if -1, no filter';
  echo '   ### DIVERSITY ###';
  echo '   -s <int>     number of genomes to select for core genome construction (default : -1); if -1, no selection (all genomes are kept)';
  echo '   -g <inFile>  reference genome fasta file -- annotation and core genome construction steps';
  echo '   -a <inFile>  reference genome genbank file -- annotation step';
  echo "   -e <inFile>  prefix string of contig ids found in the reference genome fasta and genbank files (eg id in fields LOCUS and ACCESSION of genbank file) (ex : 'NC_'; 'NZ_' ; 'AKAC' ; etc)";
  echo '   -t <int>     number of threads used during diversity and annotation steps (default : 1)';
  echo '   -h           print help';
  echo '';
  echo 'EXAMPLES :';
  echo '     #provided reference genome and related genbank annotation, the functional annotation of reference genome will be transferred to the other genomes:';
  echo "     $0 -d klpn5refannot -n klpn -g MGH78578_NC.fasta -a MGH78578_NC.gb -e NC_ -t 4";
  echo '';
  echo '     #provided reference genome but not provided related genbank annotation, the reference genome will be de novo annotated:';
  echo "     $0 -d klpn5refannot -n klpn -g MGH78578_NC.fasta -t 4";
  echo ''; 
  echo '     #not provided reference genome and related genbank annotation, the reference genome will be the first on genome list sorted alphabetically:';
  echo "     $0 -d klpn5refannot -n klpn -t 4";  
} 



# replace_names_in_matrix
# This function replaces the short names of fasta files by their original name in the ouputed andi distance matrix.
# Parameters
# - 1) Input directory - A string containing the directory path where the fasta files are.
# - 2) Input file - A string containing the file path where the short and long names of fasta files are stored.
# - 3) Input matrix - A string containing the matrix file path outputed by andi.
# - 4) Output file - A string containing the file path where the long and original names of fasta files are stored.
# - 5) Output matrix - A string containing the output matrix file path where the names of fasta files have been replaced by their old names.
# Outputs 
# - The output matrix where short names of fasta files have been replaced by their original names.
function replace_names_in_matrix {

  local dir_ inlist_ inmatrix_ outlist_ outmatrix_
  local genome oldname linkname linkprefix genomeprefix
  if [[ -z "$1" ]]; then echo '[ERROR] the function replace_names_in_matrix expects an input directory to run' 2>>"${LOG}"; exit 1; fi
  if [[ -z "$2" ]]; then echo '[ERROR] the function replace_names_in_matrix expects an input file to run' 2>>"${LOG}"; exit 1; fi
  if [[ -z "$3" ]]; then echo '[ERROR] the function replace_names_in_matrix expects an input matrix file to run' 2>>"${LOG}"; exit 1; fi
  if [[ -z "$4" ]]; then echo '[ERROR] the function replace_names_in_matrix expects an output file name to run' 2>>"${LOG}"; exit 1; fi
  if [[ -z "$5" ]]; then echo '[ERROR] the function replace_names_in_matrix expects an output matrix name to run' 2>>"${LOG}"; exit 1; fi
  if [[ $# -ne 5 ]]; then echo "[ERROR] the function replace_names_in_matrix expects 5 parameters, $# parameters given" 2>>"${LOG}"; exit 1; fi

  dir_="$1"
  inlist_="$2"  ## renamed_genomes_list.txt
  inmatrix_="$3"
  outlist_="$4"
  outmatrix_="$5"

  if [[ ! -e "${dir_}" ]]; then echo "[ERROR] the input directory ${dir_} does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${dir_}" ]]; then echo "[ERROR] the input directory ${dir_} is empty" 2>>"${LOG}"; exit 1; fi
  if [[ ! -e "${inlist_}" ]]; then echo "[ERROR] the input file ${inlist_} does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${inlist_}" ]]; then echo "[ERROR] the input file ${inlist_} is empty" 2>>"${LOG}"; exit 1; fi
  if [[ ! -e "${inmatrix_}" ]]; then echo "[ERROR] the input file ${inmatrix_} does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${inmatrix_}" ]]; then echo "[ERROR] the input file ${inmatrix_} is empty" 2>>"${LOG}"; exit 1; fi
  if [[ -e "${outlist_}" ]]; then rm "${outlist_}"; fi
  # sedcmd="cat $inmatrix_ |sed '"

  for fasta in ${dir_}/diversity/genomes/*; do
    if [[ -L "${fasta}" ]]; then   #if it is a symbolic link
      genome="$(readlink -f "${fasta}")"
      oldname="$(get_old_name_from_list "${inlist_}" "${genome}")"
      linkname="$(basename "${fasta}")"
      linkprefix="${linkname%%.*}"
      genomeprefix="$(basename $oldname .fasta)"
      # sedcmd=$sedcmd"s/\\b$linkprefix\\b/$genomeprefix/; "
      echo -e "${linkprefix}\t${genomeprefix}" >> "${outlist_}"
    fi
  done
  awk 'BEGIN{FS="\t"}NR==FNR{table[$1]=$2}NR>FNR && FNR==1{print $0;FS=" "}NR>FNR && FNR>1{$1=table[$1];print $0}' "${outlist_}" "${inmatrix_}" 1>"${outmatrix_}" 2>>"${LOG}"
  # sedcmd=$sedcmd"' 1> $outmatrix_ 2>>"${LOG}""
  # $(eval $sedcmd)  #error when a big number of genomes : 'error : /bin/sed: Argument list too long' 
  if [[ ! -e "${outlist_}" ]]; then echo "[ERROR] the output file ${outlist_} does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${outlist_}" ]]; then echo "[ERROR] the output file ${outlist_} is empty" 2>>"${LOG}"; exit 1; fi
  if [[ ! -e "${outmatrix_}" ]]; then echo "[ERROR] the output file ${outmatrix_} does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${outmatrix_}" ]]; then echo "[ERROR] the output file ${outmatrix_} is empty" 2>>"${LOG}"; exit 1; fi 
}



# remove_ref_from_matrix
# This function removes reference genome data from distance matrix.
# Parameters
# - 1) Input matrix - A string containing the matrix file path.
# - 2) A string containing the name of the reference genome.
# - 3) Output matrix - A string containing the output matrix file path where the reference genome data have been removed.
# Outputs
# - The output distance matrix where reference genome data have been removed.
function remove_ref_from_matrix {

  local matrix searchstr outmatrix linenb 
  if [[ -z "$1" ]]; then echo '[ERROR] the function remove_ref_from_matrix expects an input matrix file to run' 2>>"${LOG}"; exit 1; fi
  if [[ -z "$2" ]]; then echo '[ERROR] the function remove_ref_from_matrix expects an input string to run' 2>>"${LOG}"; exit 1; fi
  if [[ -z "$3" ]]; then echo '[ERROR] the function remove_ref_from_matrix expects an output matrix file path to run' 2>>"${LOG}"; exit 1; fi
  if [[ $# -ne 3 ]]; then echo "[ERROR] the function remove_ref_from_matrix expects 3 parameters, $# parameters given" 2>>"${LOG}"; exit 1; fi
  matrix="$1"
  searchstr="$2"
  outmatrix="$3"
  if [[ ! -e "${matrix}" ]]; then echo "[ERROR] the input file ${matrix} does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${matrix}" ]]; then echo "[ERROR] the input file ${matrix} is empty" 2>>"${LOG}"; exit 1; fi 

  awk -F 'FS' 'BEGIN{FS=" "} { print $1 }' "${matrix}" 1> "${matrix}.rr" 2>>"${LOG}"
  if grep -Fxq "${searchstr}" "${matrix}.rr" ; then
    # get line number of the search string
    linenb="$(grep -Fxn "${searchstr}" "${matrix}.rr" |cut -f1 -d:)"

    #delete the 'linenb'th column
    awk -v nb="${linenb}" -F 'FS' 'BEGIN{FS=" "}{for (i=1; i<=NF; i++) if(i!=nb)  {printf $i FS} };{printf "\n"}' "${matrix}" 1> "${outmatrix}" 2>>"${LOG}"

    #delete the 'linenb'th raw
    sed "${linenb}d" "${outmatrix}" 1> "${outmatrix}.tmp" 2>>"${LOG}"
    mv "${outmatrix}.tmp" "${outmatrix}" >> "${LOG}" 2>&1
  fi 
  if [[ -e "${matrix}.rr" ]]; then rm "${matrix}.rr" >> "${LOG}" 2>&1 ; fi
  if [[ ! -e "${outmatrix}" ]]; then echo "[ERROR] the output file ${outmatrix} does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${outmatrix}" ]]; then echo "[ERROR] the output file ${outmatrix} is empty" 2>>"${LOG}"; exit 1; fi  
}



# get_phylotree_and_clustering
# This function computes a hierarchical clustering (average linkage = UPGMA) from the input distance matrix and cuts it into clusters using the given threshold cluster number.
# Parameters :
# - 1) Input file - A string containing the distance matrix file path.
# - 2) Ouput file - A string containing the pdf file path where the hierarchical clustering and the cutted hierarchical clustering are plotted.
# - 3) An integer - The cluster number threshold.
# - 4) Output file -  A string containing the file path where the list of cluster numbers and associated fasta files is recorded.
# Outputs
# - The output pdf file
# - The output cluster list file.
function get_phylotree_and_clustering {

  local distance_matrix_ plotpdf_ clusternb_ genomeclusterlist_ Rcommand
  if [[ -z "$1" ]]; then echo '[ERROR] the function get_phylotree_and_clustering expects an input matrix file to run' 2>>"${LOG}"; exit 1; fi
  if [[ -z "$2" ]]; then echo '[ERROR] the function get_phylotree_and_clustering expects an output file name (pdf) to run' 2>>"${LOG}"; exit 1; fi
  if [[ -z "$3" ]]; then echo '[ERROR] the function get_phylotree_and_clustering expects an integer to run' 2>>"${LOG}"; exit 1; fi
  if [[ -z "$4" ]]; then echo '[ERROR] the function get_phylotree_and_clustering expects an output file name to run' 2>>"${LOG}"; exit 1; fi
  if [[ $# -ne 4 ]]; then echo "[ERROR] the function <functionName> expects 4 parameters, $# parameters given" 2>>"${LOG}"; exit 1; fi
  distance_matrix_="$1"
  plotpdf_="$2"
  clusternb_="$3"
  genomeclusterlist_="$4"
  if [[ ! -e "${distance_matrix_}" ]]; then echo "[ERROR] the input file ${distance_matrix_} does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${distance_matrix_}" ]]; then echo "[ERROR] the input file ${distance_matrix_} is empty" 2>>"${LOG}"; exit 1; fi

  Rcommand='';
  ### require R package 'stats'
    
  # function 'dropNA' that removes matrix rows with NaN values 
  Rcommand="${Rcommand} dropNA <- function(m) {";
  #count the number of NaN values in each raw (it's the same result for columns because it's a squared matrix)
  Rcommand="${Rcommand} na_rows <- apply(m,1,function(x){sum(x=='NaN')});";
  # while some NAN values exist, remove them
  Rcommand="${Rcommand} while (max(na_rows)>0) {";
  # get column index where number of NaN values is equal to the maximal number of NaN values
  Rcommand="${Rcommand} na_rows.drop <- c(1:length(na_rows))[na_rows==max(na_rows)];";
  # remove the rows related to these indexes
  Rcommand="${Rcommand} m <- m[-na_rows.drop,-na_rows.drop];";
  Rcommand="${Rcommand} na_rows <- apply(m,1,function(x){sum(x=='NaN')});";
  Rcommand="${Rcommand} };";
  Rcommand="${Rcommand} m;";  #return matrix m
  Rcommand="${Rcommand} };";
    
  Rcommand="${Rcommand} d<-read.table('${distance_matrix_}', header=F, skip=1, row.names=1);";
  Rcommand="${Rcommand} colnames(d) <- rownames(d);";
  Rcommand="${Rcommand} d <- dropNA(d);"; ## call of function dropNA
  Rcommand="${Rcommand} distances <- as.dist(d);";                      ## get distance matrix
  Rcommand="${Rcommand} ngenomes <- nrow(d);";
  Rcommand="${Rcommand} scale <- ceiling(ngenomes * 0.25);";
  Rcommand="${Rcommand} if (scale < 7) scale <- 7;";
  Rcommand="${Rcommand} hc =  hclust(distances, method='average');";    ## average linkage clustering (UPGMA)
  Rcommand="${Rcommand} pdf(file='${plotpdf_}', height=scale, width=scale);";                        ## open pdf file
  Rcommand="${Rcommand} plot(hc, main='Average linkage clustering (UPGMA)', hang=-1);";					   ## first plot : tree
  Rcommand="${Rcommand} plot(hc, main='Average linkage clustering (UPGMA)', hang=-1);";					   ## to ensure this will be plot, it is repeated...
  ## WARNING !!!! check that (genome_nb > cluster_nb) AND (cluster_nb >= 2 ); if not : R function 'cutree' and 'rect.hclust' fail.
  Rcommand="${Rcommand} groups <- cutree(hc, k=${clusternb_});";          ## cut tree into n clusters (n=number of genomes to select) 
  Rcommand="${Rcommand} rect.hclust(hc, k=${clusternb_}, border='red');"; ## 2nd plot : clustering result
                                                                      ## write (genome cluster) list to file
  Rcommand="${Rcommand} write.table(groups,file='${genomeclusterlist_}', sep='\t', col.names=FALSE, row.names=TRUE);";
  Rcommand="${Rcommand} graphics.off();";                               ## close pdf file

  echo -e "\n${Rcommand}" >> "${LOG}" 2>&1 ;
  echo "${Rcommand}" | $Rexe --slave >> "${LOG}" 2>&1 ;

  if [[ ! -e "${genomeclusterlist_}" ]]; then echo "[ERROR] the input file ${genomeclusterlist_} does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${genomeclusterlist_}" ]]; then echo "[ERROR] the input file ${genomeclusterlist_} is empty" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${plotpdf_}" ]]; then echo "[ERROR] the output file ${plotpdf_} is empty" 2>>"${LOG}"; fi

  ### example of output '$genomeclusterlist'
  ### <genome_name> <cluster_no>
  # "lm001" 1
  # "lm006" 2
  # "lm010" 2
}



# get_phylotree
# This function computes a hierarchical clustering (average linkage = UPGMA) from the input distance matrix .
# Parameters :
# - 1) Input file - A string containing the distance matrix file path.
# - 2) Ouput file - A string containing the pdf file path where the hierarchical clustering and the cutted hierarchical clustering are plotted.
# Outputs
# - The output pdf file.
function get_phylotree {

  local distance_matrix_ plotpdf_ Rcommand
  if [[ -z "$1" ]]; then echo '[ERROR] the function get_phylotree_and_clustering expects an input matrix file to run' 2>>"${LOG}"; exit 1; fi
  if [[ -z "$2" ]]; then echo '[ERROR] the function get_phylotree_and_clustering expects an output file name (pdf) to run' 2>>"${LOG}"; exit 1; fi
  if [[ $# -ne 2 ]]; then echo "[ERROR] the function <functionName> expects 2 parameters, $# parameters given" 2>>"${LOG}"; exit 1; fi
  distance_matrix_="$1"
  plotpdf_="$2"
  if [[ ! -e "${distance_matrix_}" ]]; then echo "[ERROR] the input file ${distance_matrix_} does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${distance_matrix_}" ]]; then echo "[ERROR] the input file ${distance_matrix_} is empty" 2>>"${LOG}"; exit 1; fi

  Rcommand='';
  ### require R package 'stats'
   
  # function 'dropNA' that removes matrix rows with NaN values 
  Rcommand="${Rcommand} dropNA <- function(m) {";
  #count the number of NaN values in each raw (it's the same result for columns because it's a squared matrix)
  Rcommand="${Rcommand} na_rows <- apply(m,1,function(x){sum(x=='NaN')});";
  # while some NAN values exist, remove them
  Rcommand="${Rcommand} while (max(na_rows)>0) {";
  # get column index where number of NaN values is equal to the maximal number of NaN values
  Rcommand="${Rcommand} na_rows.drop <- c(1:length(na_rows))[na_rows==max(na_rows)];";
  # remove the rows related to these indexes
  Rcommand="${Rcommand} m <- m[-na_rows.drop,-na_rows.drop];";
  Rcommand="${Rcommand} na_rows <- apply(m,1,function(x){sum(x=='NaN')});";
  Rcommand="${Rcommand} };";
  Rcommand="${Rcommand} m;";  #return matrix m
  Rcommand="${Rcommand} };";
   
  Rcommand="${Rcommand} d<-read.table('${distance_matrix_}', header=F, skip=1, row.names=1);";
  Rcommand="${Rcommand} colnames(d) <- rownames(d);";
  Rcommand="${Rcommand} d <- dropNA(d);"; ## call of function dropNA
  Rcommand="${Rcommand} distances <- as.dist(d);";                      ## get distance matrix
  Rcommand="${Rcommand} ngenomes <- nrow(d);";
  Rcommand="${Rcommand} scale <- ceiling(ngenomes * 0.25);";
  Rcommand="${Rcommand} if (scale < 7) scale <- 7;";
  Rcommand="${Rcommand} hc =  hclust(distances, method='average');";    ## average linkage clustering (UPGMA)
  Rcommand="${Rcommand} pdf(file='${plotpdf_}', height=scale, width=scale);";                        ## open pdf file
  Rcommand="${Rcommand} plot(hc, main='Average linkage clustering (UPGMA)', hang=-1);";					   ## first plot : tree
  Rcommand="${Rcommand} graphics.off();";                               ## close pdf file

  echo -e "\n${Rcommand}" >> "${LOG}" 2>&1 ;
  echo "${Rcommand}" | $Rexe --slave >> "${LOG}" 2>&1 ;
  if [[ ! -s "${plotpdf_}" ]]; then echo "[ERROR] the output file ${plotpdf_} is empty" 2>>"${LOG}"; fi
}



# get_assembly_stats
# This function gets assembly metrics (N50 contig length, genome size) of each fasta file in the list of clusters (file \t cluster).
# Parameters
# - 1) Input file - A string containing the file path where the list of clusters and associated files is recorded.
# Outputs
# - A file containing informations of the input file AND the N50 contig length and genome size values.
function get_assembly_stats {
 
  local genomeclusterlist__ toout f fasta cluster contiginfos n50 genomeSize nbcontigs toout
  if [[ -z "$1" ]]; then echo '[ERROR] the function get_assembly_stats expects an input file to run' 2>>"${LOG}"; exit 1; fi
  if [[ $# -ne 1 ]]; then echo "[ERROR] the function get_assembly_stats expects 1 parameter, $# parameters given" 2>>"${LOG}"; exit 1; fi
  genomeclusterlist__="$1"
  if [[ ! -e "${genomeclusterlist__}" ]]; then echo "[ERROR] the input file ${genomeclusterlist__} does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${genomeclusterlist__}" ]]; then echo "[ERROR] the input file ${genomeclusterlist__} is empty" 2>>"${LOG}"; exit 1; fi

  while read line; do
    toout="${line}";
    f="$(echo "${line}" | sed 's/"//g' | awk '{print $1}')";
    fasta="$(ls "${GENOME_DIR}/${f}."*.fasta)"
    cluster="$(echo "${line}" | sed 's/"//g' | awk  '{print $2}')";
    contiginfos="$(${CGB_BIN}contig_info.sh "${fasta}")";
    exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] the command '${CGB_BIN}contig_info.sh ${fasta}' exited with an error : ${contiginfos} - exit code ${exit_status}" 2>>"${LOG}"; exit 1; fi
    n50="$(echo "${contiginfos}" | grep 'N50' | awk '{print $2}')";
    genomeSize="$(echo "${contiginfos}" | grep 'Total' | awk '{print $2}')";
    nbcontigs="$(echo "${contiginfos}" | grep 'Number' | awk '{print $4}')";
    toout="${fasta} ${cluster} ${n50} ${genomeSize} ${nbcontigs}";
    echo "${toout}" >> "${genomeclusterlist__}.infos" ;   ## add n50 and genome size stats for each (genome, cluster_no)
  done < ${genomeclusterlist__}

  if [[ ! -e "${genomeclusterlist__}.infos" ]]; then echo "[ERROR] the output file ${genomeclusterlist__}.infos does not exist" 2>>"${LOG}"; fi
  if [[ ! -s "${genomeclusterlist__}.infos" ]]; then echo "[ERROR] the output file ${genomeclusterlist__}.infos is empty" 2>>"${LOG}"; fi 
}



# extract_selected_genomes 
# This function creates the list of selected genomes used as input for the next steps (annotation and core genome building).
# Parameters
# - 1) Input file - A string containing the file path where the list of clusters is stored.
# - 2) The number of genomes to select (number of clusters) - An integer.
# Outputs
# - A file containing the list of selected genomes.
function extract_selected_genomes {

  local genomeclusterlist_ nbofgenometoselect_ clusters 
  if [[ -z "$1" ]]; then echo '[ERROR] the function extract_selected_genomes expects an input file (cluster list) to run' 2>>"${LOG}"; exit 1; fi
  if [[ -z "$2" ]]; then echo '[ERROR] the function extract_selected_genomes expects an integer to run' 2>>"${LOG}"; exit 1; fi
  if [[ $# -ne 2 ]]; then echo "[ERROR] the function extract_selected_genomes expects 2 parameters, $# parameters given" 2>>"${LOG}"; exit 1; fi
  genomeclusterlist_="$1"
  nbofgenometoselect_="$2"
  if [[ ! -e "${genomeclusterlist_}.infos" ]]; then echo "[ERROR] the input file ${genomeclusterlist_}.infos does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${genomeclusterlist_}.infos" ]]; then echo "[ERROR] the input file ${genomeclusterlist_}.infos is empty" 2>>"${LOG}"; exit 1; fi

  ## sort genomes by cluster (2nd field) - increasing order
  ## and then sort genomes of each cluster by N50 (3rd field) - decreasing order
  sort -k2,2n -k3,3nr "${genomeclusterlist_}.infos" > "${genomeclusterlist_}.infos.sort"
 
  if [[ ! -e "${genomeclusterlist_}.infos.sort" ]]; then echo "[ERROR] the input file ${genomeclusterlist_}.infos.sort does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${genomeclusterlist_}.infos.sort" ]]; then echo "[ERROR] the input file ${genomeclusterlist_}.infos.sort is empty" 2>>"${LOG}"; exit 1; fi
  clusters="$(cat "${genomeclusterlist_}.infos.sort" |awk '{print $2}' |sort -nu)" ;
 
  # keep one genome per cluster
  # and then if number_genomes > number_clusters, catch one genome of each cluster until value equal to number_genomes
  # if  number_genomes < number_clusters, we conserve all genomes (skip steps N50, only do plot with andi)
  # => principle of script : 'extract_genomes_from_clusters.sh'

  # extract_genomes_from_clusters.sh <list_genomes_clusters_N50_sortedInputFile> <nb_genomes_to_extract> <fileOut>
  echo "SELECTED_GENOMES_NUMBER: ${nbofgenometoselect_}" >> "${LOG}" 2>&1
  echo "SELECTED_GENOMES_NUMBER: ${nbofgenometoselect_}"
  echo "${CGB_BIN}extract_genomes_from_clusters.sh ${genomeclusterlist_}.infos.sort ${nbofgenometoselect_} ${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt.tmp" >> "${LOG}" 2>&1
  ${CGB_BIN}extract_genomes_from_clusters.sh "${genomeclusterlist_}.infos.sort" "${nbofgenometoselect_}" "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt.tmp"
  exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] the command '${CGB_BIN}extract_genomes_from_clusters.sh ${genomeclusterlist_}.infos.sort ${nbofgenometoselect_} ${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt.tmp' exited with an error! - exit code ${exit_status}" 2>>"${LOG}"; exit 1; fi 
  awk '{print $1}' "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt.tmp" > "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt"
  if [[ ! -e "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt" ]]; then echo "[ERROR] the output file ${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt" ]]; then echo "[ERROR] the output file ${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt is empty" 2>>"${LOG}"; exit 1; fi   
  if [[ -e "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt.tmp" ]]; then rm "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt.tmp"; fi
}



########################################################################################
# MAIN FUNCTION                                                                        #
########################################################################################
# main
# Parameters : See 'getopts' part.
function main {

  # check whether user had supplied -h or --help . If yes display usage 
  if [[ "$1" = "-?" ]] || [[ "$1" = "-h" ]] || [[ "$1" = "--help" ]]; then 
    display_usage
    exit 0;
  fi  

  # if less than two arguments supplied, display usage 
  if [[  $# -le 1 ]]; then 
    display_usage
    exit 1;
  fi  
 

  ## catch option values
  while getopts :C:d:n:g:a:e:s:f:y:N:c:t: option
  do
    if [[ -z "${OPTARG}" ]]; then echo "[ERROR] empty argument for option -${option}"; exit 1; fi
    case "${option}" in
      C)  
        CGB_CONFIG="${OPTARG}";    
        if [[ ! -d "${CGB_CONFIG}" ]]; then echo "[ERROR] input directory '${CGB_CONFIG}' does not exist (option -C)." ; exit 1 ; else source ${CGB_CONFIG}/config_env.txt; fi  
        ;; # -C <inConfigDirectory>
      d)  
        DIRECTORY="${OPTARG}";
        if [[ ! -d "${DATA}/${DIRECTORY}" ]]; then echo "[ERROR] input directory '${DATA}/${DIRECTORY}' does not exist (option -d)." ; exit 1 ; fi  
        ;; # -d <inDirectory>
      n)  
        NAME="${OPTARG}";
        size=${#NAME}; 
        if [[ $size -ne 4 ]]; then echo '[ERROR] name of 4 letters is required (option -n).' ; exit 1 ; fi  
        ;; # -n <name (4 letters)>
      g)    
        REFGENOME="${OPTARG}";
        ;;     # -g <reference genome FASTA infile>
      a)    
        REFANNOTATION="${OPTARG}";
        ;; # -a <annotation linked to reference genome GENBANK infile>
      e)    
        REFIDPATTERN="${OPTARG}";
        ;;  # -e <prefix of reference contig id>
      s)  
        SELECTEDGENOMENB="${OPTARG}";
        if ! [[ "${SELECTEDGENOMENB}" =~ ^[0-9]+$ ]] && [[ "${SELECTEDGENOMENB}" != '-1' ]]; then echo '[ERROR] Bad value for option -s' ; exit 1 ; fi; 
        if [[ "${SELECTEDGENOMENB}" != '-1' ]] && [[ $SELECTEDGENOMENB -lt 2 ]]; then
          echo '[ERROR] the number of genomes to select for core genome construction must be greater than 1 (option -s).' ; exit 1 ; 
        fi  
        ;; # -s <number of selected genomes threshold>
      f)  
        CONTIGLENGTH="${OPTARG}";
        if ! [[ "${CONTIGLENGTH}" =~ ^[0-9]+$ ]] || [[ $CONTIGLENGTH -lt 1 ]]; then echo '[ERROR] the sequence length threshold  must be greater than 0 (option -f).' ; exit 1 ; fi
        ;; # -f <contig length threshold> 
      y)  
        N50THRESHOLD="${OPTARG}";
        if ! [[ "${N50THRESHOLD}" =~ ^[0-9]+$ ]]; then echo '[ERROR] Bad value for option -y'; exit 1; fi
        ;; # -y <N50 length threshold> 
      N)  
        CUTN="${OPTARG}";
        if ! [[ "${CUTN}" =~ ^[0-9]+$ ]] && [[ "${CUTN}" != "-1" ]]; then echo '[ERROR] Bad value for option -N'; exit 1; fi; 
        if [[ $CUTN -eq 0 ]]; then echo '[ERROR] Bad value for option -N'; exit 1; fi
        ;;  
      c)  
        NBCTGTHRESHOLD="${OPTARG}";
        if ! [[ "${NBCTGTHRESHOLD}" =~ ^[0-9]+$ ]] && [[ "${NBCTGTHRESHOLD}" != "-1" ]]; then echo '[ERROR] Bad value for option -c'; exit 1; fi
        if [[ "${NBCTGTHRESHOLD}" != "-1" ]] && [[ $NBCTGTHRESHOLD -le 0 ]]; then echo '[ERROR] the sequence number threshold must be greater than 0 (option -c).' ; exit 1 ; fi 
        ;; # -c <nb_contig threshold> 
      t)  
        THREADS="${OPTARG}";
        nproc="$(nproc)";
        if ! [[ "${THREADS}" =~ ^[0-9]+$ ]] || [[ $THREADS -lt 1 ]]; then echo '[ERROR] the number of threads must be greater than 0 (option -t).' ; exit 1 ; fi; 
        if [[ $THREADS -gt $nproc ]]; then echo "[ERROR] too much threads requested, use ${nproc} thread(s) instead"; exit 1; fi
        ;; # -t <number of threads> 
      :)  
        echo "[ERROR] option ${OPTARG} : missing argument" ; exit 1   
        ;;  
      \?) 
        echo "[ERROR] ${OPTARG} : option invalide" ; exit 1   
        ;;  
    esac
  done

  readonly CGB_CONFIG
  readonly DIRECTORY NAME REFGENOME SELECTEDGENOMENB CONTIGLENGTH N50THRESHOLD CUTN NBCTGTHRESHOLD
  readonly THREADS REFANNOTATION REFIDPATTERN



  ### checking input directory
  if [[ "${DIRECTORY}" = 'N.O.D.I.R' ]]; then echo '[ERROR] no input directory supplied (mandatory option -d)' ; exit 1 ; fi
  if [[ ! -e "${DATA}/${DIRECTORY}" ]]; then echo "[ERROR] directory ${DATA}/${DIRECTORY} does not exist, please create it." ; exit 1 ; fi


  ### create log
  mkdir -p "${DATA}/${DIRECTORY}/logs"
  LOG="${DATA}/${DIRECTORY}/logs/${DIRECTORY}.diversity.2.log";
  readonly LOG;
  echo '' > "${LOG}";
  

  ### checking name 
  if [[ "${NAME}" = 'N.O.N.A.M.E' ]]; then echo '[ERROR] no name supplied (mandatory option -n)' ; exit 1 ; fi

  ### checking if directory containing genomes (assemblies) exists and is not empty
  if [[ ! -e "${DATA}/${DIRECTORY}/assemblies" ]]; then echo "[ERROR] directory ${DATA}/${DIRECTORY}/assemblies does not exist, please create it." ; exit 1 ; fi
  if [[ ! "$(ls -A ${DATA}/${DIRECTORY}/assemblies)" ]]; then echo "[ERROR] directory ${DATA}/${DIRECTORY}/assemblies is empty, please add some genome fasta files into it" ; exit 1 ; fi

  ### checking if reference genome is supplied
  if [[ "${REFGENOME}" = 'N.O.R.E.F' ]]; then echo "Reference genome will be the first fasta file appearing in directory ${DIRECTORY}." ; fi
  if [[ "${REFGENOME}" != 'N.O.R.E.F' ]] && [[ ! -e "${DATA}/${DIRECTORY}/assemblies/${REFGENOME}" ]]; then echo "[ERROR] input fasta file '${DATA}/${DIRECTORY}/assemblies/${REFGENOME}' does not exist (option -g)." ; exit 1 ; fi
  if [[ "${REFGENOME}" != 'N.O.R.E.F' ]] && [[ ! -s "${DATA}/${DIRECTORY}/assemblies/${REFGENOME}" ]]; then echo "[ERROR] input fasta file '${DATA}/${DIRECTORY}/assemblies/${REFGENOME}' is empty (option -g)." ; exit 1 ; fi

  ## checking if reference genbank file is supplied (annotation)
  if [[ "${REFANNOTATION}" = 'N.O.R.E.F' ]]; then echo 'No reference genome annotation provided (genbank file).' ; fi
  if [[ "${REFANNOTATION}" != 'N.O.R.E.F' ]] && [[ ! -e "${DATA}/${DIRECTORY}/ref_gbk_annotation" ]]; then echo "[ERROR] directory containing input genbank file, '${DATA}/${DIRECTORY}/ref_gbk_annotation' does not exist." ; exit 1 ; fi
  if [[ "${REFANNOTATION}" != 'N.O.R.E.F' ]] && [[ ! -e "${DATA}/${DIRECTORY}/ref_gbk_annotation/${REFANNOTATION}" ]]; then echo "[ERROR] input genbank file '${DATA}/${DIRECTORY}/ref_gbk_annotation/${REFANNOTATION}' does not exist (option -a)." ; exit 1 ; fi
  if [[ "${REFANNOTATION}" != 'N.O.R.E.F' ]] && [[ ! -s "${DATA}/${DIRECTORY}/ref_gbk_annotation/${REFANNOTATION}" ]]; then echo "[ERROR] input genbank file '${DATA}/${DIRECTORY}/ref_gbk_annotation/${REFANNOTATION}' is empty (option -a)." ; exit 1 ; fi
  if [[ "${REFANNOTATION}" != 'N.O.R.E.F' ]] && [[ "${REFIDPATTERN}" = 'N.O.P.A.T.T.E.R.N' ]]; then 
    echo "[ERROR] prefix (only letters) of genome sequence ids found in reference fasta and genbank files REQUIRED (option -e) if reference genbank file is supplied as reference annotation (option -a). Examples : 'NC_', 'NZ_', 'AKAC'." ; exit 1 ; 
  fi  

  ## checking if the step of selection of a subset of genome sequences is activated
  if [[ "${SELECTEDGENOMENB}" = '-1' ]]; then echo 'Skipping step of genome selection.' ; fi

  ## checking if filter of sequences based on their length is not activated
  if [[ $CONTIGLENGTH -eq 0 ]]; then echo 'Skipping step of filtering out sequences based on their length.' ; fi
  if [[ $N50THRESHOLD -lt 1 ]]; then echo 'Skipping step of filtering out sequences based on the genome N50 length.' ; fi
  if [[ $NBCTGTHRESHOLD -lt 0 ]]; then echo 'Skipping step of filtering out sequences based on the genome sequence number.' ; fi
  if [[ $CUTN -gt 0 ]]; then 
    echo 'We will cut scaffolds into contigs (except for the provided reference genome sequence).'; 
  elif [[ "${CUTN}" == "-1" ]]; then 
    echo 'Skipping step of cutting scaffolds into contigs.'; 
  else 
    echo '[ERROR] Bad value for option -N'; exit 1;  
  fi  




  ## prepare input and output directories
  GENOME_DIR=$(real_path_dir "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/genomes")
  readonly GENOME_DIR


  ## 1 - run andi on all genomes
  echo 'run andi to get genome distance matrix' >> "${LOG}" 2>&1 ;
  echo 'run andi to get genome distance matrix';
  ## 1-a
  genome_number=0

  genome_list=$(for f in $(ls "${GENOME_DIR}/"*.fasta); do \
    if [[ -h "${f}" ]]; then \
      fasta="$(basename "${f}")"; \
      suffix="${fasta#$NAME.}"; \
      echo -e "${GENOME_DIR}/${suffix} ";\
    fi
  done)

  for f in $(ls "${GENOME_DIR}/"*.fasta); do
    if [[ -h "${f}" ]]; then 
      genome_number=$((${genome_number} + 1))
    fi
  done

  distance_matrix="${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/${NAME}.${genome_number}.mat" 
  andi_version="$($ANDI --version |head -n 1)"
  echo "andi version :${andi_version}"
  echo "andi version :${andi_version}" >> "${LOG}" 2>&1 ;
  echo "${ANDI} -t ${THREADS} -j ${genome_list} > ${distance_matrix}"   ## resulting matrix
  $ANDI -t "${THREADS}" -j ${genome_list} 1> "${distance_matrix}" 2>>"${LOG}"  ## resulting matrix
  exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] the program '${ANDI}' exited with an error! - exit code ${exit_status} -- see previous line to get the full command line" 2>>"${LOG}"; exit 1; fi

  echo "andi output genome distance matrix : ${distance_matrix}" >> "${LOG}" 2>&1 ;
  echo "andi output genome distance matrix : ${distance_matrix}";

  replace_names_in_matrix "${DATA}/${DIRECTORY}" "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/renamed_genomes_list.txt" "${distance_matrix}" "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/renamed_genomes_list.links.txt" "${distance_matrix}.renamed"



  ## 1-b - the same but if reference fasta name is provided

  genome_number_wo_ref=0
  genome_to_select_wo_ref=$((${SELECTEDGENOMENB} - 1))
  genome_list_wo_ref="$(echo "${genome_list}")"      # list without reference
  distance_matrix_wo_ref=''

  if [[ "${REFGENOME}" != 'N.O.R.E.F' ]]; then

    genome_number_wo_ref=$((${genome_number} - 1))

    # get link name of reference genome
    refgenomeprefix="${REFGENOME%.*}"  #file without extension 
    newname="$(get_new_name_from_list "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/renamed_genomes_list.txt" "${refgenomeprefix}.fasta")"
    # echo -e "NEWNAME_REFNAME:${newname}" >> "${LOG}" 2>&1
    fasta="$(basename "${newname}")"
    suffix="${fasta#$NAME.}"
    rename="${NAME:0:1}${NAME:2:1}"
    reflinkname="$(real_path_dir "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/genomes")"
    reflinkname="${reflinkname}/${rename}${suffix}"    ## ex: klpn2/diversity/genomes/kp002.c01.00.fasta
    reflink="${rename}${suffix}"
    reflink="${reflink%%.*}"                          ## ex: kp002
    # echo "REFLINKNAME:${reflinkname}" >> "${LOG}" 2>&1

    genome_list_wo_ref=$(echo "${genome_list/$reflinkname/}")

    distance_matrix_wo_ref="${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/${NAME}.${genome_number_wo_ref}.woref.mat" 
     
    # compute matrix without reference
    ## $ANDI -t $THREADS -j $genome_list_wo_ref > $distance_matrix_wo_ref   ## resulting matrix
    ## echo "$ANDI -t $THREADS -j $genome_list_wo_ref > $distance_matrix_wo_ref"   ## resulting matrix
    echo 'get matrix without the reference' >> "${LOG}" 2>&1
    echo "remove_ref_from_matrix ${distance_matrix} ${reflink} ${distance_matrix_wo_ref}" >> "${LOG}" 2>&1 ;
    remove_ref_from_matrix "${distance_matrix}" "${reflink}" "${distance_matrix_wo_ref}"
    echo 'get matrix without the reference but with original isolate names' >> "${LOG}" 2>&1     
    echo "replace_names_in_matrix ${DATA}/${DIRECTORY} ${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/renamed_genomes_list.txt ${distance_matrix_wo_ref} ${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/renamed_genomes_list.woref.links.txt ${distance_matrix_wo_ref}.renamed" >> "${LOG}" 2>&1 ;
    replace_names_in_matrix "${DATA}/${DIRECTORY}" "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/renamed_genomes_list.txt" "${distance_matrix_wo_ref}" "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/renamed_genomes_list.woref.links.txt" "${distance_matrix_wo_ref}.renamed"

  fi




  ## 2 - Computes hierarchical clustering (average linkage = UPGMA) from distance matrix and optionally cuts this hierarchical clustering


  ### VARIABLE USED FOR ALL THE NEXT STEPS
  SELECT_ALL=0

  if [[ "${SELECTEDGENOMENB}" = '-1' ]]; then SELECT_ALL=1; fi
  if [[ "${SELECTEDGENOMENB}" != '-1' ]] && [[ $SELECTEDGENOMENB -eq $genome_number ]]; then SELECT_ALL=1; fi
  if [[ "${SELECTEDGENOMENB}" != '-1' ]] && [[ $SELECTEDGENOMENB -gt $genome_number ]]; then SELECT_ALL=1; fi
  if [[ "${SELECTEDGENOMENB}" != '-1' ]] && [[ $SELECTEDGENOMENB -lt 2 ]]; then SELECT_ALL=1; fi

  if [[ "${SELECTEDGENOMENB}" != '-1' ]] && [[ $SELECTEDGENOMENB -ge 2 ]] && [[ $SELECTEDGENOMENB -lt $genome_number ]] ; then SELECT_ALL=0; fi


  # genome diversity plots
  plotpdf="${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/diversity_${NAME}_plots.pdf"
  touch "${plotpdf}"
  plotpdf="$(real_path_file "${plotpdf}")"
  if [[ ! -s "${plotpdf}" ]]; then rm "${plotpdf}"; fi;

  plotpdfrenamed="${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/diversity_${NAME}_plots_renamed.pdf"
  touch "${plotpdfrenamed}"
  plotpdfrenamed="$(real_path_file "${plotpdfrenamed}")"
  if [[ ! -s "${plotpdfrenamed}" ]]; then rm "${plotpdfrenamed}"; fi;

  genomeclusterlist="${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/diversity_${NAME}_genome_cluster_list.txt"
  touch "${genomeclusterlist}"
  genomeclusterlist="$(real_path_file "${genomeclusterlist}")"
  if [[ ! -s "${genomeclusterlist}" ]]; then rm "${genomeclusterlist}"; fi;

  genomeclusterlistrenamed="${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/diversity_${NAME}_genome_cluster_list_renamed.txt"
  touch "${genomeclusterlistrenamed}"
  genomeclusterlistrenamed="$(real_path_file "${genomeclusterlistrenamed}")"
  if [[ ! -s "${genomeclusterlistrenamed}" ]]; then rm "${genomeclusterlistrenamed}"; fi;

  ## if reference provided
  plotpdf_wo_ref="${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/diversity_${NAME}_plots_wo_ref.pdf"
  plotpdf_wo_ref_renamed="${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/diversity_${NAME}_plots_wo_ref_renamed.pdf"
  genomeclusterlist_wo_ref="${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/diversity_${NAME}_genome_cluster_list_wo_ref.txt"
  genomeclusterlist_wo_ref_renamed="${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/diversity_${NAME}_genome_cluster_list_wo_ref_renamed.txt"



  ## if reference genome provided or not

  #CASE 1 - REF or NO_REF : if select_all and nb_genomes > 2
  if [[ $SELECT_ALL -eq 1 ]] && [[ $genome_number -gt 2 ]]; then 
    echo 'CASE 1 - REF or NO_REF : if select_all and nb_genomes > 2'
    echo 'phylotree'
    echo "generate phylogenetic tree into ${plotpdf}" >> "${LOG}" 2>&1
    echo "generate phylogenetic tree into ${plotpdf}";
    echo "get_phylotree ${distance_matrix} ${plotpdf}" >> "${LOG}" 2>&1
    get_phylotree "${distance_matrix}" "${plotpdf}"
    echo "generate phylogenetic tree with original isolate names into ${plotpdfrenamed}" >> "${LOG}" 2>&1 ;
    echo "generate phylogenetic tree with original isolate names into ${plotpdfrenamed}";
    echo "get_phylotree ${distance_matrix}.renamed ${plotpdfrenamed}" >> "${LOG}" 2>&1 
    get_phylotree "${distance_matrix}.renamed" "${plotpdfrenamed}"     
  fi

  #CASE 2 - REF or NO_REF : if select_all and nb_genomes = 2
  if [[ $SELECT_ALL -eq 1 ]] && [[ $genome_number -eq 2 ]]; then 
    echo 'CASE 2 - REF or NO_REF : if select_all and nb_genomes = 2'
    echo 'no phylotree, only andi matrix'
  fi


  ## if no reference genome provided only

  #CASE 3 - NO_REF : if nb_selected_genomes >= 2 and nb_selected_genomes < nb_genomes 
  if [[ "${REFGENOME}" = 'N.O.R.E.F' ]] && [[ $SELECT_ALL -eq 0 ]] && [[ $SELECTEDGENOMENB -ge 2 ]]; then 
    echo 'CASE 3 - NO_REF : if nb_selected_genomes >= 2 and nb_selected_genomes < nb_genomes'
    echo 'phylotree AND clustering'
    clusternb="${SELECTEDGENOMENB}"
    echo "generate phylogenetic tree and hierarchical clustering into ${plotpdf}" >> "${LOG}" 2>&1
    echo "generate phylogenetic tree and hierarchical clustering into ${plotpdf}";
    echo "get_phylotree_and_clustering ${distance_matrix} ${plotpdf} ${clusternb} ${genomeclusterlist}"
    get_phylotree_and_clustering "${distance_matrix}" "${plotpdf}" "${clusternb}" "${genomeclusterlist}"
    echo "generate phylogenetic tree and hierarchical clustering with original isolate names into ${plotpdfrenamed}" >> "${LOG}" 2>&1 ;
    echo "generate phylogenetic tree and hierarchical clustering with original isolate names into ${plotpdfrenamed}";
    echo "get_phylotree_and_clustering ${distance_matrix}.renamed ${plotpdfrenamed} ${clusternb} ${genomeclusterlistrenamed}"
    get_phylotree_and_clustering "${distance_matrix}.renamed" "${plotpdfrenamed}" "${clusternb}" "${genomeclusterlistrenamed}"
  fi


  ## if reference genome provided only

  #CASE 4 - REF : if nb_selected_genomes = 2 and nb_selected_genomes < nb_genomes 
  if [[ "${REFGENOME}" != 'N.O.R.E.F' ]] && [[ $SELECT_ALL -eq 0 ]] && [[ $SELECTEDGENOMENB -eq 2 ]]; then
    echo 'CASE 4 - REF : if nb_selected_genomes = 2 and nb_selected_genomes < nb_genomes' 
    echo 'phylotree AND selectOneAmongGenomesInClusterDifferentOfRefGenome'
    clusternb="${SELECTEDGENOMENB}"
    echo "generate phylogenetic tree and hierarchical clustering into ${plotpdf}" >> "${LOG}" 2>&1
    echo "generate phylogenetic tree and hierarchical clustering into ${plotpdf}";
    echo "get_phylotree_and_clustering ${distance_matrix} ${plotpdf} ${clusternb} ${genomeclusterlist}"
    get_phylotree_and_clustering "${distance_matrix}" "${plotpdf}" "${clusternb}" "${genomeclusterlist}"
    echo "generate phylogenetic tree and hierarchical clustering with original isolate names into ${plotpdfrenamed}" >> "${LOG}" 2>&1 ;
    echo "generate phylogenetic tree and hierarchical clustering with original isolate names into ${plotpdfrenamed}";
    echo "get_phylotree_and_clustering ${distance_matrix}.renamed ${plotpdfrenamed} ${clusternb} ${genomeclusterlistrenamed}"
    get_phylotree_and_clustering "${distance_matrix}.renamed" "${plotpdfrenamed}" "${clusternb}" "${genomeclusterlistrenamed}"
  fi

  #CASE 5 - REF : if nb_selected_genomes > 2 and nb_selected_genomes < nb_genomes
  if [[ "${REFGENOME}" != 'N.O.R.E.F' ]] && [[ $SELECT_ALL -eq 0 ]] && [[ $SELECTEDGENOMENB -gt 2 ]]; then
    echo 'CASE 5 - REF : if nb_selected_genomes > 2 and nb_selected_genomes < nb_genomes'
    echo 'phylotree AND clustering'
    echo "generate phylogenetic tree into ${plotpdf}" >> "${LOG}" 2>&1
    echo "generate phylogenetic tree into ${plotpdf}";
    echo "get_phylotree ${distance_matrix} ${plotpdf}" >> "${LOG}" 2>&1
    get_phylotree "${distance_matrix}" "${plotpdf}"
    echo "generate phylogenetic tree with original isolate names into ${plotpdfrenamed}" >> "${LOG}" 2>&1 ;
    echo "generate phylogenetic tree with original isolate names into ${plotpdfrenamed}";
    echo "get_phylotree ${distance_matrix}.renamed ${plotpdfrenamed}" >> "${LOG}" 2>&1
    get_phylotree "${distance_matrix}.renamed" "${plotpdfrenamed}"
        
    touch "${plotpdf_wo_ref}" >> "${LOG}" 2>&1
    plotpdf_wo_ref="$(real_path_file "${plotpdf_wo_ref}")"
    if [[ ! -s "${plotpdf_wo_ref}" ]]; then rm "${plotpdf_wo_ref}"; fi;
    touch "${genomeclusterlist_wo_ref}" >> "${LOG}" 2>&1        
    genomeclusterlist_wo_ref="$(real_path_file "${genomeclusterlist_wo_ref}")"
    if [[ ! -s "${genomeclusterlist_wo_ref}" ]]; then rm "${genomeclusterlist_wo_ref}"; fi;

    touch "${plotpdf_wo_ref_renamed}" >> "${LOG}" 2>&1
    plotpdf_wo_ref_renamed="$(real_path_file "${plotpdf_wo_ref_renamed}")"
    if [[ ! -s "${plotpdf_wo_ref_renamed}" ]]; then rm "${plotpdf_wo_ref_renamed}"; fi;
    touch "${genomeclusterlist_wo_ref_renamed}" >> "${LOG}" 2>&1        
    genomeclusterlist_wo_ref_renamed="$(real_path_file "${genomeclusterlist_wo_ref_renamed}")"
    if [[ ! -s "${genomeclusterlist_wo_ref_renamed}" ]]; then rm "${genomeclusterlist_wo_ref_renamed}"; fi;
  
    clusternb_wo_ref=$((${SELECTEDGENOMENB} - 1))           
    echo "generate phylogenetic tree and hierarchical clustering into ${plotpdf_wo_ref}" >> "${LOG}" 2>&1
    echo "generate phylogenetic tree and hierarchical clustering into ${plotpdf_wo_ref}";
    echo "get_phylotree_and_clustering ${distance_matrix_wo_ref} ${plotpdf_wo_ref} ${clusternb_wo_ref} ${genomeclusterlist_wo_ref}" >> "${LOG}" 2>&1
    get_phylotree_and_clustering "${distance_matrix_wo_ref}" "${plotpdf_wo_ref}" "${clusternb_wo_ref}" "${genomeclusterlist_wo_ref}"
    echo "generate phylogenetic tree and hierarchical clustering with original isolate names into ${plotpdf_wo_ref_renamed}" >> "${LOG}" 2>&1
    echo "generate phylogenetic tree and hierarchical clustering with original isolate names into ${plotpdf_wo_ref_renamed}";
    echo "get_phylotree_and_clustering ${distance_matrix_wo_ref}.renamed ${plotpdf_wo_ref_renamed} ${clusternb_wo_ref} ${genomeclusterlist_wo_ref_renamed}" >> "${LOG}" 2>&1
    get_phylotree_and_clustering "${distance_matrix_wo_ref}.renamed" "${plotpdf_wo_ref_renamed}" "${clusternb_wo_ref}" "${genomeclusterlist_wo_ref_renamed}"
  fi




  ## 3 - get N50 and Genome Size of each fasta file in the list of (file \t cluster)
  ####### TODO : get NG50 (comparable values (N50 not comparables; (get max genome size)))


  ## if no reference provided
  #CASE 3 - NO_REF : if nb_selected_genomes >= 2 and nb_selected_genomes < nb_genomes 
  if [[ "${REFGENOME}" = 'N.O.R.E.F' ]] && [[ $SELECT_ALL -eq 0 ]] && [[ $SELECTEDGENOMENB -ge 2 ]]; then 
    echo "get_assembly_stats ${genomeclusterlist}" >> "${LOG}" 2>&1
    get_assembly_stats "${genomeclusterlist}"
  fi


  ## if reference provided
  #CASE 4 - REF : if nb_selected_genomes = 2 and nb_selected_genomes < nb_genomes 
  if [[ "${REFGENOME}" != 'N.O.R.E.F' ]] && [[ $SELECT_ALL -eq 0 ]] && [[ $SELECTEDGENOMENB -eq 2 ]]; then
    echo "get_assembly_stats ${genomeclusterlist}" >> "${LOG}" 2>&1
    get_assembly_stats "${genomeclusterlist}"
  fi

  #CASE 5 - REF : if nb_selected_genomes > 2 and nb_selected_genomes < nb_genomes
  if [[ "${REFGENOME}" != 'N.O.R.E.F' ]] && [[ $SELECT_ALL -eq 0 ]] && [[ $SELECTEDGENOMENB -gt 2 ]]; then
    echo "get_assembly_stats ${genomeclusterlist_wo_ref}" >> "${LOG}" 2>&1
    get_assembly_stats "${genomeclusterlist_wo_ref}"
  fi




  ## 4 - select input genomes for core genome construction

  echo "creation of the list of selected genomes used as input for the next steps (annotation and core genome building) : ${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt" >> "${LOG}" 2>&1 ;
  echo "creation of the list of selected genomes for the next steps (annotation and core genome building) : ${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt";



  ## if no reference provided
  if [[ "${REFGENOME}" == 'N.O.R.E.F' ]]; then

    #CASE 1 AND CASE 2 - NO_REF : if select_all and nb_genomes >= 2
    if [[ $SELECT_ALL -eq 1 ]] && [[ $genome_number -ge 2 ]]; then
      if [[ -e "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt" ]]; then rm "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt"; fi
      for f in ${genome_list}; do
        echo "${f}" >> "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt" 
      done
      if [[ ! -e "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt" ]]; then echo "[ERROR] the output file ${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt does not exist" 2>>"${LOG}"; exit 1; fi
      if [[ ! -s "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt" ]]; then echo "[ERROR] the output file ${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt is empty" 2>>"${LOG}"; exit 1; fi   
    fi

    #CASE 3 - NO_REF : if nb_selected_genomes >= 2 and nb_selected_genomes < nb_genomes 
    if [[ $SELECT_ALL -eq 0 ]] && [[ $SELECTEDGENOMENB -ge 2 ]]; then
      echo "extract_selected_genomes ${genomeclusterlist} ${SELECTEDGENOMENB}" >> "${LOG}" 2>&1
      extract_selected_genomes "${genomeclusterlist}" "${SELECTEDGENOMENB}"
    fi
  fi


  ## if reference provided
  if [[ "${REFGENOME}" != 'N.O.R.E.F' ]] && [[ -e "${DATA}/${DIRECTORY}/assemblies/${REFGENOME}" ]]; then   ##if fasta file exists

    # get link name of reference genome
    refgenomeprefix="${REFGENOME%.*}";  #file without extension
    newname="$(get_new_name_from_list "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/renamed_genomes_list.txt" "${refgenomeprefix}.fasta")"
    # echo -e "NEWNAME_REFNAME:${newname}" >> "${LOG}" 2>&1
    fasta="$(basename "${newname}")";
    suffix="${fasta#$NAME.}";
    rename="${NAME:0:1}${NAME:2:1}"
    reflinkname="$(real_path_dir "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/genomes")"
    reflinkname="${reflinkname}/${rename}${suffix}"
    # echo -e "NEWNAME_REFLINKNAME:${reflinkname}" >> "${LOG}" 2>&1 
   
    #CASE 1 AND CASE 2 - REF : if select_all and nb_genomes >= 2
    if [[ $SELECT_ALL -eq 1 ]] && [[ $genome_number -ge 2 ]]; then
      echo "${reflinkname}" 1> "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt" 2>>"${LOG}"
      for f in ${genome_list_wo_ref}; do
        echo "${f}" >> "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt" 
      done
      if [[ ! -e "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt" ]]; then echo "[ERROR] the output file ${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt does not exist" 2>>"${LOG}"; exit 1; fi
      if [[ ! -s "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt" ]]; then echo "[ERROR] the output file ${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt is empty" 2>>"${LOG}"; exit 1; fi   
    fi

    ## if reference provided
    #CASE 4 - REF : if nb_selected_genomes = 2 and nb_selected_genomes < nb_genomes 
    if [[ $SELECT_ALL -eq 0 ]] && [[ $SELECTEDGENOMENB -eq 2 ]]; then
      echo "extract_selected_genomes ${genomeclusterlist} ${SELECTEDGENOMENB}" >> "${LOG}" 2>&1
      extract_selected_genomes "${genomeclusterlist}" "${SELECTEDGENOMENB}"
      og1="$(head -n 1 "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt")";
      og2="$(tail -n 1 "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt")";
      cl_ref="$(grep "${reflinkname}" "${genomeclusterlist}.infos" | awk '{print $2}')" # cluster of reference genome 
      grep "${og1}" "${genomeclusterlist}.infos" 1> "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.1.txt" #name and cluster of selected genome 1
      grep "${og2}" "${genomeclusterlist}.infos" 1>> "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.1.txt" #name and cluster of selected genome 2
      awk -v cluster="${cl_ref}" '{if($2 != cluster){print $1} }' "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.1.txt" 1>"${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.2.txt" 2>>"${LOG}"
      echo "${reflinkname}" 1> "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.1.txt" 2>>"${LOG}"
      cat "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.1.txt" "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.2.txt" 1> "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt" 2>>"${LOG}"
      if [[ ! -e "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt" ]]; then echo "[ERROR] the output file ${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt does not exist" 2>>"${LOG}"; exit 1; fi
      if [[ ! -s "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt" ]]; then echo "[ERROR] the output file ${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt is empty" 2>>"${LOG}"; exit 1; fi   
      if [[ -e "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.1.txt" ]]; then rm "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.1.txt"; fi
      if [[ -e "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.2.txt" ]]; then rm "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.2.txt"; fi
    fi

    #CASE 5 - REF : if nb_selected_genomes > 2 and nb_selected_genomes < nb_genomes
    if [[ $SELECT_ALL -eq 0 ]] && [[ $SELECTEDGENOMENB -gt 2 ]]; then
      echo "extract_selected_genomes ${genomeclusterlist_wo_ref} ${genome_to_select_wo_ref}" >> "${LOG}" 2>&1
      extract_selected_genomes "${genomeclusterlist_wo_ref}" "${genome_to_select_wo_ref}"
      mv "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt" "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.2.txt" >> "${LOG}" 2>&1
      echo "${reflinkname}" 1> "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.1.txt" 2>>"${LOG}"
      cat "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.1.txt" "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.2.txt" 1> "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt" 2>>"${LOG}"
      if [[ ! -e "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt" ]]; then echo "[ERROR] the output file ${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt does not exist" 2>>"${LOG}"; exit 1; fi
      if [[ ! -s "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt" ]]; then echo "[ERROR] the output file ${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt is empty" 2>>"${LOG}"; exit 1; fi   
      if [[ -e "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.1.txt" ]]; then rm "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.1.txt"; fi
      if [[ -e "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.2.txt" ]]; then rm "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.2.txt"; fi
    fi

  fi
}

main "$@"

