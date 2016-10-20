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
# run_ANNOTATION.sh: computes genome annotation.                                                                    #
# Authors : Elise Larsonneur, Damien Mornico                                                                        #
#####################################################################################################################


set -o pipefail


########################################################################################
# GLOBAL VARIABLES & SOURCE LIBRARY AND CONFIGURATION FILES                            #
########################################################################################

source "${CGB_BIN}utils.sh"

## parameters initialization
CGB_CONFIG='../config/'    ## -C
DIRECTORY='N.O.D.I.R'      ## -d
NAME='N.O.N.A.M.E'         ## -n
REFGENOME='N.O.R.E.F'      ## -g
THREADS=1                  ## -t
REFANNOTATION='N.O.R.E.F'  ## -a
REFIDPATTERN='N.O.P.A.T.T.E.R.N'   ## -e

#ecamber specific parameter
AS='prodigal'

LOG=''



########################################################################################
# FUNCTIONS                                                                            #
########################################################################################

# display_usage 
# This function displays the usage of this program.
# No parameters
function display_usage {
  echo '' ;
  echo 'USAGE :' ;
  echo "     $0 [options] -d <input_directory> -n <name_of_4_chars>" ;
  echo '   -d <inDirectory>  directory where are sequence fasta files';
  echo '   -n <string>  four letter name (ex : esco (es:escherichia; co:coli))';
  echo "  where 'options' are :" ;
  echo '   ### ANNOTATION ###';
  echo '   -g <inFile>  reference genome fasta file -- annotation and core genome construction steps';
  echo '   -a <inFile>  reference genome genbank file -- annotation step';
  echo "   -e <inFile>  prefix string of sequence ids found in the reference genome fasta and genbank files (ex : 'NC_'; 'NZ_' ; 'AKAC' ; etc)";
  echo '   ### GENERAL OPTIONS ###';
  echo '   -t <int>     number of threads used during diversity and annotation steps (default : 1)';
  echo '' ;
  echo ' EXAMPLE :'
  echo "     $0 -d esco_dir -n esco -t 8" ;
}



# correctEcamberAnnsParsed
# This functions corrects errors of an ecamber annotation file, errors in start, stop, strand and contig id values.
# Parameters
# - 1) Input file - A string containing the file path of an output ecamber annotation file (ecamber/datasets/directory_name/anns_parsed/*/*.txt).
# Outputs
# - The corrected annotation file (same file name as input annotation file).
function correctEcamberAnnsParsed {

        if [[ -z "$1" ]]; then echo '[ERROR] the function correctEcamberAnnsParsed expects an input file path to run' 2>>"${LOG}"; exit 1; fi
        if [[ $# -ne 1 ]]; then echo "[ERROR] the function correctEcamberAnnsParsed expects 1 parameter, $# parameters given" 2>>"${LOG}"; exit 1; fi 
        local file_=$1   ## ecamber/datasets/directory_name/anns_parsed/*/*.txt       
        if [[ ! -e "${file_}" ]]; then echo "[ERROR] the input file ${file_} does not exist" 2>>"${LOG}"; exit 1; fi
        if [[ ! -s "${file_}" ]]; then echo "[ERROR] the input file ${file_} is empty" 2>>"${LOG}"; exit 1; fi
 
        gawk '  BEGIN{FS="\t";};
        {
          if(NR == 1){
              contig = $6; 
              borne_min_p = $2;
              if(strand == "-") {borne_min_p = $3};
              gsub(">","",borne_min_p);
              gsub("<","",borne_min_p); 
              gsub("\\(","",borne_min_p); 
              gsub("\\)","",borne_min_p); 
          }
          
          borne_min = $2;
          if(strand == "-") {borne_min = $3};         
          gsub(">","",borne_min);
          gsub("<","",borne_min); 
          gsub("\\(","",borne_min); 
          gsub("\\)","",borne_min); 
 
          if(contig != $6){
               if(borne_min > borne_min_p){ 
                        print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" contig;
               }else{
                        print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6;
               }
                
          }else{
                print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6;
          }
          contig = $6;
          borne_min_p = borne_min;
                
        }' "${file_}" 1> "${file_}.tmp" 2>>"${LOG}"
        mv "${file_}.tmp" "${file_}" >> "${LOG}" 2>&1
        if [[ ! -e "${file_}" ]]; then echo "[ERROR] the output file ${file_} does not exist" 2>>"${LOG}"; exit 1; fi
        if [[ ! -s "${file_}" ]]; then echo "[ERROR] the output file ${file_} is empty" 2>>"${LOG}"; exit 1; fi
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
    exit 0
  fi

  # if less than two arguments supplied, display usage 
  if [[  $# -le 1 ]]; then
    display_usage
    exit 1
  fi


  ## catch option values
  while getopts :C:d:n:g:a:e:t: option
  do
    if [[ -z "$OPTARG" ]]; then echo "[ERROR] empty argument for option -${option}"; exit 1; fi
    case $option in
      C)
        CGB_CONFIG="$OPTARG";
        if [[ ! -d "${CGB_CONFIG}" ]]; then echo "[ERROR] input directory '$CGB_CONFIG' does not exist (option -C)." ; exit 1 ; else source ${CGB_CONFIG}/config_env.txt; fi
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
  readonly DIRECTORY NAME REFGENOME THREADS REFANNOTATION REFIDPATTERN
 


  ### checking input directory
  if [[ "${DIRECTORY}" = 'N.O.D.I.R' ]]; then echo '[ERROR] no input directory supplied (mandatory option -d)' ; exit 1 ; fi
  if [[ ! -e "${DATA}/${DIRECTORY}" ]]; then echo "[ERROR] directory ${DATA}/${DIRECTORY} does not exist, please create it." ; exit 1 ; fi
   
    
  ### create log
  mkdir -p "${DATA}/${DIRECTORY}/logs"
  LOG="${DATA}/${DIRECTORY}/logs/${DIRECTORY}.annotation.2.log";
  readonly LOG;
  echo '' > "${LOG}";


  ### checking name 
  if [[ "${NAME}" = 'N.O.N.A.M.E' ]]; then echo '[ERROR] no name supplied (mandatory option -n)' ; exit 1 ; fi

  ### checking if directory containing genomes (assemblies) exists and is not empty
  if [[ ! -e "${DATA}/${DIRECTORY}/assemblies" ]]; then echo "[ERROR] directory ${DATA}/${DIRECTORY}/assemblies does not exist, please create it." ; exit 1 ; fi
  if [[ ! "$(ls -A ${DATA}/${DIRECTORY}/assemblies)" ]]; then echo "[ERROR] directory ${DATA}/${DIRECTORY}/assemblies is empty, please add some genome fasta files into it" ; exit 1; fi

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
    


if [[ "${REFGENOME}" != 'N.O.R.E.F' ]] && [[ "${REFANNOTATION}" != 'N.O.R.E.F' ]] && [[ "${REFIDPATTERN}" != 'N.O.P.A.T.T.E.R.N' ]]; then

   base="$(basename "${REFGENOME}")";
   prefix="${base%.*}";  #file without extension
   refname="${prefix}.fasta";
   newrefname="$(grep "${refname}" "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/renamed_genomes_list.txt" |awk '{ print $2 }' )";
   newprefixrefgenome="${newrefname%.*}";
   newrefannotname="${newprefixrefgenome}.txt";  #genbank file will have the same prefix than the reference fasta file

   annotation_ecamber="${DATA_ECAMBER_OUT}/${DIRECTORY}/anns_parsed/refseq_gbk/${newrefannotname}";

   # prepare annotations (format to ecamber annotation format) and genomes (format blast db)
   echo "ecamber format reference annotation and genomes : python ${ECAMBER}/ecamber.py -a f -d ${DIRECTORY} -w ${THREADS}"
   echo "=====================================================" >> "${LOG}"
   echo "STEP 0 : python ${ECAMBER}/ecamber.py -a f -d ${DIRECTORY} -w ${THREADS}" >> "${LOG}"
   python "${ECAMBER}/ecamber.py" -a f -d "${DIRECTORY}" -w "${THREADS}" >> "${LOG}" 2>&1
   exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] the command 'python $ECAMBER/ecamber.py -a f -d $DIRECTORY -w $THREADS' exited with an error! - exit code ${exit_status}" 2>>"${LOG}"; exit 1; fi

   # add locus_tag as gene_name when no gene_name is defined
   # some annotation lines (in ecamber format) are preceded by '#'
   awk 'BEGIN {FS="\t"}; { if ($5 == "") {locustag=$1 ; gsub(/^#/, "", locustag); print $1"\t"$2"\t"$3"\t"$4"\t"locustag"\t"$6} else {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6} }' "${annotation_ecamber}" 1> "${annotation_ecamber}.tmp" 2>>"${LOG}"
   mv "${annotation_ecamber}.tmp" "${annotation_ecamber}" >> "${LOG}" 2>&1
fi


echo "ecamber prodigal annotation : python ${ECAMBER}/ecamber.py -a pr -d ${DIRECTORY} -w ${THREADS}"
echo "=====================================================" >> "${LOG}"
echo "STEP 1 : python ${ECAMBER}/ecamber.py -a pr -d ${DIRECTORY} -w ${THREADS}" >> "${LOG}"
python "${ECAMBER}/ecamber.py" -a pr -d "${DIRECTORY}" -w "${THREADS}" >> "${LOG}" 2>&1
exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] the command 'python ${ECAMBER}/ecamber.py -a pr -d ${DIRECTORY} -w ${THREADS}' exited with an error! - exit code ${exit_status}" 2>>"${LOG}"; exit 1; fi

echo "prodigal parsed annotations correction" >> "${LOG}" 2>&1
echo "prodigal parsed annotations correction"

for f in $(ls -1 "${DATA_ECAMBER_OUT}/${DIRECTORY}/anns_parsed/prodigal/"*.txt); do
        correctEcamberAnnsParsed "${f}" >> "${LOG}" 2>&1
done

# if reference annotation
if [[ "${REFGENOME}" != 'N.O.R.E.F' ]] && [[ "${REFANNOTATION}" != 'N.O.R.E.F' ]] && [[ "${REFIDPATTERN}" != 'N.O.P.A.T.T.E.R.N' ]]; then
   
   base="$(basename "${REFGENOME}")";
   prefix="${base%.*}";  #file without extension
   refname="${prefix}.fasta";
   newrefname="$(grep "${refname}" "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/renamed_genomes_list.txt" | awk '{ print $2 }')";
   newprefixrefgenome="${newrefname%.*}";
   newrefannotname="${newprefixrefgenome}.txt";  #genbank file will have the same prefix than the reference fasta file

   annot_genbank_file="${DATA_ECAMBER_OUT}/${DIRECTORY}/anns_parsed/refseq_gbk/${newrefannotname}";


   #remove de novo annotation (prodigal) of reference genome
   rm "${DATA_ECAMBER_OUT}/${DIRECTORY}/anns_input/prodigal/${newprefixrefgenome}.txt" >> "${LOG}" 2>&1  # prodigal output annotation
   rm "${DATA_ECAMBER_OUT}/${DIRECTORY}/anns_parsed/prodigal/${newprefixrefgenome}.txt" >> "${LOG}" 2>&1 # prodigal output annotation formated as ecamber annotation
   rm "${DATA_ECAMBER_OUT}/${DIRECTORY}/prodigal/train/${newprefixrefgenome}" >> "${LOG}" 2>&1
   rm "${DATA_ECAMBER_OUT}/${DIRECTORY}/prodigal/scores/${newprefixrefgenome}.txt" >> "${LOG}" 2>&1



   # the same correction to reference annotation .txt
   echo 'prodigal parsed reference annotation correction' >> "${LOG}" 2>&1
   echo 'prodigal parsed reference annotation correction'
   correctEcamberAnnsParsed "${annot_genbank_file}";
   
   #copy genbank annotation to prodigal directory (ecamber annotation format)  #is this really neccessary ???
   cp -p "${DATA_ECAMBER_OUT}/${DIRECTORY}/anns_parsed/refseq_gbk/${newprefixrefgenome}.txt" "${DATA_ECAMBER_OUT}/${DIRECTORY}/anns_parsed/prodigal/." >> "${LOG}" 2>&1

fi


echo "ecamber Blast alignment : python ${ECAMBER}/ecamber.py -a ph1 -d ${DIRECTORY} -w ${THREADS} -as ${AS}"
echo "=====================================================" >> "${LOG}"
echo "STEP 2 : python ${ECAMBER}/ecamber.py -a ph1 -d ${DIRECTORY} -w ${THREADS} -as ${AS}" >> "${LOG}"
python "${ECAMBER}/ecamber.py" -a ph1 -d "${DIRECTORY}" -w "${THREADS}" -as "${AS}" >> "${LOG}" 2>&1
exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] the command 'python ${ECAMBER}/ecamber.py -a ph1 -d ${DIRECTORY} -w ${THREADS} -as ${AS}' exited with an error! - exit code ${exit_status}" 2>>"${LOG}"; exit 1; fi 

echo "ecamber TIS voting : python ${ECAMBER}/ecamber.py -a ph2 -d ${DIRECTORY} -w ${THREADS} -as ${AS}"
echo "=====================================================" >> "${LOG}"
echo "STEP 3 : python ${ECAMBER}/ecamber.py -a ph2 -d ${DIRECTORY} -w ${THREADS} -as ${AS}" >> "${LOG}"
python "${ECAMBER}/ecamber.py" -a ph2 -d "${DIRECTORY}" -w "${THREADS}" -as "${AS}" >> "${LOG}" 2>&1
exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] the command 'python ${ECAMBER}/ecamber.py -a ph2 -d ${DIRECTORY} -w ${THREADS} -as ${AS}' exited with an error! - exit code ${exit_status}" 2>>"${LOG}"; exit 1; fi

echo "ecamber output : python ${ECAMBER}/ecamber.py -a out -d ${DIRECTORY} -w ${THREADS} -as ${AS}"
echo "=====================================================" >> "${LOG}"
echo "STEP 4 : python ${ECAMBER}/ecamber.py -a out -d ${DIRECTORY} -w ${THREADS} -as ${AS}" >> "${LOG}"
python "${ECAMBER}/ecamber.py" -a out -d "${DIRECTORY}" -w "${THREADS}" -as "${AS}" >> "${LOG}" 2>&1
exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] the command 'python ${ECAMBER}/ecamber.py -a out -d ${DIRECTORY} -w ${THREADS} -as ${AS}' exited with an error! - exit code ${exit_status}" 2>>"${LOG}"; exit 1; fi



# if reference annotation
# parse ecamber output files TO TRANSFER REFERENCE GENE ANNOTATIONS TO OTHERS (DE NOVO) ANNOTATIONS
if [[ "${REFGENOME}" != 'N.O.R.E.F' ]] && [[ "${REFANNOTATION}" != 'N.O.R.E.F' ]] && [[ "${REFIDPATTERN}" != 'N.O.P.A.T.T.E.R.N' ]]; then

    echo 'ecamber - transfer reference annotation to other annotations' >> "${LOG}" 2>&1    
    echo 'ecamber - transfer reference annotation to other annotations';    

    base="$(basename "${REFGENOME}")";
    prefix="${base%.*}";  #file without extension
    refname="${prefix}.fasta";
    newrefname="$(grep "${refname}" "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/renamed_genomes_list.txt" | awk '{ print $2 }')";
    newprefixrefgenome="${newrefname%.*}";
    newrefannotname="${newprefixrefgenome}.txt";  #genbank file will have the same prefix than the reference fasta file
    strainname="${prefix}";
 
    # NEXT STEPS IN ORDER TO TRANSFER REFERENCE GENE ANNOTATIONS TO OTHERS (DE NOVO) ANNOTATIONS:
    # 1- modification of ecamber output reference annotation - rename gene_name and store clusters into 'cluster_gene_names.modif.txt'
    # 2- rename gene in de novo prodigal annotations output by ecamber using the info from file 'cluster_gene_names.modif.txt'
    # 3- modify 'output/clusters_table.txt' and output it to 'output/clusters_table.modif.txt'

    echo "${CGB_BIN}get_cds_info_from_genbank.py -i ${DATA}/${DIRECTORY}/ref_gbk_annotation/${REFANNOTATION}" >> "${LOG}" 2>&1
    ${CGB_BIN}get_cds_info_from_genbank.py -i "${DATA}/${DIRECTORY}/ref_gbk_annotation/${REFANNOTATION}"    
    exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] the command '${CGB_BIN}get_cds_info_from_genbank.py -i ${DATA}/${DIRECTORY}/ref_gbk_annotation/${REFANNOTATION}' exited with an error! - exit code ${exit_status}" 2>>"${LOG}"; exit 1; fi    
 
    # 1- modification of ecamber output reference annotation - rename gene_name and store clusters into 'cluster_gene_names.modif.txt'
    ## a- genes_aa
    fastain="${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_aa/${newprefixrefgenome}.fasta"
    fastaout="${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_aa/${newprefixrefgenome}.modif.fasta"
    clustersout="${DATA_ECAMBER_OUT}/${DIRECTORY}/output/cluster_gene_names.ecamber_aa.txt"
    annotgblst="${DATA}/${DIRECTORY}/ref_gbk_annotation/${REFANNOTATION}.txt"
    echo "${CGB_BIN}fasta_rename_reference_genes.py -i ${fastain} -o ${fastaout} -c ${clustersout} -a ${annotgblst} -n ${strainname}" >> "${LOG}" 2>&1
    ${CGB_BIN}fasta_rename_reference_genes.py -i "${fastain}" -o "${fastaout}" -c "${clustersout}" -a "${annotgblst}" -n "${strainname}" >> "${LOG}" 2>&1
    exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] the command '${CGB_BIN}fasta_rename_reference_genes.py -i ${fastain} -o ${fastaout} -c ${clustersout} -a ${annotgblst} -n ${strainname}' exited with an error! - exit code ${exit_status}" 2>>"${LOG}"; exit 1; fi   
 
    ## b- genes_dna
    fastain="${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_dna/${newprefixrefgenome}.fasta"
    fastaout="${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_dna/${newprefixrefgenome}.modif.fasta"
    clustersout="${DATA_ECAMBER_OUT}/${DIRECTORY}/output/cluster_gene_names.ecamber_dna.txt"
    echo "${CGB_BIN}fasta_rename_reference_genes.py -i ${fastain} -o ${fastaout} -c ${clustersout} -a ${annotgblst} -n ${strainname}" >> "${LOG}" 2>&1
    ${CGB_BIN}fasta_rename_reference_genes.py -i "${fastain}" -o "${fastaout}" -c "${clustersout}" -a "${annotgblst}" -n "${strainname}" >> "${LOG}" 2>&1
    exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] the command '${CGB_BIN}fasta_rename_reference_genes.py -i ${fastain} -o ${fastaout} -c ${clustersout} -a ${annotgblst} -n ${strainname}' exited with an error! - exit code ${exit_status}" 2>>"${LOG}"; exit 1; fi



    # 2- rename gene_name in de novo prodigal annotations output by ecamber using the info from file 'cluster_gene_names.modif.txt'

    ## a- genes_aa
    for fasta in $(ls -1 "${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_aa/"*.fasta |grep -v "${newprefixrefgenome}.fasta" |grep -v "${newprefixrefgenome}.modif.fasta"); do
       prefix="$(basename "${fasta}" '.fasta')"; 
       fastain="${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_aa/${prefix}.fasta"
       fastaout="${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_aa/${prefix}.modif.fasta"
       clustersin="${DATA_ECAMBER_OUT}/${DIRECTORY}/output/cluster_gene_names.ecamber_aa.txt"
       strainname=''
       strainname="$(get_old_name_from_list "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/renamed_genomes_list.txt" "${prefix}.fasta")"
       if [[ -n "${strainname}" ]]; then strainname="$(basename "${strainname}" '.fasta')"; fi
       ${CGB_BIN}fasta_rename_slave_genes.py -i "${fastain}" -o "${fastaout}" -c "${clustersin}" -n "${strainname}" >> "${LOG}" 2>&1
       exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] the command '${CGB_BIN}fasta_rename_slave_genes.py -i ${fastain} -o ${fastaout} -c ${clustersin} -n ${strainname}' exited with an error! - exit code ${exit_status}" 2>>"${LOG}"; exit 1; fi
    done

    ## b- genes_dna
    for fasta in $(ls -1 "${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_dna/"*.fasta |grep -v "${newprefixrefgenome}.fasta" |grep -v "${newprefixrefgenome}.modif.fasta"); do
       prefix="$(basename "${fasta}" '.fasta')";
       fastain="${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_dna/${prefix}.fasta"
       fastaout="${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_dna/${prefix}.modif.fasta"
       clustersin="${DATA_ECAMBER_OUT}/${DIRECTORY}/output/cluster_gene_names.ecamber_dna.txt"
       strainname=""
       strainname="$(get_old_name_from_list "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/renamed_genomes_list.txt" "${prefix}.fasta")"
       if [[ -n "${strainname}" ]]; then strainname="$(basename "${strainname}" '.fasta')"; fi
       ${CGB_BIN}fasta_rename_slave_genes.py -i "${fastain}" -o "${fastaout}" -c "${clustersin}" -n "${strainname}" >> "${LOG}" 2>&1
       exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] the command '${CGB_BIN}fasta_rename_slave_genes.py -i ${fastain} -o ${fastaout} -c ${clustersin} -n ${strainname}' exited with an error! - exit code ${exit_status}" 2>>"${LOG}"; exit 1; fi
    done



    #rename '*.modif.fasta' files into '*.fasta'
    ##  a- genes_aa
    for fasta in $(ls -1 "${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_aa/"*.modif.fasta); do
       prefix="$(basename "${fasta}" '.modif.fasta')";
       mv "${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_aa/${prefix}.modif.fasta" "${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_aa/${prefix}.fasta" >> "${LOG}" 2>&1
    done

    ## b- genes_dna
    for fasta in $(ls -1 "${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_dna/"*.modif.fasta); do
       prefix="$(basename "${fasta}" '.modif.fasta')";
       mv "${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_dna/${prefix}.modif.fasta" "${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_dna/${prefix}.fasta" >> "${LOG}" 2>&1
    done



    # 3- modify 'output/clusters_table.txt' and output it to 'output/clusters_table.modif.txt'
    clustersin="${DATA_ECAMBER_OUT}/${DIRECTORY}/output/cluster_gene_names.ecamber_dna.txt"
    clusterstablein="${DATA_ECAMBER_OUT}/${DIRECTORY}/output/clusters_table.txt"
    clusterstableout="${DATA_ECAMBER_OUT}/${DIRECTORY}/output/clusters_table.ecamber.txt"
    echo "${CGB_BIN}update_cluster_table.py -i ${clustersin} -c ${clusterstablein} -o ${clusterstableout}" >> "${LOG}" 2>&1
    echo "${CGB_BIN}update_cluster_table.py -i ${clustersin} -c ${clusterstablein} -o ${clusterstableout}";
    ${CGB_BIN}update_cluster_table.py -i "${clustersin}" -c "${clusterstablein}" -o "${clusterstableout}" >> "${LOG}" 2>&1
    exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] the command '${CGB_BIN}update_cluster_table.py -i ${clustersin} -c ${clusterstablein} -o ${clusterstableout}' exited with an error! - exit code ${exit_status}" 2>>"${LOG}"; exit 1; fi

fi




#if no reference annotation provided - add additionnal fields ($8 locus_tag ; $9  gene_product ; $10 strain) in fasta headers
if [[ "${REFANNOTATION}" == 'N.O.R.E.F' ]]; then

## a- genes_aa
    for fasta in $(ls -1 "${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_aa/"*.fasta |grep -v '.modif.fasta'); do  
       prefix="$(basename "${fasta}" '.fasta')"; 
       fastain="${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_aa/${prefix}.fasta"
       fastaout="${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_aa/${prefix}.modif.fasta"
       strainname=''
       strainname="$(get_old_name_from_list "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/renamed_genomes_list.txt" "${prefix}.fasta")"
       if [[ -n "${strainname}" ]]; then strainname="$(basename "${strainname}" '.fasta')"; fi
       strtoadd="|x|x|${strainname}"
       tr -d '\15\32' < "${fastain}" |awk -v stradd="${strtoadd}" 'BEGIN{seq=""} !/^>/ {seq=seq$0} /^>/ {if(seq!=""){print seq;seq=""}{id=$0 stradd; print id}} END{print seq}' |grep -v "^$" > "${fastaout}" #add additionnal fields in fasta headers
   done

     ## b- genes_dna
    for fasta in $(ls -1 "${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_dna/"*.fasta |grep -v '.modif.fasta'); do
       prefix="$(basename "${fasta}" '.fasta')";
       fastain="${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_dna/${prefix}.fasta"
       fastaout="${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_dna/${prefix}.modif.fasta"
       strainname=''
       strainname="$(get_old_name_from_list "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/renamed_genomes_list.txt" "${prefix}.fasta")"
       if [[ -n "${strainname}" ]]; then strainname="$(basename "${strainname}" '.fasta')"; fi
       strtoadd="|x|x|${strainname}"
       tr -d '\15\32' < "${fastain}" |awk -v stradd="${strtoadd}" 'BEGIN{seq=""} !/^>/ {seq=seq$0} /^>/ {if(seq!=""){print seq;seq=""}{id=$0 stradd; print id}} END{print seq}' |grep -v "^$" > "${fastaout}" #add additionnal fields in fasta headers
   done

    #rename '*.modif.fasta' files into '*.fasta'
    ##  a- genes_aa
    for fasta in $(ls -1 "${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_aa/"*.modif.fasta); do
       prefix="$(basename "${fasta}" '.modif.fasta')";
       mv "${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_aa/${prefix}.modif.fasta" "${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_aa/${prefix}.fasta" >> "${LOG}" 2>&1
    done

    ## b- genes_dna
    for fasta in $(ls -1 "${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_dna/"*.modif.fasta); do
       prefix="$(basename "${fasta}" '.modif.fasta')";
       mv "${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_dna/${prefix}.modif.fasta" "${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_dna/${prefix}.fasta" >> "${LOG}" 2>&1
    done

fi



# 4- output number of orthoLOG clusters detected by ECAMBER where at most one gene per strain is found
# number of genomes annotated by ecamber 
nb_genomes="$(ls -1 "${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_dna/"*.fasta |wc -l)"
# nombre maximal de souches retrouvees dans chaque cluster
nb_strains="$(tail -n+2 "${DATA_ECAMBER_OUT}/${DIRECTORY}/output/clusters_table.txt" |awk 'BEGIN{FS="\t"} ; {if(max==""){max=$6}; if($6>max) {max=$6} } END {print max}')"
# nombre de clusters pour lesquels il existe un gene de chaque souche parmi les 'nb_strains':
clusters="$(tail -n+2 "${DATA_ECAMBER_OUT}/${DIRECTORY}/output/clusters_table.txt" | awk -v nbstrains="${nb_strains}" 'BEGIN{FS="\t"} ; {if(nb_clusters==""){nb_clusters=0}; if($6==nbstrains) {nb_clusters+=1} } END {print nb_clusters}')";
clusters_all="$(tail -n+2 "${DATA_ECAMBER_OUT}/${DIRECTORY}/output/clusters_table.txt" |wc -l)"
echo -e "ECAMBER --> homolog clusters = \t${clusters_all}" >> "${LOG}" 2>&1
echo -e "ECAMBER --> homolog clusters = \t${clusters_all}"
echo -e "ECAMBER --> core ortholog clusters = clusters shared by at least and at most ${nb_strains}, from ${nb_genomes} annotated strains" >> "${LOG}" 2>&1
echo -e "ECAMBER --> core ortholog clusters = clusters shared by at least and at most ${nb_strains}, from ${nb_genomes} annotated strains"
echo -e "ECAMBER --> CORE ORTHOLOG CLUSTERS = \t${clusters}" >> "${LOG}" 2>&1
echo -e "ECAMBER --> CORE ORTHOLOG CLUSTERS = \t${clusters}";


## ouput/genes_aa/*.fasta example :
## grep ">" ../ext-tools/ecamber/datasets/klpn5refannot/output/genes_aa/klpn.0000001.c001.fasta |head
## >KPN_RS00005|7944|KPN_RS00005|+|21|137|klpn_0000001_c001|KPN_RS00005|hypothetical_protein|MGH78578_NC
## >KPN_RS00010|2848|thrA|+|340|2802|klpn_0000001_c001|KPN_RS00010|bifunctional_aspartokinase_I/homoserine_dehydrogenase_I|MGH78578_NC
## >KPN_RS00015|6009|KPN_RS00015|+|2804|3733|klpn_0000001_c001|KPN_RS00015|homoserine_kinase|MGH78578_NC
## >KPN_RS00020|2784|KPN_RS00020|+|3737|5017|klpn_0000001_c001|KPN_RS00020|threonine_synthase|MGH78578_NC
## >KPN_RS00025|1898|KPN_RS00025|+|5354|5710|klpn_0000001_c001|KPN_RS00025|hypothetical_protein|MGH78578_NC
## >KPN_RS00030|91|KPN_RS00030|-|6552|5779|klpn_0000001_c001|KPN_RS00030|hypothetical_protein|MGH78578_NC
## >KPN_RS00035|3068|KPN_RS00035|-|8060|6630|klpn_0000001_c001|KPN_RS00035|sodium:alanine_symporter|MGH78578_NC
## >KPN_RS00040|311|KPN_RS00040|+|8263|9216|klpn_0000001_c001|KPN_RS00040|transaldolase|MGH78578_NC
## >KPN_RS00045|2878|mogA|+|9316|9903|klpn_0000001_c001|KPN_RS00045|molybdopterin_adenylyltransferase|MGH78578_NC
## >KPN_RS00050|3851|KPN_RS00050|+|10023|11327|klpn_0000001_c001|KPN_RS00050|MFS_transporter|MGH78578_NC

## legend of the fields separated by '|':
## $1  gene_id (=locus_tag if reference)
## $2  cluster_multigene_id
## $3  gene_name
## $4  strand (+/-)
## $5  start/stop
## $6  stop/start
## $7  contig_id (=fastaname if genome have only ONE contig)
## $8  locus_tag (one or several transferred locus_tag(s) separated by a ';' OR 'x') 
## $9  gene_product(s)
## $10 strain
}

main "$@"

