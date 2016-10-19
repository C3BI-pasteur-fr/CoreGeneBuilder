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
# post_DIVERSITY.sh: steps after genomic sequence classification based on estimated genomic evolutionary            #
#                    distances.                                                                                     #
# Author: Elise Larsonneur                                                                                          #
#####################################################################################################################

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
  LOG="${DATA}/${DIRECTORY}/logs/${DIRECTORY}.diversity.3.log";
  readonly LOG;
  echo '' > "${LOG}";
  

  ### checking name 
  if [[ "${NAME}" = 'N.O.N.A.M.E' ]]; then echo '[ERROR] no name supplied (mandatory option -n)' ; exit 1 ; fi

  ### checking if directory containing genomes (assemblies) exists and is not empty
  if [[ ! -e "${DATA}/${DIRECTORY}/assemblies" ]]; then echo "[ERROR] directory ${DATA}/${DIRECTORY}/assemblies does not exist, please create it." ; exit 1 ; fi
  if [[ ! -s "${DATA}/${DIRECTORY}/assemblies" ]]; then echo "[ERROR] directory ${DATA}/${DIRECTORY}/assemblies is empty, please add some genome fasta files into it" ; exit 1 ; fi

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




  #files transfers of diversity output to ecamber input directory (fasta)
  echo "copy genomes to ecamber input directory (ext-tools/ecamber/datasets/${DIRECTORY}) from the list 'selected_genomes_list.txt'." >> "${LOG}" 2>&1
  echo "copy genomes to ecamber input directory (ext-tools/ecamber/datasets/${DIRECTORY}) from the list 'selected_genomes_list.txt'.";
  mkdir -p "${DATA_ECAMBER_OUT}/${DIRECTORY}" >> "${LOG}" 2>&1
  mkdir -p "${DATA_ECAMBER_OUT}/${DIRECTORY}/genomes" >> "${LOG}" 2>&1
  genomelist="$(real_path_file "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt")"
  if [[ ! -s "${genomelist}" ]]; then echo "[ERROR] Genome list ${genomelist} is empty, exiting program..."; exit 1; fi;
  while read line; do
    fasta="$(readlink -f "${line}")"
    cp "${fasta}" "${DATA_ECAMBER_OUT}/${DIRECTORY}/genomes/" >> "${LOG}" 2>&1
  done < ${genomelist}

  #merge logs
  cat "${DATA}/${DIRECTORY}/logs/${DIRECTORY}.diversity."?.log > "${DATA}/${DIRECTORY}/logs/${DIRECTORY}.diversity.log"
  rm "${DATA}/${DIRECTORY}/logs/${DIRECTORY}.diversity."?.log;
}

main "$@"

