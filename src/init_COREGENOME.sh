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
# init_COREGENOME.sh: initial steps before building a core or persistent gene set.                                  #
# Authors: Elise Larsonneur, Damien Mornico                                                                         #
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
IDPRCT=80                  ## -i
PROTLENGTHRATIO=1.3        ## -l
SYNTGENESUM=4              ## -S
SYNTRADIUSSIZE=5           ## -R
REFGENOME='N.O.R.E.F'      ## -g
SELECTEDGENOMENB='-1'      ## -s
CGGENOMEPRCT=95            ## -p



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
  echo '   ### DIVERSITY ###';
  echo '   -s <int>     number of genomes to select for core genome construction (default : -1); if -1, no selection (all genomes are kept)';
  echo '   ### CORE-GENOME ###';
  echo '   -i <int>     similarity percent (default : 80) -- core genome construction step';
  echo '   -l <float>   protein length ratio (default : 1.2) -- core genome construction step';
  echo '   -S <int>     synteny - synteny threshold : minimal number of syntenic genes found in the neighborhood of a given homologous gene, ';
  echo '                      the boundaries of explored neighborood are defined by window size (option -R); (default : 4); if 0, no synteny criteria is applied';
  echo '   -R <int>     synteny - radius size around each analyzed homologous gene (a number of genes) (default : 5) -- window_size=(radius_size * 2 + 1)';
  echo "   -p <int>     core genes are present at least in p% of the genomes (default : 95) -- if '-p 100' is supplied, core gene set is output; else if 'p' is lower than 100, persistent gene set is output";
  echo '   ### GENERAL OPTIONS ###';
  echo "   -z <STEPS>   steps to run ; ex : '-z DAC' or '-z DA' or '-z D' ; D : diversity - A : annotation - C : coregenome (default: DAC)";
  echo '   -g <inFile>  reference genome fasta file -- annotation and core genome construction steps';
  echo '   -h           print help';
  echo '';
  echo 'EXAMPLES :';
  echo '     #provided reference genome, the functional annotation of reference genome will be transferred to the other genomes:';
  echo "     $0 -d klpn5refannot -n klpn -g MGH78578_NC.fasta -p 95 -t 4";
  echo '';
  echo '     #not provided reference genome, the reference genome will be the first on the genome list sorted alphabetically:';
  echo "     $0 -d klpn5refannot -n klpn -p 100 -t 4";  
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
  while getopts :C:d:n:i:l:S:R:g:s:p: option
  do
    if [[ -z "${OPTARG}" ]]; then echo "[ERROR] empty argument for option -${option}"; exit 1; fi
    case "${option}" in
      C)  
        CGB_CONFIG="${OPTARG}";    
        if [[ ! -d "${CGB_CONFIG}" ]]; then echo "[ERROR] input directory '${CGB_CONFIG}' does not exist (option -C)." ; exit 1 ; else source ${CGB_CONFIG}/config_env.txt; source ${CGB_CONFIG}/config_coregenome_param.txt; fi  
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
      i)  
        IDPRCT="${OPTARG}";
        if ! [[ "${IDPRCT}" =~ ^[0-9]+$ ]] || [[ $IDPRCT -lt 0 ]] || [[ $IDPRCT -gt 100 ]]; then echo '[ERROR] the similarity percent threshold must range from 0 to 100 (option -q).' ; exit 1 ; fi  
        ;; # -i <similarity percent threshold>
      l)    
        PROTLENGTHRATIO="${OPTARG}";
        num1="$(echo "$PROTLENGTHRATIO <= 0.0" |bc -l)";
        num2="$(echo "$PROTLENGTHRATIO >= 10.0" |bc -l)"; 
        if [[ $num1 -eq 1 ]] || [[ $num2 -eq 1 ]] ; then echo '[ERROR] the protein length ratio threshold must be greater than 0 (option -l).' ; exit 1 ; fi
        ;; # -l < protein length ratio threshold>
      S)    
        SYNTGENESUM="${OPTARG}";
        if ! [[ "${SYNTGENESUM}" =~ ^[0-9]+$ ]] || [[ $SYNTGENESUM -lt 0 ]]; then echo '[ERROR] the syntenic gene number threshold must be greater or equal to 0 (option -S).' ; exit 1 ; fi
        ;; # -S <syntenic gene number threshold> 
      R)    
        SYNTRADIUSSIZE="${OPTARG}";
        if ! [[ "${SYNTRADIUSSIZE}" =~ ^[0-9]+$ ]] || [[ $SYNTRADIUSSIZE -lt 1 ]]; then echo '[ERROR] the syntenic radius size threshold must be greater than 0 (option -R).' ; exit 1 ; fi
        ;; # -R <syntenic radius size threshold>
      g)    
        REFGENOME="${OPTARG}";
        ;;     # -g <reference genome FASTA infile>
      s)
        SELECTEDGENOMENB="${OPTARG}";
        if ! [[ "${SELECTEDGENOMENB}" =~ ^[0-9]+$ ]] && [[ "${SELECTEDGENOMENB}" != '-1' ]]; then echo '[ERROR] Bad value for option -s' ; exit 1 ; fi;
        if [[ "${SELECTEDGENOMENB}" != '-1' ]] && [[ $SELECTEDGENOMENB -lt 2 ]]; then
          echo '[ERROR] the number of genomes to select for core genome construction must be greater than 1 (option -s).' ; exit 1 ;
        fi
        ;; # -s <number of selected genomes threshold>
      p)  
        CGGENOMEPRCT="${OPTARG}";
        if ! [[ "${CGGENOMEPRCT}" =~ ^[0-9]+$ ]] || [[ $CGGENOMEPRCT -lt 0 ]] || [[ $CGGENOMEPRCT -gt 100 ]]; then 
          echo "[ERROR] the core gene set will be present in '-p <int>' genomes over 100 genomes'. The percentage threshold 'p' must range from 0 to 100 (option -p)." ; exit 1 ; 
        fi  
        ;; # -p <homologs present in <p>% of the proteomes>
      :)  
        echo "[ERROR] option ${OPTARG} : missing argument" ; exit 1   
        ;;  
      \?) 
        echo "[ERROR] ${OPTARG} : option invalide" ; exit 1   
        ;;  
    esac
  done

  readonly CGB_CONFIG
  readonly DIRECTORY NAME IDPRCT PROTLENGTHRATIO SYNTGENESUM SYNTRADIUSSIZE REFGENOME SELECTEDGENOMENB CGGENOMEPRCT



  ### checking input directory
  if [[ "${DIRECTORY}" = 'N.O.D.I.R' ]]; then echo '[ERROR] no input directory supplied (mandatory option -d)' ; exit 1 ; fi
  if [[ ! -e "${DATA}/${DIRECTORY}" ]]; then echo "[ERROR] directory ${DATA}/${DIRECTORY} does not exist, please create it." ; exit 1 ; fi


  ### create log
  mkdir -p "${DATA}/${DIRECTORY}/logs"
  LOG="${DATA}/${DIRECTORY}/logs/${DIRECTORY}.coregenome.1.log";
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

  ## checking if the step of selection of a subset of genome sequences is activated
  if [[ "${SELECTEDGENOMENB}" = '-1' ]]; then echo 'Skipping step of genome selection.' ; fi




  echo 'starting core genome building step' >> "${LOG}" 2>&1 ;
  echo 'starting core genome building step';
  echo 'prepare coregenome directories and files' >> "${LOG}" 2>&1
  echo 'prepare coregenome directories and files';

  #CoreGenome directory creation
  if [[ -d "${DATA}/${DIRECTORY}/${DATA_CORE_DIR}" ]]; then rm -r "${DATA}/${DIRECTORY}/${DATA_CORE_DIR}" >> "${LOG}" 2>&1; fi
  mkdir -p "${DATA}/${DIRECTORY}/${DATA_CORE_DIR}" >> "${LOG}" 2>&1
  mkdir -p "${DATA}/${DIRECTORY}/${DATA_CORE_DIR}/${temp_dir}" >> "${LOG}" 2>&1   ## $DATA_CORE_DIR/temp

  #create input dir
  mkdir -p "${DATA}/${DIRECTORY}/${DATA_CORE_DIR}/${input_dir}" >> "${LOG}" 2>&1  ## $DATA_CORE_DIR/temp/Proteins
  for file in ${DATA}/${DIRECTORY}/${DATA_CORE_DIR}/${input_dir}/*; do if [[ -f "${file}" ]]; then rm "${file}" >> "${LOG}" 2>&1; fi; done
  cp "${DATA}/${DIRECTORY}/${DATA_PROT_DIR}/"*.prt "${DATA}/${DIRECTORY}/${DATA_CORE_DIR}/${input_dir}" >> "${LOG}" 2>&1

  #create output directory for this analysis
  mkdir -p "${DATA}/${DIRECTORY}/${DATA_CORE_DIR}/${output_dir}" >> "${LOG}" 2>&1 ## $DATA_CORE_DIR/temp/output

  #create TMP dir
  mkdir -p "${DATA}/${DIRECTORY}/${DATA_CORE_DIR}/${tmp_dir}" >> "${LOG}" 2>&1    ## $DATA_CORE_DIR/temp/TMP


  #Create proteome list 'Genomes-${DIRECTORY}.lst'       
  echo "create proteome list : ${DATA}/${DIRECTORY}/${DATA_CORE_DIR}/Genomes-${DIRECTORY}.lst" >> "${LOG}" 2>&1
  echo "create proteome list : ${DATA}/${DIRECTORY}/${DATA_CORE_DIR}/Genomes-${DIRECTORY}.lst";

  ##version 1 - only succesful annotations from ecamber :
  ls -1 "${DATA}/${DIRECTORY}/${DATA_CORE_DIR}/${input_dir}" | gawk '{ gsub(".prt","",$0); print $0}' 1> "${DATA}/${DIRECTORY}/${DATA_CORE_DIR}/Genomes-${DIRECTORY}.lst" 2>>"${LOG}"

  ##version 2 - warning : successfull and unsuccesfull annotations from ecamber :
  ## the next code cancels the previous line
  ## to have the reference file in first position in file 'Genomes-${DIRECTORY}.lst'  (the list selected_genomes_list.txt was previously sorted !)
  if [[ -e "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt" ]]; then
    cp -p "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt" "${DATA}/${DIRECTORY}/${DATA_CORE_DIR}/Genomes-${DIRECTORY}.lst.tmp" >> "${LOG}" 2>&1
    genomelist="${DATA}/${DIRECTORY}/${DATA_CORE_DIR}/Genomes-${DIRECTORY}.lst.tmp"

    while read line; do
      fasta="$(readlink -f "${line}")"   ## get real path of file pointed by the symbolic link
      base="$(basename "${fasta}")"
      echo "${base%.*}"  >> "${DATA}/${DIRECTORY}/${DATA_CORE_DIR}/Genomes-${DIRECTORY}.lst.tmptmp"
    done < ${genomelist}
    
    if [[ -e "${genomelist}" ]]; then rm "${genomelist}" >> "${LOG}" 2>&1; fi
    mv "${DATA}/${DIRECTORY}/${DATA_CORE_DIR}/Genomes-${DIRECTORY}.lst.tmptmp" "${DATA}/${DIRECTORY}/${DATA_CORE_DIR}/Genomes-${DIRECTORY}.lst" >> "${LOG}" 2>&1
  fi
}

main "$@"

