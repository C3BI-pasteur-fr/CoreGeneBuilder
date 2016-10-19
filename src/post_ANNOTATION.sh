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
# post_ANNOTATION.sh: steps after genome annotation.                                                                #
# Authors: Elise Larsonneur, Damien Mornico                                                                         #
#####################################################################################################################

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
  LOG="${DATA}/${DIRECTORY}/logs/${DIRECTORY}.annotation.3.log";
  readonly LOG;
  echo '' > "${LOG}";


  ### checking name 
  if [[ "${NAME}" = 'N.O.N.A.M.E' ]]; then echo '[ERROR] no name supplied (mandatory option -n)' >> "${LOG}" 2>&1; exit 1 ; fi

  ### checking if directory containing genomes (assemblies) exists and is not empty
  if [[ ! -e "${DATA}/${DIRECTORY}/assemblies" ]]; then echo "[ERROR] directory ${DATA}/${DIRECTORY}/assemblies does not exist, please create it." >> "${LOG}" 2>&1 ; exit 1 ; fi
  if [[ ! -s "${DATA}/${DIRECTORY}/assemblies" ]]; then echo "[ERROR] directory ${DATA}/${DIRECTORY}/assemblies is empty, please add some genome fasta files into it" >> "${LOG}" 2>&1 ; exit 1 ; fi

  ### checking if reference genome is supplied
  if [[ "${REFGENOME}" = 'N.O.R.E.F' ]]; then echo "Reference genome will be the first fasta file appearing in directory ${DIRECTORY}." >> "${LOG}" 2>&1 ; fi
  if [[ "${REFGENOME}" != 'N.O.R.E.F' ]] && [[ ! -e "${DATA}/${DIRECTORY}/assemblies/${REFGENOME}" ]]; then echo "[ERROR] input fasta file '${DATA}/${DIRECTORY}/assemblies/${REFGENOME}' does not exist (option -g)." >> "${LOG}" 2>&1 ; exit 1 ; fi
  if [[ "${REFGENOME}" != 'N.O.R.E.F' ]] && [[ ! -s "${DATA}/${DIRECTORY}/assemblies/${REFGENOME}" ]]; then echo "[ERROR] input fasta file '${DATA}/${DIRECTORY}/assemblies/${REFGENOME}' is empty (option -g)." >> "${LOG}" 2>&1 ; exit 1 ; fi

  ## checking if reference genbank file is supplied (annotation)
  if [[ "${REFANNOTATION}" = 'N.O.R.E.F' ]]; then echo 'No reference genome annotation provided (genbank file).' >> "${LOG}" 2>&1 ; fi
  if [[ "${REFANNOTATION}" != 'N.O.R.E.F' ]] && [[ ! -e "${DATA}/${DIRECTORY}/ref_gbk_annotation" ]]; then echo "[ERROR] directory containing input genbank file, '${DATA}/${DIRECTORY}/ref_gbk_annotation' does not exist." >> "${LOG}" 2>&1 ; exit 1 ; fi
  if [[ "${REFANNOTATION}" != 'N.O.R.E.F' ]] && [[ ! -e "${DATA}/${DIRECTORY}/ref_gbk_annotation/${REFANNOTATION}" ]]; then echo "[ERROR] input genbank file '${DATA}/${DIRECTORY}/ref_gbk_annotation/${REFANNOTATION}' does not exist (option -a)." >> "${LOG}" 2>&1 ; exit 1 ; fi
  if [[ "${REFANNOTATION}" != 'N.O.R.E.F' ]] && [[ ! -s "${DATA}/${DIRECTORY}/ref_gbk_annotation/${REFANNOTATION}" ]]; then echo "[ERROR] input genbank file '${DATA}/${DIRECTORY}/ref_gbk_annotation/${REFANNOTATION}' is empty (option -a)." >> "${LOG}" 2>&1 ; exit 1 ; fi
  if [[ "${REFANNOTATION}" != 'N.O.R.E.F' ]] && [[ "${REFIDPATTERN}" = 'N.O.P.A.T.T.E.R.N' ]]; then
    echo "[ERROR] prefix (only letters) of genome sequence ids found in reference fasta and genbank files REQUIRED (option -e) if reference genbank file is supplied as reference annotation (option -a). Examples : 'NC_', 'NZ_', 'AKAC'." >> "${LOG}" 2>&1 ; exit 1 ;
  fi



  #Gene and protein directory creation
  mkdir -p "${DATA}/${DIRECTORY}/${DATA_GENE_DIR}" >> "$LOG" 2>&1
  mkdir -p "${DATA}/${DIRECTORY}/${DATA_PROT_DIR}" >> "$LOG" 2>&1
  mkdir -p "${DATA}/${DIRECTORY}/${DATA_ECAMBER_OUTPUT}" >> "$LOG" 2>&1

  #Copy of eCamber outputs
  echo "copy of ecamber outputs to directory ${DIRECTORY}" >> "$LOG" 2>&1
  echo "copy of ecamber outputs to directory ${DIRECTORY}";
  cp "${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_dna/"* "${DATA}/${DIRECTORY}/${DATA_GENE_DIR}" >> "$LOG" 2>&1
  cp "${DATA_ECAMBER_OUT}/${DIRECTORY}/output/genes_aa/"* "${DATA}/${DIRECTORY}/${DATA_PROT_DIR}" >> "$LOG" 2>&1
  cp "${DATA_ECAMBER_OUT}/${DIRECTORY}/output/clusters_table."* "${DATA}/${DIRECTORY}/${DATA_ECAMBER_OUTPUT}" >> "$LOG" 2>&1
  cp "${DATA_ECAMBER_OUT}/${DIRECTORY}/output/cluster_gene_names."* "${DATA}/${DIRECTORY}/${DATA_ECAMBER_OUTPUT}" >> "$LOG" 2>&1

  #Edit fasta format for Gembase
  echo 'format ecamber annotations to coregenome annotation input format' >> "$LOG" 2>&1
  echo 'format ecamber annotations to coregenome annotation input format'
  for f in $(ls -1 "${DATA}/${DIRECTORY}/${DATA_GENE_DIR}"); do
    filename="$(basename "${f}" '.fasta')"
    ## echo -e "${CGB_BIN}ecamber2gembase.sh ${filename} ${DATA}/${DIRECTORY}/${DATA_GENE_DIR} ${DATA}/${DIRECTORY}/${DATA_PROT_DIR}" >> "$LOG" 2>&1
    ${CGB_BIN}ecamber2gembase.sh "${filename}" "${DATA}/${DIRECTORY}/${DATA_GENE_DIR}" "${DATA}/${DIRECTORY}/${DATA_PROT_DIR}" >> "$LOG" 2>&1
    exit_status=$? ; if [ $exit_status -ne 0 ]; then echo "[ERROR] the command '${CGB_BIN}ecamber2gembase.sh $filename $DATA/$DIRECTORY/$DATA_GENE_DIR $DATA/$DIRECTORY/$DATA_PROT_DIR' exited with an error! - exit code ${exit_status}" 2>>"$LOG"; exit 1; fi
  done

  #move ecamber outputs to directory 'DATA_ECAMBER_OUTPUT'
  mv "${DATA}/${DIRECTORY}/anns_input" "${DATA}/${DIRECTORY}/anns_parsed" "${DATA}/${DIRECTORY}/blast" "${DATA}/${DIRECTORY}/cambervis" "${DATA}/${DIRECTORY}/formatdb.log" "${DATA}/${DIRECTORY}/genomes" "${DATA}/${DIRECTORY}/genomes_dbs" "${DATA}/${DIRECTORY}/output" "${DATA}/${DIRECTORY}/prodigal" "${DATA}/${DIRECTORY}/results" "${DATA}/${DIRECTORY}/strains.txt" "${DATA}/${DIRECTORY}/${DATA_ECAMBER_OUTPUT}/."

  #Delete eCamber outputs (previously copied)
  if [[ -d "${DATA}/${DIRECTORY}/${DATA_ECAMBER_OUTPUT}/anns_input" ]]; then rm -rf "${DATA}/${DIRECTORY}/${DATA_ECAMBER_OUTPUT}/anns_input"; fi
  if [[ -d "${DATA}/${DIRECTORY}/${DATA_ECAMBER_OUTPUT}/anns_parsed" ]]; then rm -rf "${DATA}/${DIRECTORY}/${DATA_ECAMBER_OUTPUT}/anns_parsed"; fi
  if [[ -d "${DATA}/${DIRECTORY}/${DATA_ECAMBER_OUTPUT}/blast" ]]; then rm -rf "${DATA}/${DIRECTORY}/${DATA_ECAMBER_OUTPUT}/blast"; fi
  if [[ -d "${DATA}/${DIRECTORY}/${DATA_ECAMBER_OUTPUT}/cambervis" ]]; then rm -rf "${DATA}/${DIRECTORY}/${DATA_ECAMBER_OUTPUT}/cambervis"; fi
  if [[ -d "${DATA}/${DIRECTORY}/${DATA_ECAMBER_OUTPUT}/genomes" ]]; then rm -rf "${DATA}/${DIRECTORY}/${DATA_ECAMBER_OUTPUT}/genomes"; fi
  if [[ -d "${DATA}/${DIRECTORY}/${DATA_ECAMBER_OUTPUT}/genomes_dbs" ]]; then rm -rf "${DATA}/${DIRECTORY}/${DATA_ECAMBER_OUTPUT}/genomes_dbs"; fi
  if [[ -d "${DATA}/${DIRECTORY}/${DATA_ECAMBER_OUTPUT}/output" ]]; then rm -rf "${DATA}/${DIRECTORY}/${DATA_ECAMBER_OUTPUT}/output"; fi
  if [[ -d "${DATA}/${DIRECTORY}/${DATA_ECAMBER_OUTPUT}/prodigal" ]]; then rm -rf "${DATA}/${DIRECTORY}/${DATA_ECAMBER_OUTPUT}/prodigal"; fi
  if [[ -d "${DATA}/${DIRECTORY}/${DATA_ECAMBER_OUTPUT}/results" ]]; then rm -rf "${DATA}/${DIRECTORY}/${DATA_ECAMBER_OUTPUT}/results"; fi
  if [[ -f "${DATA}/${DIRECTORY}/${DATA_ECAMBER_OUTPUT}/formatdb.log" ]]; then rm "${DATA}/${DIRECTORY}/${DATA_ECAMBER_OUTPUT}/formatdb.log"; fi

  #merge logs
  cat "${DATA}/${DIRECTORY}/logs/${DIRECTORY}.annotation."?.log > "${DATA}/${DIRECTORY}/logs/${DIRECTORY}.annotation.log";
  rm "${DATA}/${DIRECTORY}/logs/${DIRECTORY}.annotation."?.log;
}

main "$@"

