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
# run_COREGENOME.sh: computes a core or persistent gene set from input proteomes.                                   #
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



# createCoreGeneFasta
# This function extracts the core gene sequences of each genome and writes them into fasta files (.gen and .prt), from a input list of core gene ids.
#   There are one nucleic fasta file and one amino-acid fasta file for each genome.
# Parameters
# - 1) Input file - A string containing the file path where the core gene ids are stored.
# - 2) Directory - A string containing the directory path where the core gene sequences will be outputed.
# - 3) Input directory 'genes' - A string containing the directory path where the annotation files are stored (nucleic acid format '.gen').
# - 4) Input directory 'proteins' - A string containing the directory path where the annotation files are stored (amino-acid format '.prt').
# - 5) Input file - A string containing the file path where the list of genomes is stored.
# Outputs
# - Fasta files containing amino-acid or nucleic acid sequences of core genes for each genome of the input list of genomes.
function createCoreGeneFasta {

  local cglist coredir genedir protdir genomeslist;
  local prefix pattern listcg allgenfasta cggenfastaout allprtfasta cgprtfastaout;
  if [[ -z "$1" ]]; then echo '[ERROR] the function createCoreGeneFasta expects an input file (core gene list) to run' 2>>"${LOG}"; exit 1; fi
  if [[ -z "$2" ]]; then echo '[ERROR] the function createCoreGeneFasta expects a directory (core_genome) to run' 2>>"${LOG}"; exit 1; fi
  if [[ -z "$3" ]]; then echo '[ERROR] the function createCoreGeneFasta expects a directory (genes) to run' 2>>"${LOG}"; exit 1; fi
  if [[ -z "$4" ]]; then echo '[ERROR] the function createCoreGeneFasta expects a directory (proteins) to run' 2>>"${LOG}"; exit 1; fi
  if [[ -z "$5" ]]; then echo '[ERROR] the function createCoreGeneFasta expects an input file (genome list) to run' 2>>"${LOG}"; exit 1; fi
  if [[ $# -ne 5 ]]; then echo "[ERROR] the function createCoreGeneFasta expects 5 parameters, $# parameters given" 2>>"${LOG}"; exit 1; fi
  cglist="$1"
  coredir="$2"
  genedir="$3"
  protdir="$4"
  genomeslist="$5"
  if [[ ! -e "${cglist}" ]]; then echo "[ERROR] the input file ${cglist} does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${cglist}" ]]; then echo "[ERROR] the input file ${cglist} is empty" 2>>"${LOG}"; exit 1; fi
  if [[ ! -e "${coredir}" ]]; then echo "[ERROR] the input directory ${coredir} does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${coredir}" ]]; then echo "[ERROR] the input directory ${coredir} is empty" 2>>"${LOG}"; exit 1; fi
  if [[ ! -e "${genedir}" ]]; then echo "[ERROR] the input directory ${genedir} does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${genedir}" ]]; then echo "[ERROR] the input directory ${genedir} is empty" 2>>"${LOG}"; exit 1; fi
  if [[ ! -e "${protdir}" ]]; then echo "[ERROR] the input directory ${protdir} does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${protdir}" ]]; then echo "[ERROR] the input directory ${protdir} is empty" 2>>"${LOG}"; exit 1; fi
  if [[ ! -e "${genomeslist}" ]]; then echo "[ERROR] the input file ${genomeslist} does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${genomeslist}" ]]; then echo "[ERROR] the input file ${genomeslist} is empty" 2>>"${LOG}"; exit 1; fi

  mkdir -p "${coredir}/${core_genes_dir}" >> "${LOG}" 2>&1

  #without ref 
  for fna in $(tail -n+2 "${genomeslist}"); do 
    prefix="${fna%.*}"
    pattern="$(echo "${prefix}" |awk '{ gsub("\\.","",$1); print toupper($1) }')"

    listcg="${coredir}/${core_genes_dir}/${fna}.coregenes.${z}${prct}.lst"
    allgenfasta="${genedir}/${fna}.gen"
    cggenfastaout="${coredir}/${core_genes_dir}/${fna}.coregenes.${z}${prct}.gen"
    allprtfasta="${protdir}/${fna}.prt"
    cgprtfastaout="${coredir}/${core_genes_dir}/${fna}.coregenes.${z}${prct}.prt"

    while read -r line; do
      echo "${line}" |grep -P "${pattern}\D+_\d+" -o  >> "${listcg}"
    done < ${cglist}

    echo "${CGB_BIN}fasta_keep_seq_from_list.py ${allgenfasta} ${listcg} ${cggenfastaout}" >> "${LOG}" 2>&1
    ${CGB_BIN}fasta_keep_seq_from_list.py "${allgenfasta}" "${listcg}" "${cggenfastaout}" >> "${LOG}" 2>&1
    exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] the command '${CGB_BIN}fasta_keep_seq_from_list.py ${allgenfasta} ${listcg} ${cggenfastaout}' exited with an error! - exit code ${exit_status}" 2>>"${LOG}"; exit 1; fi
    echo "${CGB_BIN}fasta_keep_seq_from_list.py ${allprtfasta} ${listcg} ${cgprtfastaout}" >> "${LOG}" 2>&1
    ${CGB_BIN}fasta_keep_seq_from_list.py "${allprtfasta}" "${listcg}" "${cgprtfastaout}" >> "${LOG}" 2>&1
    exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] the command '${CGB_BIN}fasta_keep_seq_from_list.py ${allprtfasta} ${listcg} ${cgprtfastaout}' exited with an error! - exit code ${exit_status}" 2>>"${LOG}"; exit 1; fi
  
    #regenerate '.lst' with complete header of fasta files (.prt and .gen)
    grep -P '^>' "${cgprtfastaout}" |awk '{gsub(">","",$0); print $0}' > "${listcg}"

  done
}



# extractProteinFromFasta
# This function extracts the fasta sequence of one protein id.
# Parameters
# - 1) Protein id - A string.
# - 2) Input file - A string containing a fasta file path.
# - 3) Output file - A string containing a fasta file path.
# Appends
# - The fasta sequence associated to the input protein id into the given output fasta file.
function extractProteinFromFasta {
  local prot fasta outfasta
  if [[ -z "$1" ]]; then echo '[ERROR] the function extractProteinFromFasta expects a string (protein name) to run' 2>>"${LOG}"; exit 1; fi
  if [[ -z "$2" ]]; then echo '[ERROR] the function extractProteinFromFasta expects an input fasta path to run' 2>>"${LOG}"; exit 1; fi
  if [[ -z "$3" ]]; then echo '[ERROR] the function extractProteinFromFasta expects an output fasta path to run' 2>>"${LOG}"; exit 1; fi
  prot="$1"
  fasta="$2"
  outfasta="$3"
  if [[ ! -e "${fasta}" ]]; then echo "[ERROR] the input file ${fasta} does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${fasta}" ]]; then echo "[ERROR] the input file ${fasta} is empty" 2>>"${LOG}"; exit 1; fi
  
  grep -A 1 -m 1 "^>${prot}" "${fasta}" 1>> "${outfasta}";
  
  if [[ ! -e "${outfasta}" ]]; then echo "[ERROR] the input file ${outfasta} does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${outfasta}" ]]; then echo "[ERROR] the input file ${outfasta} is empty" 2>>"${LOG}"; exit 1; fi
}



# createCoreGeneByGeneFasta
# This function extracts the core gene sequences of each genome and writes them into fasta files (.gen and .prt), from a input list of core gene ids.
#    There are one nucleic fasta file and one amino-acid fasta file for each core gene. Each of these output fasta file contains the allelic sequences of the core gene.
# Parameters
# - 1) Input file - A string containing the file path where the core gene ids are stored.
# - 2) Prefix name of the output fasta files - A string of four letters. 
# - 3) Directory - A string containing the directory path where the core gene sequences will be outputed.
# - 4) Input directory 'genes' - A string containing the directory path where the annotation files are stored (nucleic acid format '.gen').
# - 5) Input directory 'proteins' - A string containing the directory path where the annotation files are stored (amino-acid format '.prt').
# Outputs
# - Fasta files containing amino-acid or nucleic acid allelic sequences for each core gene.
function createCoreGeneByGeneFasta {

  local cglist name cgdir genedir protdir;
  local PATTERN cgbycgdir;
  local genomes ref_genome proteins ref_prot prefix cggenfastaout cgprtfastaout prot_array_len;
  local prot genome allgenfastain allprtfastain;
  if [[ -z "$1" ]]; then echo '[ERROR] the function createCoreGeneByGeneFasta expects an input file (core gene list) to run' 2>>"${LOG}"; exit 1; fi
  if [[ -z "$2" ]]; then echo '[ERROR] the function createCoreGeneByGeneFasta expects a string (name) to run' 2>>"${LOG}"; exit 1; fi
  if [[ -z "$3" ]]; then echo '[ERROR] the function createCoreGeneByGeneFasta expects a directory (core_genome) to run' 2>>"${LOG}"; exit 1; fi
  if [[ -z "$4" ]]; then echo '[ERROR] the function createCoreGeneByGeneFasta expects a directory (genes) to run' 2>>"${LOG}"; exit 1; fi
  if [[ -z "$5" ]]; then echo '[ERROR] the function createCoreGeneByGeneFasta expects a directory (proteins) to run' 2>>"${LOG}"; exit 1; fi
  if [[ $# -ne 5 ]]; then echo "[ERROR] the function createCoreGeneByGeneFasta expects 5 parameters, $# parameters given" 2>>"${LOG}"; exit 1; fi
  cglist="$1"
  name="$2"
  cgdir="$3"
  genedir="$4"
  protdir="$5"
  if [[ ! -e "${cglist}" ]]; then echo "[ERROR] the input file ${cglist} does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${cglist}" ]]; then echo "[ERROR] the input file ${cglist} is empty" 2>>"${LOG}"; exit 1; fi
  if [[ ! -e "${cgdir}" ]]; then echo "[ERROR] the input directory ${cgdir} does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${cgdir}" ]]; then echo "[ERROR] the input directory ${cgdir} is empty" 2>>"${LOG}"; exit 1; fi
  if [[ ! -e "${genedir}" ]]; then echo "[ERROR] the input directory ${genedir} does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${genedir}" ]]; then echo "[ERROR] the input directory ${genedir} is empty" 2>>"${LOG}"; exit 1; fi
  if [[ ! -e "${protdir}" ]]; then echo "[ERROR] the input directory ${protdir} does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${protdir}" ]]; then echo "[ERROR] the input directory ${protdir} is empty" 2>>"${LOG}"; exit 1; fi

  PATTERN="$(echo "${name}" |awk '{print toupper($1)}')";
  cgbycgdir="${cgdir}/core_genes_by_gene";
  mkdir -p "${cgbycgdir}";

  while read -r line; do
    genomes=( $(echo "${line}" |grep -P "${PATTERN}\d+" -o |sort -u |awk -v pattern="${PATTERN}" '{gsub(pattern,pattern".",$1); print tolower($1)}') );
    ref_genome="${genomes[0]}";
    proteins=( $(echo "${line}" |grep -P "${PATTERN}\d+\D+_\d+" -o |sort -u) );
    ref_prot="${proteins[0]}";
    prefix="${ref_genome}.${ref_prot}";
    cggenfastaout="${cgbycgdir}/${prefix}.gen";
    cgprtfastaout="${cgbycgdir}/${prefix}.prt";
    prot_array_len="${#proteins[@]}";
    prot_array_len=$((${prot_array_len} - 1));
 
    for i in $(seq 0 $prot_array_len); do 
      prot="${proteins[$i]}";
      genome="${genomes[$i]}"; 
      allgenfastain="${genedir}/${genome}.c001.gen";
      allprtfastain="${protdir}/${genome}.c001.prt";
      extractProteinFromFasta "${prot}" "${allgenfastain}" "${cggenfastaout}";
      extractProteinFromFasta "${prot}" "${allprtfastain}" "${cgprtfastaout}";
    done
  done < ${cglist}
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
  LOG="${DATA}/${DIRECTORY}/logs/${DIRECTORY}.coregenome.2.log";
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




  SYNTMAXDIST="${SYNTRADIUSSIZE}"
  readonly SYNTMAXDIST


  echo 'building core genome...' >> "${LOG}" 2>&1
  echo 'building core genome...';

  DATATEMPRLK="$(readlink -f "${DATA}/${DIRECTORY}/${DATA_CORE_DIR}/${temp_dir}")"

  echo '=====================================================' >> "${LOG}" 2>&1
  echo "${CGB_BIN}do_core_genome ${DIRECTORY} ${IDPRCT} ${PROTLENGTHRATIO} ${DATATEMPRLK} ${SYNTGENESUM} ${SYNTRADIUSSIZE} ${SYNTMAXDIST} ${OPSCAN_MATRIX}" >> "${LOG}" 2>&1
  echo "${CGB_BIN}do_core_genome ${DIRECTORY} ${IDPRCT} ${PROTLENGTHRATIO} ${DATATEMPRLK} ${SYNTGENESUM} ${SYNTRADIUSSIZE} ${SYNTMAXDIST} ${OPSCAN_MATRIX}";
  ${CGB_BIN}do_core_genome "${DIRECTORY}" "${IDPRCT}" "${PROTLENGTHRATIO}" "${DATATEMPRLK}" "${SYNTGENESUM}" "${SYNTRADIUSSIZE}" "${SYNTMAXDIST}" "${OPSCAN_MATRIX}" >> "${LOG}" 2>&1
  exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] the command '${CGB_BIN}do_core_genome ${DIRECTORY} ${IDPRCT} ${PROTLENGTHRATIO} ${DATATEMPRLK} ${SYNTGENESUM} ${SYNTRADIUSSIZE} ${SYNTMAXDIST} ${OPSCAN_MATRIX}' exited with an error! - exit code ${exit_status}" 2>>"${LOG}"; exit 1; fi

  # STEP-1 - 
  #get_ortholog_distribution.sh
  dir="$(readlink -f "${DATA}/${DIRECTORY}/${DATA_CORE_DIR}/${tmp_dir}")"  ## readlink => absolute path of this directory
  genomeslist="${DATA}/${DIRECTORY}/${DATA_CORE_DIR}/Genomes-${DIRECTORY}.lst"
  echo "${CGB_BIN}get_ortholog_distribution.sh ${dir} ${genomeslist} ${DATA}/${DIRECTORY}/${DATA_CORE_DIR}" >> "${LOG}" 2>&1
  ${CGB_BIN}get_ortholog_distribution.sh "${dir}" "${genomeslist}" "${DATA}/${DIRECTORY}/${DATA_CORE_DIR}" >> "${LOG}" 2>&1 ### output file '$refname.histo'
  exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] the command '${CGB_BIN}get_ortholog_distribution.sh ${dir} ${genomeslist} ${DATA}/${DIRECTORY}/${DATA_CORE_DIR}' exited with an error! - exit code ${exit_status}" 2>>"${LOG}"; exit 1; fi


  # STEP-2
  ### get number of genomes
  nb="$(wc -l "${genomeslist}" | awk '{ print $1 }')"
  ### orthologs will be present in $CUTOFF % of the proteomes
  CUTOFF="$(echo "${nb}*${CGGENOMEPRCT}/100" | bc -l)"  # float value


  z=''; boolean="$(echo "${CGGENOMEPRCT} < 100.0" |bc -l)"; if [[ $boolean -eq 1 ]]; then z=0; fi

  ref="$(head -n 1 "${genomeslist}")"  						              # reference genome is always at the first line in this file  ; ex: "limo.001.c001"
  refname="${ref}"					  			              # remove the last 3 characters ; ex: "limo.001.c01"
  listprot="${DATA}/${DIRECTORY}/${DATA_CORE_DIR}/${refname}.${z}${CGGENOMEPRCT}.listprot"      # list of ortholog proteins ; ex: "limo.001.c01.095.listprot"
  species="${NAME}"                                                                             #=> name ; ex: "limo"
  outputdir="${DATA}/${DIRECTORY}/${DATA_CORE_DIR}/${output_dir}"                               #=> relative path for finding $refname*.synt
  prct="${CGGENOMEPRCT}"

  echo "awk -v cutoff=${CUTOFF} '{ if (\$2 >= cutoff) print \$1; }' ${DATA}/${DIRECTORY}/${DATA_CORE_DIR}/${refname}.histo" 1> "${listprot}" 2>>"${LOG}";
  awk -v cutoff="${CUTOFF}" '{ if ($2 >= cutoff) print $1; }' "${DATA}/${DIRECTORY}/${DATA_CORE_DIR}/${refname}.histo" 1> "${listprot}" 2>>"${LOG}" ;
  exit_status=$?
  if [[ ! -e "${listprot}" ]]; then echo "[ERROR] the input file ${listprot} does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${listprot}" ]]; then echo "[ERROR] the input file ${listprot} is empty" 2>>"${LOG}"; exit 1; fi
  if [[ $exit_status -ne 0 ]]; then 
    echo 'No core genes found, restart program with new parameters !' 2>>"${LOG}";
    echo 'No core genes found, restart program with new parameters !';
    exit 1
  fi

  # STEP-3  - get list of ortholog genes shared by >= $CGGENOMEPRCT of the genomes
  echo "${CGB_BIN}get_orthologs_from_prot_list.sh ${refname} ${DIRECTORY} ${outputdir} ${listprot} ${prct}" >> "${LOG}" 2>&1
  ${CGB_BIN}get_orthologs_from_prot_list.sh "${refname}" "${DIRECTORY}" "${outputdir}" "${listprot}" "${prct}" >> "${LOG}" 2>&1
  exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] the command '${CGB_BIN}get_orthologs_from_prot_list.sh ${refname} ${DIRECTORY} ${outputdir} ${listprot} ${prct}' exited with an error! - exit code ${exit_status}" 2>>"${LOG}"; exit 1; fi
  ## output file : $DATA/$DIRECTORY/$DATA_CORE_DIR/$output_dir/CoreGenome-<name>-<similarity_prct>.lst

  echo "cp -p ${DATA}/${DIRECTORY}/${DATA_CORE_DIR}/${output_dir}/CoreGenome-${DIRECTORY}*.lst" "${DATA}/${DIRECTORY}/${DATA_CORE_DIR}/." >> "${LOG}" 2>&1
  cp -p "${DATA}/${DIRECTORY}/${DATA_CORE_DIR}/${output_dir}/CoreGenome-${DIRECTORY}"*.lst "${DATA}/${DIRECTORY}/${DATA_CORE_DIR}/." >> "${LOG}" 2>&1;



  coredir="${DATA}/${DIRECTORY}/${DATA_CORE_DIR}"
  genedir="${DATA}/${DIRECTORY}/${DATA_GENE_DIR}"
  protdir="${DATA}/${DIRECTORY}/${DATA_PROT_DIR}"



  # format fasta (.prt and .gen) - DOS2UNIX then print one sequence per line of FASTA (instead of multi-line per sequence) 
  for prt in $(ls "${protdir}/"*.prt); do
    tr -d '\15\32' < "${prt}" |awk 'BEGIN{seq=""} !/^>/ {seq=seq$0} /^>/ {if(seq!=""){print seq;seq=""}{print $0}} END{print seq}' |grep -v "^$" > "${prt}.tmp"
    mv "${prt}.tmp" "${prt}"
  done

  for gen in $(ls "${genedir}/"*.gen); do
    tr -d '\15\32' < "${gen}" |awk 'BEGIN{seq=""} !/^>/ {seq=seq$0} /^>/ {if(seq!=""){print seq;seq=""}{print $0}} END{print seq}' |grep -v "^$" > "${gen}.tmp"
    mv "${gen}.tmp" "${gen}"
  done



  # STEP-4 - extract core genes into fasta file
  echo 'extract core genes of reference isolate' >> "${LOG}" 2>&1
  echo 'extract core genes of reference isolate';

  core_genes_dir='core_genes_by_genome'
  mkdir -p "${coredir}/${core_genes_dir}" >> "${LOG}" 2>&1

  listcoregenes="${coredir}/CoreGenome-${DIRECTORY}-${z}${prct}.lst"          # listcoregenes="$DATA/$DIRECTORY/$DATA_CORE_DIR/CoreGenome-limo-095.lst"
  listcg="${coredir}/${core_genes_dir}/${refname}.coregenes.${z}${prct}.lst"
  awk '{ print $1 }' "${listcoregenes}" 1> "${listcg}" 2>>"${LOG}" ;


  allgenfasta="${genedir}/${refname}.gen"
  cggenfastaout="${coredir}/${core_genes_dir}/${refname}.coregenes.${z}${prct}.gen"
  allprtfasta="${protdir}/${refname}.prt"
  cgprtfastaout="${coredir}/${core_genes_dir}/${refname}.coregenes.${z}${prct}.prt"

  # allgenesfasta="${ref}.gen"
  # allgenesfasta="limo.001.c001.gen"

  echo "${CGB_BIN}fasta_keep_seq_from_list.py ${allgenfasta} ${listcg} ${cggenfastaout}" >> "${LOG}" 2>&1
  ${CGB_BIN}fasta_keep_seq_from_list.py "${allgenfasta}" "${listcg}" "${cggenfastaout}" >> "${LOG}" 2>&1
  exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] the command '${CGB_BIN}fasta_keep_seq_from_list.py ${allgenfasta} ${listcg} ${cggenfastaout}' exited with an error! - exit code $exit_status" 2>>"${LOG}"; exit 1; fi
  echo "${CGB_BIN}fasta_keep_seq_from_list.py ${allprtfasta} ${listcg} ${cgprtfastaout}" >> "${LOG}" 2>&1
  ${CGB_BIN}fasta_keep_seq_from_list.py "${allprtfasta}" "${listcg}" "${cgprtfastaout}" >> "${LOG}" 2>&1
  exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] the command '${CGB_BIN}fasta_keep_seq_from_list.py ${allprtfasta} ${listcg} ${cgprtfastaout}' exited with an error! - exit code ${exit_status}" 2>>"${LOG}"; exit 1; fi

  #regenerate '.lst' with complete header of fasta files (.prt and .gen)
  grep -P '^>' "${cgprtfastaout}" |awk '{gsub(">","",$0); print $0}' > "${listcg}"


  # extract core gene sequences for each genome to fasta files (.gen and .prt)
  echo 'extract core genes for each other isolate' >> "${LOG}" 2>&1
  echo 'extract core genes for each other isolate';

  cglist="${listcoregenes}"
  echo "createCoreGeneFasta ${cglist} ${coredir} ${genedir} ${protdir} ${genomeslist}" >> "${LOG}" 2>&1
  createCoreGeneFasta "${cglist}" "${coredir}" "${genedir}" "${protdir}" "${genomeslist}"



  # STEP 5 - create core genes fasta by gene

  cgbycgdir="${coredir}/core_genes_by_gene"
  if [[ -d "${cgbycgdir}" ]]; then rm -r "${cgbycgdir}"; fi

  echo 'extract core genes by gene' >> "${LOG}" 2>&1
  echo 'extract core genes by gene';

  echo "createCoreGeneByGeneFasta ${cglist} ${NAME} ${coredir} ${genedir} ${protdir}" >> "${LOG}" 2>&1;
  createCoreGeneByGeneFasta "${cglist}" "${NAME}" "${coredir}" "${genedir}" "${protdir}";
}

main "$@"

