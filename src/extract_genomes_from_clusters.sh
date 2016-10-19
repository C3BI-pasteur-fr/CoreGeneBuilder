#!/bin/bash

#####################################################################################################################
# CoreGeneBuilder - extracts a core genome or a persistent genome from a set of bacterial genomes.                  #
# Authors : Elise Larsonneur, Marie Touchon, Damien Mornico, Alexis Criscuolo, Sylvain Brisse, Eduardo P. C. Rocha  #
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
# extract_genomes_from_clusters.sh: extracts one genome per cluster of distance-based hierarchical clustering.      #              
# Author: Elise Larsonneur                                                                                          #
#####################################################################################################################


set -o pipefail


########################################################################################
# GLOBAL VARIABLES                                                                     #
########################################################################################

## storage of genomes added to the output list
declare -A ADDED_GENOMES_ARRAY


########################################################################################
# FUNCTIONS                                                                            #
########################################################################################
    
# recursive_extraction_of_genomes 
# This function recursively extracts one genome per cluster of distance-based hierarchical clustering
# Parameters
# - Input file - ordered list of the genomes and their clusters. A string containing the file path.
# - Output file - list of selected genomes from the input list. A string containing the file path.
# - An integer.
function recursive_extraction_of_genomes {
    
  local file_i file_o diff 
  local last_cluster cluster genome_id 
  if [[ -z "$1" ]]; then echo "[ERROR] the function recursive_extraction_of_genomes expects an input file to run" >&2; exit 1; fi
  if [[ -z "$2" ]]; then echo "[ERROR] the function recursive_extraction_of_genomes expects an output file to run" >&2; exit 1; fi
  if [[ $# -eq 0 ]] || [ $# -gt 2 ]; then echo "[ERROR] the function recursive_extraction_of_genomes expects 3 parameters, $# parameters given" >&2; exit 1; fi
  file_i=$1
  file_o=$2
  diff=$3

  if [[ ! -e "${file_i}" ]]; then echo "[ERROR] the input file ${file_i} does not exist" >&2; exit 1; fi
  if [[ ! -s "${file_i}" ]]; then echo "[ERROR] the input file ${file_i} is empty" >&2; exit 1; fi
    
  last_cluster="0"

  while read line 
  do
    cluster=$(echo "${line}" | awk '{print $2}')    ## cluster no
    genome_id=$(echo "${line}" | awk '{print $1}')  ## genome id        

    if [[ "${diff}" -eq "0" ]]; then
      break
    elif [[ "${diff}" -ne "0" ]]; then
      if [[ "${last_cluster}" -ne "${cluster}" ]] && [[ "${ADDED_GENOMES_ARRAY[${genome_id}]}" -eq "0" ]] ; then
        echo "${line}" >> "${file_o}"  
        ADDED_GENOMES_ARRAY[${genome_id}]=3 
        #echo "keep line : $line"
        diff=$((${diff} - 1))
        last_cluster="${cluster}"
       fi
     fi
  done < "${file_i}"
}



# get_organism_by_cluster
# This function extract one genome per cluster of distance-based hierarchical clustering
# Parameters
# - Input file - ordered list of the genomes and their clusters. A string containing the file path.
# - Number of genomes to select among those of the list. An integer.
# - Output file - list of selected genomes from the input list. A string containing the file path.
function get_organism_by_cluster {
 
  local filein selected_nb fileout
  local no_cluster orga_to_keep_by_cluster last_cluster keep_nb cluster genome_id real_nb_to_keep
  if [[ -z "$1" ]]; then echo "[ERROR] the function get_organism_by_cluster expects an input file to run" >&2; exit 1; fi
  if [[ -z "$2" ]]; then echo "[ERROR] the function get_organism_by_cluster expects an integer to run" >&2; exit 1; fi
  if [[ -z "$3" ]]; then echo "[ERROR] the function get_organism_by_cluster expects an output file to run" >&2; exit 1; fi
  if [[ $# -eq 0 ]] || [[ $# -gt 3 ]]; then echo "[ERROR] the function get_organism_by_cluster expects 3 parameters, $# parameters given" >&2; exit 1; fi
  filein="$1"
  selected_nb="$2"
  fileout="$3"

  if [[ ! -e "${filein}" ]]; then echo "[ERROR] the input file ${filein} does not exist" >&2; exit 1; fi
  if [[ ! -s "${filein}" ]]; then echo "[ERROR] the input file ${filein} is empty" >&2; exit 1; fi

  no_cluster=$(awk '{print $2}' "${filein}" | sort -u | wc -l)
  echo "Number of different clusters : ${no_cluster}"

  orga_to_keep_by_cluster=$(echo "${selected_nb}/${no_cluster}" | bc)
  echo "Number of organism to select by cluster : ${orga_to_keep_by_cluster}"




  ## FIRST PART  - get the first n="$orga_to_keep_by_cluster" genomes for each cluster

  last_cluster=$(head -n 1 "${filein}" | awk '{print $2}')
  keep_nb=0

  while read line; do
    cluster=$(echo "${line}" | awk '{print $2}')    ## cluster no
    genome_id=$(echo "${line}" | awk '{print $1}')  ## genome id
   
    if [[ "${cluster}" -eq "${last_cluster}" ]] && [[ "${keep_nb}" -lt "${orga_to_keep_by_cluster}" ]]; then
      echo "${line}" >> "${fileout}"
      ADDED_GENOMES_ARRAY[${genome_id}]=1
      #echo "keep line : ${line}"
      keep_nb=$((${keep_nb} + 1))
       
    elif [[ "${cluster}" -ne "${last_cluster}" ]]; then
      echo "${line}" >> "${fileout}"
      ADDED_GENOMES_ARRAY[${genome_id}]=1
      #echo "keep line : ${line}"
      last_cluster=$(echo "${line}" | awk '{print $2}')
      keep_nb=1
    fi   
  done < "${filein}"

  real_nb_to_keep=$(cat "${fileout}" |sort -u | wc -l)
  echo "Number of organism really kept : ${real_nb_to_keep}"

  diff=$((${selected_nb} - ${real_nb_to_keep}))

  if [[ "${real_nb_to_keep}" -lt "${selected_nb}" ]]; then
    echo "organisms missed"
    echo "${diff} more organims needed"
  fi




  ### SECOND PART - get the next genomes, one genome per cluster at a time

  #call the function 'recursive_extraction_of_genomes'
  i=0

  while [[ "${diff}" -gt "0" ]]; do
    recursive_extraction_of_genomes "${filein}" "${fileout}"   ## version with $diff and $ADDED_GENOMES_ARRAY as global variables
    i=$(echo "${i}+1" | bc)
    n=0
    # display the content of associative array
    for k in "${!ADDED_GENOMES_ARRAY[@]}"; do
      n=$(echo "${n}+1" | bc)
      #echo "${k}"
    done
    #echo "len dic : ${n}"
  done

  #echo "Number of iterations : ${i}"
  #echo "Dictionnary :"
  #for k in "${!ADDED_GENOMES_ARRAY[@]}"; do
  #     echo "${k}: ${ADDED_GENOMES_ARRAY[$k]}"
  #done
}



# display_usage 
# This function displays the usage of this program.
# No parameters
function display_usage {
  echo -ne "Usage: $0 <sorted genome cluster file> <nb selection>\n" 1>&2
}



########################################################################################
# MAIN FUNCTION                                                                        #
########################################################################################
# main
# Parameters
# - 1) Input file - ordered list of the genomes and their clusters. A string containing the file path.
# - 2) An integer.
# - 3) Output file - list of selected genomes from the input list. A string containing the file path.
# Outputs 
# - a file containing the list of selected genomes from the input list.
function main {
  if [[ $# -ne 3 ]]; then
    display_usage
    exit 1
  fi

  # check whether user had supplied -h or --help . If yes display usage 
  if [[ $# == "--help" ]] || [[ $# == "-h" ]]; then
    display_usage
    exit 0
  fi

  infile="$1"
  echo "Input file : ${infile}"
  nb_selection="$2"
  echo "Number of organism to select : ${nb_selection}"
  outfile="$3"
  echo "Output file : ${outfile}"
  get_organism_by_cluster "${infile}" "${nb_selection}" "${outfile}"
}

main "$@"
