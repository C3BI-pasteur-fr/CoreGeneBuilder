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
# get_orthologs_from_prot_list.sh: computes the list of persistent genes from BBH lists.                            #
# Author: Elise Larsonneur                                                                                          #
#####################################################################################################################


set -o pipefail


########################################################################################
# FUNCTIONS                                                                            #
########################################################################################

# display_usage 
# This function displays the usage of this program.
# No parameters
function display_usage { 
  echo -e "Usage: $0 <homolog_dir> <four_letter_prefix> <dir_path> <protein_list> <cutoff>" 
}   



########################################################################################
# MAIN FUNCTION                                                                        #
########################################################################################
# main
# Computes the list of persistent genes from homolog lists.
# Parameters:
#  - 1) a string as prefix of the *.synt files (homologs).
#  - 2) a string as name of the analysis directory.
#  - 3) a directory where are input and output files. 
#  - 4) listprot : a file with extension '.listprot' containing for each remaining homologous gene (after applying 
#                  the threshold given below), the number of strains having this gene.
#  - 5) cutoff : a threshold used to output a subset of genes (persistent genes). 
# Returns:
#  - a file with the extension '.lst' containing the list of genes found in <cutoff> % of the genomes.  
#       if 'cutoff=100' the list of strict core genes is output; else if 'cutoff<100', persistent gene set is output  
function main {
  # if less than five arguments supplied, display usage 
  if [[  $# -ne 5 ]]; then 
    display_usage
    exit 1
  fi  
 
  # check whether user had supplied -h or --help . If yes display usage 
  if [[ $# == "--help" ]] || [[ $# == "-h" ]]; then 
    display_usage
    exit 0
  fi  


  echo "ls $3/../${1}*.synt"
  data=( $(ls $3/../${1}*.synt) )
  name="$2"
  listprot="$4"
  cutoff=$5

  if [[ ! -e "${listprot}" ]]; then echo "[ERROR] the input file ${listprot} does not exist" >&2; exit 1; fi
  if [[ ! -s "${listprot}" ]]; then echo "[ERROR] the input file ${listprot} is empty" >&2; exit 1; fi

  z=""
  if [[ $cutoff -lt 100 ]]; then z="0"; fi 

  if [[ "${#data[@]}" -eq "1" ]]; then
    awk '{print $1,$3, $2, $5}' ${data[0]} 1> "$3/CoreGenome-${name}-${z}${cutoff}.lst"
    exit $?
  fi

  res=""; for (( i=0; i<${#data[@]}; i++ )); do res="${res} "${data[i]}; done
  echo "${#data[@]} ${res}"

  join "${listprot}" "${data[0]}" |awk '{print $1,$3,$2,$5}' 1> "$3/firsta"  || { echo "[ERROR] joining error" >&2; exit 1; }
  i=1
  last_indice=$((${#data[@]} - 1))

  echo "Next join"
  while [[ $i -le $last_indice ]]; do
    join -a1 "$3/firsta" ${data[$i]} 1> "$3/joina";
    gawk '
      #storage of number_of_fields (NF) and each line ($0)
      FILENAME == ARGV[1] { one[FNR]=NF; }; 
      FILENAME == ARGV[2] { two[FNR]=NF; rec2[FNR]=$0 };
      END { 
        #check empty lines
        for (i=1; i<=length(one); i++) { 
            if (length(one[i]) == 0){ one[i]=0 }; 
            if (length(two[i]) == 0){ two[i]=0 };
        } 
        #and print
        for (i=1; i<=length(one); i++) { 
          if (one[i]==two[i]) { 
            print rec2[i];
          } else {
            n=split(rec2[i],a," "); 
            res = a[1]; for (j=2;j<=(n-5); j++){ res = res " " a[j]; }; print res, a[n-2]; 
          }
        }
      }' "$3/firsta" "$3/joina" 1> "$3/secondd";

    cp -f "$3/secondd" "$3/firsta";
    i=$((${i} + 1));
  done

  mv "$3/firsta" "$3/CoreGenome-${name}-${z}${cutoff}.lst" || { echo "[ERROR] unable to change name of the file '$3/firsta'"; }
  if [[ -e "$3/joina" ]]; then rm "$3/joina"; fi 
  if [[ -e "$3/secondd" ]]; then rm "$3/secondd"; fi

  #output file
  file="$3/CoreGenome-${name}-${z}${cutoff}.lst"
  if [[ ! -e $file ]]; then echo "[ERROR] the output file $file (list of core genes with threshold ${cutoff}) does not exist" >&2; exit 1; fi
  if [[ ! -s $file ]]; then echo "[ERROR] the output file $file (list of core genes with threshold ${cutoff}) is empty" >&2; exit 1; fi 
  cp -p "$3/CoreGenome-${name}-${z}${cutoff}.lst" "$3/../../." || { echo "[ERROR] unable to copy $3/CoreGenome-${name}-${z}${cutoff}.lst to directory $3/../../." >&2; exit 1; }
}

main "$@"
 
