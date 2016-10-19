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
# rename_genome_fasta_files.sh: change names of genomic fasta files.                                                #
# Authors: Elise Larsonneur, Damien Mornico                                                                         #
#####################################################################################################################


set -o pipefail


########################################################################################
# SOURCE LIBRARY FILE                                                                  #
########################################################################################

source ${CGB_BIN}utils.sh 


########################################################################################
# FUNCTIONS                                                                            #
########################################################################################

# display_usage 
# This function displays the usage of this program.
# No parameters
function display_usage { 
  echo -e "Usage:\n$0 [fasta_directory] [new_name] [ref_name] [outfile_list] \n ex : $0 fasta/ esco NOREF renamedlist.txt\n $0 fasta/ esco ref.fasta renamedlist.txt" 
} 



# renameFiles
# This function changes names of input fasta files.
# Parameters
# - 1) A list of fasta file paths - A string.
# - 2) An integer.
# - 3) Output directory - A string containing the directory path where the fasta files are copied after renaming.
# - 4) A prefix used in the name of new fasta files - A string.
# - 5) An output file where are stored the old and new names of renamed fasta files - A string containing the file path. 
function renameFiles {

  local i_ dir_ prefixname_ outlist_ ind_ base 
  if [[ -z "$1" ]]; then echo "[ERROR] the function renameFiles expects a list of files to run" >&2; exit 1; fi 
  if [[ -z "$2" ]]; then echo "[ERROR] the function renameFiles expects an integer to run" >&2; exit 1; fi 
  if [[ -z "$3" ]]; then echo "[ERROR] the function renameFiles expects a directory to run" >&2; exit 1; fi 
  if [[ -z "$4" ]]; then echo "[ERROR] the function renameFiles expects a string to run" >&2; exit 1; fi 
  if [[ -z "$5" ]]; then echo "[ERROR] the function renameFiles expects an output file name to run" >&2; exit 1; fi
  if [[ $# -eq 0 ]] || [[ $# -gt 5 ]]; then echo "[ERROR] the function renameFiles expects 5 parameters, $# parameters given" >&2; exit 1; fi
 
  declare -a cmd=("${!1}")
  #echo "${cmd[@]}"
  i_=$2    #indice
  dir_="$3"
  prefixname_="$4"
  outlist_="$5"

  if [[ ! -e "${dir_}" ]]; then echo "[ERROR] the directory ${dir_} does not exist" >&2; exit 1; fi
  if [[ ! -s "${dir_}" ]]; then echo "[ERROR] the directory ${dir_} is empty" >&2; exit 1; fi

  for f in "${cmd[@]}"; do
    i_=$((${i_} + 1));
    ind_="$(printf "%07d" "${i_}")";	 
    base="$(basename "${f}")";	
    echo "${base}	${prefixname_}.${ind_}.c001.fasta" >> $outlist_;
    mv "$f" "${dir_}/${prefixname_}.${ind_}.c001.fasta" || { echo "[ERROR] error when changing the name of file ${f}" >&2; exit 1; }
  done
   
  if [[ ! -e "${outlist_}" ]]; then echo "[ERROR] the output file ${outlist_} does not exist" >&2; exit 1; fi
  if [[ ! -s "${outlist_}" ]]; then echo "[ERROR] the output file ${outlist_} is empty" >&2; exit 1; fi  
}



########################################################################################
# MAIN FUNCTION                                                                        #
########################################################################################
# main
# Change names of genomic fasta files.
# Parameters:
#  - 1) dir : a directory containing the input genomic fasta files to rename. A path.
#  - 2) prefixname : a four-letter string used as prefix for the new names of fasta files. A string.
#  - 3) reffastafilename : the name of the reference genomic fasta file. A string. 
#  - 4) outlist : the name of the file output by this program, the file containing the list of old and new names of each genomic fasta file. A path. 
# Ouputs :
# - a file containing the list of old and new names of each genomic fasta file.
function main {

  # if less than four arguments supplied, display usage 
  if [[  $# -ne 4 ]]; then 
    display_usage
    exit 1
  fi 
 
  # check whether user had supplied -h or --help . If yes display usage 
  if [[ "$#" == "--help" ]] || [[ "$#" == "-h" ]]; then 
    display_usage
    exit 0
  fi


  if [[ -z "$1" ]]; then echo "[ERROR] $0 expects a directory to run" >&2; exit 1; fi 
  if [[ -z "$2" ]]; then echo "[ERROR] $0 expects a string (prefix name) to run" >&2; exit 1; fi 
  if [[ -z "$3" ]]; then echo "[ERROR] $0 expects an string (fasta name) to run" >&2; exit 1; fi 
  if [[ -z "$4" ]]; then echo "[ERROR] $0 expects an output file name to run" >&2; exit 1; fi 
  dir="$1"
  prefixname="$2"         ## 4 characters 
  reffastafilename="$3"
  outlist="$4"


  if [[ ! -e "${dir}" ]]; then echo "[ERROR] the directory ${dir} does not exist" >&2; exit 1; fi
  if [[ ! -s "${dir}" ]]; then echo "[ERROR] the directory ${dir} is empty" >&2; exit 1; fi
  if [[ -e "${outlist}" ]]; then rm "${outlist}"; fi



  if [[ "${reffastafilename}" = "NOREF" ]]; then
    i=0
    cmd1=( $(find "${dir}/" -maxdepth 1 ! -type l -name "*.fasta" |sort -u) )
    renameFiles cmd1[@] $i "${dir}" "${prefixname}" "${outlist}"
  else
    i=0
    cmd1=( $(find "${dir}/" -maxdepth 1 ! -type l -name "*.fasta" |sort -u |grep "${reffastafilename}") )
    renameFiles cmd1[@] $i "${dir}" "${prefixname}" "${outlist}"

    newrefname="$(get_new_name_from_list "${outlist}" "${dir}/${reffastafilename}")"
    newrefname="$(basename "${newrefname}")";

    i=1
    cmd2=( $(find "${dir}/" -maxdepth 1 ! -type l -name "*.fasta" |sort -u |grep -v "${newrefname}") )
    renameFiles cmd2[@] $i "${dir}" "${prefixname}" "${outlist}"
  fi
}

main "$@"

