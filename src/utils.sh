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
# utils.sh: contains utility functions.                                                                             #
# Author: Elise Larsonneur                                                                                          #
#####################################################################################################################


set -o pipefail


########################################################################################
# FUNCTIONS                                                                            #
########################################################################################

# real_path_file - this function returns the absolute path of a file.
# Parameter:
#  - a relative path of a file
# Returns:
#  - the absolute path of a file.
function real_path_file {
  echo $(cd $(dirname "$1"); pwd)/$(basename "$1");
}


# real_path_dir - this function returns the absolute path of a directory.
# Parameter:
#  - a relative path of a directory
# Returns:
#  - the absolute path of a directory.
function real_path_dir {
  echo $(cd "$1"; pwd);
}


# get_new_name_from_list - this function returns the new name of a file from a list of old and new names of files.
# Parameters:
#  - list : list of old and new names of files. A path.
#  - oldname : old name of the file
# Returns:
#  - the new name of the file.
function get_new_name_from_list {
  local list oldname oldfilename dirlist tmpfile newfilename linenb
  if [[ -z "$1" ]]; then echo "[ERROR] the function get_new_name_from_list expects a file path to run" >&2; exit 1; fi
  if [[ -z "$2" ]]; then echo "[ERROR] the function get_new_name_from_list expects a string (old fasta name) to run" >&2; exit 1; fi
  if [[ $# -ne 2 ]]; then echo "[ERROR] the function get_new_name_from_list expects 2 parameters, $# parameters given" >&2; exit 1; fi
  list="$1"		# path of the list
  oldname="$2"          # old name of fasta file 
  if [[ ! -e "${list}" ]]; then echo "[ERROR] the input file ${list} does not exist" >&2; exit 1; fi
  if [[ ! -s "${list}" ]]; then echo "[ERROR] the input file ${list} is empty" >&2; exit 1; fi

  oldfilename=$(basename "${oldname}")
  dirlist=$(dirname "${list}")
  tmpfile="${dirlist}/newtmplist.txt"
  newfilename="";
  awk '{ print $1 }' "${list}" > "${tmpfile}"
  if grep -Fxq "${oldfilename}" "${tmpfile}" ; then
    linenb=$(grep -Fxn "${oldfilename}" "${tmpfile}" |cut -f1 -d:)
    #get new filename of reference fasta file
    newfilename=$( sed "${linenb}q;d" "${list}" | awk '{ print $2 }')
  fi
  
  if [[ -e "${tmpfile}" ]]; then rm "${tmpfile}"; fi

  echo "${newfilename}";
}


# get_old_name_from_list - this function returns the old name of a file from a list of old and new names of files.
# Parameters:
#  - list : list of old and new names of files. A path.
#  - newname : new name of the file
# Returns:
#  - the old name of the file.
function get_old_name_from_list {
  local list newname newfilename dirlist tmpfile oldfilename linenb
  if [[ -z "$1" ]]; then echo "[ERROR] the function get_old_name_from_list expects a file path to run" >&2; exit 1; fi
  if [[ -z "$2" ]]; then echo "[ERROR] the function get_old_name_from_list expects a string (new fasta name) to run" >&2; exit 1; fi
  if [[ $# -ne 2 ]]; then echo "[ERROR] the function get_old_name_from_list expects 2 parameters, $# parameters given" >&2; exit 1; fi
  list="$1"		# path of the list
  newname="$2"          # new name of fasta file
  if [[ ! -e "${list}" ]]; then echo "[ERROR] the input file ${list} does not exist" >&2; exit 1; fi
  if [[ ! -s "${list}" ]]; then echo "[ERROR] the input file ${list} is empty" >&2; exit 1; fi

  newfilename=$(basename "${newname}")
  dirlist=$(dirname "${list}")
  tmpfile="${dirlist}/oldtmplist.txt"
  oldfilename="";
  awk '{ print $2 }' "${list}" > "${tmpfile}"
  if grep -Fxq "${newfilename}" "${tmpfile}" ; then
    linenb=$(grep -Fxn "${newfilename}" "${tmpfile}" |cut -f1 -d:)
    #get old filename of reference fasta file
    oldfilename=$( sed "${linenb}q;d" "${list}" | awk '{ print $1 }')
  fi
   
  if [[ -e "${tmpfile}" ]]; then rm "${tmpfile}"; fi

  echo "${oldfilename}";
}

#test 
#get_new_name_from_list $CGPIPELINE/data/klpn3/diversity/renamed_genomes_list.txt klpn.006.c01.10.fasta
#get_old_name_from_list $CGPIPELINE/data/klpn3/diversity/renamed_genomes_list.txt klpn.001.c01.00.fasta
 
