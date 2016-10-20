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
# init_ANNOTATION.sh: initial steps before genome annotation.                                                       #
# Authors: Elise Larsonneur, Damien Mornico                                                                         #
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

LOG='';



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



# removeGenbankIdFromList
# This function removes genbank entries of contig to remove, using contig ids.
# Parameters :
# - 1) Input file - A string containing the file path where the contig ids to remove are stored.
# - 2) Input file - A string containing a genbank file path.
# - 3) Output file -  A string containing an output genbank file path.
# Outputs
# - A genbank file where some contig entries have been removed.
function removeGenbankIdFromList {

  local contig_id_list_file genbank_file genbank_file_out
  local nb_lines line_nb_stop line_nb_start locus_precedent start_precedent oldid search searchnext searchnblines search_start search_next_start search_stop sedcmd
  if [[ -z "$1" ]]; then echo '[ERROR] the function removeGenbankIdFromList expects an input file containing contig ids to run' >> "${LOG}" 2>&1; exit 1; fi
  if [[ -z "$2" ]]; then echo '[ERROR] the function removeGenbankIdFromList expects an input genbank file to run' >> "${LOG}" 2>&1; exit 1; fi
  if [[ -z "$3" ]]; then echo '[ERROR] the function removeGenbankIdFromList expects an output genbank file name to run' >> "${LOG}" 2>&1; exit 1; fi
  if [[ $# -eq 0 ]] || [[ $# -gt 3 ]]; then echo "[ERROR] the function removeGenbankIdFromList expects 3 parameters, $# parameters given" >> "${LOG}" 2>&1; exit 1; fi
  contig_id_list_file="$1"   # path of the list
  genbank_file="$2"          # path of input genbank file
  genbank_file_out="$3"      # path of output modified genbank file

  if [[ ! -e "${contig_id_list_file}" ]]; then echo "[ERROR] the input file ${contig_id_list_file} does not exist" >> "${LOG}" 2>&1; exit 1; fi
  if [[ ! -s "${contig_id_list_file}" ]]; then echo "[ERROR] the input file ${contig_id_list_file} is empty" >> "${LOG}" 2>&1; exit 1; fi
  if [[ ! -e "${genbank_file}" ]]; then echo "[ERROR] the input file ${genbank_file} does not exist" >> "${LOG}" 2>&1; exit 1; fi
  if [[ ! -s "${genbank_file}" ]]; then echo "[ERROR] the input file ${genbank_file} is empty" >> "${LOG}" 2>&1; exit 1; fi


  cp -p "${genbank_file}" "${genbank_file}.tmp" >> "${LOG}" 2>&1
  genbank_file="${genbank_file}.tmp"

  while read line; do
    nb_lines="$(wc -l "${genbank_file}" |awk '{print $1}')"
    line_nb_stop=0
    line_nb_start=0
    locus_precedent=""
    start_precedent=0

    oldid="$(echo "${line}" |awk '{print $1}')"  #ex : NC_009648
 
    search="$(grep -Fn 'LOCUS' "${genbank_file}" |grep -A 1 "${oldid}" | head -n 1)"
    searchnext="$(grep -Fn 'LOCUS' "${genbank_file}" |grep -A 1 "${oldid}" |tail -n 1)"
    searchnblines="$(grep -Fn 'LOCUS' "${genbank_file}" |grep -A 1 "${oldid}" |wc -l |awk '{print $1}')" 

    search_start="$(echo "${search}" |awk 'BEGIN{FS=":"}; {print $1}')"
    search_next_start="$(echo "${searchnext}" |awk 'BEGIN{FS=":"}; {print $1}')"
    search_stop=$((${search_next_start} - 1))
    if [[ $searchnblines -eq 1 ]]; then search_stop="${nb_lines}"; fi

    # old_id start_line stop_line
    #echo $oldid $search_start $search_stop;
    #remove the LOCUS from genbank file starting from line 'search_start' and ending at line 'search_stop'
    sedcmd="sed '${search_start},${search_stop}d' ${genbank_file} 1> ${genbank_file}.removed 2>>${LOG}"
    echo ${sedcmd} |sh
    if [[ "$?" -ne 0 ]]; then echo "[ERROR] Failure to remove filtered out contigs from the annotation file ${genbank_file}"; exit 1; fi
    mv "${genbank_file}.removed" "${genbank_file}" >> "${LOG}" 2>&1
  done < ${contig_id_list_file}

  mv "${genbank_file}" "${genbank_file_out}" >> "${LOG}" 2>&1
  if [[ ! -e "${genbank_file_out}" ]]; then echo "[ERROR] the output file ${genbank_file_out} does not exist" >> "${LOG}" 2>&1; exit 1; fi
  if [[ ! -s "${genbank_file_out}" ]]; then echo "[ERROR] the output file ${genbank_file_out} is empty" >> "${LOG}" 2>&1; exit 1; fi
}



# replaceGenbankId
# This function replaces refseq contig ids by ecamber contig ids in the reference annotation genbank file.
# Parameters 
# - 1) Input file - A string containing the file path where the list of contig ids is recorded.
# - 2) Input file - A string containing the input genbank file path with contig id genbank nomenclature.
# - 3) Output file - A string containing the output genbank file path with contig id ecamber nomenclature.
# Outputs
# - The output genbank file where the genbank contig ids have been replaced by ecamber contig ids. 
function replaceGenbankId {

  local contig_id_list_file genbank_file genbank_file_out
  local oldid newid result linenumber contigid
  if [[ -z "$1" ]]; then echo '[ERROR] the function replaceGenbankId expects an input file containing contig ids to run' >> "${LOG}" 2>&1; exit 1; fi
  if [[ -z "$2" ]]; then echo '[ERROR] the function replaceGenbankId expects an input genbank file to run' >> "${LOG}" 2>&1; exit 1; fi
  if [[ -z "$3" ]]; then echo '[ERROR] the function replaceGenbankId expects an output genbank file name to run' >> "${LOG}" 2>&1; exit 1; fi
  if [[ $# -eq 0 ]] || [[ $# -gt 3 ]]; then echo "[ERROR] the function replaceGenbankId expects 3 parameters, $# parameters given" >> "${LOG}" 2>&1; exit 1; fi
  contig_id_list_file="$1"  #contains old and new contigs id (old : refseq id ; new : ecamber id)
  genbank_file="$2"
  genbank_file_out="$3"
  
  if [[ ! -e "${contig_id_list_file}" ]]; then echo "[ERROR] the input file ${contig_id_list_file} does not exist" >> "${LOG}" 2>&1; exit 1; fi
  if [[ ! -s "${contig_id_list_file}" ]]; then echo "[ERROR] the input file ${contig_id_list_file} is empty" >> "${LOG}" 2>&1; exit 1; fi
  if [[ ! -e "${genbank_file}" ]]; then echo "[ERROR] the input file ${genbank_file} does not exist" >> "${LOG}" 2>&1; exit 1; fi
  if [[ ! -s "${genbank_file}" ]]; then echo "[ERROR] the input file ${genbank_file} is empty" >> "${LOG}" 2>&1; exit 1; fi

  cp -p "${genbank_file}" "${genbank_file_out}" >> "${LOG}" 2>&1

  while read line; do
    oldid="$(echo "${line}" |awk '{print $1}')"
    newid="$(echo "${line}" |awk '{print $NF}')"
    #echo "oldid:$oldid;newid:$newid"
    result="$(grep -Hn 'ACCESSION' "${genbank_file_out}" |grep "${oldid}")"
    if [[ -n "${result}" ]]; then
      linenumber="$(echo "${result}" | awk 'BEGIN{FS=":"}; {print $2}')"
      contigid="$(echo "${result}" | awk 'BEGIN{FS=":"}; {print $3}'| awk '{print $2}')"
      #echo "linenumber:$linenumber; contigid:$contigid"
      awk -v startt="${contigid}" -v stopp="${newid}" -v lines="${linenumber}" 'NR==lines { sub(startt, stopp) } { print }' "${genbank_file_out}" 1> "${genbank_file_out}.test" 2>>"${LOG}"
      mv "${genbank_file_out}.test" "${genbank_file_out}" >> "${LOG}" 2>&1 ;
    #else
    #  echo "old_id:$oldid does not exist";
    fi  
  done < ${contig_id_list_file}

  if [[ ! -e "${genbank_file_out}" ]]; then echo "[ERROR] the output file ${genbank_file_out} does not exist" >> "${LOG}" 2>&1; exit 1; fi
  if [[ ! -s "${genbank_file_out}" ]]; then echo "[ERROR] the output file ${genbank_file_out} is empty" >> "${LOG}" 2>&1; exit 1; fi
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
  LOG="${DATA}/${DIRECTORY}/logs/${DIRECTORY}.annotation.1.log";
  readonly LOG;
  echo '' > "${LOG}";

  GENOME_DIR="${DATA_ECAMBER_OUT}/${DIRECTORY}/genomes"
  readonly GENOME_DIR;


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



  echo 'starting annotation step...' >> "${LOG}" 2>&1 ;
  echo 'starting annotation step...';
  echo 'transfer of genomes to ecamber directory (ext-tools/ecamber/datasets) and format fasta headers' >> "${LOG}" 2>&1
  echo 'transfer of genomes to ecamber directory (ext-tools/ecamber/datasets) and format fasta headers';


  #files transfers of diversity output to ecamber input directory (fasta)
  echo "copy genomes to ecamber input directory (ext-tools/ecamber/datasets/${DIRECTORY}) from the list 'selected_genomes_list.txt'." >> "${LOG}" 2>&1
  echo "copy genomes to ecamber input directory (ext-tools/ecamber/datasets/${DIRECTORY}) from the list 'selected_genomes_list.txt'.";
  mkdir -p "${DATA_ECAMBER_OUT}/${DIRECTORY}" >> "${LOG}" 2>&1
  if [[ -d "${DATA_ECAMBER_OUT}/${DIRECTORY}/genomes" ]]; then rm -r "${DATA_ECAMBER_OUT}/${DIRECTORY}/genomes" >> "${LOG}" 2>&1; fi
  mkdir -p "${DATA_ECAMBER_OUT}/${DIRECTORY}/genomes" >> "${LOG}" 2>&1
  genomelist="$(real_path_file "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/selected_genomes_list.txt")"
  if [[ ! -s "${genomelist}" ]]; then echo "[ERROR] Genome list ${genomelist} is empty, exiting program..."; exit 1; fi; 
  while read line; do
    fasta="$(readlink -f "${line}")"
    cp "${fasta}" "${DATA_ECAMBER_OUT}/${DIRECTORY}/genomes/" >> "${LOG}" 2>&1
  done < ${genomelist}


  #fasta header modification
  #2 cases :
  #-1 sequence -> header = name of the strain
  #-multifasta -> header = contig_name strain_name 

  for f in $(ls "${GENOME_DIR}"); do
    #count number of fasta sequence in the file
    seq_nb="$(grep -c '>' "${GENOME_DIR}/${f}")"
	
    gawk -v file="$(basename "${f}" ".fasta")" -v seq_nb="${seq_nb}" '{

      #if header of fasta file
      if($0 ~ "^>.*"){
        gsub(">","",$1); 
	gsub("\\.","_",$1);
	gsub("\\|","",$1);
	gsub("\\(","",$1);
	gsub(")","",$1);    

	if(seq_nb==1){
          print ">" file;
	}else{
	  print ">" $1 " " file;	
	}
    	}else{	
       	  print $0
	}

    }' "${GENOME_DIR}/${f}" 1> "${GENOME_DIR}/${f}.tmp" 2>>"${LOG}"

    mv "${GENOME_DIR}/${f}.tmp" "${GENOME_DIR}/${f}" >> "${LOG}" 2>&1
  done


  #catches each genome fasta file name to create strains.txt (required by ecamber as input)
  ls -1 "${GENOME_DIR}" | sed 's/.fasta//g' 1> "${DATA_ECAMBER_OUT}/${DIRECTORY}/strains.txt" 2>>"${LOG}"

  if [[ ! -s "${DATA_ECAMBER_OUT}/${DIRECTORY}/strains.txt" ]]; then 
    echo "[ERROR] the file ${DATA_ECAMBER_OUT}/${DIRECTORY}/strains.txt is empty, but it is required for ecamber." 2>>"${LOG}"; exit 1; 
  fi


  #get again modified contig id of the reference genome (ecamber contig id) 
  if [[ "${REFGENOME}" != 'N.O.R.E.F' ]] && [[ "${REFANNOTATION}" != 'N.O.R.E.F' ]] && [[ "${REFIDPATTERN}" != 'N.O.P.A.T.T.E.R.N' ]]; then

    base="$(basename "${REFGENOME}")";
    extension="${base##*.}";
    prefix="${base%.*}";  #file without extension
    refname="${prefix}.fasta";
    newrefname="$(grep "${refname}" "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/renamed_genomes_list.txt" |awk '{ print $2 }')";
    grep '>' "${GENOME_DIR}/${newrefname}" |awk '{print $1}' |sed 's/>//g' 1> "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_filtered_contigs.new.txt" 2>>"${LOG}";
    paste "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_filtered_contigs.txt" "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_filtered_contigs.new.txt" \
        1> "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_filtered_contigs.txt.tmp" 2>>"${LOG}"
    mv "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_filtered_contigs.txt.tmp" "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_filtered_contigs.txt" >> "${LOG}" 2>&1
    if [[ -e "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_filtered_contigs.new.txt" ]]; then rm "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_filtered_contigs.new.txt"; fi

   
    basea="$(basename "${REFANNOTATION}")";
    extensiona="${base##*.}";
    prefixa="${base%.*}";  #file without extension
    refannotname="${basea}"
    newprefixrefgenome="${newrefname%.*}";
    newrefannotname="${newprefixrefgenome}.txt";  #genbank file will have the same prefix than the reference fasta file


    #rename annotation file as .txt
    if [[ -d "${DATA}/${DIRECTORY}/ref_gbk_annotation/parsed" ]]; then rm -r "${DATA}/${DIRECTORY}/ref_gbk_annotation/parsed" >> "${LOG}" 2>&1; fi
    mkdir -p "${DATA}/${DIRECTORY}/ref_gbk_annotation/parsed" >> "${LOG}" 2>&1
    cp -p "${DATA}/${DIRECTORY}/ref_gbk_annotation/${refannotname}" "${DATA}/${DIRECTORY}/ref_gbk_annotation/parsed/${newrefannotname}" >> "${LOG}" 2>&1 #ecamber requires .txt extension for annotation inputs


    #remove sequences from genbank file if they have been filtered out (N50, number_of_contigs, genome_size filters in script init_DIVERSITY.sh)
    if [[ -e "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_contigs.removed.txt" ]] && [[ -s "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_contigs.removed.txt" ]]; then
      echo -e 'INIT_ECAMBER - remove filtered contigs in genbank reference file';
      genbank_file_="${DATA}/${DIRECTORY}/ref_gbk_annotation/parsed/${newrefannotname}";
      genbank_file_out_="${DATA}/${DIRECTORY}/ref_gbk_annotation/parsed/${newrefannotname}.filtered";
      contig_id_list_file_="${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_contigs.removed.txt";
      removeGenbankIdFromList "${contig_id_list_file_}" "${genbank_file_}" "${genbank_file_out_}";
      mv "${DATA}/${DIRECTORY}/ref_gbk_annotation/parsed/${newrefannotname}.filtered" "${DATA}/${DIRECTORY}/ref_gbk_annotation/parsed/${newrefannotname}" >> "${LOG}" 2>&1
    fi


    #rename ACCESSION in reference genbank file
    contig_id_list_file_="${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_filtered_contigs.txt";
    genbank_file_="${DATA}/${DIRECTORY}/ref_gbk_annotation/parsed/${newrefannotname}";
    genbank_file_out_="${DATA}/${DIRECTORY}/ref_gbk_annotation/parsed/${newrefannotname}.modifiedaccessions";
    if [[ -e "${contig_id_list_file_}" ]] && [[ -e "${genbank_file_}" ]]; then 
      replaceGenbankId "${contig_id_list_file_}" "${genbank_file_}" "${genbank_file_out_}";
      mv "${genbank_file_out_}" "${genbank_file_}" >> "${LOG}" 2>&1
    else
      echo -e "Missing files : \n\t ${contig_id_list_file_} \t OR : ${genbank_file_}." >> "${LOG}" 2>&1
    fi

    ## copy genbank annotation to ecamber directory
    if [[ -d "${DATA_ECAMBER_OUT}/${DIRECTORY}/anns_input" ]]; then rm -r "${DATA_ECAMBER_OUT}/${DIRECTORY}/anns_input" >> "${LOG}" 2>&1; fi
    mkdir -p "${DATA_ECAMBER_OUT}/${DIRECTORY}/anns_input" >> "${LOG}" 2>&1
    mkdir -p "${DATA_ECAMBER_OUT}/${DIRECTORY}/anns_input/refseq_gbk" >> "${LOG}" 2>&1
    cp -p "${genbank_file_}" "${DATA_ECAMBER_OUT}/${DIRECTORY}/anns_input/refseq_gbk/." >> "${LOG}" 2>&1

  fi


  if [[ -d "${DATA}/${DIRECTORY}/${DATA_ECAMBER_OUTPUT}" ]]; then rm -r "${DATA}/${DIRECTORY}/${DATA_ECAMBER_OUTPUT}"; fi
  if [[ -d "${DATA_ECAMBER_OUT}/${DIRECTORY}/anns_parsed" ]]; then rm -r "${DATA_ECAMBER_OUT}/${DIRECTORY}/anns_parsed"; fi
  if [[ -d "${DATA_ECAMBER_OUT}/${DIRECTORY}/blast" ]]; then rm -r "${DATA_ECAMBER_OUT}/${DIRECTORY}/blast"; fi
  if [[ -d "${DATA_ECAMBER_OUT}/${DIRECTORY}/cambervis" ]]; then rm -r "${DATA_ECAMBER_OUT}/${DIRECTORY}/cambervis"; fi
  if [[ -d "${DATA_ECAMBER_OUT}/${DIRECTORY}/genomes_dbs" ]]; then rm -r "${DATA_ECAMBER_OUT}/${DIRECTORY}/genomes_dbs"; fi
  if [[ -d "${DATA_ECAMBER_OUT}/${DIRECTORY}/output" ]]; then rm -r "${DATA_ECAMBER_OUT}/${DIRECTORY}/output"; fi
  if [[ -d "${DATA_ECAMBER_OUT}/${DIRECTORY}/prodigal" ]]; then rm -r "${DATA_ECAMBER_OUT}/${DIRECTORY}/prodigal"; fi
  if [[ -d "${DATA_ECAMBER_OUT}/${DIRECTORY}/results" ]]; then rm -r "${DATA_ECAMBER_OUT}/${DIRECTORY}/results"; fi
  if [[ -f "${DATA_ECAMBER_OUT}/${DIRECTORY}/formatdb.log" ]]; then rm "${DATA_ECAMBER_OUT}/${DIRECTORY}/formatdb.log"; fi
}

main "$@"

