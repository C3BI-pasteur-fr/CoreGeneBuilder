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
# init_DIVERSITY.sh: initial steps before genomic sequence classification based on estimated genomic evolutionary   # 
#                    distances.                                                                                     #
# Author: Elise Larsonneur                                                                                          #
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



# filterGenomesBasedOnSize
# This function gets assembly metrics of each fasta file of the input directory and add suffix to files that don't pass the defined 
#    thresholds (median genome size, contig number, N50 contig length). 
# Parameters 
# - 1) Input directory - A string containing the directory path where the fasta files are present.
# - 2) Input file - A string containing a reference genome fasta file path.
# - 3) Output file - A string containing a file path that records the assembly metrics of each fasta file of the input directory. 
function filterGenomesBasedOnSize {
  ## get Genome Size of each fasta of the directory genome_dir_ 
  local genome_dir_ outfile_ ref_
  local base contiginfos genomeSize genomeN50 genomeNbCtgs toout median outfile_n50_ baseref baserefwoext baseremovedg
  if [[ -z "$1" ]]; then echo '[ERROR] the function filterGenomesBasedOnSize expects a directory to run' 2>>"${LOG}"; exit 1; fi
  if [[ -z "$3" ]]; then echo '[ERROR] the function filterGenomesBasedOnSize expects a reference fasta file to run' 2>>"${LOG}"; exit 1; fi
  if [[ -z "$2" ]]; then echo '[ERROR] the function filterGenomesBasedOnSize expects a file to run' 2>>"${LOG}"; exit 1; fi
  if [[ $# -eq 0 ]] || [[ $# -gt 3 ]]; then echo "[ERROR] the function filterGenomesBasedOnSize expects 3 parameters, $# parameters given" 2>>"${LOG}"; exit 1; fi
  genome_dir_="$1"
  outfile_="$2"
  ref_="$3"
  if [[ ! -e "${genome_dir_}" ]]; then echo "[ERROR] the directory ${genome_dir_} does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${genome_dir_}" ]]; then echo "[ERROR] the directory ${genome_dir_} is empty" 2>>"${LOG}"; exit 1; fi

  for fasta in $(ls "${genome_dir_}/"*.fasta); do
    base="$(basename "${fasta}" '.fasta')"
    contiginfos="$(${CGB_BIN}contig_info.sh -m 500 "${fasta}")";
    exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] the command '${CGB_BIN}contig_info.sh -m 500 ${fasta}' exited with an error : ${contiginfos} - exit code ${exit_status}" 2>>"${LOG}"; exit 1; fi
    genomeSize="$(echo "${contiginfos}" |grep 'Total' |awk '{print $2}')";
    genomeN50="$(echo "${contiginfos}" |grep 'N50' |awk '{print $2}')";
    genomeNbCtgs="$(echo "${contiginfos}" |grep 'Number' |awk '{print $4}')";

    toout="${fasta}\t${genomeSize}\t${genomeN50}\t${genomeNbCtgs}";
    echo -e "${toout}";   ## fasta_name genome_size for each genome
  done  1>"${outfile_}" 2>>"${LOG}"

  if [[ ! -e "${outfile_}" ]]; then echo "[ERROR] the file ${outfile_} does not exist" 2>>"${LOG}"; exit 1; fi
  if [[ ! -s "${outfile_}" ]]; then echo "[ERROR] the file ${outfile_} is empty" 2>>"${LOG}"; exit 1; fi

  # compute the median genome size
  median="$(awk 'BEGIN {FS="\t"}; {print $2}' "${outfile_}" |sort -n |awk '{ a[i++]=$1; } END { x=int(i/2); if (x < (i/2)) print (a[x-1]+a[x])/2; else print a[x-1]; }')"
  echo "Median Genome Size :${median}" >> "${LOG}" 2>&1
  echo "Median Genome Size :${median}";

  # filter out the genomes the size of which is inferior to 0.5 x median or superior to 1.5 x median
  echo "filter out the genomes the size of which is inferior to 0.5 x median or superior to 1.5 x median" >> "${LOG}" 2>&1
  echo "filter out the genomes the size of which is inferior to 0.5 x median or superior to 1.5 x median";

  awk -v med="${median}" 'BEGIN {FS="\t"}; {if($2>(med*1.5) || $2<(med*0.5)) {print $1"\t"$2}} ' "${outfile_}" 1> "${outfile_}.removedgenomes" 2>>"${LOG}"

  # filter out the genomes the N50 sequence length of which is inferior to N50 threshold OR the genomes the sequence number is greater than the sequence number threshold 
  echo 'filter out the genomes the N50 sequence length of which is inferior to N50 threshold OR the genomes the sequence number is greater than the sequence number threshold' >> "${LOG}" 2>&1
  echo 'filter out the genomes the N50 sequence length of which is inferior to N50 threshold OR the genomes the sequence number is greater than the sequence number threshold';
  outfile_n50_="${outfile_}.n50nbctg.removedgenomes"
  if [[ "${N50THRESHOLD}" -ge 1 ]] && [[ "${NBCTGTHRESHOLD}" -gt 0 ]]; then 
       awk -v n50length="${N50THRESHOLD}" -v nbctgth="${NBCTGTHRESHOLD}" 'BEGIN {FS="\t"}; {if($3<(n50length) || $4>(nbctgth)) {print $1"\t"$2"\t"$3"\t"$4}} ' "${outfile_}" 1> "${outfile_n50_}" 2>>"${LOG}"
  fi   #N50 + nb_seq filters
  if [[ "${N50THRESHOLD}" -ge 1 ]] && [[ "${NBCTGTHRESHOLD}" -lt 0 ]]; then 
       awk -v n50length="${N50THRESHOLD}" 'BEGIN {FS="\t"}; {if($3<(n50length)) {print $1"\t"$2"\t"$3"\t"$4}} ' "${outfile_}" 1> "${outfile_n50_}" 2>>"${LOG}"
  fi   #N50 filter only
  if [[ "${N50THRESHOLD}" -lt 1 ]] && [[ "${NBCTGTHRESHOLD}" -gt 0 ]]; then 
       awk -v nbctgth="${NBCTGTHRESHOLD}" 'BEGIN {FS="\t"}; {if($4>(nbctgth)) {print $1"\t"$2"\t"$3"\t"$4}} ' "${outfile_}" 1> "${outfile_n50_}" 2>>"${LOG}"
  fi   #nb_seq filter only

  # if some problematic genomes, print them
  if [[ -s "${outfile_}.removedgenomes" ]]; then   #if is not zero size 
    echo 'genomes with a problematic genome size (inferior to 0.5 * median genome size OR superior to 1.5 * median genome size) :' >> "${LOG}" 2>&1
    echo 'genomes with a problematic genome size (inferior to 0.5 * median genome size OR superior to 1.5 * median genome size) :'
    # rename them 
    for fasta in $(awk 'BEGIN {FS="\t"}; {print $1}' "${outfile_}.removedgenomes"); do
      if [[ -e "${fasta}" ]]; then
         echo "     ${fasta}" >> "${LOG}" 2>&1
         echo "     ${fasta}";
         mv "${fasta}" "${fasta}.sizeproblem" >> "${LOG}" 2>&1 ;
      fi
    done
  else  
    if [[ -e "${outfile_}.removedgenomes" ]]; then rm "${outfile_}.removedgenomes"; fi
  fi

  if [[ -s "${outfile_n50_}" ]]; then   #if is not zero size 
    if [[ "${N50THRESHOLD}" -ge 1 ]] && [[ "${NBCTGTHRESHOLD}" -gt 0 ]]; then 
        echo "genomes with a problematic N50 and/or contigs_number value(s), eg ( N50 < ${N50THRESHOLD} bp ) AND/OR ( nb_contigs > ${NBCTGTHRESHOLD} ) :" >> "${LOG}" 2>&1 ; 
        echo "genomes with a problematic N50 and/or contigs_number value(s), eg ( N50 < ${N50THRESHOLD} bp ) AND/OR ( nb_contigs > ${NBCTGTHRESHOLD} ) :"; 
    fi
    if [[ "${N50THRESHOLD}" -ge 1 ]] && [[ "${NBCTGTHRESHOLD}" -lt 0 ]]; then 
        echo "genomes with a problematic N50 and/or contigs_number value(s), eg ( N50 < ${N50THRESHOLD} bp ) :" >> "${LOG}" 2>&1 ;  
        echo "genomes with a problematic N50 and/or contigs_number value(s), eg ( N50 < ${N50THRESHOLD} bp ) :";  
    fi
    if [[ "${N50THRESHOLD}" -lt 1 ]] && [[ "${NBCTGTHRESHOLD}" -gt 0 ]]; then 
        echo "genomes with a problematic N50 and/or contigs_number value(s), eg ( nb_contigs > ${NBCTGTHRESHOLD} ) :" >> "${LOG}" 2>&1 ;
        echo "genomes with a problematic N50 and/or contigs_number value(s), eg ( nb_contigs > ${NBCTGTHRESHOLD} ) :";
    fi

    # rename them 
    for fasta in $(awk 'BEGIN {FS="\t"}; {print $1}' "${outfile_n50_}"); do
      if [[ -e "${fasta}" ]]; then
         echo "     ${fasta}" >> "${LOG}" 2>&1
         echo "     ${fasta}";
         mv "${fasta}" "${fasta}.n50_nbctg_problem" >> "${LOG}" 2>&1 ; 
      fi
    done
  else
    if [[ -e "${outfile_n50_}" ]]; then rm "${outfile_n50_}"; fi
  fi

  # check that reference is not in the problematic genomes
  if [[ "${ref_}" != 'N.O.R.E.F' ]]; then 
       baseref="$(basename "${ref_}")";
       baserefwoext="${baseref%.*}";  
       if [[ -s "${outfile_}.removedgenomes" ]]; then
         for fasta in $(awk 'BEGIN {FS="\t"}; {print $1}' "${outfile_}.removedgenomes"); do
             baseremovedg="$(basename "${fasta}" '.fasta')"; 
             if [[ "${baserefwoext}" = "${baseremovedg}" ]]; then 
                  echo '[ERROR] reference genome has a problematic genome size! Exiting program...' >> "${LOG}" 2>&1 ; 
                  echo '[ERROR] reference genome has a problematic genome size! Exiting program...';
                  exit 1;
             fi
         done
       fi

       if [[ -s "${outfile_n50_}" ]]; then 
         for fasta in $(awk 'BEGIN {FS="\t"}; {print $1}' "${outfile_n50_}"); do
             baseremovedg="$(basename "${fasta}" '.fasta')";
             if [[ "${baserefwoext}" = "${baseremovedg}" ]]; then 
                  echo '[ERROR] reference genome has problematic N50 and/or contigs_number value(s)! Exiting program...' >> "${LOG}" 2>&1 ;
                  echo '[ERROR] reference genome has problematic N50 and/or contigs_number value(s)! Exiting program...';
                  exit 1;
             fi
         done

       fi
  fi
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
  LOG="${DATA}/${DIRECTORY}/logs/${DIRECTORY}.diversity.1.log";
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

  ## checking if reference genome fasta file and reference genbank file have the same contig id
  if [[ "${REFGENOME}" != "N.O.R.E.F" ]] && [[ "${REFANNOTATION}" != "N.O.R.E.F" ]] && [[ "${REFIDPATTERN}" != "N.O.P.A.T.T.E.R.N" ]]; then
    grep ">" "${DATA}/${DIRECTORY}/assemblies/${REFGENOME}" |grep -qP "${REFIDPATTERN}\d+" -o
    retvalfasta=$?
    grep "ACCESSION" "${DATA}/${DIRECTORY}/ref_gbk_annotation/${REFANNOTATION}" |grep -qP "${REFIDPATTERN}\d+" -o
    retvalgb=$?
    if [[ $retvalfasta -eq 1 ]] || [[ $retvalgb -eq 1 ]]; then echo "[ERROR] Files ${DATA}/${DIRECTORY}/assemblies/${REFGENOME} and ${DATA}/${DIRECTORY}/ref_gbk_annotation/${REFANNOTATION} don't have the same contig_id, their contig id must start both with provided prefix : '${REFIDPATTERN}'."; exit 1; fi
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



  echo 'starting diversity step...' >> "${LOG}" 2>&1
  echo 'starting diversity step...';
  echo 'prepare files and directories before starting diversity step' >> "${LOG}" 2>&1
  echo 'prepare files and directories before starting diversity step';
  #files transfers from assembly output to diversity input (fasta)
  mkdir -p "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}" >> "${LOG}" 2>&1
  mkdir -p "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/genomes" >> "${LOG}" 2>&1   ## directory where fasta files will be pre-processed
  touch "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/renamed_genomes_list.txt" >> "${LOG}" 2>&1
  if [[ "${REFANNOTATION}" != 'N.O.R.E.F' ]]; then touch "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_filtered_contigs.txt" >> "${LOG}" 2>&1 ; fi  ## we keep the original contig id of the reference genome
  cp -p "${DATA}/${DIRECTORY}/${DATA_ASSEMBLY_DIR}/"* "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/genomes/" >> "${LOG}" 2>&1

  GENOME_DIR="$(real_path_dir "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/genomes")"
  readonly GENOME_DIR

  #remove space in filenames
  ls ${GENOME_DIR}/*\ * 2>/dev/null && for f in ${GENOME_DIR}/*\ *; do echo "mv ${f} ${f// /}"; mv "${f}" "${f// /}"; done || echo 'no space characters in fasta files'


  #rename fasta files with .fasta extension
  echo 'rename fasta files with .fasta extension' ;  
  echo 'rename fasta files with .fasta extension' >> "${LOG}" 2>&1 ;  
  #rename fasta files with .fasta extension 
  for file in $(ls ${GENOME_DIR} | sort); do
    if ! [[ -h "${file}" ]]; then      ### !!! ignoring existing symbolic links  !!!
      base="$(basename "${file}")"
      prefix="${base%.*}";
      extension="${base##*.}";
      #check extension (.fa ; .fasta ; .fna ; .fas)
      if [[ "${extension}" = "fa" ]] || [[ "${extension}" = "fna" ]] || [[ "${extension}" = "fas" ]]; then 
        mv "${GENOME_DIR}/${file}" "${GENOME_DIR}/${prefix}.fasta" >> "${LOG}" 2>&1;  ## if extension = fa or fna or fas , replace it by fasta  
      fi
    fi
  done


  ## 1a-a- cut scaffolds into contigs (if strechtes of 'N' in sequence) - except for reference genome
  if [[ "${CUTN}" != "-1" ]] && [[ "${CUTN}" -gt "0" ]]; then
   
    #cut N-mers
    echo "cut scaffolds into contigs (if 10-mers of 'N' found), except for the provided reference genome" ;
    echo "cut scaffolds into contigs (if 10-mers of 'N' found), except for the provided reference genome" >> "${LOG}" 2>&1 ;
    list_N_cut='';
    
    if [[ "${REFGENOME}" != 'N.O.R.E.F' ]]; then 
       baseref="$(basename "${REFGENOME}")";
       baserefwoext="${baseref%.*}"; 
       list_N_cut="$(ls "${GENOME_DIR}/"*.fasta |grep -v "${baserefwoext}")";
    else
       list_N_cut="$(ls "${GENOME_DIR}/"*.fasta)";
    fi

    for fasta in ${list_N_cut}; do
      if [[ -s "${fasta}" ]]; then 
        echo "${CGB_BIN}fasta2agp -n ${CUTN} -i ${fasta} -o ${fasta}.tmp;" >> "${LOG}" 2>&1;
        ${CGB_BIN}fasta2agp -n "${CUTN}" -i "${fasta}" -o "${fasta}.tmp" >> "${LOG}" 2>&1;
        exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] the command '${CGB_BIN}fasta2agp -n ${CUTN} -i ${fasta} -o ${fasta}.tmp' exited with an error! - exit code ${exit_status}" 2>>"${LOG}"; exit 1; fi
        prefix="$(basename "${fasta}" '.fasta')";

        if [[ -e "${fasta}.tmp" ]]; then 
          nb_scaff="$(grep -c '>' "${fasta}")";
          nb_ctgs="$(grep -c '>' "${fasta}.tmp")";
          if [[ "${nb_scaff}" != "${nb_ctgs}" ]]; then echo -e "n_cut:\t${prefix}.fasta"; fi
          mv "${fasta}.tmp" "${fasta}"; 
        fi
      fi
    done

    rm "${GENOME_DIR}/"*.fasta.agp 
  fi


  ## 1a-b- keep contigs above the threshold size
  #Usage: fasta_keep_seq_above_threshold_size.py input.fasta threshold_length output.fasta

  echo "Keeping sequences of length above ${CONTIGLENGTH} bp." >> "${LOG}" 2>&1 ; 
  echo "Keeping sequences of length above ${CONTIGLENGTH} bp." ; 
  echo "running script ${CGB_BIN}fasta_keep_seq_above_threshold_size.py" >> "${LOG}" 2>&1
  
  for file in $(ls "${GENOME_DIR}" | sort); do
    if ! [[ -h "${file}" ]]; then      ### !!! ignoring existing symbolic links  !!!
      base="$(basename "${file}")"
      prefix="${base%.*}";
      extension="${base##*.}";

      #check extension (.fa ; .fasta ; .fna ; .fas)
      if [[ "${extension}" = 'fa' ]] || [[ "${extension}" = 'fasta' ]] || [[ "${extension}" = 'fna' ]] || [[ "${extension}" = 'fas' ]]; then 
        ${CGB_BIN}fasta_keep_seq_above_threshold_size.py "${GENOME_DIR}/${file}" "${CONTIGLENGTH}" "${GENOME_DIR}/${prefix}.tmp.fasta" >> "${LOG}" 2>&1 ;
        exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] the command '${CGB_BIN}fasta_keep_seq_above_threshold_size.py ${GENOME_DIR}/${file} ${CONTIGLENGTH} ${GENOME_DIR}/${prefix}.tmp.fasta' exited with an error! - exit code ${exit_status}" 2>>"${LOG}"; exit 1; fi
        if [[ "${GENOME_DIR}/${file}" != "${GENOME_DIR}/${prefix}.fasta" ]]; then mv "${GENOME_DIR}/${file}" "${GENOME_DIR}/${prefix}.fasta" >> "${LOG}" 2>&1; fi  ## if extension = fa or fna or fas , replace it by fasta - added the 2016-01-06  
        mv "${GENOME_DIR}/${prefix}.tmp.fasta" "${GENOME_DIR}/${prefix}.fasta" >> "${LOG}" 2>&1 ;
        seq_nb="$(grep -c '>' "${GENOME_DIR}/${prefix}.fasta")" ## empty fasta file ?
        if [[ $seq_nb -eq 0 ]]; then mv "${GENOME_DIR}/${prefix}.fasta" "${GENOME_DIR}/${prefix}.fasta.empty" >> "${LOG}" 2>&1 ; fi; 
      fi
    fi
  done



  ## 1b- check that the genome size of each fasta is around the median genome size  
  ###### (if not, fasta is discarded, because ecamber will not annotate them)

  echo 'check that the genome size of each fasta is around the median genome size' >> "${LOG}" 2>&1
  echo 'check that the genome size of each fasta is around the median genome size';


  filterGenomesBasedOnSize "${GENOME_DIR}" "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/genomes_and_genomesize_list" "${REFGENOME}" 
  echo "filterGenomesBasedOnSize ${GENOME_DIR} ${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/genomes_and_genomesize_list ${REFGENOME}" >> "${LOG}" 2>&1 ;



  #number of fasta files :
  nb_fasta="$(find "${GENOME_DIR}/." -mindepth 1 ! -type l -name '*.fasta' | wc -l)"
  nb_empty_fasta="$(find "${GENOME_DIR}/." -mindepth 1 ! -type l -name '*.fasta.empty' | wc -l)"
  nb_sizeproblem_fasta="$(find "${GENOME_DIR}/." -mindepth 1 ! -type l -name '*.fasta.sizeproblem' | wc -l)"
  nb_n50_nbctg_problem_fasta="$(find "${GENOME_DIR}/." -mindepth 1 ! -type l -name '*.fasta.n50_nbctg_problem' | wc -l)"
  echo "number of not-empty fasta files : ${nb_fasta}" >> "${LOG}" 2>&1 ;
  echo "number of not-empty fasta files : ${nb_fasta}";
  echo "number of empty fasta files : ${nb_empty_fasta}" >> "${LOG}" 2>&1 ;
  echo "number of empty fasta files : ${nb_empty_fasta}";
  echo "number of fasta files with genome size problem : ${nb_sizeproblem_fasta}" >> "${LOG}" 2>&1 ;
  echo "number of fasta files with genome size problem : ${nb_sizeproblem_fasta}";
  echo "number of fasta files with N5O length and/or contigs number problem : ${nb_n50_nbctg_problem_fasta}" >> "${LOG}" 2>&1 ;
  echo "number of fasta files with N5O length and/or contigs number problem : ${nb_n50_nbctg_problem_fasta}";

  if [[ $nb_fasta -lt 2 ]]; then 
    echo '[ERROR] The number of remaining genomes is lower than 2 ! Exiting program...' >> "${LOG}" 2>&1 ;
    echo '[ERROR] The number of remaining genomes is lower than 2 ! Exiting program...';
    exit 1;
  fi




  ## 2- rename sequence fasta files
  echo 'rename sequence fasta files' >> "${LOG}" 2>&1 ;
  echo 'rename sequence fasta files';

  if [[ "${REFGENOME}" = 'N.O.R.E.F' ]]; then
    #rename_genome_fasta_files.sh <directory of fasta files> <2 first letters of Genus;2 first letters of species> <ref_fasta_name> <outlist>
    #rename_genome_fasta_files.sh data/limo/assemblies/ limo NOREF renamed_list.txt    ##### listeria (li) monocytogenes (mo)
    # !!! WARNING : rename_genome_fasta_files.sh script ignores symbolic links !!!!
    echo "${CGB_BIN}rename_genome_fasta_files.sh ${GENOME_DIR} ${NAME} NOREF ${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/renamed_genomes_list.txt" >> "${LOG}" 2>&1 ;
    ${CGB_BIN}rename_genome_fasta_files.sh "${GENOME_DIR}" "${NAME}" "NOREF" "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/renamed_genomes_list.txt" >> "${LOG}" 2>&1 ;
    exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] the command '${CGB_BIN}rename_genome_fasta_files.sh ${GENOME_DIR} ${NAME} NOREF ${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/renamed_genomes_list.txt' exited with an error! - exit code ${exit_status}" 2>>"${LOG}"; exit 1; fi
  fi

  if [[ "${REFGENOME}" != 'N.O.R.E.F' ]]; then 
    # ref genome name without extension :
    base="$(basename "${REFGENOME}")";
    extension="${base##*.}";
    prefix="${base%.*}";  #file without extension
    refname="${prefix}.fasta";

    if [[ "${REFANNOTATION}" != 'N.O.R.E.F' ]] && [[ "${REFIDPATTERN}" != 'N.O.P.A.T.T.E.R.N' ]]; then
      echo 'remove filtered out contigs from the genbank reference annotation file' >> "${LOG}" 2>&1 ;
      echo 'remove filtered out contigs from the genbank reference annotation file';
      grep '>' "${GENOME_DIR}/${refname}" | grep -P "${REFIDPATTERN}\d+" -o 1> "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_filtered_contigs.txt" 2>>"${LOG}" ; #get refseq accession id of each contig from filtered genome
      grep '>' "${DATA}/${DIRECTORY}/assemblies/${REFGENOME}" | grep -P "${REFIDPATTERN}\d+" -o 1> "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_contigs.txt" 2>>"${LOG}" ; #get refseq accession id of each contig from not filtered genome
      #get filtered out contigs
      awk '{print $1}' "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_filtered_contigs.txt" 1> "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_filtered_contigs.txt.tmp" 2>>"${LOG}"
      awk '{print $1}' "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_contigs.txt" 1> "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_contigs.txt.tmp" 2>>"${LOG}"
      grep -v -F -x -f "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_filtered_contigs.txt.tmp" "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_contigs.txt.tmp" 1> "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_contigs.removed.txt" 2>>"${LOG}"
      if [[ -e "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_filtered_contigs.txt.tmp" ]]; then rm "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_filtered_contigs.txt.tmp" >> "${LOG}" 2>&1 ; fi
      if [[ -e "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_contigs.txt.tmp" ]]; then rm "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_contigs.txt.tmp" >> "${LOG}" 2>&1 ; fi
      #if empty remove the file (eg no contig has been removed) :
      if [[ ! -s "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_contigs.removed.txt" ]]; then rm "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_contigs.removed.txt" >> "${LOG}" 2>&1 ; fi
    fi

    echo "${CGB_BIN}rename_genome_fasta_files.sh ${GENOME_DIR} ${NAME} ${refname} ${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/renamed_genomes_list.txt" >> "${LOG}" 2>&1 ;
    echo "${CGB_BIN}rename_genome_fasta_files.sh ${GENOME_DIR} ${NAME} ${refname} ${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/renamed_genomes_list.txt";
    ${CGB_BIN}rename_genome_fasta_files.sh "${GENOME_DIR}" "${NAME}" "${refname}" "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/renamed_genomes_list.txt" >> "${LOG}" 2>&1 ;   
    exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] the command '${CGB_BIN}rename_genome_fasta_files.sh ${GENOME_DIR} ${NAME} ${refname} ${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/renamed_genomes_list.txt' exited with an error! - exit code ${exit_status}" 2>>"${LOG}"; exit 1; fi 
  
  fi

  echo "list containing old and new name of each isolate : ${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/renamed_genomes_list.txt" >> "${LOG}" 2>&1 ;
  echo "list containing old and new name of each isolate : ${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/renamed_genomes_list.txt";



  ## 3-preprocess fasta before running ecamber (preprocessing : renaming of fasta sequence ids, uppercase )
  echo 'preprocess fasta before running ecamber (preprocessing : renaming of fasta sequence ids, uppercase)' >> "${LOG}" 2>&1 ;
  echo 'preprocess fasta before running ecamber (preprocessing : renaming of fasta sequence ids, uppercase)';
  echo "running script '${CGB_BIN}preprocess_fasta_seq_id.py' for each isolate" >> "${LOG}" 2>&1 ;

  for file in $(ls "${GENOME_DIR}/${NAME}"*.fasta); do 
    if ! [[ -h "${file}" ]]; then ${CGB_BIN}preprocess_fasta_seq_id.py -i "${file}" -o "${file}.renamed" >> "${LOG}" 2>&1 ; # ignoring existing symbolic links
    exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] the command '${CGB_BIN}preprocess_fasta_seq_id.py -i ${file} -o ${file}.renamed' exited with an error! - exit code ${exit_status}" 2>>"${LOG}"; exit 1; fi
    mv "${file}.renamed" "${file}" >> "${LOG}" 2>&1 ; fi; 
  done

  # get modified contig id of the reference genome
  if [[ "${REFANNOTATION}" != 'N.O.R.E.F' ]] && [[ "${REFIDPATTERN}" != 'N.O.P.A.T.T.E.R.N' ]]; then
    base="$(basename "${REFGENOME}")";
    prefix="${base%.*}";  #file without extension
    refname="${prefix}.fasta";
    newrefname="$(get_new_name_from_list "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/renamed_genomes_list.txt" "${refname}")"
    grep '>' "${GENOME_DIR}/${newrefname}" |awk '{print $1}' |sed 's/>//g' 1> "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_filtered_contigs.new.txt" 2>>"${LOG}" ;
    paste "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_filtered_contigs.txt" "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_filtered_contigs.new.txt" 1> "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_filtered_contigs.id.txt" 2>>"${LOG}"
    mv "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_filtered_contigs.id.txt" "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_filtered_contigs.txt" >> "${LOG}" 2>&1 ;
    if [[ -e "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_filtered_contigs.new.txt" ]]; then rm "${DATA}/${DIRECTORY}/${DATA_DIVERSITY_DIR}/ref_filtered_contigs.new.txt" >> "${LOG}" 2>&1 ; fi
  fi


  ## 4 - create symbolic links for andi input files (andi : compute distances between genomes)
  echo "create symbolic links as andi input files" >> "${LOG}" 2>&1 ;
  echo "create symbolic links as andi input files";
  for f in $(find "${GENOME_DIR}/." -mindepth 1 ! -type l -name '*.fasta'); do   # ignoring existing symbolic links 
    fasta="$(basename "${f}")";
    suffix="${fasta#$NAME.}";
    rename="${NAME:0:1}${NAME:2:1}"
    ln -s "${GENOME_DIR}/${fasta}" "${GENOME_DIR}/${rename}${suffix}" >> "${LOG}" 2>&1
  done
}

main "$@"

