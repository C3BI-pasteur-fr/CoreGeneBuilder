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
# get_ortholog_distribution.sh: computes for each homologous gene, the number of strains having this gene.          # 
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
  echo -ne "Usage: $0 <directory_of_homologs> <file_of_genome_list> <output_directory>\n" 1>&2
}



########################################################################################
# MAIN FUNCTION                                                                        #
########################################################################################
# main
# Computes for each homologous gene, the number of strains having this gene.
# Parameters:
#  - 1) directory_ortho : A directory path where are *.synt files (homologs).
#  - 2) list_all_genomes : A file path containing a list of fasta files (genomes).
#  - 3) dir_out : A directory where output files (*.histo and *.histo.pdf). 
# Returns:
#  - a file with extension '.histo' containing for each homologous gene, the number of strains having this gene.  
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

  #input parameters
  ## ref_genome_annot=$1       # limo.001.c01
  directory_ortho="$1"  # /pasteur/projets/gem/programs/PhyloPipeline/ext-tools/CoreGenome/script/TMP
  list_all_genomes="$2"  # /pasteur/projets/gem/programs/PhyloPipeline/ext-tools/CoreGenome/script/Genomes-limo.lst  
  dir_out="$3"


  # $ head Genomes-limo.lst 
  #limo.001.c001
  #limo.002.c001

  # files containing homologs
  #limo.001.c0Ø1.limo.001.cØ01.opsc.synt

  ref_genome_annot=$(head -n 1 "${list_all_genomes}")  #suppose the first line contains the name of the reference annotation
  # ref_genome_annot="${ref_genome_annot%.*}"

  echo -e "Reference Genome Annotation : ${ref_genome_annot}"
  echo -e "Genome list: ${list_all_genomes}"

  ref_genome_annot_synt="${directory_ortho}/${ref_genome_annot}.${ref_genome_annot}.opsc.synt"
  echo -e "Reference Genome Orthologs : ${ref_genome_annot_synt}"
  orthologs_ref=$(awk '{print $1}' "${ref_genome_annot_synt}")
  #echo -e "$orthologs_ref"

  echo -n "" > "${dir_out}/${ref_genome_annot}.histo"
  #echo -e "#ORTHOLOG #NUMBER" > $dir_out/$ref_genome_annot.histo
  echo -e "SEARCH FOR REFERENCE ORTHOLOGS IN OTHER GENOMES:"

  occ_protein=0

  for ortholog in $(echo "${orthologs_ref}"); do
    occ_protein=0;

    for elt in $(cat "${list_all_genomes}"); do 
      currentfile=$(echo "${directory_ortho}/${ref_genome_annot}.${elt}.opsc.synt");
      occurrence=$(grep -c "${ortholog}" "${currentfile}");
      occ_protein=$((${occ_protein} + ${occurrence}));
    done 
  
    echo -e "${ortholog} ${occ_protein}" >> "${dir_out}/${ref_genome_annot}.histo"

  done

  file="${dir_out}/${ref_genome_annot}.histo"
  if [[ ! -e $file ]]; then echo "[ERROR] the output file $file does not exist" >&2; exit 1; fi
  if [[ ! -s $file ]]; then echo "[ERROR] the output file $file is empty" >&2; exit 1; fi


  #plot histogram with R (ortholog frequency in the genomes)
  histodata="${dir_out}/${ref_genome_annot}.histo"   #input
  plotpdf="${dir_out}/${ref_genome_annot}.histo.pdf" #output
  Rcommand="";
  Rcommand="${Rcommand} orthologs<-read.table(file='$histodata', header=FALSE, sep=\"\", col.name=c('orthologId','occurrence'));"; 
  Rcommand="${Rcommand} max=max(orthologs\$occurrence, na.rm=TRUE);";
  # percentage instead of number
  Rcommand="${Rcommand} prct=orthologs\$occurrence/max*100;";
  # orthologs contains now the vector 'prct' as a third column 
  Rcommand="${Rcommand} orthologs = cbind(orthologs,prct);";
  Rcommand="${Rcommand} scale <- 10;";
  Rcommand="${Rcommand} pdf(file='$plotpdf', height=scale, width=scale);"; 
  Rcommand="${Rcommand} h=hist(orthologs\$prct, breaks=100, plot=FALSE);"; 
  Rcommand="${Rcommand} max_counts=max(h\$counts)+100;";
  Rcommand="${Rcommand} plot(h, col='grey', xlab='Percentage of genomes', ylab='Orthologs frequency', main='Frequency of orthologs into genomes', axes=FALSE);";
  Rcommand="${Rcommand} abline(v=95);";  #add x=95
  Rcommand="${Rcommand} axis(side=1, at=seq(0, 100, by=10));";
  Rcommand="${Rcommand} axis(side=2, at=seq(0, max_counts, by=100));";
  Rcommand="${Rcommand} graphics.off();";

  Rexe=R;
  echo "${Rcommand}" |$Rexe --slave;
  exit_status=$? ; if [[ $exit_status -ne 0 ]]; then echo "[ERROR] pdf file '$plotpdf' was not output! - exit code $exit_status" >&2; fi
}

main "$@"

