#!/bin/bash
#$ -S /bin/bash

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
# contig_info.sh: a BASH script to estimate standard statistics from FASTA contig file.                             #
# Author:  Alexis Criscuolo                                                                                         #
#####################################################################################################################

  VERSION=0.2.2

########################################################################################
#                                                                                      #
# ================                                                                     #
# = INSTALLATION =                                                                     #
# ================                                                                     #
#                                                                                      #
# Just give  the execute  permission on  the script  contig_info.sh by  using the      #
# following command:                                                                   #
#                                                                                      #
#   chmod +x contig_info.sh                                                            #
#                                                                                      #
########################################################################################

########################################################################################
#                                                                                      #
# =============                                                                        #
# = EXECUTION =                                                                        #
# =============                                                                        #
#                                                                                      #
# You can launch contig_info with the following command:                               #
#                                                                                      #
#   ./contig_info.sh [options] <contig_file>                                           #
#                                                                                      #
# where <contig_file> if the name of a FASTA formatted contig sequence file.           #
#                                                                                      #
# Available options:                                                                   #
#                                                                                      #
# -m <int>              minimum  contig  length;  all contigs  shorter than  this      #
#                       value will be discarded (default: 0)                           #
#                                                                                      #
# -g <int>              expected genome size;  this value will  be only used when      #
#                       computing the N50, N75 and N90  statistics  (default: the      #
#                       total sum of all considered contig sequence lengths)           #
#                                                                                      #
# -d                    print contig sequence length distribution                      #
#                                                                                      #
# -l                    print length of each contig sequence                           #
#                                                                                      #
# -r                    print residue counts                                           #
#                                                                                      #
# -t                    tab-delimited output                                           #
#                                                                                      #
########################################################################################


########################################################################################
# Help displaying                                                                      #
########################################################################################

if [ "$1" = "-?" ] || [ $# -lt 1 ]
then
  echo "" ;
  echo " contig_info v.$VERSION" ;
  echo "" ;
  echo " USAGE :" ;
  echo "" ;
  echo "    contig_info.sh [options] <contig_file> " ;
  echo "" ;
  echo "  where 'options' are :" ;
  echo "" ;
  echo "   -m <int>    minimum  contig  length;  all contigs  shorter than  this      "
  echo "               value will be discarded (default: 0)                           "
  echo "" ;
  echo "   -g <int>    expected genome size;  this value will  be only used when      "
  echo "               computing the N50, N75 and N90  statistics  (default: the      "
  echo "               total sum of all considered contig sequence lengths)           "
  echo "" ;
  echo "   -d          print contig sequence length distribution                      "
  echo "" ;
  echo "   -l          print length of each contig sequence                           "
  echo "" ;
  echo "   -r          print residue counts                                           "
  echo "" ;
  echo "   -t          tab-delimited output                                           "
  echo "" ;
  exit
fi
 

########################################################################################
# Functions                                                                            #
########################################################################################

# randomfile  $INFILE
#
#   INFILE: file name, this will returns the randomly generated file name INFILE.RANDOM
#
randomfile() {
  rdmf=$1.$RANDOM; while [ -e $rdmf ]; do rdmf=$1.$RANDOM ; done
  echo $rdmf ;
}



########################################################################################
# Beginning contig_info                                                                #
########################################################################################

###############################
#  Reading options            #
###############################
INFILE="XXX"; for param in "$@"; do INFILE="$param"; done
if [ ! -e $INFILE ]; then echo "   no input file" ; exit 1; fi
if [ ! -s $INFILE ]; then echo "   empty input file" ; exit 1; fi

MIN_CONTIG_LGT=0
GENOME_SIZE=0
PRINT_LGT=false
PRINT_RESIDUE=false
PRINT_DIST=false
CSVOUT=false
while getopts m:g:ldrt option
do
  case $option in
  m)
    MIN_CONTIG_LGT=$OPTARG
    if [ $MIN_CONTIG_LGT -lt 0 ]; then echo "   the min contig length threshold must be a positive integer (option -m)" ; exit 1; fi
   ;;
  g)
    GENOME_SIZE=$OPTARG
    if [ $GENOME_SIZE -lt 0 ]; then echo "   the expected genome size must be a positive integer (option -g)" ; exit 1; fi
   ;;
  r)
    PRINT_RESIDUE=true
   ;;
  l)
    PRINT_LGT=true
   ;;
  d)
    PRINT_DIST=true
   ;;
  t)
    CSVOUT=true
   ;;
  esac
done

###############################
#  Conversion (for dos files) #
###############################
INFILE_UNIX=$(randomfile $INFILE)
tr -d '\15\32' < $INFILE > $INFILE_UNIX ;

###############################
#  Computing info             #
###############################
FH=$(randomfile $INFILE);          ##### FH : txt file containing fhs
SEQS=$(randomfile $INFILE);        ##### SEQS : txt file with one sequence per line
SUBSEQS=$(randomfile $INFILE);     ##### SUBSEQS : txt file with sequence per line larger than $MIN_CONTIG_LGT
LGTFILE=$(randomfile $INFILE);     ##### LGTFILE : length of each sequence
INFOFILE=$(randomfile $INFILE);    ##### INFOFILE : output if $PRINT_LGT="true"

grep "^>" $INFILE_UNIX | cut -c2- > $FH ;
awk 'BEGIN {seq=""} !/^>/ {seq=seq$0} /^>/ {if(seq!=""){print seq ; seq=""} {print $0}} END {print seq}' $INFILE_UNIX | grep -v "^>" | tr '[:lower:]' '[:upper:]' > $SEQS ;

(echo -e -n "\r     \r0.0%" >&2);  awk -v threshold=$MIN_CONTIG_LGT '( length($0) >= threshold ) {print length($0)}' $SEQS > $LGTFILE ;
(echo -e -n "\r     \r20.0%" >&2); NUMBER_OF_SEQUENCE=$(cat $LGTFILE | wc -l);
(echo -e -n "\r     \r40.0%" >&2); R=$(awk '{ sum += $1 } END { print sum }' $LGTFILE);
(echo -e -n "\r     \r60.0%" >&2); if $PRINT_RESIDUE ; then awk -v threshold=$MIN_CONTIG_LGT '( length($0) >= threshold ) {print $0}' $SEQS > $SUBSEQS ; fi
(echo -e -n "\r     \r80.0%" >&2); if $PRINT_LGT ; then lgt=$(cat $SEQS | wc -l); i=0; while read seq; do i=$((i+1)); (echo -e -n "\r     \r$(( 80 + 19 * $i / $lgt )).0%" >&2); r=${#seq}; if [ $r -ge $MIN_CONTIG_LGT ]; then echo "  $(sed -n "$i p" $FH)        $r" >> $INFOFILE ; fi ; done < $SEQS ; fi

if $PRINT_RESIDUE
then
  (echo -e -n "\r     \r99.0%" >&2); A=$(cat $SUBSEQS | tr -d '[:cntrl:]' | awk '{ x=0; x+=gsub("A",""); print x }')
  (echo -e -n "\r     \r99.2%" >&2); C=$(cat $SUBSEQS | tr -d '[:cntrl:]' | awk '{ x=0; x+=gsub("C",""); print x }')
  (echo -e -n "\r     \r99.4%" >&2); G=$(cat $SUBSEQS | tr -d '[:cntrl:]' | awk '{ x=0; x+=gsub("G",""); print x }')
  (echo -e -n "\r     \r99.6%" >&2); T=$(cat $SUBSEQS | tr -d '[:cntrl:]' | awk '{ x=0; x+=gsub("T",""); print x }')
  (echo -e -n "\r     \r99.8%" >&2); N=$(cat $SUBSEQS | tr -d '[:cntrl:]' | awk '{ x=0; x+=gsub("N",""); print x }')
  rm $SUBSEQS ;
fi
(echo -e -n "\r     \r" >&2);

if ! $CSVOUT
then
  echo ;
  echo "File                           $(basename $INFILE)" ;
  echo ;
  echo "Number of sequences            $NUMBER_OF_SEQUENCE" ;
  echo ;
  echo "Residue counts:" ;
  if $PRINT_RESIDUE
  then
    echo "  Number of A's                $A  $(echo "scale=2;100*$A/$R" | bc -l) %" ;
    echo "  Number of C's                $C  $(echo "scale=2;100*$C/$R" | bc -l) %" ;
    echo "  Number of G's                $G  $(echo "scale=2;100*$G/$R" | bc -l) %" ;
    echo "  Number of T's                $T  $(echo "scale=2;100*$T/$R" | bc -l) %" ;
    echo "  Number of N's                $N  $(echo "scale=2;100*$N/$R" | bc -l) %" ;
  fi
  echo "  Total                        $R" ;
  if [ $GENOME_SIZE -ne 0 ]; then echo "  Expected genome size         $GENOME_SIZE"; fi
else
  CSVCAPT="#File\tNseq"; CSVLINE="$(basename $INFILE)\t$NUMBER_OF_SEQUENCE";
  if $PRINT_RESIDUE
  then
    CSVCAPT="$CSVCAPT\tA\tC\tG\tT\tN\t%A\t%C\t%G\t%T\t%N";
    CSVLINE="$CSVLINE\t$A\t$C\t$G\t$T\t$N";
    CSVLINE="$CSVLINE\t$(echo "scale=2;100*$A/$R" | bc -l)%";
    CSVLINE="$CSVLINE\t$(echo "scale=2;100*$C/$R" | bc -l)%";
    CSVLINE="$CSVLINE\t$(echo "scale=2;100*$G/$R" | bc -l)%";
    CSVLINE="$CSVLINE\t$(echo "scale=2;100*$T/$R" | bc -l)%";
    CSVLINE="$CSVLINE\t$(echo "scale=2;100*$N/$R" | bc -l)%";
  fi
  CSVCAPT="$CSVCAPT\tNres"; CSVLINE="$CSVLINE\t$R";
  if [ $GENOME_SIZE -ne 0 ]; then CSVCAPT="$CSVCAPT\tExpSize"; CSVLINE="$CSVLINE\t$GENOME_SIZE"; fi
fi

tmp=$(randomfile $INFILE); cat $LGTFILE | sort -n > $tmp; mv $tmp $LGTFILE ;
MIN=$(head -1 $LGTFILE);
MAX=$(tail -1 $LGTFILE);
Q25=$(head -$(( $NUMBER_OF_SEQUENCE / 4 )) $LGTFILE | tail -1);
Q50=$(head -$(( $NUMBER_OF_SEQUENCE / 2 )) $LGTFILE | tail -1);
Q75=$(head -$(( 3 * $NUMBER_OF_SEQUENCE / 4 )) $LGTFILE | tail -1);
AVG=$(echo "scale=2;$R/$NUMBER_OF_SEQUENCE" | bc -l);

tmp=$(randomfile $INFILE); cat $LGTFILE | sort -r -n > $tmp; mv $tmp $LGTFILE ;
N50=0
N75=0
N90=0
if [ $GENOME_SIZE -eq 0 ]; then GENOME_SIZE=$R; fi
size=0
while read lgt
do
  size=$(( $size + $lgt ))
  if [ $N50 -eq 0 ] && [ $size -ge $(( $GENOME_SIZE / 2 )) ];      then N50=$lgt; fi
  if [ $N75 -eq 0 ] && [ $size -ge $(( 3 * $GENOME_SIZE / 4 )) ];  then N75=$lgt; fi
  if [ $N90 -eq 0 ] && [ $size -ge $(( 9 * $GENOME_SIZE / 10 )) ]; then N90=$lgt; break; fi
done < $LGTFILE

if ! $CSVOUT
then
  echo ;
  echo "Sequence lengths:" ;
  echo "  Minimum                      $MIN" ;
  echo "  Quartile 25%                 $Q25" ;
  echo "  Median                       $Q50" ;
  echo "  Quartile 75%                 $Q75" ;
  echo "  Maximum                      $MAX" ;
  echo "  Average                      $AVG" ;
  echo "  N50                          $N50" ;
  echo "  N75                          $N75" ;
  echo "  N90                          $N90" ;
  echo ;
else
  CSVCAPT="$CSVCAPT\tMin\tQ25\tMed\tQ75\tMax\tAvg\tN50\tN75\tN90";
  CSVLINE="$CSVLINE\t$MIN\t$Q25\t$Q50\t$Q75\t$MAX\t$AVG\t$N50\t$N75\t$N90";
  echo -e "$CSVCAPT" ; 
  echo -e "$CSVLINE" ;
fi



if $PRINT_LGT ; then cat $INFOFILE; echo; fi

if $PRINT_DIST
then 
  tmp=$(randomfile $INFILE); cat $LGTFILE | sort -n > tmp; mv tmp $LGTFILE ;
  echo "Sequence length distribution:"
  min=1
  max=100
  while [ $max -le 1000 ]
  do
    nb=0; while read lgt; do if [ $lgt -lt $min ]; then continue; elif [ $lgt -le $max ]; then nb=$(( $nb + 1 )); else break; fi; done < $LGTFILE;
    echo "  $min $max   $nb"; min=$(( $min + 100 )); max=$(( $max + 100 )); 
  done
  max=2000
  while [ $max -le 10000 ]
  do
    nb=0; while read lgt; do if [ $lgt -lt $min ]; then continue; elif [ $lgt -le $max ]; then nb=$(( $nb + 1 )); else break; fi; done < $LGTFILE;
    echo "  $min $max   $nb"; min=$(( $max + 1 )); max=$(( $max + 1000 )); 
  done
  up=$(( $MAX + 5000 ))
  max=15000
  while [ $max -le $up ]
  do
    nb=0; while read lgt; do if [ $lgt -lt $min ]; then continue; elif [ $lgt -le $max ]; then nb=$(( $nb + 1 )); else break; fi; done < $LGTFILE;
    echo "  $min $max   $nb"; min=$(( $max + 1 )); max=$(( $max + 5000 )); 
  done
  echo
fi

rm $LGTFILE $FH $INFILE_UNIX $SEQS;
if [ -e $INFOFILE ]; then rm $INFOFILE; fi








