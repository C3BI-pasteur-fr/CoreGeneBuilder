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
# ecamber2gembase.sh: formats sequence headers of fasta files outputed by ecamber into specific header format.      #
# Authors: Damien Mornico, Elise Larsonneur                                                                         #
#####################################################################################################################


set -o pipefail


##########################################################
# Inputs : genes and proteins outputed by ecamber (fasta)
# Outputs : genes, proteins in specific format (fasta)
# 
# Script in two parts : gene parsing then protein parsing
# For each part, we parse a first time the fasta file
# (awk BEGIN) to get informations (start and stop codons, etc) 
# We parse a second time to format header informations and write them; and to write sequence linked to each header
#########################################################


########################################################################################
# FUNCTIONS                                                                            #
########################################################################################

# display_usage 
# This function displays the usage of this program.
# No parameters
function display_usage { 
  echo -e "Usage:\n$0 genome_name gene_directory prot_directory \n ex : $0 esco.001.c01.60 ../data/esco/genes ../data/esco/proteins" 
} 


########################################################################################
# MAIN FUNCTION                                                                        #
########################################################################################
# main
# Parameters
# - 1) Input fasta file basename - A string containing the file path.
# - 2) Input directory 'genes' - A string containing the directory path where is the fasta file.
# - 3) Input directory 'proteins'  - A string containing the directory path where is the fasta file.
# Outputs 
# -  a fasta file (modified input fasta). 
function main {

  # if less than two arguments supplied, display usage 
  if [[ $# -ne 3 ]]; then 
    display_usage
    exit 1
  fi 
 
  # check whether user had supplied -h or --help . If yes display usage 
  if [[ $# == "--help" ]] || [[ $# == "-h" ]]; then
    display_usage
    exit 0
  fi

  #name of intputs
  gene_file="$2/$1.fasta"
  prot_file="$3/$1.fasta"

  #delete existing output files
  rm "$2/$1.lst" > /dev/null 2>&1
  rm "$3/$1.lst" > /dev/null 2>&1
  rm "$2/$1.gen" > /dev/null 2>&1
  rm "$3/$1.prt" > /dev/null 2>&1

  # 4-letters code followed by the strain number (ex : klpn.0000001)
  genome_name=$(echo "$1" | gawk -v FS="." '{print toupper($1) $2}')



  # fasta headers to parse :

  #case 1 : de novo prodigal annotation
  ## grep ">" ../ext-tools/ecamber/datasets/klpn5/output/genes_aa/klpn.0000001.c001.fasta |head
  ## >1_1|2837|x|+|340|2802|klpn_0000001_c001|x|x|MGH_78578_25
  ## >1_2|5929|x|+|2804|3733|klpn_0000001_c001|x|x|MGH_78578_25
  ## >1_3|2779|x|+|3737|5017|klpn_0000001_c001|x|x|MGH_78578_25
  ## >1_4|1924|x|+|5306|5710|klpn_0000001_c001|x|x|MGH_78578_25
  ## >1_5|90|x|-|6552|5779|klpn_0000001_c001|x|x|MGH_78578_25
  ## >1_6|3057|x|-|8060|6630|klpn_0000001_c001|x|x|MGH_78578_25
  ## >1_7|308|x|+|8263|9216|klpn_0000001_c001|x|x|MGH_78578_25
  ## >1_8|2869|x|+|9316|9903|klpn_0000001_c001|x|x|MGH_78578_25
  ## >1_9|3842|x|+|10023|11327|klpn_0000001_c001|x|x|MGH_78578_25
  ## >1_10|84|x|-|11956|11390|klpn_0000001_c001|x|x|MGH_78578_25

  #case 2 : de novo prodigal annotation + reference genbank annotation 
  #         (transfer of gene names of the reference annotation (genbank) to the others annotations)
  ## grep ">" ../ext-tools/ecamber/datasets/klpn5refannot/output/genes_aa/klpn.0000001.c001.fasta |head
  ## >KPN_RS00005|7944|KPN_RS00005|+|21|137|klpn_0000001_c001|KPN_RS00005|hypothetical_protein|MGH78578_NC
  ## >KPN_RS00010|2848|thrA|+|340|2802|klpn_0000001_c001|KPN_RS00010|bifunctional_aspartokinase_I/homoserine_dehydrogenase_I|MGH78578_NC
  ## >KPN_RS00015|6009|KPN_RS00015|+|2804|3733|klpn_0000001_c001|KPN_RS00015|homoserine_kinase|MGH78578_NC
  ## >KPN_RS00020|2784|KPN_RS00020|+|3737|5017|klpn_0000001_c001|KPN_RS00020|threonine_synthase|MGH78578_NC
  ## >KPN_RS00025|1898|KPN_RS00025|+|5354|5710|klpn_0000001_c001|KPN_RS00025|hypothetical_protein|MGH78578_NC
  ## >KPN_RS00030|91|KPN_RS00030|-|6552|5779|klpn_0000001_c001|KPN_RS00030|hypothetical_protein|MGH78578_NC
  ## >KPN_RS00035|3068|KPN_RS00035|-|8060|6630|klpn_0000001_c001|KPN_RS00035|sodium:alanine_symporter|MGH78578_NC
  ## >KPN_RS00040|311|KPN_RS00040|+|8263|9216|klpn_0000001_c001|KPN_RS00040|transaldolase|MGH78578_NC
  ## >KPN_RS00045|2878|mogA|+|9316|9903|klpn_0000001_c001|KPN_RS00045|molybdopterin_adenylyltransferase|MGH78578_NC
  ## >KPN_RS00050|3851|KPN_RS00050|+|10023|11327|klpn_0000001_c001|KPN_RS00050|MFS_transporter|MGH78578_NC
 
  ## legend of the fields separated by '|' to parse :
  ## $1  gene_id (=locus_tag if reference)
  ## $2  cluster_multigene_id
  ## $3  gene_name
  ## $4  strand (+/-)
  ## $5  start/stop
  ## $6  stop/start
  ## $7  contig_id (=fastaname if genome have only ONE contig)
  ## $8  locus_tag (one or several transferred locus_tag(s) separated by a ';' OR 'x') 
  ## $9  gene_product(s)
  ## $10 strain



  #Parsing DNA fasta (GENES)
  tr -d '\15\32' < "${gene_file}" | awk 'BEGIN{seq=""} !/^>/ {seq=seq$0} /^>/ {if(seq!=""){print seq;seq=""}{print $0}} END{print seq}' | grep -v "^$" > "${gene_file}.tmp" || { echo "[ERROR] error when formatting the file ${gene_file}" >&2; exit 1; } 
  mv "${gene_file}.tmp" "${gene_file}" || { echo "[ERROR] error when changing the name of file ${gene_file}.tmp" >&2; exit 1; }

  gawk -v genome_file="$1" -v genes="${gene_file}" -v geneoutput="$2" -v genome="${genome_name}" '

    BEGIN {
		
	Tdir[1]="";
	Tbeg[1]="";
	Tend[1]="";
	Tcontig[1]="";
	TStart[1]="";
	LineNb=0;
	HeaderLineNb=0;
	LastLineSeq="";
	Tfasta[1]="";
	GeneNb=0;
	Contig="";
	ContigNb=0;
        LastTcontigId="";
        LastTcontigName="";

	#First parsing of the fasta file, to store information		
	 while ((getline line < genes) > 0){
	
		#As NR is unavailable in BEGIN, Line nb is stored in a variable
	 	LineNb++
		#File line keep in a table, used for Stop codon
		Tfasta[LineNb]=line;
		
		
		
		#Header case
	 	if(line ~ "^>.*"){
		
			#New gene/prot
			GeneNb++
			#Keep header line nb to identify first sequence line
			HeaderLineNb=LineNb
			
			#header split by "|"
			nb=split(line,T,"|")
			
			id=T[1] T[5];
			Tid[id]=GeneNb	#GeneNb stored in a table
			TGene[GeneNb]=id	#Id stored in a table
			
			#Direction of gene stored in table (with Direct and Complement, instead of + and -)
			if(T[4]=="+"){
				Tdir[id]="D"
				
			}else{
				Tdir[id]="C"
			}
			
			#Begin and ENd positions stored in table
			if(T[5]<T[6]){
			
				Tbeg[id]=T[5] 		
				Tend[id]=T[6] 
			
			}else{
			
				Tbeg[id]=T[6] 		
				Tend[id]=T[5]
			
			}
			 		
			
			#Format contig id (AA, AB,AC...)
			if(T[7]!=Contig){
				
				Tcontig[TGene[GeneNb-1]]="g" contig
				
				ContigNb++
                                
                                #Format contig id (AA, AB,AC...,ZZ)
				#p=97+int(ContigNb/26); 
				#q=(ContigNb%26)+96; 
                                #if (q == (0 + 96)){
                                #    q = q + 26
                                #    p = p - 1
                                #} 
				#contig=toupper(sprintf("%c%c", p, q));  ## ASCII number formatting

				#Format contig id (AAA, AAB, AAC...ZZZ)
                                p=97+int(ContigNb/26);
                                q=(ContigNb%26)+96;
                                if (q == (0 + 96)) {
                                     q = q + 26
                                     p = p - 1
                                };
                                o=97+int(ContigNb/676);
                                oi=int(ContigNb/676); 
                                oii=(ContigNb%676);
                                p = p - (26*oi)      
                                if (oii == 0) {
                                     o = o - 1
                                     p = 26 + 96 
                                };
                                contig=toupper(sprintf("%c%c%c", o, p, q));  ## ASCII number formatting

                                #Id of the contig stored in table
				Tcontig[id]="g" contig
                                #Id and Name of the last contig stored in the table 
                                LastTcontigId=id;
                                LastTcontigName="g" contig			      
    	
				Contig=T[7]		
			}else{
				#Id of the contig stored in table
				Tcontig[id]="i" contig
                                #Id and Name of the last contig stored in the table
                                LastTcontigId=id; 
                                LastTcontigName="g" contig
			}
			
			
			
			
			
		#Sequence case			
		}else{
			#start codon
			if(LineNb==HeaderLineNb+1){
				
				TStart[id]=substr(line,1,3)
                                TStop[id]=substr(line,length(line)-2,length(line))
			}
			
		}
	 }
    }

    #Second parsing
    {
	
	
	#Header case
	 if($0 ~ "^>.*"){
	 
		#header split by "|"
		nb=split($0,T,"|")
		
		id=T[1] T[5]
		
		#Nb of 0 in gene_id
		if(length(Tid[id])==4){
			zero="0";
		}else if(length(Tid[id])==3){
			zero="00";		
		}else if(length(Tid[id])==2){
			zero="000";		
		}else if(length(Tid[id])==1){
			zero="0000";		
		}
	        
                if(id==LastTcontigId){
                    Tcontig[id] = LastTcontigName
                }
	
		printf genome  Tcontig[id] "_" zero Tid[id] "0 " Tdir[id] " " TStart[id] " " TStop[id] " " Tbeg[id] " " Tend[id] " CDS " T[3] " " Tend[id]-Tbeg[id]+1 " " T[7] " " T[8] " " T[9] " " T[10] "\n" >> geneoutput "/" genome_file".lst"
		printf ">"genome  Tcontig[id] "_" zero Tid[id] "0 " Tdir[id] " " TStart[id] " " TStop[id] " " Tbeg[id] " " Tend[id] " CDS " T[3] " " Tend[id]-Tbeg[id]+1 " " T[7] " " T[8] " " T[9] " " T[10] "\n" >> geneoutput "/"genome_file".gen"
                #debug #if(genome == "KLPN002" && Tcontig[id] == "iAA" && Tid[id] == "20"){ printf "line214-"TStop[id] "\n" }
	}else{
		printf $0 "\n" >> geneoutput "/"genome_file".gen"
	
	}
	
    } ' "${gene_file}" || { echo "[ERROR] error when parsing the file ${gene_file}" >&2; exit 1; }




  #parsing AA fasta (PROTEINS)
  tr -d '\15\32' < "${prot_file}" | awk 'BEGIN{seq=""} !/^>/ {seq=seq$0} /^>/ {if(seq!=""){print seq;seq=""}{print $0}} END{print seq}' | grep -v "^$" > $prot_file.tmp || { echo "[ERROR] error when formatting the file ${prot_file}" >&2; exit 1; }
  mv $prot_file.tmp $prot_file || { echo "[ERROR] error when changing the name of file ${prot_file}.tmp" >&2; exit 1; }

  gawk -v genome_file="$1" -v genes="${gene_file}" -v prot="${prot_file}" -v protoutput="$3" -v genome="${genome_name}" '

    BEGIN {
	
	Tdir[1]="";
	Tbeg[1]="";
	Tend[1]="";
	Tcontig[1]="";
	TStart[1]="";
	LineNb=0;
	HeaderLineNb=0;
	LastLineSeq="";
	Tfasta[1]="";
	GeneNb=0
	Contig=""
	ContigNb=0
        LastTcontigId="";
        LastTcontigName="";

	#First parsing of the fasta file, to store information		
	while ((getline line < genes) > 0){
	
		#As NR is unavailable in BEGIN, Line nb is stored in a variable
	 	LineNb++
		#File line keep in a table, used for Stop codon
		Tfasta[LineNb]=line;
		
		
		
		#Header case
	 	if(line ~ "^>.*"){
		
			#New gene/prot
			GeneNb++
			#Keep header line nb to identify first sequence line
			HeaderLineNb=LineNb
			
			#header split by "|"
			nb=split(line,T,"|")
			
			id=T[1] T[5];
			Tid[id]=GeneNb	#GeneNb stored in a table
			TGene[GeneNb]=id	#Id stored in a table
			
			#Direction of gene stored in table (with Direct and Complement, instead of + and -)
			if(T[4]=="+"){
				Tdir[id]="D"
				
			}else{
				Tdir[id]="C"
			}
			
			#Begin and ENd positions stored in table
			if(T[5]<T[6]){
			
				Tbeg[id]=T[5] 		
				Tend[id]=T[6] 
			
			}else{
			
				Tbeg[id]=T[6] 		
				Tend[id]=T[5]
			
			}
			 		
			
			#Format contig id (AA, AB,AC...)
			if(T[7]!=Contig){
				
				Tcontig[TGene[GeneNb-1]]="g" contig
				
				ContigNb++
                                
				#Format contig id (AAA, AAB, AAC...ZZZ)
                                p=97+int(ContigNb/26);
                                q=(ContigNb%26)+96;
                                if (q == (0 + 96)) {
                                     q = q + 26
                                     p = p - 1
                                };
                                o=97+int(ContigNb/676);
                                oi=int(ContigNb/676); 
                                oii=(ContigNb%676);
                                p = p - (26*oi)      
                                if (oii == 0) {
                                     o = o - 1
                                     p = 26 + 96 
                                };
                                contig=toupper(sprintf("%c%c%c", o, p, q));  ## ASCII number formatting

				#Id of the contig stored in table
				Tcontig[id]="g" contig 
                                #Id and Name of the last contig stored in the table 
                                LastTcontigId=id;
                                LastTcontigName="g" contig
				
				Contig=T[7]		
			}else{
				#Id of the contig stored in table
				Tcontig[id]="i" contig
                                #Id and Name of the last contig stored in the table 
                                LastTcontigId=id;
                                LastTcontigName="g" contig
			}
			
			
			
			
			
		#Sequence case			
		}else{
			#start codon
			if(LineNb==HeaderLineNb+1){
				
				TStart[id]=substr(line,1,3)
                                TStop[id]=substr(line,length(line)-2,length(line))
			}
			
		}
	 }
    }

    #Second parsing
    {
	
	
	#Header case
	 if($0 ~ "^>.*"){
	 
		#header split by "|"
		nb=split($0,T,"|")
		
		id=T[1] T[5]
		
		#Nb of 0 in gene_id
		if(length(Tid[id])==4){
			zero="0";
		}else if(length(Tid[id])==3){
			zero="00";		
		}else if(length(Tid[id])==2){
			zero="000";		
		}else if(length(Tid[id])==1){
			zero="0000";		
		}

                if(id==LastTcontigId){
                    Tcontig[id] = LastTcontigName
                }

	        #debug #if(genome == "KLPN002" && Tcontig[id] == "iAA" && Tid[id] == "20"){ printf "line385-"TStop[id]"\n" }	
		printf genome  Tcontig[id] "_" zero Tid[id] "0 " Tdir[id] " " TStart[id] " " TStop[id] " " Tbeg[id] " " Tend[id] " CDS " T[3] " " Tend[id]-Tbeg[id]+1 " " T[7] " " T[8] " " T[9] " " T[10] "\n" >> protoutput "/" genome_file".lst"
                printf ">"genome  Tcontig[id] "_" zero Tid[id] "0 " Tdir[id] " " TStart[id] " " TStop[id] " " Tbeg[id] " " Tend[id] " CDS " T[3] " " Tend[id]-Tbeg[id]+1 " " T[7] " " T[8] " " T[9] " " T[10] "\n" >> protoutput "/" genome_file".prt"
	}else{
		printf $0 "\n" >> protoutput "/" genome_file".prt"
	
	}
	


    } ' "${prot_file}" || { echo "[ERROR] error when parsing the file ${prot_file}" >&2; exit 1; } 
}

main "$@"

