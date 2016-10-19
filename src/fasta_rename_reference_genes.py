#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
# fasta_rename_reference_genes.py: updates description of ecamber homolog clusters by adding gene features          # 
#                                  (gene name, locus_tag, product) from the reference annotation file, in the       #
#                                  reference proteome fasta file.                                                   #
# Author: Elise Larsonneur                                                                                          #
#####################################################################################################################


from Bio import SeqIO
import sys
import os.path
import getopt


def usage():
      print 'Usage: ' + str(sys.argv[0]) + \
      ' -i ref_ecamber_annotation_input.fasta -o ref_annotation_output.fasta -c clusters.txt -a ref_gb_annotation_in.lst -n strain_name\n' + \
      'Updates description of ecamber homolog clusters by adding gene features (gene name, locus_tag, product) ' + \
      'from the reference annotation file, in the reference proteome fasta file.'
  

def main(argv):

 fasta_in = ''      # -i 
 gb_annot_in = ''   # -a 
 strain_in = ''     # -n
 fasta_out = ''     # -o
 clusters_out = ''  # -c

 if (len(sys.argv) <= 1):
      usage()
      sys.exit(2)

 try:
      opts, args = getopt.getopt(argv,"hi:o:c:a:n:")
 except getopt.GetoptError:
      usage()
      sys.exit(2)

 if not opts:
     print "Bad argument"
     usage()
     sys.exit(2)

 options = []

 for opt, arg in opts:
      options.append(opt)
      if (opt == "-h"):
         usage()
         sys.exit()
      elif (opt == "-i"):
         fasta_in = arg
         if not os.path.exists(fasta_in):
             print '[ERROR] file \'' + str(fasta_in) + '\' does not exist'
             sys.exit(1)
      elif (opt == "-a"):
         gb_annot_in = arg
         if not os.path.exists(gb_annot_in):
             print '[ERROR] file \'' + str(gb_annot_in) + '\' does not exist'
             sys.exit(1)
      elif (opt == "-n"):
         strain_in = arg
      elif (opt == "-o"):
         fasta_out = arg
         if not os.path.exists(os.path.dirname(fasta_out)):
             print '[ERROR] directory of output file \'' + str(fasta_out) + '\' does not exist'
             sys.exit(1)
      elif (opt == "-c"):
         clusters_out = arg
         if not os.path.exists(os.path.dirname(clusters_out)):
             print '[ERROR] directory of output file \'' + str(clusters_out) + '\' does not exist'
             sys.exit(1)

 if ("-i" in options) and ("-a" in options) and ("-n" in options) and ("-o" in options) and ("-c" in options):
    pass
 else:
    usage()
    sys.exit(2)

 #print "input_fasta:\t\t", fasta_in
 #print "input_gbk_annotation_lst:\t\t", gb_annot_in
 #print "strain_name:\t\t", strain_in
 #print "output_fasta:\t\t", fasta_out
 #print "output_clusters:\t\t", clusters_out


 #store CDS features from formatted genbank file (.gb.txt) to get gene products
 gbk_cds = {}
 fh_gbk = open(gb_annot_in,"r")
 for line in fh_gbk.readlines():
     info = line.replace("\n","").split("\t")
     a_accession = info[0]
     a_locus_tag = info[1]
     a_gene_name = info[2]
     a_start = info[3]
     a_stop = info[4]
     a_strand = info[5]
     a_product = info[6] # space characters replaced by '_'
     tuple_info = (a_accession, a_gene_name, a_start, a_stop, a_strand, a_product)
     if a_locus_tag not in gbk_cds :
         gbk_cds[a_locus_tag] = tuple_info


 # fastain example :
 # grep '>' datasets/testlimoecamber/output/genes_aa/NC_003210_limo_EGDe.fasta|less
 # >lmo0001|2414|dnaA|+|318|1673|NC_003210_limo_EGDe                          # case 1 : reference annotation with provided gene_name
 # >lmo0583|1282|_synonym="azi|+|620805|623135|NC_003210_limo_EGDe            # case 1
 # >lmo0003|1305|x|+|3121|4464|NC_003210_limo_EGDe                            # case 2 : reference annotation with no gene_name provided but locus_tag
 # >x|2203|x|-|1275778|1275584|NC_003210_limo_EGDe                            # case 3 : new annotation by ecamber
 # > gene_id|multigene_cluster_id|gene_name|strand|start|stop|contig_id       #legend


 #store clusters
 clusters={}

 # first parsing of fasta file : we only store the clusters
 output_fasta_h = open(fasta_out,"w")

 for seq_record in SeqIO.parse(fasta_in, "fasta"):

    id=seq_record.id
    id_list = id.split("|")
    id_geneid = id_list[0]
    id_clusterid =  id_list[1]
    id_genename = id_list[2]
    id_strand = id_list[3] 
    id_start = id_list[4]
    id_stop = id_list[5]
    id_contigid = id_list[6]
    product = 'x'
 
    if id_geneid in gbk_cds:
             product = gbk_cds[id_geneid][5]

    # case 1 (remove pattern '_synonym="') 
    if ( '_synonym="' in str(id_genename) ):
        id_genename = str(id_genename).replace('_synonym="','')   
 
    # if (case 2)
    if ( str(id_geneid) != "x" ) and ( str(id_genename) == "x" ):
        id_genename = str(id_geneid)       

    # if (case 1)
    if ( str(id_geneid) != "x" ) and ( str(id_genename) != "x" ) :
         if id_clusterid not in clusters:
             clusters[id_clusterid] = (str(id_genename), str(id_geneid), str(product))
         else :
             cl_genename = str(clusters[id_clusterid][0])
             if (str(id_genename) not in cl_genename.split(";")):    
                 cl_genename = str(clusters[id_clusterid][0]) + ";" + str(id_genename)
             cl_geneid = str(clusters[id_clusterid][1])
             if (str(id_geneid) not in cl_geneid.split(";")): 
                 cl_geneid = str(clusters[id_clusterid][1]) + ";" + str(id_geneid)
             cl_product = str(clusters[id_clusterid][2])
             if (str(product) not in cl_product.split(";")):
                 cl_product = str(clusters[id_clusterid][2]) + ";" + str(product)
             clusters[id_clusterid] = (str(cl_genename), str(cl_geneid), str(cl_product))
 
  
 print "ecamber_cluster_number :", len(clusters)

 output_clusters_h = open(clusters_out,"w")
 output_clusters_h.write("cluster_id\tcluster_gene_name\tcluster_gene_id\tcluster_gene_product\n")

 for cl_id in clusters:
     cl_genename = clusters[cl_id][0]
     cl_geneid = clusters[cl_id][1]
     cl_product = clusters[cl_id][2]
     output_clusters_h.write(str(cl_id) + "\t" + str(cl_genename) + "\t" + str(cl_geneid) + "\t" + str(cl_product)  + "\n")

 output_fasta_h.close()
 output_clusters_h.close()


 #parse again the fasta file : we generate the output fasta file using clusters dictionnary
 output_fasta_h = open(fasta_out,"w")

 for seq_record in SeqIO.parse(fasta_in, "fasta"):

    sequence=str(seq_record.seq).upper()

    id=seq_record.id
    id_list = id.split("|")
    id_geneid = id_list[0]
    id_clusterid =  id_list[1]
    id_genename = id_list[2]
    id_strand = id_list[3] 
    id_start = id_list[4]
    id_stop = id_list[5]
    id_contigid = id_list[6]
    product = 'x'
    
    #get gene product
    if id_geneid in gbk_cds:
         product = gbk_cds[id_geneid][5]
    else:
         if id_clusterid in clusters :
             product = str(clusters[id_clusterid][2])

    # case 1 (remove pattern '_synonym="') 
    if ( '_synonym="' in str(id_genename) ):
        id_genename = str(id_genename).replace('_synonym="','')   
    
    # if (case 2)
    if ( str(id_geneid) != "x" ) and ( str(id_genename) == "x" ):
        #if id_clusterid in clusters :
        #     cl_genename = str(clusters[id_clusterid][0])
        #     cl_geneid = str(clusters[id_clusterid][1])
        #     id_genename = str(cl_genename)   #we keep one or several locus_tag values as gene_name
        #else :
        #     id_genename = str(id_geneid) 
        id_genename = str(id_geneid)           #we keep only one locus_tag value as gene_name

    # if (case 1)
    if ( str(id_geneid) != "x" ) and ( str(id_genename) != "x" ) :
        #if id_clusterid in clusters :
        #     cl_genename = str(clusters[id_clusterid][0])
        #     cl_geneid = str(clusters[id_clusterid][1])
        #     id_genename = str(cl_genename)    #we keep one or several gene_name values as gene_name
        id_genename = str(id_genename)          #we keep only one gene_name value as gene_name

    # id - eg output fasta headers
    # $1  gene_id (locus_tag or 'x')
    # $2  cluster_id
    # $3  gene_name (gene_name or locus_tag or 'x')
    # $4  strand
    # $5  start
    # $6  stop
    # $7  contig_id
    # $8  gene_id (locus_tag or 'x')
    # $9  gene_product (one product or several products separated by a ';')
    # $10 strain
    id = str(id_geneid) + "|" + str(id_clusterid) + "|" + str(id_genename) + "|" + str(id_strand) + "|" + str(id_start) + "|" + str(id_stop) + "|" + str(id_contigid) + "|" + str(id_geneid) + "|"  + str(product) + "|" + str(strain_in)
    output_fasta_h.write(">" + id + "\n" + sequence  + "\n" )

 output_fasta_h.close()


# MAIN
if __name__ == "__main__":
   main(sys.argv[1:])

