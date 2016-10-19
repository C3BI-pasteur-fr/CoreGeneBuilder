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
# fasta_rename_slave_genes.py: transfers gene products of a reference annotation file to the current proteome fasta #
#                              file using ecamber homolog cluster features.                                         #
# Author: Elise Larsonneur                                                                                          #
#####################################################################################################################


from Bio import SeqIO
import sys
import os.path
import getopt


def usage():
      print 'Usage: ' + str(sys.argv[0]) + \
      ' -i ecamber_annotation_input.fasta -o ecamber_annotation_output.fasta -c clusters_in.txt -n strain_name \n' + \
      'Transfers gene products of a reference annotation file to the current proteome ' + \
      'fasta file using ecamber homolog cluster features.'


def main(argv):

 fasta_in = ''      # -i 
 strain_in = ''     # -n
 clusters_in = ''   # -c
 fasta_out = ''     # -o

 if (len(sys.argv) <= 1):
      usage()
      sys.exit(2)

 try:
      opts, args = getopt.getopt(argv,"hi:o:c:n:")
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
      elif (opt == "-n"):
         strain_in = arg
      elif (opt == "-c"):
         clusters_in = arg
         if not os.path.exists(clusters_in):
             print '[ERROR] file \'' + str(clusters_in) + '\' does not exist'
             sys.exit(1)
      elif (opt == "-o"):
         fasta_out = arg
         if not os.path.exists(os.path.dirname(fasta_out)):
             print '[ERROR] directory of output file \'' + str(fasta_out) + '\' does not exist'
             sys.exit(1)

 if ("-i" in options) and ("-n" in options) and ("-c" in options) and ("-o" in options):
    pass
 else:
    usage()
    sys.exit(2)
 
 #print "input_fasta:\t\t", fasta_in
 #print "strain_name:\t\t", strain_in
 #print "input_clusters:\t\t", clusters_in
 #print "output_fasta:\t\t", fasta_out

 output_fasta_h = open(fasta_out,"w")
 input_clusters_h = open(clusters_in,"r")

 # fastain example :
 #grep '>' datasets/testlimoecamber/output/genes_aa/NC_003210_limo_EGDe.fasta|less
 #>x|2203|x|-|1275778|1275584|NC_003210_limo_EGDe                            # case 3 : new annotation by ecamber
 #> gene_id|multigene_cluster_id|gene_name|strand|start|stop|contig_id       #legend

 # input_clusters example (new genes annotated by ecamber are not in this file) :
 #cluster_id  cluster_gene_name   cluster_gene_id    cluster_gene_product            # header
 #2414    dnaA    lmo0001    chromosome_replication_initiator_DnaA 
 #1320    dnaN    lmo0002    DNA_polymerase_III_subunit_beta 
 #1305    lmo0003 lmo0003    hypothetical_protein     
 #2951    lmo0174;lmo0329;lmo0827 lmo0174;lmo0329;lmo0827 transposase 


 # store input_clusters data in dictionnary
 clusters={}
 for cluster_line in input_clusters_h.readlines():
    cluster_line = cluster_line.replace("\n","")
    if "cluster_id" in cluster_line:  #skip header line
      continue
    
    cluster_info = cluster_line.split("\t")
    cluster_id = cluster_info[0]
    cluster_gene_name = cluster_info[1]
    cluster_gene_id = cluster_info[2]
    cluster_gene_product = cluster_info[3]

    clusters[cluster_id] = (str(cluster_gene_name), str(cluster_gene_id), str(cluster_gene_product))
 input_clusters_h.close()
 # print "number_of_input_clusters: ", len(clusters)

 for seq_record in SeqIO.parse(fasta_in, "fasta"):

    sequence=str(seq_record.seq).upper()
    #breaks sequence every <freq> characters  
    #freq=100
    #sequence=('\n'.join(sequence[i:i+freq] for i in xrange(0, len(sequence), freq)))

    id=seq_record.id
    id_list = id.split("|")
    id_geneid = id_list[0]
    id_clusterid =  id_list[1]
    id_genename = id_list[2]
    id_strand = id_list[3] 
    id_start = id_list[4]
    id_stop = id_list[5]
    id_contigid = id_list[6]
    locustag_t = 'x'
    product_t = 'x'

    if (str(id_clusterid) in clusters):
       id_genename_t = clusters[id_clusterid][0]
       id_geneid_t = clusters[id_clusterid][1]
       id_genename = id_genename_t  
       locustag_t = id_geneid_t     
       product_t = clusters[id_clusterid][2]

    # id - eg output fasta headers
    # $1  gene_id (prodigal_id or 'x')
    # $2  cluster_id
    # $3  gene_name (tranferred gene_name(s) OR locus_tag(s) OR 'x')
    # $4  strand
    # $5  start
    # $6  stop
    # $7  contig_id
    # $8  locus_tag (one or several transferred locus_tag(s) separated by a ';' OR 'x')
    # $9  gene_product (one product or several transferred product(s) separated by a ';' OR 'x')
    # $10 strain
    id = str(id_geneid) + "|" + str(id_clusterid) + "|" + str(id_genename) + "|" + str(id_strand) + "|" + str(id_start) + "|" + str(id_stop) + "|" + str(id_contigid) + "|" + str(locustag_t) + "|"  + str(product_t) + "|" + str(strain_in)
    output_fasta_h.write(">" + id + "\n" + sequence  + "\n" )

 output_fasta_h.close()


# MAIN
if __name__ == "__main__":
   main(sys.argv[1:])

