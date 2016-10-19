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
# update_cluster_table.py: extends the list of homolog clusters in ecamber table 'clusters_table.txt'               # 
#                          (adding clusters with locus_tag but no gene name, and clusters without                   #
#                          gene_name or locus_tag) and outputs it to a new table.                                   #
# Author: Elise Larsonneur                                                                                          #
#####################################################################################################################


from Bio import SeqIO
import sys
import os.path
import getopt


def usage():
   print 'Usage: ' + str(sys.argv[0]) + ' ecamber_annotation_input.fasta annotation_output.fasta clusters_input.txt\n' + \
   'Extends the list of homolog clusters in ecamber table \'clusters_table.txt\' ' + \
   '(adding clusters with locus_tag but no gene name, and clusters without ' + \
   'gene_name or locus_tag) and outputs it to a new table.'
    

def main(argv):

   clusters_in = ''          # -i 
   clusters_table_in = ''    # -c
   clusters_table_out = ''   # -o

   if (len(sys.argv) <= 1):
      usage()
      sys.exit(2)

   try:
      opts, args = getopt.getopt(argv,"hi:c:o:")
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
         clusters_in = arg
         if not os.path.exists(clusters_in):
             print '[ERROR] file ' + str(clusters_in) + ' does not exist'
             sys.exit(1)
      elif (opt == "-c"):
         clusters_table_in = arg
         if not os.path.exists(clusters_table_in):
             print '[ERROR] file ' + str(clusters_table_in) + ' does not exist'
             sys.exit(1)
      elif (opt == "-o"):
         clusters_table_out = arg
         if not os.path.exists(os.path.dirname(clusters_table_out)):
             print '[ERROR] directory of output file' + str(clusters_table_out) + ' does not exist'
             sys.exit(1)

   if ("-i" in options) and ("-c" in options) and ("-o" in options):
      pass
   else:
      usage()
      sys.exit(2)

   #print "input_clusters:\t\t", clusters_in
   #print "input_clusters_table:\t\t", clusters_table_in
   #print "output_clusters_table:\t\t", clusters_table_out

   clusters_in_h = open(clusters_in,"r")
   clusters_table_in_h = open(clusters_table_in,"r")
   clusters_table_out_h = open(clusters_table_out,"w")

   # input_clusters example (new genes annotated by ecamber are not in this file) :
   #cluster_id  cluster_gene_name   cluster_gene_id            # header
   #2414    dnaA    lmo0001
   #1320    dnaN    lmo0002
   #1305    lmo0003 lmo0003

   # head datasets/testlimoecamber/output/clusters_table.txt
   #cluster_id  cluster_gene_name   cluster_type    cluster_mg_count    cluster_mg_ann  cluster_strain_count    cluster_strain_ann  NC_003210_limo_EGDe NC_012488_limo_Clip80459    NZ_LMXJ00000000_limo_SLCC208
   #1200    x   ANCHOR  1   1   1   1           2_435.407998.409521.+.contig_2
   #1615    flgG    ANCHOR  3   3   3   3   lmo0682.717670.718449.+ 1_693.726153.726932.+   10_9.8219.7440.-.contig_10
   #347 x   ANCHOR  2   2   2   2   lmo2300.2385332.2384001.-       2_406.384206.385537.+.contig_2
   #340 _synonym="ychF  ANCHOR  3   3   3   3   lmo2779.2864661.2863561.-   1_2728.2837207.2836107.-    2_206.198344.197244.-.contig_2


   # store input_clusters data in dictionnary
   clusters={}
   for cluster_line in clusters_in_h.readlines():
       cluster_line = cluster_line.replace("\n","")
       if "cluster_id" in cluster_line:
         continue
       else : 
         cluster_info = cluster_line.split("\t")
         cluster_id = cluster_info[0]
         cluster_gene_name = cluster_info[1]
         cluster_gene_id = cluster_info[2]
         if cluster_id not in clusters :
            clusters[cluster_id] = (str(cluster_gene_name), str(cluster_gene_id))
#         else :
#            print "dic_cluster_id", cluster_id,  "\t", clusters[cluster_id][0], "\t", clusters[cluster_id][1]
#            print "file_cluster_id", cluster_id, "\t", cluster_gene_name, "\t", cluster_gene_id
   clusters_in_h.close()
   print "number_of_input_clusters: ", len(clusters)

   # parse input clusters_table.txt
   for cluster_tab_line in clusters_table_in_h.readlines():
       if "cluster_id" in cluster_tab_line:
          clusters_table_out_h.write(cluster_tab_line)
        
       else :
          cl_infos = cluster_tab_line.split("\t")
          cl_cluster_id = cl_infos[0]
          #print "cl_cluster_id"+str(cl_cluster_id)+"___END" 
          cl_cluster_gene_name = cl_infos[1]

          # write into new clusters_table.txt
          if (str(cl_cluster_id) in clusters):
             cl_gene_name_t = clusters[cl_cluster_id][0]
             cl_gene_id_t = clusters[cl_cluster_id][1]

             strout = str(cl_cluster_id) + "\t" + str(cl_gene_name_t)
             for i in range(2,len(cl_infos)-1):
                strout = strout + "\t" + str(cl_infos[i])
             clusters_table_out_h.write( strout + "\n")
          else :
             clusters_table_out_h.write( str(cluster_tab_line) )

   clusters_table_in_h.close()
   clusters_table_out_h.close()


if __name__ == '__main__':
   main(sys.argv[1:])
