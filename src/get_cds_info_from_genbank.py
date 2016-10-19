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
# get_cds_info_from_genbank.py: outputs the list of CDS and their features from a genbank annotation file           #
#                               to a tab-delimited file.                                                            #
# Author: Elise Larsonneur                                                                                          #
#####################################################################################################################


from Bio import SeqIO
import sys
import os.path
import getopt


def usage():
      print 'Usage: ' + str(sys.argv[0]) + ' -i annotation.gb\n' + \
      'Outputs the list of CDS and their features from a genbank annotation file ' + \
      'to a tab-delimited file.'

     
def main(argv):

   gb_file = ''      # -i 

   if (len(sys.argv) <= 1):
      usage()
      sys.exit(2)
 
   try:
      opts, args = getopt.getopt(argv,"hi:")
   except getopt.GetoptError:
      usage() 
      sys.exit(2)

   if not opts:
     print "Bad argument"
     usage()
     sys.exit(2)

   for opt, arg in opts:
      if (opt == "-h"):
         usage()
         sys.exit()
      elif (opt == "-i"):
         gb_file = arg
         if not os.path.exists(gb_file):
             print '[ERROR] file \'' + str(gb_file) + '\' does not exist'
             sys.exit(1)

   out_file = open(str(gb_file) + ".txt", 'w')
   #print "infile:", gb_file
   #print "outfile:", str(gb_file) + ".txt" 

   out_file.write("accession\tlocus_tag\tgene\tstart\tend\tstrand\tproduct\n")

   #case multiple records
   for gb_record in SeqIO.parse(gb_file, "genbank"):
      accession=gb_record.name
      id=gb_record.id

      for feature in gb_record.features:
         type=feature.type
         if (type ==  "CDS"):
            locus_tag=feature.qualifiers['locus_tag'][0]
            location=feature.location
            strand=location.strand #1 for the forward (plus) strand, -1 for the reverse (negative) strand, 0 for stranded but strand unknown
            if(strand == 1): strand='+'
            if(strand == -1): strand='-'
            start=location.start
            end=location.end
            gene=""
            if('gene' in feature.qualifiers): gene=feature.qualifiers['gene'][0]
            product=""
            if('product' in feature.qualifiers): product=feature.qualifiers['product'][0].replace(" ","_")
            out_file.write(str(accession) + "\t" + str(locus_tag) + "\t" + str(gene) + "\t" + str(start) + "\t" + str(end) + "\t" + str(strand) + "\t" + str(product) + "\n")
 
   out_file.close()


if __name__ == '__main__':
   main(sys.argv[1:]) 
