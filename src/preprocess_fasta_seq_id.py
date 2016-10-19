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
# preprocess_fasta_seq_id.py: renames sequence IDs in a fasta file.                                                 #
# Author: Elise Larsonneur                                                                                          #
#####################################################################################################################


from Bio import SeqIO
import sys
import os.path
import getopt
import datetime


def usage():
   print 'Usage: ' + str(sys.argv[0]) + ' input.fasta outputfile\n' + \
   '\t\tFasta filename must start by <two_first_letters_of_genus><two_first_letters_of_species>.<genomeNumberId> : \n' + \
   '\t\texample: limo.102.*.fasta\n' + \
   'Renames sequence IDs in a fasta file.'


def main(argv):

   fasta_file = ''      # -i 
   output_file = ''     # -o

   if (len(sys.argv) <= 1):
      usage()
      sys.exit(2)

   try:
      opts, args = getopt.getopt(argv,"hi:o:")
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
         fasta_file = arg
         if not os.path.exists(fasta_file):
             print '[ERROR] file \'' + str(fasta_file) + '\' does not exist'
             sys.exit(1)
      elif (opt == "-o"):
         output_file = arg
         if not os.path.exists(os.path.dirname(output_file)):
             print '[ERROR] directory of output file \'' + str(output_file) + '\' does not exist'
             sys.exit(1)
   
   if ("-i" in options) and ("-o" in options):
      pass
   else:
      usage()
      sys.exit(2)

   basename_fasta=os.path.basename(fasta_file)

   prefix=basename_fasta.split('.')[:2]
   prefix=('.'.join(prefix)).lower() 
   #print "PREFIX",prefix   #example : limo.102

   out_file_handler = open(output_file,"w")

   now = datetime.datetime.now()
   day=now.strftime("%d")
   month=(now.strftime("%B"))[:3]
   year=now.strftime("%Y")
   datestr= day + "-" + month + "-" + year

   #nb_records
   nb_records=0
   for seq_record in SeqIO.parse(fasta_file, "fasta"):
      nb_records = nb_records + 1

   i=0
   ind=""
   digit_nb = len(str(nb_records))
   if (digit_nb < 3) :  digit_nb = 3  ## minimum number of digits = 3

   for seq_record in SeqIO.parse(fasta_file, "fasta"):
      i=i+1
      ind=str(i).zfill(digit_nb)  # add some zero as a prefix of id i  to have 'digit_nb' digits
      id=seq_record.id
      id=prefix + ".c" + str(ind)   
      header=(">" + id + " " + datestr + " " + str(len(seq_record))) + " bp" #Example of header: >LIMO.102.C01_00 27-Aug-2015 4774 bp
      sequence=str(seq_record.seq).upper()
      freq=60
      #breaks sequence every <freq> characters  
      sequence=('\n'.join(sequence[i:i+freq] for i in xrange(0, len(sequence), freq)))
      out_file_handler.write(header + "\n" + sequence  + "\n" )

   out_file_handler.close()


if __name__ == '__main__':
    main(sys.argv[1:])
