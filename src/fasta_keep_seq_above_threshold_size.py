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
# fasta_keep_seq_above_threshold_size.py: removes sequences of length lower than the threshold from a fasta file.   #   
# Author: Sylvain Brisse                                                                                            #
#####################################################################################################################

 
from string import *
import sys

if len(sys.argv) == 4:
	filename = sys.argv[1]
	length = sys.argv[2]
	fileout = sys.argv[3]
else:
	print 'Usage: ' + str(sys.argv[0]) + ' input.fasta threshold_length output.fasta\n' + \
        'Removes sequences of length lower than the given threshold from a fasta file.'
	sys.exit(2)

def fasta_read1(fh_in):
    """
    Returns information of one sequence (sequence id, sequence comments and sequence) from the input fasta file handler. 
    """
    seq_holder = {}
    flag = 'header'
    name = ''
    comm = ''
    seq  = ''
    pos = fh_in.tell()
    line = fh_in.readline()
    while line != '':
        if line[0] == '>':
            if flag == 'header':
                p = find(line, ' ')
                name = line [1:p]
                comm = line [p+1:]
                flag = 'seq'
            else:
                fh_in.seek(pos, 0)
                break
        else:
            seq = seq +line
        pos = fh_in.tell()
        line = fh_in.readline()

    if seq != '':
        seq_holder['comm'] = comm
        seq_holder['seq'] = seq
        if name == '':
            name = 'anonymous'
        seq_holder['name'] = name
    return seq_holder


#_______________________
if __name__ == '__main__':
	
	ifh = open(filename,"r")
	out_file_handler = open(fileout,"w")

	seq_holder = fasta_read1(ifh)
	while seq_holder:
		name = '>'+seq_holder['name']+'\n'
		seq = replace(seq_holder['seq'],'\n','') + '\n'
	
		if len(seq) > int(length):
			out_file_handler.write(name)
			out_file_handler.write(seq)
		
		seq_holder = fasta_read1(ifh)
	
	
	ifh.close()
	out_file_handler.close()
		
