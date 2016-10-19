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
# get_syntenic_genes.awk: computes gene neiborhood and filter out genes for which the syntenic threshold is not     #
#                         reached.                                                                                  #
# Authors: Marie Touchon, Eduardo P. C. Rocha                                                                       #
#####################################################################################################################


#it searches in among the memory last ones how many are closer than the threshold
BEGIN{
	if (! memory)
		memory = 5
	if (! threshold)
		threshold = 10
}

(NR <= memory){
	counter = 0 
	current = $NF
	for (i=1;i< NR; i++){
		dist = current - last[i] 
		if (dist <= threshold && dist >= -threshold)
			counter ++
#		print dist, " op", current, " op", last[i]
		}
	print $0, counter 
	last[NR] = $NF
	next
}

{
	counter = 0 
	current = $NF
	for (i=1;i<= memory; i++){
		dist = current - last[i] 
		if (dist <= threshold && dist >= -threshold)
			counter ++
		}

	for (i=1;i< memory; i++)
		last[i] = last[i+1]

	last[memory] = $NF

	print $0, counter 
} 
