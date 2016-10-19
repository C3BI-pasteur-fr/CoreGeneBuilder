/*
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
#  FASTA2AGP: creating agp file from FASTA-formatted scaffold sequence file.                                        #
#  Author: Alexis Criscuolo                                                                                         #
#####################################################################################################################
*/

//  [Version 1.1]                                                                                                    

import java.io.*;
import java.util.*;

public class FASTA2AGP {

    // constants
    static String NOFILE = "N.O.F.I.L.E";

    // options
    static File infile;  // -i
    static File ctgFile; // -o
    static File agpFile; // -a
    static int polyNlgt; // -n

    // io
    static BufferedReader in;
    static BufferedWriter outc, outa;

    // data
    static int size;
    static ArrayList<String> fh, sq;
    
    // stuffs
    static int o, i, s, sN, eN, l, cptctg, part_number, base;
    static String line, hdr, polyN;
    static StringBuilder sb;

    public static void main(String[] args) throws IOException {

	// ############################
	// ### man              #######
	// ############################
	if ( args.length < 2 ) {
	    System.out.println("USAGE: FASTA2AGP -i <scaffolds.fasta> [-o <contigs.fasta>] [-a <info.agp>]");
	    System.out.println("where options are:");
	    System.out.println(" -i <infile>   FASTA-formatted scaffold sequence file (mandatory)");
	    System.out.println(" -o <outfile>  output FASTA contig sequence file name (default: <infile>.fna)");
	    System.out.println(" -a <outfile>  output AGP file name (default: <infile>.agp)");
	    System.out.println(" -n <integer>  minimum length of scaffolding stretch of Ns (default: 10)");
	    System.exit(0);
	}

	// ############################
	// ### parsing options  #######
	// ############################
	infile = new File(NOFILE); 
	ctgFile = new File(NOFILE);
	agpFile = new File(NOFILE);
	polyNlgt = 10;
	o = -1;
	while ( ++o < args.length ) {
	    if ( args[o].equals("-i") ) { infile = new File(args[++o]); continue; }
	    if ( args[o].equals("-o") ) { ctgFile = new File(args[++o]); continue; }
	    if ( args[o].equals("-a") ) { agpFile = new File(args[++o]); continue; }
	    if ( args[o].equals("-n") ) {
		try { polyNlgt = Integer.parseInt(args[++o]); }
		catch ( NumberFormatException e ) { System.out.println("incorrect integer value: " + args[o] + " (option -n)"); System.exit(1); }
		if ( polyNlgt < 1 ) { System.out.println("incorrect integer value: " + polyNlgt + " (option -n)"); System.exit(1); }
		continue;
	    }
	}
	if ( ! infile.exists() ) { System.out.println(infile.toString() + " not found (options -i)"); System.exit(1); }
	if ( ctgFile.toString().equals(NOFILE) ) ctgFile = new File(infile.toString() + ".fna");	
	if ( agpFile.toString().equals(NOFILE) ) agpFile = new File(infile.toString() + ".agp");
	sb = new StringBuilder(""); s = polyNlgt; while ( --s >= 0 ) sb = sb.append('N');
	polyN = sb.toString();	

	// ############################
	// ### reading infile  ########
	// ############################
	fh = new ArrayList<String>(); sq = new ArrayList<String>(); sb = new StringBuilder("");
	in = new BufferedReader(new FileReader(infile));
	while ( true ) {
	    try { line = in.readLine().trim(); } catch ( NullPointerException e ) { in.close(); break; }
	    if ( line.startsWith(">") ) { if ( sb.length() != 0 ) { sq.add(sb.toString()); sb = new StringBuilder(""); } fh.add(line); continue; }
	    sb = sb.append(line);
	}
	if ( sb.length() != 0 ) sq.add(sb.toString());
	size = fh.size();
	o = 10; while ( o < size ) o *= 10; base = ("" + (o *= 10)).length();

	// ##############################
	// ### writing outfiles  ########
	// ##############################
	outc = new BufferedWriter(new FileWriter(ctgFile)); outa = new BufferedWriter(new FileWriter(agpFile)); 
	cptctg = 0; i = -1;
	while ( ++i < size ) {
	    line = sq.get(i).toUpperCase(); l = line.length(); s = 0; sN = line.indexOf(polyN);
	    if ( s == -1 ) {
		hdr = "contig_" + frmt(++cptctg,base); outc.write(">" + hdr); outc.newLine(); outc.write(line); outc.newLine();
		// =======================================================================================================================================================
		//          object                           object_beg  object_end  part_number  component_type  component_id   component_beg  component_end  orientation
		//          |                                |           |           |            |               |              |              |              |
		outa.write("scaffold_" + frmt(i+1,base) + "\t1\t" +      l +      "\t1\t" +      "W\t" +          hdr +       "\t1\t" +         l +         "\t+"); 
		outa.newLine();
		continue;
	    }
	    part_number = 0;
	    while ( sN != -1 ) {
		eN = sN; while ( (++eN < l) && (line.charAt(eN) == 'N') ) {}
		// ========================================================
		//        s                  sN               eN
		//        |                  |                |
		//  ...NNNAACTGTCACTACGAATGCTNNNNNNNNNNNNNNNNNACGTACGT
		hdr = "contig_" + frmt(++cptctg,base); outc.write(">" + hdr); outc.newLine(); outc.write(line.substring(s, sN)); outc.newLine();
		// =====================================================================================================================================================================
		//          object                               object_beg     object_end  part_number          component_type  component_id  component_beg  component_end  orientation
		//          |                                    |              |           |                    |               |             |              |              |
		outa.write("scaffold_" + frmt(i+1,base) + "\t" + (s+1) + "\t" + sN + "\t" + (++part_number) + "\tW\t" +          hdr +      "\t1\t" +         (sN-s) +    "\t+"); outa.newLine();
		// ============================================================================================================================================================
		//          object                               object_beg      object_end  part_number          component_type gap_length   gap_type  linkage   link.evidence
		//          |                                    |               |           |                    |              |            |         |         |
		outa.write("scaffold_" + frmt(i+1,base) + "\t" + (sN+1) + "\t" + eN + "\t" + (++part_number) + "\tN\t" +         (eN-sN) + "\tscaffold\tyes\t" + "paired-ends"); outa.newLine();
		s = eN; sN = line.indexOf(polyN, s);
	    }
	    // =====================================
	    //        s                  
	    //        |                  
	    //  ...NNNAACTGTCACTACGAATGCTACGTACGT
	    hdr = "contig_" + frmt(++cptctg,base); outc.write(">" + hdr); outc.newLine(); outc.write(line.substring(s)); outc.newLine();
	    // ======================================================================================================================================================================
	    //          object                               object_beg     object_end  part_number          component_type  component_id   component_beg  component_end  orientation
	    //          |                                    |              |           |                    |               |              |              |              |
	    outa.write("scaffold_" + frmt(i+1,base) + "\t" + (s+1) + "\t" + l + "\t" +  (++part_number) + "\tW\t" +          hdr +       "\t1\t" +         (l-s) +     "\t+"); outa.newLine();
	}
	outc.close(); outa.close();


    }

    static String frmt( int x, int base ) {
	return Integer.toString((int)(Math.pow(10, base) + x)).substring(1);
    }

}

