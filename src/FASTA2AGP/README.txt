FASTA2AGP: creating agp file from FASTA-formatted scaffold sequence file

1. COMPILATION
To obtain a executable jar, launch the following command line in the src/
directory:
  ./JarMaker.sh FASTA2AGP.java
To obtain a JAVA class, launch the following command line in the src/
directory:
  javac FASTA2AGP.java
To obtain a native machine code, launch the following command line in
the src/ directory:
  gcj -O3 --main=FASTA2AGP FASTA2AGP.java -o fasta2agp

2. EXECUTION
To execute the executable jar, launch the following command line:
  java -jar FASTA2AGP.jar <options>
To execute the JAVA class, launch the following command line:
  java FASTA2AGP <option>
To execute the native machine code, launch the following command line:
  ./fasta2agp <options>
A short documentation is available when executed without option.
