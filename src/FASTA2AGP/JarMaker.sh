#!/bin/bash
#$ -S /bin/bash 

# USAGE : ./JarMaker.sh Program.java

p=$1
c=${p%.*};
javac $p ;
if [ ! -e $c.class ]; then exit; fi;

mkdir JARBUILDER/
cp *.class JARBUILDER/ ;

mkdir JARBUILDER/META-INF/ ;
echo "Main-Class: $c" > JARBUILDER/META-INF/MANIFEST.MF ;
echo "Class-Path: " >> JARBUILDER/META-INF/MANIFEST.MF ;

cd JARBUILDER/ ;
jar cvfm $c.jar ./META-INF/MANIFEST.MF . ;
cd ../ ;

cp JARBUILDER/$c.jar . ;

rm -r JARBUILDER/ ; 
rm *.class ;
