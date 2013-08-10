#!/bin/sh

BinDir="../../../../bin/LiDIA/sparc8-sun-solaris2.6/c++"

DatDir="appl_data"
DatPrefix="data"

OutDir="new_appl_data"
OutPrefix="out"

programs="complex_periods_appl minimal_model_appl elliptic_curve_bigint_appl point_bigint_appl"
programs2="curve_isomorphism_appl"
data="1 2 3 4 5 6 7" 


echo "This test tool executes the application files"
echo "$programs with the input data appl_data/data*"
echo "and writes the results to new_appl_data/"
echo "You may use diff to check the correctness of the output."
echo " "


for p in $programs; do

 for d in $data; do
 
 InFile=$DatDir/$DatPrefix$d
 OutFile=$OutDir/$p.$OutPrefix$d

 Prog=$BinDir/$p

 echo "Testing $Prog with $InFile"

 $Prog < $InFile > $OutFile

 done 

 echo " "
done


for p in $programs2; do

 OutFile=$OutDir/$p.$OutPrefix

 Prog=$BinDir/$p

 echo "Testing $Prog"

 $Prog > $OutFile

 echo " "
done


echo DONE.
