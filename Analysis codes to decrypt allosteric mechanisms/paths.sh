#!/bin/bash

file=paths.txt

if [ -f paths_$1.txt ]
then
rm paths_$1.txt
fi

if [ -f paths_init.txt ]
then
rm paths_init.txt
fi


npair=`wc $file | awk '{print $1}'`
for ((j=1; j<= $npair; j++)) ; do
st=`awk '{if (NR == k) print $1}' k=$j $file` ; echo $st 
en=`awk '{if (NR == k) print $2}' k=$j $file` ; echo $en 


sed -e "s/k_paths_wanted 10/k_paths_wanted $1/g" ksp.c > ksp_0.c
sed -e "s/int start_node=X/int start_node=$st/g" ksp_0.c > ksp_1.c
sed -e "s/int destination_node=X/int destination_node=$en/g" ksp_1.c > ksp_"$st"_"$en".c 
rm ksp_1.c ksp_0.c

gcc ksp_"$st"_"$en".c -lm -o ksp_"$st"_"$en".x

echo "$3 paths for $st $en">> paths_$1.txt 
echo "                    ">> paths_$1.txt 
./ksp_"$st"_"$en".x >> paths_$1.txt 
echo "                    ">> paths_$1.txt 

rm ksp_"$st"_"$en".*

#awk '{if ($1 == k) print $0 }' k=$st ../initialpath.txt  | awk '{ if (NR == e) print $0}' e=$en | awk '{for (i = 1; i <= NF; i=i+1) $i=$i+1} {print $0}' >> paths_init.txt
done




