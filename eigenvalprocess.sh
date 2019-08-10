#/usr/bin/bash

#reading the q-points from output of QE eigenvalues

name=$1
awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' $1 > rearranged.tmp
grep 'q =' rearranged.tmp > listqpoints.tmp
readarray qpoints < listqpoints.tmp

# splitting the eigenvectors in individual files

numq=$(wc -l < listqpoints.tmp) 

for (( c=1; c<$numq; c++))
do
   eval sed -n \'/${qpoints[$c-1]}/\,/${qpoints[$c]}/p\' rearranged.tmp > q-point-$c.tmp 
   awk '{print $2,$3,$4,$5,$6,$7}' q-point-$c.tmp > q--point-$c.tmp    
   head -n -3 q--point-$c.tmp > qpointla-$c.tmp
   awk 'NR > 2' qpointla-$c.tmp > q-pointlala-$c.tmp
done

for (( c=1; c<$numq; c++))
do
   sed '/the dynamical matrix/q' q-pointlala-$c.tmp > q-point---$c.tmp
   head -n -1 q-point---$c.tmp > q-point-$c.dat
done

rm -r *.tmp





  

    





    
