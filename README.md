# AKPE
ALGORITHM:
AKPE -> AKPEnum
AKPE+ -> AKPEnum*
AKPE+nocolor -> A*-nocolor


EXAMPLE:
compute maximal antagonistic k-plex in test.txt, set t=4 k=1
g++ AKPE.cpp -o AKPE -O3
./AKPE test.txt 4 1  
g++ AKPE+.cpp -o AKPE+ -O3
./AKPE+ test.txt 4 1
g++ AKPE+nocolor.cpp -o AKPE+nocolor -O3
./AKPE+nocolor test.txt 4 1
