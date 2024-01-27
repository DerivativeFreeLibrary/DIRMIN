#!/bin/sh

echo `date` >> ris.tex

EXEC=./dirmin
FILES='fort.2'

echo "         INTEGER, PARAMETER :: nn=2  " > parametri.fi
make clean
make FUNCT=Schubert
echo " Schubert  " > nomefun
echo "                        -186.7309D0         " > PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
rm -f fort.*

echo $NUMPROC
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=2  " > parametri.fi
make clean
make FUNCT=Schubertpen1
echo " Schub. pen. 1 " > nomefun
echo "                        -186.7309D0         " > PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=2  " > parametri.fi
make clean
make FUNCT=Schubertpen2
echo " Schub. pen. 2  " > nomefun
echo "                        -186.7309D0         " > PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=2  " > parametri.fi
make clean
make FUNCT=fun6humpcb
echo " S-H. Camel B.  " > nomefun
echo "                           -1.0316        " > PROBDAT.d
echo "       -2.5                 2.5          " >> PROBDAT.d
echo "       -1.5                 1.5          " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=2  " > parametri.fi
make clean
make FUNCT=fungoldprice
echo " Goldstein-Price  " > nomefun
echo "                            3.0          " > PROBDAT.d
echo "       -2.0                 2.0          " >> PROBDAT.d
echo "       -2.0                 2.0          " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=2  " > parametri.fi
make clean
make FUNCT=Treccanimod
echo " Treccani mod. " > nomefun
echo "                            0.0D0           " > PROBDAT.d
echo "      -25.D-1              25.D-1           " >> PROBDAT.d
echo "      -15.D-1              15.D-1           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=2  " > parametri.fi
make clean
make FUNCT=Quartic
echo " Quartic  " > nomefun
echo "                           -0.352            " > PROBDAT.d
echo "       -1.D+1               1.D+1           " >> PROBDAT.d
echo "       -1.D+1               1.D+1           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=4  " > parametri.fi
make clean
make FUNCT=funshekel5
echo " Shekel \$m=5\$  " > nomefun
echo "                          -10.1532D0       " > PROBDAT.d
echo "        0.0                10.0           " >> PROBDAT.d
echo "        0.0                10.0           " >> PROBDAT.d
echo "        0.0                10.0           " >> PROBDAT.d
echo "        0.0                10.0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=4  " > parametri.fi
make clean
make FUNCT=funshekel7
echo " Shekel \$m=7\$  " > nomefun
echo "                          -10.4029       " > PROBDAT.d
echo "        0.0                10.0           " >> PROBDAT.d
echo "        0.0                10.0           " >> PROBDAT.d
echo "        0.0                10.0           " >> PROBDAT.d
echo "        0.0                10.0           " >> PROBDAT.d
rm -f fort.*
echo $NUMPROC
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=4  " > parametri.fi
make clean
make FUNCT=funshekel10
echo " Shekel \$m=10\$  " > nomefun
echo "                          -10.5364d0       " > PROBDAT.d    
echo "        0.0                10.0           " >> PROBDAT.d
echo "        0.0                10.0           " >> PROBDAT.d
echo "        0.0                10.0           " >> PROBDAT.d
echo "        0.0                10.0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=2  " > parametri.fi
make clean
make FUNCT=funexpmod
echo " Espon. mod. " > nomefun
echo "                           -1.0d0         " > PROBDAT.d    
echo "       -1.0                 1.0           " >> PROBDAT.d
echo "       -1.0                 1.0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=4  " > parametri.fi
make clean
make FUNCT=funexpmod
echo " Espon. mod. " > nomefun
echo "                           -1.0d0         " > PROBDAT.d    
echo "       -1.0                 1.0           " >> PROBDAT.d
echo "       -1.0                 1.0           " >> PROBDAT.d
echo "       -1.0                 1.0           " >> PROBDAT.d
echo "       -1.0                 1.0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=2  " > parametri.fi
make FUNCT=funcosmixmod
echo " Cos-mix mod. " > nomefun
echo "                           -0.2d0         " > PROBDAT.d    
echo "       -1.0                 1.0           " >> PROBDAT.d
echo "       -1.0                 1.0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=4  " > parametri.fi
make clean
make FUNCT=funcosmixmod
echo " Cos-mix mod. " > nomefun
echo "                           -0.4d0         " > PROBDAT.d    
echo "       -1.0                 1.0           " >> PROBDAT.d
echo "       -1.0                 1.0           " >> PROBDAT.d
echo "       -1.0                 1.0           " >> PROBDAT.d
echo "       -1.0                 1.0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=3  " > parametri.fi
make clean
make FUNCT=funhart3
echo " Hartman  " > nomefun
echo "                           -3.8627d0      " > PROBDAT.d 
echo "        0.0                 1.0           " >> PROBDAT.d
echo "        0.0                 1.0           " >> PROBDAT.d
echo "        0.0                 1.0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=6  " > parametri.fi
make clean
make FUNCT=funhart6
echo " Hartman  " > nomefun
echo "                           -3.3223d0       " > PROBDAT.d 
echo "        0.0                 1.0           " >> PROBDAT.d
echo "        0.0                 1.0           " >> PROBDAT.d
echo "        0.0                 1.0           " >> PROBDAT.d
echo "        0.0                 1.0           " >> PROBDAT.d
echo "        0.0                 1.0           " >> PROBDAT.d
echo "        0.0                 1.0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=2  " > parametri.fi
make clean
make FUNCT=fun5nmin
echo " \$5^n\$ loc-min  " > nomefun
echo "                            0.0d0        " > PROBDAT.d 
echo "      -10.0                10.0          " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=5  " > parametri.fi
make clean
make FUNCT=fun5nmin
echo " \$5^n\$ loc-min  " > nomefun
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=10  " > parametri.fi
make clean
make FUNCT=fun5nmin
echo " \$5^n\$ loc-min  " > nomefun
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=20  " > parametri.fi
make clean
make FUNCT=fun5nmin
echo " \$5^n\$ loc-min  " > nomefun
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=2  " > parametri.fi
make clean
make FUNCT=fun10nmin
echo " \$10^n\$ loc-min  " > nomefun
echo "                            0.0d0          " > PROBDAT.d 
echo "      -10.0D0              10.0D0          " >> PROBDAT.d
echo "      -10.0D0              10.0D0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=5  " > parametri.fi
make clean
make FUNCT=fun10nmin
echo " \$10^n$ loc-min  " > nomefun
echo "      -10.0D0              10.0D0           " >> PROBDAT.d
echo "      -10.0D0              10.0D0           " >> PROBDAT.d
echo "      -10.0D0              10.0D0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=10  " > parametri.fi
make clean
make FUNCT=fun10nmin
echo " \$10^n\$ loc-min  " > nomefun
echo "      -10.0D0              10.0D0           " >> PROBDAT.d
echo "      -10.0D0              10.0D0           " >> PROBDAT.d
echo "      -10.0D0              10.0D0           " >> PROBDAT.d
echo "      -10.0D0              10.0D0           " >> PROBDAT.d
echo "      -10.0D0              10.0D0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=20  " > parametri.fi
make clean
make FUNCT=fun10nmin
echo " \$10^n\$ loc-min  " > nomefun
echo "      -10.0D0              10.0D0           " >> PROBDAT.d
echo "      -10.0D0              10.0D0           " >> PROBDAT.d
echo "      -10.0D0              10.0D0           " >> PROBDAT.d
echo "      -10.0D0              10.0D0           " >> PROBDAT.d
echo "      -10.0D0              10.0D0           " >> PROBDAT.d
echo "      -10.0D0              10.0D0           " >> PROBDAT.d
echo "      -10.0D0              10.0D0           " >> PROBDAT.d
echo "      -10.0D0              10.0D0           " >> PROBDAT.d
echo "      -10.0D0              10.0D0           " >> PROBDAT.d
echo "      -10.0D0              10.0D0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=2  " > parametri.fi
make clean
make FUNCT=fun15nmin
echo " \$15^n\$ loc-min  " > nomefun
echo "                            0.0d0       " > PROBDAT.d 
echo "      -10.0                10.0         " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=5  " > parametri.fi
make clean
make FUNCT=fun15nmin
echo " \$15^n\$ loc-min  " > nomefun
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=10  " > parametri.fi
make clean
make FUNCT=fun15nmin
echo " \$15^n\$ loc-min  " > nomefun
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=20  " > parametri.fi
make clean
make FUNCT=fun15nmin
echo " \$15^n\$ loc-min  " > nomefun
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
echo "      -10.0                10.0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=2  " > parametri.fi
make clean
make FUNCT=Griewankmod
echo " Griewank mod. " > nomefun
echo "                            0.0d0         " > PROBDAT.d 
echo "       -6.D+2              6.D+2          " >> PROBDAT.d
echo "       -6.D+2              6.D+2           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=5  " > parametri.fi
make clean
make FUNCT=Griewankmod
echo " Griewank mod. " > nomefun
echo "       -6.D+2              6.D+2          " >> PROBDAT.d
echo "       -6.D+2              6.D+2           " >> PROBDAT.d
echo "       -6.D+2              6.D+2          " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=10  " > parametri.fi
make clean
make FUNCT=Griewankmod
echo " Griewank mod. " > nomefun
echo "       -6.D+2              6.D+2          " >> PROBDAT.d
echo "       -6.D+2              6.D+2           " >> PROBDAT.d
echo "       -6.D+2              6.D+2          " >> PROBDAT.d
echo "       -6.D+2              6.D+2           " >> PROBDAT.d
echo "       -6.D+2              6.D+2          " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=20  " > parametri.fi
make clean
make FUNCT=Griewankmod
echo " Griewank mod. " > nomefun
echo "       -6.D+2              6.D+2          " >> PROBDAT.d
echo "       -6.D+2              6.D+2           " >> PROBDAT.d
echo "       -6.D+2              6.D+2          " >> PROBDAT.d
echo "       -6.D+2              6.D+2           " >> PROBDAT.d
echo "       -6.D+2              6.D+2          " >> PROBDAT.d
echo "       -6.D+2              6.D+2          " >> PROBDAT.d
echo "       -6.D+2              6.D+2           " >> PROBDAT.d
echo "       -6.D+2              6.D+2          " >> PROBDAT.d
echo "       -6.D+2              6.D+2           " >> PROBDAT.d
echo "       -6.D+2              6.D+2          " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=2  " > parametri.fi
make clean
make FUNCT=pinter1
echo " Pinter  " > nomefun
echo "                            0.0d0          " > PROBDAT.d 
echo "       -5.D+0              5.D+0           " >> PROBDAT.d
echo "       -5.D+0              5.D+0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=5  " > parametri.fi
make clean
make FUNCT=pinter1
echo " Pinter  " > nomefun
echo "       -5.D+0              5.D+0           " >> PROBDAT.d
echo "       -5.D+0              5.D+0           " >> PROBDAT.d
echo "       -5.D+0              5.D+0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=10  " > parametri.fi
make clean
make FUNCT=pinter1
echo " Pinter  " > nomefun
echo "                            0.0d0          " > PROBDAT.d 
echo "       -5.D+0              5.D+0           " >> PROBDAT.d
echo "       -5.D+0              5.D+0           " >> PROBDAT.d
echo "       -5.D+0              5.D+0           " >> PROBDAT.d
echo "       -5.D+0              5.D+0           " >> PROBDAT.d
echo "       -5.D+0              5.D+0           " >> PROBDAT.d
echo "       -5.D+0              5.D+0           " >> PROBDAT.d
echo "       -5.D+0              5.D+0           " >> PROBDAT.d
echo "       -5.D+0              5.D+0           " >> PROBDAT.d
echo "       -5.D+0              5.D+0           " >> PROBDAT.d
echo "       -5.D+0              5.D+0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=20  " > parametri.fi
make clean
make FUNCT=pinter1
echo " Pinter  " > nomefun
echo "       -5.D+0              5.D+0           " >> PROBDAT.d
echo "       -5.D+0              5.D+0           " >> PROBDAT.d
echo "       -5.D+0              5.D+0           " >> PROBDAT.d
echo "       -5.D+0              5.D+0           " >> PROBDAT.d
echo "       -5.D+0              5.D+0           " >> PROBDAT.d
echo "       -5.D+0              5.D+0           " >> PROBDAT.d
echo "       -5.D+0              5.D+0           " >> PROBDAT.d
echo "       -5.D+0              5.D+0           " >> PROBDAT.d
echo "       -5.D+0              5.D+0           " >> PROBDAT.d
echo "       -5.D+0              5.D+0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done

echo "         INTEGER, PARAMETER :: nn=2  " > parametri.fi
make clean
make FUNCT=ShiftedGriewank
echo " Griewrot2 " > nomefun
echo "                         -180.0d0          " > PROBDAT.d 
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done

echo "         INTEGER, PARAMETER :: nn=10  " > parametri.fi
make clean
make FUNCT=ShiftedGriewank
echo " Griewrot2 " > nomefun
echo "                         -180.0d0          " > PROBDAT.d 
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done

echo "         INTEGER, PARAMETER :: nn=30  " > parametri.fi
make clean
make FUNCT=ShiftedGriewank
echo " Griewrot2 " > nomefun
echo "                         -180.0d0          " > PROBDAT.d 
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done

echo "         INTEGER, PARAMETER :: nn=50  " > parametri.fi
make clean
make FUNCT=ShiftedGriewank
echo " Griewrot2 " > nomefun
echo "                         -180.0d0          " > PROBDAT.d 
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
echo "       -1.D+1              1.D+1           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=2  " > parametri.fi
make clean
make FUNCT=funAckley
echo " Ackley  " > nomefun
echo "                             0.0d0          " > PROBDAT.d 
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=10  " > parametri.fi
make clean
make FUNCT=funAckley
echo " Ackley  " > nomefun
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=30  " > parametri.fi
make clean
make FUNCT=funAckley
echo " Ackley  " > nomefun
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=50  " > parametri.fi
make clean
make FUNCT=funAckley
echo " Ackley  " > nomefun
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
echo "       -15.D+0             30.D+0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=2  " > parametri.fi
make clean
make FUNCT=Dixon_Price
echo " Dixon Price  " > nomefun
echo "                             0.0d0          " > PROBDAT.d 
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=10  " > parametri.fi
make clean
make FUNCT=Dixon_Price
echo " Dixon Price  " > nomefun
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=25  " > parametri.fi
make clean
make FUNCT=Dixon_Price
echo " Dixon Price  " > nomefun
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=50  " > parametri.fi
make clean
make FUNCT=Dixon_Price
echo " Dixon Price  " > nomefun
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
echo "       -10.D+0             10.D+0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=2  " > parametri.fi
make clean
make FUNCT=Easom
echo " Easom  " > nomefun
echo "                             -1.0d0          " > PROBDAT.d 
echo "       -100.D+0             100.D+0           " >> PROBDAT.d
echo "       -100.D+0             100.D+0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=2  " > parametri.fi
make clean
make FUNCT=Michalewics
echo " Michalewics  " > nomefun
echo "                            -1.8013d0          " > PROBDAT.d 
echo "         0.D+0             3.1415923D0           " >> PROBDAT.d
echo "         0.D+0             3.1415923D0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=5  " > parametri.fi
make clean
make FUNCT=Michalewics
echo " Michalewics  " > nomefun
echo "                            -4.687658d0          " > PROBDAT.d 
echo "         0.D+0             3.1415923D0           " >> PROBDAT.d
echo "         0.D+0             3.1415923D0           " >> PROBDAT.d
echo "         0.D+0             3.1415923D0           " >> PROBDAT.d
echo "         0.D+0             3.1415923D0           " >> PROBDAT.d
echo "         0.D+0             3.1415923D0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=10  " > parametri.fi
make clean
make FUNCT=Michalewics
echo " Michalewics  " > nomefun
echo "                            -9.66015d0          " > PROBDAT.d 
echo "         0.D+0             3.1415923D0           " >> PROBDAT.d
echo "         0.D+0             3.1415923D0           " >> PROBDAT.d
echo "         0.D+0             3.1415923D0           " >> PROBDAT.d
echo "         0.D+0             3.1415923D0           " >> PROBDAT.d
echo "         0.D+0             3.1415923D0           " >> PROBDAT.d
echo "         0.D+0             3.1415923D0           " >> PROBDAT.d
echo "         0.D+0             3.1415923D0           " >> PROBDAT.d
echo "         0.D+0             3.1415923D0           " >> PROBDAT.d
echo "         0.D+0             3.1415923D0           " >> PROBDAT.d
echo "         0.D+0             3.1415923D0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=2  " > parametri.fi
make clean
make FUNCT=funRastrigin
echo " Rastrigin  " > nomefun
echo "                                   0.d0          " > PROBDAT.d 
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=10  " > parametri.fi
make clean
make FUNCT=funRastrigin
echo " Rastrigin  " > nomefun
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=30  " > parametri.fi
make clean
make FUNCT=funRastrigin
echo " Rastrigin  " > nomefun
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=50  " > parametri.fi
make clean
make FUNCT=funRastrigin
echo " Rastrigin  " > nomefun
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
echo "         -5.12D+0             5.12D0           " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=2   " > parametri.fi
make clean
make FUNCT=beale
echo " Beale  " > nomefun
echo "                                   0.0d0         " >  PROBDAT.d 
echo "              -4.5D+0             4.5D+0         " >> PROBDAT.d
echo "              -4.5D+0             4.5D+0         " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=2   " > parametri.fi
make clean
make FUNCT=funBohachevsky1
echo " Bohachevsky 1 " > nomefun
echo "                             0.0d0         " > PROBDAT.d 
echo "       -1.0D+2             1.0D+2          " >> PROBDAT.d
echo "       -1.0D+2             1.0D+2          " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=2   " > parametri.fi
make clean
make FUNCT=funBohachevsky2
echo " Bohachevsky 2 " > nomefun
echo "                             0.0d0         " > PROBDAT.d 
echo "       -1.0D+2             1.0D+2          " >> PROBDAT.d
echo "       -1.0D+2             1.0D+2          " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=2   " > parametri.fi
make clean
make FUNCT=funBohachevsky3
echo " Bohachevsky 3 " > nomefun
echo "                             0.0d0         " > PROBDAT.d 
echo "       -1.0D+2             1.0D+2          " >> PROBDAT.d
echo "       -1.0D+2             1.0D+2          " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done

echo "         INTEGER, PARAMETER :: nn=2   " > parametri.fi
make clean
make FUNCT=funBooth
echo " Booth  " > nomefun
echo "                             0.0d0         " > PROBDAT.d 
echo "       -1.0D+1             1.0D+1          " >> PROBDAT.d
echo "       -1.0D+1             1.0D+1          " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=4   " > parametri.fi
make clean
make FUNCT=funColville
echo " Colville  " > nomefun
echo "                             0.0d0         " > PROBDAT.d 
echo "       -1.0D+1             1.0D+1          " >> PROBDAT.d
echo "       -1.0D+1             1.0D+1          " >> PROBDAT.d
echo "       -1.0D+1             1.0D+1          " >> PROBDAT.d
echo "       -1.0D+1             1.0D+1          " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=2   " > parametri.fi
make clean
make FUNCT=perm1
echo " perm1  " > nomefun
echo "                             0.0d0        "  > PROBDAT.d 
echo "       -2.50D0             2.5D0          " >> PROBDAT.d
echo "       -2.50D0             2.5D0          " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=5   " > parametri.fi
make clean
make FUNCT=perm1
echo " perm1  " > nomefun
echo "                             0.0d0      "    > PROBDAT.d 
echo "       -5.50D0             6.D0         "  >> PROBDAT.d
echo "       -5.50D0             6.D0         "  >> PROBDAT.d
echo "       -5.50D0             6.D0         "  >> PROBDAT.d
echo "       -5.50D0             6.D0         "  >> PROBDAT.d
echo "       -5.50D0             6.D0         "  >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=2   " > parametri.fi
make clean
make FUNCT=perm2
echo " perm2  " > nomefun
echo "                             0.0d0      "    > PROBDAT.d 
echo "       -2.0D0             2.0D0         "  >> PROBDAT.d
echo "       -2.0D0             2.0D0         "  >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=5   " > parametri.fi
make clean
make FUNCT=perm2
echo " perm2  " > nomefun
echo "                             0.0d0      "    > PROBDAT.d 
echo "       -5.0D0             5.0D0         "  >> PROBDAT.d
echo "       -5.0D0             5.0D0         "  >> PROBDAT.d
echo "       -5.0D0             5.0D0         "  >> PROBDAT.d
echo "       -5.0D0             5.0D0         "  >> PROBDAT.d
echo "       -5.0D0             5.0D0         "  >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=4   " > parametri.fi
make clean
make FUNCT=powell
echo " powell  " > nomefun
echo "                             0.0d0       "   > PROBDAT.d 
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=8   " > parametri.fi
make clean
make FUNCT=powell
echo " powell  " > nomefun
echo "                             0.0d0       "   > PROBDAT.d 
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=16  " > parametri.fi
make clean
make FUNCT=powell
echo " powell  " > nomefun
echo "                             0.0d0         " > PROBDAT.d 
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=24  " > parametri.fi
make clean
make FUNCT=powell
echo " powell  " > nomefun
echo "                             0.0d0         " > PROBDAT.d 
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
echo "       -4.0D0             5.0D0          " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=4  " > parametri.fi
make clean
make FUNCT=powersum
echo " powersum  " > nomefun
echo "                             0.0d0         " > PROBDAT.d 
echo "       0.0D0             4.0D0          " >> PROBDAT.d
echo "       0.0D0             4.0D0          " >> PROBDAT.d
echo "       0.0D0             4.0D0          " >> PROBDAT.d
echo "       0.0D0             4.0D0          " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=2  " > parametri.fi
make clean
make FUNCT=schwefel
echo " schwefel  " > nomefun
echo "                            0.0d0         " > PROBDAT.d 
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=5  " > parametri.fi
make clean
make FUNCT=schwefel
echo " schwefel  " > nomefun
echo "                            0.0d0         " > PROBDAT.d 
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done


echo "         INTEGER, PARAMETER :: nn=10  " > parametri.fi
make clean
make FUNCT=schwefel
echo " schwefel  " > nomefun
echo "                            0.0d0         " > PROBDAT.d 
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done



echo "         INTEGER, PARAMETER :: nn=20  " > parametri.fi
make clean
make FUNCT=schwefel
echo " schwefel  " > nomefun
echo "                            0.0d0         " > PROBDAT.d 
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
echo "       -5.D+2              5.D+2          " >> PROBDAT.d
rm -f fort.*
$EXEC
for file in $FILES
do
	cat $file >> ris.tex
done

make clean
