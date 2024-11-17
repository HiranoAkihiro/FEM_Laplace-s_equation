#!/bin/bash

echo "***** S T A R T *****"
echo "<select analysis>"
read -p "Enter 'laplace1', 'laplace2', or, 'poisson':" case

echo "< create mesh >"
read -p "Enter '(num of divisions in x direction)' and '(num of y in y direction)' : " x y

if [ $case = laplace1 ]; then
    cd input1
    python3 mesher.py $x $y
elif [ $case = laplace2 ]; then
    cd input2
    python3 mesher.py $x $y
elif [ $case = poisson ]; then
    cd input_poisson
    python3 mesher.py $x $y
fi

echo "press any key to compile the program..."
read

cd ../
if [ $case = poisson ]; then
    make -f Makefile.poisson clean
    make -f Makefile.poisson
else
    make clean
    make
fi

echo "press any key to start the program..."
read

if [ $case = laplace1 ]; then
    cd input1
    ../bin/laplace_eq
elif [ $case = laplace2 ]; then
    cd input2
    ../bin/laplace_eq
elif [ $case = poisson ]; then
    cd input_poisson
    ../bin/poisson_eq
fi