#!/bin/bash
echo
echo "***** S T A R T *****"
echo
echo "<select analysis>"
read -p "Enter 'laplace1', 'laplace2', 'poisson', 'square', or, 'hole' :" case

echo "< create mesh >"
if [ $case != hole ]; then
    read -p "Enter num of element divisions (ex. '5 5' ) : " x y
fi

if [ $case = laplace1 ]; then
    cd input1
    python3 mesher.py $x $y
elif [ $case = laplace2 ]; then
    cd input2
    python3 mesher.py $x $y
elif [ $case = poisson ]; then
    cd input_poisson
    python3 mesher.py $x $y
elif [ $case = square ]; then
    cd input_square
    python3 mesher.py $x $y
elif [ $case = hole ]; then
    cd input_spherehole
fi

echo "press any key to compile the program..."
read

cd ../
if [ $case = poisson ]; then
    make -f M_poisson.mk clean
    make -f M_poisson.mk
elif [ $case = square ]; then
    make -f M_elastic.mk clean
    make -f M_elastic.mk
elif [ $case = hole ]; then
    make -f M_elastic.mk clean
    make -f M_elastic.mk
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
elif [ $case = square ]; then
    cd input_square
    ../bin/elastic
elif [ $case = hole ]; then
    cd input_spherehole
    ../bin/elastic
fi