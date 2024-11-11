#!/bin/bash

echo "***** S T A R T *****"

read -p "< selet case > Case 1 => 1, Case 2 => 2 : " case

read -p "< create mesh >(num of divisions in x direction) (num of y in ~) : " x y

if [ $case = 1 ]; then
    cd input1
    python3 mesher.py $x $y
elif [ $case = 2 ]; then
    cd input2
    python3 mesher.py $x $y
fi

echo "press any key to compile the program..."
read

cd ../
make clean
make

echo "press any key to start the program..."
read

if [ $case = 1 ]; then
    cd input1
elif [ $case = 2 ]; then
    cd input2
fi

../bin/laplace_eq
