#!/bin/tcsh
g++ cvrp.cpp SimpleGA.cpp Chromosome.cpp -O3 -o cvrp -march=native -mtune=native -fopenmp
./cvrp
java -jar CVRPValidater.jar
