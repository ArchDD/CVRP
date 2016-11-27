#!/bin/tcsh
#module add intel-cluster-studio/compiler/64/13.1/117
#icpc cvrp.cpp SimpleGA.cpp Chromosome.cpp -Ofast -o cvrp
g++ cvrp.cpp SimpleGA.cpp Chromosome.cpp -O3 -o cvrp -march=native -mtune=native -fopenmp
./cvrp
java -jar CVRPValidater.jar
