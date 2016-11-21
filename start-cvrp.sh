#!/bin/tcsh
#g++ -o SimpleGA SimpleGA.cpp -Wall -O3
g++ cvrp.cpp SimpleGA.cpp Chromosome.cpp -O3 -o cvrp
./cvrp
java -jar CVRPValidater.jar
