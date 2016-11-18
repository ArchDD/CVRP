#!/bin/tcsh
#g++ -o SimpleGA SimpleGA.cpp -Wall -O3
g++ cvrp.cpp SimpleGA.cpp -O3 -o cvrp
./cvrp
