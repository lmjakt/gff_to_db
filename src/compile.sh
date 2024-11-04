#!/bin/bash
## g++ -Wall -std=c++0x -o main gff_db.cpp main.cpp
## note that this gives c++11; which is required for
## some of the set iteration used.
g++ -Wall -std=gnu++0x -o main gff_db.cpp main.cpp
