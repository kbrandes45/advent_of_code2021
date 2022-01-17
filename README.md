# Advent of Code 2021

https://adventofcode.com/2021

Solved parts 1 and 2 for days 1-17 completely, except for part 2 of day 14 because my solution was too slow and I didn't want to optimize it that day. 

Written in C++ and compiled on Ubuntu 18.04 Linux machine. It can be compiled by running:
```
g++ -std=c++11 -O1 solutions.cpp -o aoc
```

And then it can be executed by running:
```
./aoc
```

Parts 1 and 2 for each day are all contained in a function named 'dayXX' (e.x. `day10()` function), which will print out the outputs (sometimes part 1 is commented/unincluded and will not run without editing the function body). To change which day is run, edit the `main` function at the bottom of the solutions.cpp file.

NOTE: Each day's data file is saved as "dayXXdata.txt" and must be in the same directory as the executable for it to work. 
