#include "linear_algebra.hpp"
#include "RLPN_solver_statistics.hpp"
#include <iostream>
#include "tools.hpp"
#include <cmath>
#include <chrono>
#include<cstdlib>
#include <stdio.h>
#include <string.h>
using namespace std;

int main(int argc, char **argv){
	int nbIt = atoi(argv[1]);
	int nbCodes = atoi(argv[2]);
	informationAverage avg1;
	informationAverage avg2;
	make_stat(nbIt, nbCodes,avg1,avg2);



	
}