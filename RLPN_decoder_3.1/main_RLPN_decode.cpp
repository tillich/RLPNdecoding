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
 void decode_codes(size_t nbIt,size_t nbCodes,size_t u,size_t t,long int treshold){
        size_t stop[4];
        for(size_t i = 0; i < 4; ++i){
            stop[i] = 0;
        }
       for(size_t i = 0; i < nbCodes; ++i){
		   auto problem = make_problem(t);
           size_t j = problem.decode_code(nbIt, u,treshold);
           stop[j] += 1;
       }
       cout << "Number of decoding : " << nbCodes << endl;
       cout << "Number of iteration for each decoding : " << nbIt << endl;
	   cout << "-----" << endl;
       cout << "Percentage of successfull decoding : " << 100*((double )stop[0])/((double)nbCodes) << "%"<<endl;
       cout << "Percentage of failed decoding : " << 100*((double )(stop[1] + stop[2] + stop[3]))/((double)nbCodes) << "%"<< endl;
        cout << "    -> Percentage of false positive : " << 100*((double )stop[1])/((double)nbCodes) << "%" <<endl;
        cout << "    -> Percentage of false negative : " << 100*((double )stop[2])/((double)nbCodes) << "%" <<endl;
        cout << "    -> Percentage of successive bad bets : " << 100*((double )stop[3])/((double)nbCodes) << "%" <<endl;
		
	}

int main(int argc, char **argv){
	int t = atoi(argv[1]);
	int u = atoi(argv[2]);
	long int treshold = atol(argv[3]);
	int nbIt = atoi(argv[4]);
	int nbCodes = atoi(argv[5]);
	decode_codes(nbIt,nbCodes, u, t, treshold);
}