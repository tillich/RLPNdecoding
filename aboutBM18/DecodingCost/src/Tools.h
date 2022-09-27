/*
 * Tools.h
 *
 *  Created on: 10 août 2022
 *      Author: Kévin Carrier
 */
#ifndef TOOLS_H_
#define TOOLS_H_

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>


#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795028841971693993751058209749445923078164
#endif /* M_PI */

#ifndef M_INF
#define M_INF numeric_limits<double>::infinity()
#endif /* M_INF */

#define NNS_PRECISION 0.0000000001
#define ENTROPY_PRECISION 0.0000000001

#define XLOGX(x)	( (x) > 0 ? (x)*log2(x) : (x) < -ENTROPY_PRECISION ? INFINITY : 0)
#define H(x)	( - XLOGX(x) - XLOGX(1 - (x)) )
#define BINOM(n, k)	( (n)*H((k)/(n)) )

using namespace std;

namespace tools{

/*
 *  @brief Golden section technique for minimizing a function
 *  @param f a function to minimize
 *  @param params some parameters
 *  @param start the lower bound of the parameter to optimize
 *  @param end the upper bound of the parameter to optimize
 *  @return {f(x), x} where f(x) is the minimum of the function f
 */
template <typename T>
double *golden_section(T& f, double start, double end, double precision){
	double golden_ratio = 2.61803398874989490252573887119; // 1 + (1 + sqrt(5))/2
	unsigned int i, min_i, max_i;
	double param_tab[4] = {start, start + (end - start)/golden_ratio , end - (end - start)/golden_ratio , end};
	double cost_tab[4];

	for (i = 0 ; i < 4 ; ++i){
		cost_tab[i] = f(param_tab[i]);
	}

	min_i = 0;
	max_i = 0;
	for(i = 1 ; i < 4 ; ++i){
		if(cost_tab[min_i] > cost_tab[i]){
			min_i = i;
		}
		if(cost_tab[max_i] < cost_tab[i]){
			max_i = i;
		}
	}
	while (param_tab[3] - param_tab[0] > precision){
		switch (min_i){
		case 0:
			cost_tab[3] = cost_tab[1];
			param_tab[3] = param_tab[1];

			param_tab[2] = param_tab[3] - (param_tab[3] - param_tab[0]) / golden_ratio;
			cost_tab[2] = f(param_tab[2]);

			param_tab[1] = param_tab[0] + (param_tab[3] - param_tab[0]) / golden_ratio;
			cost_tab[1] = f(param_tab[1]);
			break;
		case 1:
			cost_tab[3] = cost_tab[2];
			param_tab[3] = param_tab[2];

			cost_tab[2] = cost_tab[1];
			param_tab[2] = param_tab[1];

			param_tab[1] = param_tab[0] + (param_tab[3] - param_tab[0]) / golden_ratio;
			cost_tab[1] = f(param_tab[1]);
			break;
		case 2:
			cost_tab[0] = cost_tab[1];
			param_tab[0] = param_tab[1];

			cost_tab[1] = cost_tab[2];
			param_tab[1] = param_tab[2];

			param_tab[2] = param_tab[3] - (param_tab[3] - param_tab[0]) / golden_ratio;
			cost_tab[2] = f(param_tab[2]);
			break;
		case 3:
			cost_tab[0] = cost_tab[2];
			param_tab[0] = param_tab[2];

			param_tab[1] = param_tab[0] + (param_tab[3] - param_tab[0]) / golden_ratio;
			cost_tab[1] = f(param_tab[1]);

			param_tab[2] = param_tab[3] - (param_tab[3] - param_tab[0]) / golden_ratio;
			cost_tab[2] = f(param_tab[2]);
			break;
		}

		for(i= 1 ; i < 4 ; ++i){
			if(cost_tab[min_i] > cost_tab[i]){
				min_i = i;
			}
			if(cost_tab[max_i] < cost_tab[i]){
				max_i = i;
			}
		}
	}
	double *res = (double *)malloc(2*sizeof(double));
	res[0] = cost_tab[min_i];
	res[1] = param_tab[min_i];
	f(res[1]); // small extra cost for having the good other optimal parameters
	return res;
}



/*
 *  @brief Print a table containing 200 points of the function f
 *  @param f a function to plot
 *  @param params some parameters
 *  @param start the lower bound of the variable parameter
 *  @param end the upper bound of the variable parameter
 *  @return a table containing {x, f(x)} where x is in [start, end]
 */
template <typename T>
void print_values(T& f, double start, double end){
	double x = start;
	double step = (end - start)/100.0;

	fprintf(stdout, "L = [");
	for ( ; x < end ; x += step){
		fprintf(stdout, "[%.10f, %.10f],\n", x, f(x));
	}
	fprintf(stdout, "[%.10f, %.10f]]\n", end, f(end));
}


/*
 *  @brief Newton method of order 1 for inverting the Shannon entropy function
 *  @param y a double in [0,1]
 *  @return x in [0,0.5] such that H(x)=y
 */
double Hinv(double y);

/*
 *  @brief Near-collisions search complexity
 *  @param lmda the size of the lists $L := 2^{lmda n}$
 *  @param w the relative distance to achieve
 *  @return a double alpha in [1,2] such that the complexity for finding all the near-collisions at distance $w n$ in a list of size $L$ is of order $L^{alpha}$
 */
double NNS(double lmda, double w);

/*
 * @brief Number of ww-representation of a binary vector of weight w
 * @param s size of the binary vectors
 * @param w weight of the destination vector
 * @param ww weight of the representation
 * @return Number of pairs (x,y) such that |x|=|y|=ww and x+y=z where |z|=w
 */
double representations(double s, double w, double ww);
}
#endif /* TOOLS_H_ */
