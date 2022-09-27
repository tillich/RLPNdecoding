/*
 * Tools.cpp
 *
 *  Created on: 10 août 2022
 *      Author: Kévin Carrier
 */

#include "Tools.h"

namespace tools{

double Hinv(double y){
	if ( (y < -ENTROPY_PRECISION) || (y > 1 + ENTROPY_PRECISION) ){
		cerr << "Hinv(" << y << ") cannot be computed !" << endl;
		exit(EXIT_FAILURE);
	}
	if ( y <= 0 ){
		return 0.0;
	}
	if ( y >= 1.0 ){
		return 0.5;
	}
	double x_pred = 0.25;
	while (H(x_pred) - x_pred * log2((1.0 - x_pred)/x_pred) > y){
		x_pred = x_pred / 2.0;
	}
	double x = (y - H(x_pred))/log2((1.0 - x_pred) / x_pred) + x_pred;

	while (fabs(x_pred - x) > ENTROPY_PRECISION){
		x_pred = x;
		x = (y - H(x_pred))/log2((1.0 - x_pred) / x_pred) + x_pred;
	}
	return x;
}

double NNS(double lmda, double w){
    if (lmda <= NNS_PRECISION)
        return 1.0;
    if (lmda > 1)
        return 2.0;
    double ww = fmin(w, 1 - w);
    if (ww <= NNS_PRECISION)
    	return 1.0;

    double dGV = Hinv(1 - lmda);
    if ((1 - sqrt(1 - 2*ww))/2 < dGV){
        return (1 - ww)*( 1 - H( (dGV - ww/2)/(1 - ww) ) )/lmda;
    } else {
        return 2 - (1 - H(ww))/lmda;
    }
}

double representations(double s, double w, double ww){
	return BINOM(s - w, ww - w/2.0) + w;
}
}
