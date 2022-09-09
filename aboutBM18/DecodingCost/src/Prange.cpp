/*
 * Prange.cpp
 *
 *  Created on: 10 août 2022
 *      Author: Kévin Carrier
 */

#include "Prange.h"

double Prange::complexity() const{
	double Psucc;
	if (t > (1.0 - R)/2.0){
		Psucc = 0.0;
	} else {
		Psucc = BINOM(1.0 - R, t) - fmin(1.0 - R , H(t));
	}

	if (verbose){
		cout << *this << endl;

		cout << "Success probability = " << Psucc << endl;
		cout << "Complexity = " << (- Psucc) << endl;
	}

	return -Psucc;
}

double Prange::optimization(){
	// nothing to do
	return complexity();
}

ostream& operator<<(ostream& os, const Prange& c) {
	os.setf(ios::fixed);
	os.precision(12);
	os << "Prange parameters:" << endl;
    os << "R=" << c.R << ";" << endl;
    os << "t=" << c.t << ";" << endl;
    return os;
}
