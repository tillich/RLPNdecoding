/*
 * SternMayOzerov.cpp
 *
 *  Created on: 10 août 2022
 *      Author: Kévin Carrier
 */

#include "SternMayOzerov.h"



double SternMayOzerov::get_p() const {
	return p;
}

void SternMayOzerov::set_p(double p) {
	this->p = p;
}

double SternMayOzerov::complexity() const{
	double lmda, alpha, Psucc;
	if ( 1.0 - R - BINOM(1.0 - R, t - p) > BINOM(R, p) ){
		lmda = BINOM(R/2.0, p/2.0);
		alpha = NNS(lmda/(1.0-R), (t - p)/(1.0-R));
		Psucc = BINOM(R, p) + BINOM(1.0 - R, t - p) - fmin(1.0 - R, H(t));
	} else {
		lmda = (1.0 - R - BINOM(1.0 - R, t - p))/2.0;
		alpha = NNS(lmda/(1.0-R), (t - p)/(1.0-R));
		Psucc = 0.0;
	}
	if (verbose){
		cout << *this << endl;
		cout << "lmda = " << lmda << endl;
		cout << "alpha = " << lmda << endl;
		cout << "Psucc = " << lmda << endl;
		cout << "Complexity = " << alpha*lmda - Psucc << endl;
	}
	return alpha*lmda - Psucc;
}

double SternMayOzerov::optimization(){

	auto f = [this](double p){
		this->set_p(p);
		return this->complexity();
	};

	double *opt = golden_section(f, fmax(0.0, t + R - 1.0), fmin(R, t), SMO_PRECISION);
	double res = opt[0];
	free(opt);
	return res;
}

ostream& operator<<(ostream& os, const SternMayOzerov& c) {
	os.setf(ios::fixed);
	os.precision(12);
	os << "Stern-May-Ozerov parameters:" << endl;
    os << "R=" << c.R << ";" << endl;
    os << "t=" << c.t << ";" << endl;
    os << "p=" << c.p << ";" << endl;
    return os;
}
