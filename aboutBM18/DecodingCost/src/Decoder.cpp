/*
 * ISD.cpp
 *
 *  Created on: 10 août 2022
 *      Author: Kévin Carrier
 */

#include "Decoder.h"

Decoder::Decoder(double R, double t) : R(R), t((t > 0.5) ? 1.0 - t : t), large(t > 0.5){}

double Decoder::get_R() const {
	return R;
}

void Decoder::set_R(double R) {
	this->R = R;
}

double Decoder::get_t() const {
	if (large){
		return 1.0 - t;
	} else {
		return t;
	}
}

void Decoder::set_t(double t) {
	this->t = (t > 0.5) ? 1.0 - t : t;
	large = (t > 0.5);
}

void Decoder::set_verbose(bool v){
	this->verbose = v;
}

ostream& operator<<(ostream& os, const Decoder& c) {
	os.setf(ios::fixed);
	os.precision(12);
	os << "Parameters:" << endl;
    os << "R=" << c.R << ";" << endl;
    os << "t=" << c.t << ";" << endl;
    return os;
}
