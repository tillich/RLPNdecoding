/*
 * MayOzerov.cpp
 *
 *  Created on: 27 août 2022
 *      Author: kcarrier
 */

#include "MayOzerov.h"

MayOzerov::MayOzerov(double R, double t) : Decoder(R, t){
	// initialisation with Prange parameters
	for (int i = 0 ; i < 3 ; ++i){
		this->ell[i] = 0.0;
	}
	for (int i = 0 ; i < 4 ; ++i){
		this->p[i] = 0.0;
	}
}

MayOzerov::MayOzerov(double R, double t, double *ell, double *p) : Decoder(R,t){

	for (int i = 1 ; i < 3 ; ++i){
		this->ell[i] = ell[i];
	}
	for (int i = 0 ; i < 4 ; ++i){
		this->p[i] = p[i];
	}
}

MayOzerov::MayOzerov(const MayOzerov &mo) : Decoder(mo.R, mo.t){
	for (int i = 1 ; i < 3 ; ++i){
		this->ell[i] = mo.ell[i];
	}
	for (int i = 0 ; i < 4 ; ++i){
		this->p[i] = mo.p[i];
	}
}

void MayOzerov::operator=(const MayOzerov &mo) {
	this->R = mo.R;
	this->t = mo.t;
	for (int i = 1 ; i < 3 ; ++i){
		this->ell[i] = mo.ell[i];
	}
	for (int i = 0 ; i < 4 ; ++i){
		this->p[i] = mo.p[i];
	}
}


double MayOzerov::get_ell(int i) const {
	return ell[i];
}

double MayOzerov::get_p(int i) const {
	return p[i];
}

void MayOzerov::set_params(double *ell, double *p) {
	for (int i = 1 ; i < 3 ; ++i){
		this->ell[i] = ell[i];
	}
	for (int i = 0 ; i < 4 ; ++i){
		this->p[i] = p[i];
	}
}


void MayOzerov::set_optim_param(double precision) {
	this->precision = precision;
}


ostream& operator<<(ostream& os, const MayOzerov& c) {
	os.setf(ios::fixed);
	os.precision(12);
	os << "May-Ozerov parameters:" << endl;
	os << "R=" << c.R << ";" << endl;
	os << "t=" << c.t << ";" << endl;
	os << "ell[1]=" << c.ell[1] << ";" << endl;
	os << "ell[2]=" << c.ell[2] << ";" << endl;
	os << "p[0]=" << c.p[0] << ";" << endl;
	os << "p[1]=" << c.p[1] << ";" << endl;
	os << "p[2]=" << c.p[2] << ";" << endl;
	os << "p[3]=" << c.p[3] << ";" << endl;

	return os;
}

double MayOzerov::complexity() const{
	// List sizes
	double S[4];
	S[0] = BINOM((R + ell[1] + ell[2])/2.0 , p[0]);
	S[1] = BINOM(R + ell[1] + ell[2], p[1]) - ell[1];
	S[2] = BINOM(R + ell[1] + ell[2], p[2]) - ell[1] - ell[2];
	S[3] = BINOM(R + ell[1] + ell[2], p[3]) + BINOM(1.0 - R - ell[1] - ell[2], t - p[3]) - 1.0 + R;

	// Time complexity
	double T[4];
	T[0] = S[0];
	T[1] = fmax(S[0], 2*S[0] - ell[1]);
	T[2] = fmax(S[1], 2*S[1] - ell[2]);
	T[3] = S[2] * NNS(S[2]/(1.0 - R - ell[1] - ell[2]), (t - p[3])/(1.0 - R - ell[1] - ell[2]));


	// Success probability
	double Psucc = BINOM(R + ell[1] + ell[2], p[3]) + BINOM(1.0 - R - ell[1] - ell[2], t - p[3]) - H(t);


	double C = fmax(fmax(T[1], T[2]), T[3]);

	if (verbose){
		cout << *this << endl;

		cout << "List sizes:" << endl;
		for (int i = 0 ; i < 4 ; ++i){
			cout << "S("<< i << ") = " << S[i] << endl;
		}

		cout << "Time for producing the lists:" << endl;
		for (int i = 0 ; i < 4 ; ++i){
			cout << "T("<< i << ") = " << T[i] << endl;
		}

		cout << "Correctness Lemma:" << endl;
		double Correctness;
		Correctness = ell[1] - representations(R + ell[1] + ell[2], p[2], p[1]);
		cout << "Correctness 1: ell1 - nb_repr(R + ell, p2, p1) = " << Correctness << ( (Correctness >= 0) ? " is valid" : " is not valid" ) << endl;
		Correctness = ell[1] + ell[2] - representations(R + ell[1] + ell[2], p[3], p[2]);
		cout << "Correctness 2: ell - nb_repr(R + ell, p3, p2) = " << Correctness << ( (Correctness >= 0) ? " is valid" : " is not valid" ) << endl;


		cout << "Success probability = " << Psucc << endl;
		cout << "Complexity = " << (C - Psucc) << endl;
	}
	return C - Psucc;
}


double MayOzerov::optimization(){
	auto optim_ell = [this](double sum_ell){
		auto optim_p3 = [this, sum_ell](double p3){
			p[3] = p3;

			// computing p[2]
			double inf = p[3]/2;
			double sup = (R + sum_ell)/2.0;
			double cur = (sup + inf)/2.0;
			while (sup - inf > precision){
				if (sum_ell - representations(R + sum_ell, p[3], cur) > 0) {
					inf = cur;
				} else {
					sup = cur;
				}
				cur = (sup + inf)/2.0;
			}
			p[2] = inf;

			if (p[2] > sum_ell){
				return M_INF;
			}

			auto optim_ell1 = [this, sum_ell](double ell1){
				ell[1] = ell1;
				ell[2] = sum_ell - ell[1];

				// computing p[1]
				double inf = p[2]/2;
				double sup = (R + sum_ell)/2.0;
				double cur = (sup + inf)/2.0;
				while (sup - inf > precision){
					if (ell[1] - representations(R + sum_ell, p[2], cur) > 0) {
						inf = cur;
					} else {
						sup = cur;
					}
					cur = (sup + inf)/2.0;
				}
				p[1] = inf;

				// computing p[0]
				p[0] = p[1]/2.0;

				return complexity();
			};
			double *opt = golden_section(optim_ell1, p[2], sum_ell, precision);
			double res = opt[0];
			free(opt);
			return res;
		};
		double *opt = golden_section(optim_p3,0, fmin(sum_ell, t), precision);
		double res = opt[0];
		free(opt);
		return res;
	};
	double *opt = golden_section(optim_ell,0, 1-R, precision);
	double res = opt[0];
	free(opt);
	return res;
}
