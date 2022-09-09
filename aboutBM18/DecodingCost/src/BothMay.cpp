/*
 * BothMay.cpp
 *
 *  Created on: 10 août 2022
 *      Author: Kévin Carrier
 */

#include "BothMay.h"

BothMay::BothMay(double R, double t) : Decoder(R, t){
	/*
	// initialisation with Prange parameters
	this->p[0] = 0.0;
	for (int i = 1 ; i < 5 ; ++i){
		this->r[i] = 0.0;
		this->p[i] = 0.0;
		for (int j = 1 ; j <= i ; ++j){
			this->ww[j][i] = 0.0;
		}
	}
	this->r[4] = 1.0 - R;
	this->ww[4][4] = t;
	*/

	// initialisation with medium parameters
	this->r[1]=(1.0 - R)/4.0;
	this->r[2]=(1.0 - R)/4.0;
	this->r[3]=(1.0 - R)/4.0;
	this->r[4]=(1.0 - R)/4.0;
	this->p[4]=t*R;
	this->p[3]=this->p[4]/2.0;
	this->p[2]=this->p[3]/2.0;
	this->p[1]=this->p[2]/2.0;
	this->p[0]=this->p[1]/2.0;
	this->ww[1][4]=t*this->r[1];
	this->ww[1][3]=this->ww[1][4]/2.0;
	this->ww[1][2]=this->ww[1][3]/2.0;
	this->ww[1][1]=this->ww[1][2]/2.0;
	this->ww[2][4]=t*this->r[2];
	this->ww[2][3]=this->ww[2][4]/2.0;
	this->ww[2][2]=this->ww[2][3]/2.0;
	this->ww[3][4]=t*this->r[3];
	this->ww[3][3]=this->ww[3][4]/2.0;
	this->ww[4][4]=t*this->r[4];

}

BothMay::BothMay(double R, double t, double r[5], double p[5], double ww[5][5]) : Decoder(R,t){
	this->p[0] = p[0];
	for (int i = 1 ; i < 5 ; ++i){
		this->r[i] = r[i];
		this->p[i] = p[i];
		for (int j = 1 ; j <= i ; ++j){
			this->ww[j][i] = ww[j][i];
		}
	}
}

BothMay::BothMay(const BothMay &bm) : Decoder(bm.R, bm.t){
	p[0] = bm.p[0];
	for (int i = 1 ; i < 5 ; ++i){
		r[i] = bm.r[i];
		p[i] = bm.p[i];
		for (int j = 1 ; j <= i ; ++j){
			ww[j][i] = bm.ww[j][i];
		}
	}
}

void BothMay::operator=(const BothMay &bm) {
	R = bm.R;
	t = bm.t;
	p[0] = bm.p[0];
	for (int i = 1 ; i < 5 ; ++i){
		r[i] = bm.r[i];
		p[i] = bm.p[i];
		for (int j = 1 ; j <= i ; ++j){
			ww[j][i] = bm.ww[j][i];
		}
	}
}


double BothMay::get_r(int i) const {
	return r[i];
}

double BothMay::get_p(int i) const {
	return p[i];
}

double BothMay::get_ww(int j, int i) const {
	return ww[j][i];
}

void BothMay::set_params(double *r, double *p, double **ww) {
	this->p[0] = p[0];
	for (int i = 1 ; i < 5 ; ++i){
		this->r[i] = r[i];
		this->p[i] = p[i];
		for (int j = 1 ; j <= i ; ++j){
			this->ww[j][i] = ww[j][i];
		}
	}
}

void BothMay::set_optim_params(double radius, double step, double limit, double precision){
	this->radius = radius;
	this->step = step;
	this->limit = limit;
	this->precision = precision;
}


ostream& operator<<(ostream& os, const BothMay& c) {
	os.setf(ios::fixed);
	os.precision(12);
	os << "Both-May parameters:" << endl;
    os << "R=" << c.R << ";" << endl;
    os << "t=" << c.t << ";" << endl;
    os << "r[1]=" << c.r[1] << ";" << endl;
    os << "r[2]=" << c.r[2] << ";" << endl;
    os << "r[3]=" << c.r[3] << ";" << endl;
	os << "r[4]=" << c.r[4] << ";" << endl;
	os << "p[0]=" << c.p[0] << ";" << endl;
	os << "p[1]=" << c.p[1] << ";" << endl;
	os << "p[2]=" << c.p[2] << ";" << endl;
	os << "p[3]=" << c.p[3] << ";" << endl;
    os << "p[4]=" << c.p[4] << ";" << endl;
    os << "ww[1][1]=" << c.ww[1][1] << ";" << endl;
    os << "ww[1][2]=" << c.ww[1][2] << ";" << endl;
    os << "ww[1][3]=" << c.ww[1][3] << ";" << endl;
    os << "ww[1][4]=" << c.ww[1][4] << ";" << endl;
    os << "ww[2][2]=" << c.ww[2][2] << ";" << endl;
    os << "ww[2][3]=" << c.ww[2][3] << ";" << endl;
    os << "ww[2][4]=" << c.ww[2][4] << ";" << endl;
	os << "ww[3][3]=" << c.ww[3][3] << ";" << endl;
	os << "ww[3][4]=" << c.ww[3][4] << ";" << endl;
	os << "ww[4][4]=" << c.ww[4][4] << ";" << endl;
    return os;
}


double BothMay::correctness(int l) const {
	double res = representations(R, p[l], p[l - 1]);
	for (int i = 1 ; i < l ; ++i) {
		res += (representations(r[i], ww[i][l], ww[i][l - 1]) - r[i]);
	}
	return res;
}

double BothMay::complexity() const{
	int l, i;

	// Succes probability
	double Psucc = BINOM(R, p[4]);
	for ( int i = 1 ; i < 5 ; ++i ) {
		Psucc += BINOM(r[i], ww[i][4]);
	}
	Psucc -= fmin(1.0 - R, H(t));

	// List sizes
	double L[5];
	L[0] = BINOM(R/2.0, p[0]);
	L[1] = BINOM(R, p[1]) + BINOM(r[1], ww[1][1]) - r[1];
	for (l = 2 ; l < 5 ; ++l){

		// corrected version of [BM18]
		double Lbound = BINOM(R, p[l]) + BINOM(r[l], ww[l][l]) -  r[l];
		L[l] = 2*L[l-1] + BINOM(r[l], ww[l][l]) -  r[l] + representations(R, p[l], p[l-1]) + BINOM(R, p[l]) - 2.0 * BINOM(R, p[l - 1]);
		for (i = 1 ; i < l ; ++i){
			Lbound += ( BINOM(r[i], ww[i][l]) -  r[i] );
			L[l] += ( representations(r[i], ww[i][l], ww[i][l-1]) + BINOM(r[i], ww[i][l]) - 2.0 * BINOM(r[i], ww[i][l - 1]) );
		}
		L[l] = fmin(L[l], Lbound);
	}

	// Time complexity
	double T[5];
	double C = 0;
	for ( l = 1 ; l < 5 ; ++l){
		if ( (L[l - 1] > r[l]) || (L[l - 1] < 0) || (r[l] == 0) ){
			T[l] = fmax(2*L[l - 1], 0.0);
		} else {
			T[l] = L[l - 1] * NNS(L[l - 1] / r[l], ww[l][l] / r[l]);
		}
		C = fmax(C, T[l]);
	}

	if (verbose){
		cout << *this << endl;

		cout << "List sizes:" << endl;
		for (l = 0 ; l < 5 ; ++l){
			cout << "L("<< l << ") = " << L[l] << endl;
		}

		cout << "Time for producing the lists:" << endl;
		for (l = 1 ; l < 5 ; ++l){
			cout << "T("<< l << ") = " << T[l] << endl;
		}

		cout << "Correctness Lemma:" << endl;
		for (l = 2 ; l < 5 ; ++l){
			cout << "E("<< l << ") = " << correctness(l) << endl;
		}

		cout << "Success probability = " << Psucc << endl;
		cout << "Complexity = " << (C - Psucc) << endl;
	}
	return C - Psucc;
}

/*
// very large bounds
double BothMay::optimization_tmp(){
	verbose = false;

	double exp4_tmp1, exp4_tmp2, exp4_tmp3;
	list<double> list_ww33;
	double exp3_tmp1, exp3_tmp2;
	list<double> list_ww22;
	double exp2_tmp1;
	list<double> list_ww11;
	double sup, inf, cur;
	double C_cur, C_min = M_INF;
	BothMay params_C_min(*this);

	BothMay pinf(*this);
	BothMay psup(*this);

	for (int i = 0 ; i < 5 ; ++i) {
		pinf.r[i] -= radius;
		psup.r[i] += radius;
		pinf.p[i] -= radius;
		psup.p[i] += radius;
		for (int j = 1 ; j < i ; ++j) {
			pinf.ww[j][i] -= radius;
			psup.ww[j][i] += radius;
		}
	}


	for (p[4] = fmax(pinf.p[4], 0) ; p[4] <= fmin(psup.p[4], fmin(R, t)) ; p[4] += step){
		for (p[3] = fmax(pinf.p[3], p[4]/2.0) ; p[3] <= fmin(psup.p[3], R - p[4]/2.0) ; p[3] += step){
			exp4_tmp1 = representations(R, p[4], p[3]);
			for (p[2] = fmax(pinf.p[2], p[3]/2.0) ; p[2] <= fmin(psup.p[2], R - p[3]/2.0) ; p[2] += step){
				exp3_tmp1 = representations(R, p[3], p[2]);
				for (p[1] = fmax(pinf.p[1], p[2]/2.0) ; p[1] <= fmin(psup.p[1], R - p[2]/2.0) ; p[1] += step){
					p[0] = p[1]/2.0;
					exp2_tmp1 = representations(R, p[2], p[1]);
					for (r[1] = fmax(pinf.r[1], 0) ; r[1] <= fmin(psup.r[1], 1.0 - R) ; r[1] += step){
						for (ww[1][4] = fmax(pinf.ww[1][4], 0) ; ww[1][4] <= fmin(psup.ww[1][4], fmin(r[1], t - p[4])) ; ww[1][4] += step){
							for (ww[1][3] = fmax(pinf.ww[1][3], ww[1][4]/2.0) ; ww[1][3] <= fmin(psup.ww[1][3], r[1] - ww[1][4]/2.0) ; ww[1][3] += step){
								exp4_tmp2 = exp4_tmp1 + representations(r[1], ww[1][4], ww[1][3]) - r[1];
								for (ww[1][2] = fmax(pinf.ww[1][2], ww[1][3]/2.0) ; ww[1][2] <= fmin(psup.ww[1][2], r[1] - ww[1][3]/2.0) ; ww[1][2] += step){
									exp3_tmp2 = exp3_tmp1 + representations(r[1], ww[1][3], ww[1][2]) - r[1];
									list_ww11.clear();
									if ( (exp2_tmp1 + representations(r[1], ww[1][2], r[1]/2.0)) < r[1] ) {
										break;
									}
									if ( (exp2_tmp1 + representations(r[1], ww[1][2], ww[1][2]/2.0)) > r[1] ) {
										list_ww11.push_front(ww[1][2]/2.0);
									} else {
										sup = r[1]/2.0;
										inf = ww[1][2]/2.0;
										cur = (sup + inf)/2.0;
										while (sup - inf > precision){
											if ( (exp2_tmp1 + representations(r[1], ww[1][2], cur)) > r[1] ) {
												sup = cur;
											} else {
												inf = cur;
											}
											cur = (sup + inf)/2.0;
										}
										list_ww11.push_front(sup);
									}
									if ( (exp2_tmp1 + representations(r[1], ww[1][2], r[1] - ww[1][2]/2.0)) > r[1] ) {
										list_ww11.push_front(r[1] - ww[1][2]/2.0);
									} else {
										sup = r[1] - ww[1][2]/2.0;
										inf = r[1]/2.0;
										cur = (sup + inf)/2.0;
										while (sup - inf > precision){
											if (exp2_tmp1 + representations(r[1], ww[1][2], cur) > r[1]) {
												inf = cur;
											} else {
												sup = cur;
											}
											cur = (sup + inf)/2.0;
										}
										list_ww11.push_front(inf);
									}
									for (double ww11 : list_ww11) {
										ww[1][1] = ww11;
										for (r[2] = fmax(pinf.r[2], 0) ; r[2] <= fmin(psup.r[2], 1.0 - R - r[1]) ; r[2] += step){
											for (ww[2][4] = fmax(pinf.ww[2][4], 0) ; ww[2][4] <= fmin(psup.ww[2][4], fmin(r[2], t - p[4] - ww[1][4])) ; ww[2][4] += step){
												for (ww[2][3] = fmax(pinf.ww[2][3], ww[2][4]/2.0) ; ww[2][3] <= fmin(psup.ww[2][3], r[2] - ww[2][4]/2.0 ) ; ww[2][3] += step){
													exp4_tmp3 = exp4_tmp2 + representations(r[2], ww[2][4], ww[2][3]) - r[2];
													list_ww22.clear();
													if ( (exp3_tmp2 + representations(r[2], ww[2][3], r[2]/2.0)) < r[2] ) {
														break;
													}
													if ( (exp3_tmp2 + representations(r[2], ww[2][3], ww[2][3]/2.0)) > r[2] ) {
														list_ww22.push_front(ww[2][3]/2.0);
													} else {
														sup = r[2]/2.0;
														inf = ww[2][3]/2.0;
														cur = (sup + inf)/2.0;
														while (sup - inf > precision){
															if (exp3_tmp2 + representations(r[2], ww[2][3], cur) > r[2]) {
																sup = cur;
															} else {
																inf = cur;
															}
															cur = (sup + inf)/2.0;
														}
														list_ww22.push_front(sup);
													}
													if ( (exp3_tmp2 + representations(r[2], ww[2][3], r[2] - ww[2][3]/2.0)) > r[2] ) {
														list_ww22.push_front(r[2] - ww[2][3]/2.0);
													} else {
														sup = r[2] - ww[2][3]/2.0;
														inf = r[2]/2.0;
														cur = (sup + inf)/2.0;
														while (sup - inf > precision){
															if (exp3_tmp2 + representations(r[2], ww[2][3], cur) > r[2]) {
																inf = cur;
															} else {
																sup = cur;
															}
															cur = (sup + inf)/2.0;
														}
														list_ww22.push_front(inf);
													}
													for (double ww22 : list_ww22) {
														ww[2][2] = ww22;
														for (r[3] = fmax(pinf.r[3], 0) ; r[3] <= fmin(psup.r[3], 1.0 - R - r[1] - r[2]) ; r[3] += step){
															r[4] = 1.0 - R - r[1] - r[2] - r[3];
															for (ww[3][4] = fmax(pinf.ww[3][4], fmax(0, t - p[4] - ww[1][4] - ww[2][4] - r[4])) ; ww[3][4] <= fmin(psup.ww[3][4], fmin(r[3], t - p[4] - ww[1][4] - ww[2][4])) ; ww[3][4] += step){
																ww[4][4] = t - p[4] - ww[1][4] - ww[2][4] - ww[3][4];
																list_ww33.clear();
																if ( (exp4_tmp3 + representations(r[3], ww[3][4], r[3]/2.0)) < r[3] ) {
																	break;
																}
																if ( (exp4_tmp3 + representations(r[3], ww[3][4], ww[3][4]/2.0)) > r[3] ) {
																	list_ww33.push_front(ww[3][4]/2.0);
																} else {
																	sup = r[3]/2.0;
																	inf = ww[3][4]/2.0;
																	cur = (sup + inf)/2.0;
																	while (sup - inf > precision){
																		if (exp4_tmp3 + representations(r[3], ww[3][4], cur) > r[3]) {
																			sup = cur;
																		} else {
																			inf = cur;
																		}
																		cur = (sup + inf)/2.0;
																	}
																	list_ww33.push_front(sup);
																}
																if ( (exp4_tmp3 + representations(r[3], ww[3][4], r[3] - ww[3][4]/2.0)) > r[3] ) {
																	list_ww33.push_front(r[3] - ww[3][4]/2.0);
																} else {
																	sup = r[3] - ww[3][4]/2.0;
																	inf = r[3]/2.0;
																	cur = (sup + inf)/2.0;
																	while (sup - inf > precision){
																		if (exp4_tmp3 + representations(r[3], ww[3][4], cur) > r[3]) {
																			inf = cur;
																		} else {
																			sup = cur;
																		}
																		cur = (sup + inf)/2.0;
																	}
																	list_ww33.push_front(inf);
																}
																for (double ww33 : list_ww33) {
																	ww[3][3] = ww33;
																	C_cur = complexity();
																	if (C_cur < C_min){
																		C_min = C_cur;
																		params_C_min = *this;
																		cerr << C_min << endl;
																	}
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	*this = params_C_min;
	return C_min;
}
*/

/*
 // reduced bounds
double BothMay::optimization_tmp(){
	verbose = false;

	double exp4_tmp1, exp4_tmp2, exp4_tmp3;
	list<double> list_ww33;
	double exp3_tmp1, exp3_tmp2;
	list<double> list_ww22;
	double exp2_tmp1;
	list<double> list_ww11;
	double sup, inf, cur;
	double C_cur, C_min = M_INF;
	BothMay params_C_min(*this);

	BothMay pinf(*this);
	BothMay psup(*this);

	for (int i = 0 ; i < 5 ; ++i) {
		pinf.r[i] -= radius;
		psup.r[i] += radius;
		pinf.p[i] -= radius;
		psup.p[i] += radius;
		for (int j = 1 ; j < i ; ++j) {
			pinf.ww[j][i] -= radius;
			psup.ww[j][i] += radius;
		}
	}


	// binary case
	for (p[4] = fmax(pinf.p[4], 0) ; p[4] <= fmin(psup.p[4], fmin(R/2.0, t)) ; p[4] += step){
		for (p[3] = fmax(pinf.p[3], p[4]/2.0) ; p[3] <= fmin(psup.p[3], p[4]) ; p[3] += step){
			exp4_tmp1 = representations(R, p[4], p[3]);
			for (p[2] = fmax(pinf.p[2], p[3]/2.0) ; p[2] <= fmin(psup.p[2], p[3]) ; p[2] += step){
				exp3_tmp1 = representations(R, p[3], p[2]);
				for (p[1] = fmax(pinf.p[1], p[2]/2.0) ; p[1] <= fmin(psup.p[1], p[2]) ; p[1] += step){
					p[0] = p[1]/2.0;
					exp2_tmp1 = representations(R, p[2], p[1]);
					for (r[1] = fmax(pinf.r[1], 0) ; r[1] <= fmin(psup.r[1], 1.0 - R) ; r[1] += step){
						for (ww[1][4] = fmax(pinf.ww[1][4], 0) ; ww[1][4] <= fmin(psup.ww[1][4], fmin(r[1]/2.0, t - p[4])) ; ww[1][4] += step){
							for (ww[1][3] = fmax(pinf.ww[1][3], ww[1][4]/2.0) ; ww[1][3] <= fmin(psup.ww[1][3], ww[1][4]) ; ww[1][3] += step){
								exp4_tmp2 = exp4_tmp1 + representations(r[1], ww[1][4], ww[1][3]) - r[1];
								for (ww[1][2] = fmax(pinf.ww[1][2], ww[1][3]/2.0) ; ww[1][2] <= fmin(psup.ww[1][2], ww[1][3]) ; ww[1][2] += step){
									exp3_tmp2 = exp3_tmp1 + representations(r[1], ww[1][3], ww[1][2]) - r[1];
									list_ww11.clear();
									if ( (exp2_tmp1 + representations(r[1], ww[1][2], r[1]/2.0)) < r[1] ) {
										break;
									}
									if ( (exp2_tmp1 + representations(r[1], ww[1][2], ww[1][2]/2.0)) > r[1] ) {
										list_ww11.push_front(ww[1][2]/2.0);
									} else {
										sup = r[1]/2.0;
										inf = ww[1][2]/2.0;
										cur = (sup + inf)/2.0;
										while (sup - inf > precision){
											if ( (exp2_tmp1 + representations(r[1], ww[1][2], cur)) > r[1] ) {
												sup = cur;
											} else {
												inf = cur;
											}
											cur = (sup + inf)/2.0;
										}
										list_ww11.push_front(sup);
									}
									for (double ww11 : list_ww11) {
										ww[1][1] = ww11;
										for (r[2] = fmax(pinf.r[2], 0) ; r[2] <= fmin(psup.r[2], 1.0 - R - r[1]) ; r[2] += step){
											for (ww[2][4] = fmax(pinf.ww[2][4], 0) ; ww[2][4] <= fmin(psup.ww[2][4], fmin(r[2]/2.0, t - p[4] - ww[1][4])) ; ww[2][4] += step){
												for (ww[2][3] = fmax(pinf.ww[2][3], ww[2][4]/2.0) ; ww[2][3] <= fmin(psup.ww[2][3], ww[2][4] ) ; ww[2][3] += step){
													exp4_tmp3 = exp4_tmp2 + representations(r[2], ww[2][4], ww[2][3]) - r[2];
													list_ww22.clear();
													if ( (exp3_tmp2 + representations(r[2], ww[2][3], r[2]/2.0)) < r[2] ) {
														break;
													}
													if ( (exp3_tmp2 + representations(r[2], ww[2][3], ww[2][3]/2.0)) > r[2] ) {
														list_ww22.push_front(ww[2][3]/2.0);
													} else {
														sup = r[2]/2.0;
														inf = ww[2][3]/2.0;
														cur = (sup + inf)/2.0;
														while (sup - inf > precision){
															if (exp3_tmp2 + representations(r[2], ww[2][3], cur) > r[2]) {
																sup = cur;
															} else {
																inf = cur;
															}
															cur = (sup + inf)/2.0;
														}
														list_ww22.push_front(sup);
													}
													for (double ww22 : list_ww22) {
														ww[2][2] = ww22;
														for (r[3] = fmax(pinf.r[3], 0) ; r[3] <= fmin(psup.r[3], 1.0 - R - r[1] - r[2]) ; r[3] += step){
															r[4] = 1.0 - R - r[1] - r[2] - r[3];
															for (ww[3][4] = fmax(pinf.ww[3][4], fmax(0, t - p[4] - ww[1][4] - ww[2][4] - r[4])) ; ww[3][4] <= fmin(psup.ww[3][4], fmin(r[3]/2.0, t - p[4] - ww[1][4] - ww[2][4])) ; ww[3][4] += step){
																ww[4][4] = t - p[4] - ww[1][4] - ww[2][4] - ww[3][4];
																list_ww33.clear();
																if ( (exp4_tmp3 + representations(r[3], ww[3][4], r[3]/2.0)) < r[3] ) {
																	break;
																}
																if ( (exp4_tmp3 + representations(r[3], ww[3][4], ww[3][4]/2.0)) > r[3] ) {
																	list_ww33.push_front(ww[3][4]/2.0);
																} else {
																	sup = r[3]/2.0;
																	inf = ww[3][4]/2.0;
																	cur = (sup + inf)/2.0;
																	while (sup - inf > precision){
																		if (exp4_tmp3 + representations(r[3], ww[3][4], cur) > r[3]) {
																			sup = cur;
																		} else {
																			inf = cur;
																		}
																		cur = (sup + inf)/2.0;
																	}
																	list_ww33.push_front(sup);
																}
																for (double ww33 : list_ww33) {
																	ww[3][3] = ww33;
																	C_cur = complexity();
																	if (C_cur < C_min){
																		C_min = C_cur;
																		params_C_min = *this;
																		cerr << C_min << endl;
																	}
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	*this = params_C_min;
	return C_min;
}
*/

/*
// very reduced bound and use of golden section technique for optimizing parameters p4, r1, r2, r3
double BothMay::optimization_tmp(){
	verbose = false;

	BothMay pinf(*this);
	BothMay psup(*this);

	for (int i = 0 ; i < 5 ; ++i) {
		pinf.r[i] -= radius;
		psup.r[i] += radius;
		pinf.p[i] -= radius;
		psup.p[i] += radius;
		for (int j = 1 ; j < i ; ++j) {
			pinf.ww[j][i] -= radius;
			psup.ww[j][i] += radius;
		}
	}

	double lbound_p4, ubound_p4;
	if (psup.p[4] < 0.0) {
		lbound_p4 = 0.0;
		ubound_p4 = 0.0;
	} else if (pinf.p[4] > fmin(R/2.0, t)) {
		lbound_p4 = fmin(R/2.0, t);
		ubound_p4 = fmin(R/2.0, t);
	} else {
		lbound_p4 = fmax(0.0, pinf.p[4]);
		ubound_p4 = fmin(fmin(R/2.0, t), psup.p[4]);
	}
	auto optim_p4 = [this, pinf, psup](double p4){
		p[4] = p4;
		double exp4_tmp1, exp3_tmp1, exp2_tmp1;

		double C_min = 20;
		BothMay params_C_min(*this);

		double lbound_p3, ubound_p3;
		if (psup.p[3] < p[4]/2.0) {
			lbound_p3 = p[4]/2.0;
			ubound_p3 = p[4]/2.0;
		} else if (pinf.p[3] > p[4]) {
			lbound_p3 = p[4];
			ubound_p3 = p[4];
		} else {
			lbound_p3 = fmax(p[4]/2.0, pinf.p[3]);
			ubound_p3 = fmin(p[4], psup.p[3]);
		}
		for (p[3] = lbound_p3 ; p[3] <= ubound_p3 ; p[3] += step){
			exp4_tmp1 = representations(R, p[4], p[3]);
			double lbound_p2, ubound_p2;
			if (psup.p[2] < p[3]/2.0) {
				lbound_p2 = p[3]/2.0;
				ubound_p2 = p[3]/2.0;
			} else if (pinf.p[2] > p[3]) {
				lbound_p2 = p[3];
				ubound_p2 = p[3];
			} else {
				lbound_p2 = fmax(p[3]/2.0, pinf.p[2]);
				ubound_p2 = fmin(p[3], psup.p[2]);
			}
			for (p[2] = lbound_p2 ; p[2] <= ubound_p2 ; p[2] += step){
				exp3_tmp1 = representations(R, p[3], p[2]);
				double lbound_p1, ubound_p1;
				if (psup.p[1] < p[2]/2.0) {
					lbound_p1 = p[2]/2.0;
					ubound_p1 = p[2]/2.0;
				} else if (pinf.p[1] > p[2]) {
					lbound_p1 = p[2];
					ubound_p1 = p[2];
				} else {
					lbound_p1 = fmax(p[2]/2.0, pinf.p[1]);
					ubound_p1 = fmin(p[2], psup.p[1]);
				}
				for (p[1] = lbound_p1 ; p[1] <= ubound_p1 ; p[1] += step){
					p[0] = p[1]/2.0;
					exp2_tmp1 = representations(R, p[2], p[1]);

					// we must have ( exp3_tmp1 + representations(r[1], t - p[4], t - p[4]) - r[1] > 0 )
					// and we must have ( exp4_tmp1 + representations(r[1], t - p[4], t - p[4]) - r[1] > 0 )
					double inf, sup, cur, lbound_r1, ubound_r1;
					if (exp3_tmp1 < exp4_tmp1){
						inf = 2*(t - p[4]);
						sup = 1.0 - R;   //  t - p4 < t < Hinv(1 - R) < (1 - R)/2
						cur = (sup + inf)/2.0;
						while (sup - inf > precision){
							if ( exp3_tmp1 + representations(cur, t - p[4], t - p[4]) - cur > 0 ) {
								inf = cur;
							} else {
								sup = cur;
							}
							cur = (sup + inf)/2.0;
						}
					} else {
						inf = 2*(t - p[4]);
						sup = 1.0 - R;   //  t - p4 < t < Hinv(1 - R) < (1 - R)/2
						cur = (sup + inf)/2.0;
						while (sup - inf > precision){
							if ( exp4_tmp1 + representations(cur, t - p[4], t - p[4]) - cur > 0 ) {
								inf = cur;
							} else {
								sup = cur;
							}
							cur = (sup + inf)/2.0;
						}
					}
					if (psup.r[1] < 0.0) {
						lbound_r1 = 0.0;
						ubound_r1 = 0.0;
					} else if (pinf.r[1] > inf) {
						lbound_r1 = inf;
						ubound_r1 = inf;
					} else {
						lbound_r1 = fmax(0.0, pinf.r[1]);
						ubound_r1 = fmin(inf, psup.r[1]);
					}

					auto optim_r1 = [this, exp4_tmp1, exp3_tmp1, exp2_tmp1, pinf, psup](double r1){

						r[1] = r1;
						double exp4_tmp2;
						double exp3_tmp2;
						double sup, inf, cur;
						double C_min = 10;
						BothMay params_C_min(*this);

						// we must have ( exp3_tmp1 + representations(r[1], ww[1][4], ww[1][4]) - r[1] > 0 )
						// and we must have ( exp4_tmp1 + representations(r[1], ww[1][4], ww[1][4]) - r[1] > 0 )
						if (exp3_tmp1 < exp4_tmp1){
							inf = 0.0;
							sup = fmin(r[1]/2.0, t - p[4]);
							cur = (sup + inf)/2.0;
							while (sup - inf > precision){
								if ( exp3_tmp1 + representations(r[1], cur, cur) - r[1] > 0 ) {
									sup = cur;
								} else {
									inf = cur;
								}
								cur = (sup + inf)/2.0;
							}
						} else {
							inf = 0.0;
							sup = fmin(r[1]/2.0, t - p[4]);
							cur = (sup + inf)/2.0;
							while (sup - inf > precision){
								if ( exp4_tmp1 + representations(r[1], cur, cur) - r[1] > 0 ) {
									sup = cur;
								} else {
									inf = cur;
								}
								cur = (sup + inf)/2.0;
							}
						}

						double lbound_w14, ubound_w14;
						if (psup.ww[1][4] < sup) {
							lbound_w14 = sup;
							ubound_w14 = sup;
						} else if (pinf.ww[1][4] > fmin(r[1]/2.0, t - p[4])) {
							lbound_w14 = fmin(r[1]/2.0, t - p[4]);
							ubound_w14 = fmin(r[1]/2.0, t - p[4]);
						} else {
							lbound_w14 = fmax(sup, pinf.ww[1][4]);
							ubound_w14 = fmin(fmin(r[1]/2.0, t - p[4]), psup.ww[1][4]);
						}
						for (ww[1][4] = lbound_w14 ; ww[1][4] <= ubound_w14 ; ww[1][4] += step){

							double lbound_w13, ubound_w13;
							// we must have ( exp3_tmp1 + representations(r[1], ww[1][3], ww[1][3]) - r[1] > 0 )
							inf = ww[1][4]/2.0;
							sup = ww[1][4];
							cur = (sup + inf)/2.0;
							while (sup - inf > precision){
								if ( exp3_tmp1 + representations(r[1], cur, cur) - r[1] > 0 ) {
									sup = cur;
								} else {
									inf = cur;
								}
								cur = (sup + inf)/2.0;
							}

							ubound_w13 = ww[1][4];
							lbound_w13 = sup;

							// we must have exp4_tmp2 > 0  and so ( exp4_tmp1 + representations(r[1], ww[1][4], ww[1][3]) - r[1] > 0 )
							inf = lbound_w13;
							sup = ubound_w13;
							cur = (sup + inf)/2.0;
							while (sup - inf > precision){
								if ( exp4_tmp1 + representations(r[1], ww[1][4], cur) - r[1] > 0 ) {
									sup = cur;
								} else {
									inf = cur;
								}
								cur = (sup + inf)/2.0;
							}

							if (psup.ww[1][3] < sup) {
								lbound_w13 = sup;
								ubound_w13 = sup;
							} else if (pinf.ww[1][3] > ww[1][4]) {
								lbound_w13 = ww[1][4];
								ubound_w13 = ww[1][4];
							} else {
								lbound_w13 = fmax(sup, pinf.ww[1][3]);
								ubound_w13 = fmin(ww[1][4], psup.ww[1][3]);
							}
							for (ww[1][3] = lbound_w13 ; ww[1][3] <= ubound_w13 ; ww[1][3] += step){
								exp4_tmp2 = exp4_tmp1 + representations(r[1], ww[1][4], ww[1][3]) - r[1];

								// we must have ( exp3_tmp1 + representations(r[1], ww[1][3], ww[1][2]) - r[1] > 0 )
								inf = ww[1][3]/2.0;
								sup = ww[1][3];
								cur = (sup + inf)/2.0;
								while (sup - inf > precision){
									if ( exp3_tmp1 + representations(r[1], ww[1][3], cur) - r[1] > 0 ) {
										sup = cur;
									} else {
										inf = cur;
									}
									cur = (sup + inf)/2.0;
								}

								double lbound_w12, ubound_w12;
								if (psup.ww[1][2] < sup) {
									lbound_w12 = sup;
									ubound_w12 = sup;
								} else if (pinf.ww[1][2] > ww[1][3]) {
									lbound_w12 = ww[1][3];
									ubound_w12 = ww[1][3];
								} else {
									lbound_w12 = fmax(sup, pinf.ww[1][2]);
									ubound_w12 = fmin(ww[1][3], psup.ww[1][2]);
								}
								for (ww[1][2] = lbound_w12 ; ww[1][2] <= ubound_w12 ; ww[1][2] += step){

									exp3_tmp2 = exp3_tmp1 + representations(r[1], ww[1][3], ww[1][2]) - r[1]; //is greater than 0.0
									//we always have ( (exp2_tmp1 + representations(r[1], ww[1][2], r[1]/2.0)) >= r[1] )

									if ( (exp2_tmp1 + representations(r[1], ww[1][2], ww[1][2]/2.0)) > r[1] ) {
										ww[1][1] = ww[1][2]/2.0;
									} else {
										sup = r[1]/2.0;
										inf = ww[1][2]/2.0;
										cur = (sup + inf)/2.0;
										while (sup - inf > precision){
											if ( (exp2_tmp1 + representations(r[1], ww[1][2], cur)) > r[1] ) {
												sup = cur;
											} else {
												inf = cur;
											}
											cur = (sup + inf)/2.0;
										}
										ww[1][1] = sup;
									}

									// we have exp4_tmp2 > 0 and so ( exp4_tmp2 + representations(r[2], fmin(r[2]/2.0 , t - p[4] - ww[1][4]), fmin(r[2]/2.0 , t - p[4] - ww[1][4])) - r[2] > 0 )
									double lbound_r2, ubound_r2;
									if (2.0*(t - p[4] - ww[1][4]) < 1.0 - R - r[1]){
										inf = 2.0*(t - p[4] - ww[1][4]);
										sup = 1.0 - R - r[1];
										cur = (sup + inf)/2.0;
										while (sup - inf > precision){
											if ( exp4_tmp2 + representations(cur, t - p[4] - ww[1][4], t - p[4] - ww[1][4]) - cur > 0 ) {
												inf = cur;
											} else {
												sup = cur;
											}
											cur = (sup + inf)/2.0;
										}
										ubound_r2 = inf;
									} else {
										ubound_r2 = 1.0 - R - r[1];
									}
									if (psup.r[2] < 0.0) {
										lbound_r2 = 0.0;
										ubound_r2 = 0.0;
									} else if (pinf.r[2] > ubound_r2) {
										lbound_r2 = ubound_r2;
									} else {
										lbound_r2 = fmax(0.0, pinf.r[2]);
										ubound_r2 = fmin(ubound_r2, psup.r[2]);
									}
									auto optim_r2 = [this, exp4_tmp2, exp3_tmp2, pinf, psup](double r2){

										r[2] = r2;
										double exp4_tmp3;
										double sup, inf, cur;
										double C_min = 10;
										BothMay params_C_min(*this);

										// we must have ( exp4_tmp2 + representations(r[2], ww[2][4], ww[2][4]) - r[2] > 0 )
										inf = 0.0;
										sup = fmin(r[2]/2.0, t - p[4] - ww[1][4]);
										cur = (sup + inf)/2.0;
										while (sup - inf > precision){
											if ( exp4_tmp2 + representations(r[2], cur, cur) - r[2] > 0 ) {
												sup = cur;
											} else {
												inf = cur;
											}
											cur = (sup + inf)/2.0;
										}

										double lbound_w24, ubound_w24;
										if (psup.ww[2][4] < sup) {
											lbound_w24 = sup;
											ubound_w24 = sup;
										} else if (pinf.ww[2][4] > fmin(r[2]/2.0, t - p[4] - ww[1][4])) {
											lbound_w24 = fmin(r[2]/2.0, t - p[4] - ww[1][4]);
											ubound_w24 = fmin(r[2]/2.0, t - p[4] - ww[1][4]);
										} else {
											lbound_w24 = fmax(sup, pinf.ww[2][4]);
											ubound_w24 = fmin(fmin(r[2]/2.0, t - p[4] - ww[1][4]), psup.ww[2][4]);
										}
										for (ww[2][4] = lbound_w24 ; ww[2][4] <= ubound_w24 ; ww[2][4] += step){

											// we must have ( exp4_tmp2 + representations(r[2], ww[2][4], ww[2][3]) - r[2] > 0 )
											inf = ww[2][4]/2.0;
											sup = ww[2][4];
											cur = (sup + inf)/2.0;
											while (sup - inf > precision){
												if ( exp4_tmp2 + representations(r[2], ww[2][4], cur) - r[2] > 0 ) {
													sup = cur;
												} else {
													inf = cur;
												}
												cur = (sup + inf)/2.0;
											}

											double lbound_w23, ubound_w23;
											if (psup.ww[2][3] < sup) {
												lbound_w23 = sup;
												ubound_w23 = sup;
											} else if (pinf.ww[2][3] > ww[2][4]) {
												lbound_w23 = ww[2][4];
												ubound_w23 = ww[2][4];
											} else {
												lbound_w23 = fmax(sup, pinf.ww[2][3]);
												ubound_w23 = fmin(ww[2][4], psup.ww[2][3]);
											}
											for (ww[2][3] = lbound_w23 ; ww[2][3] <= ubound_w23 ; ww[2][3] += step){
												exp4_tmp3 = exp4_tmp2 + representations(r[2], ww[2][4], ww[2][3]) - r[2]; // is greater than 0.0

												// we always have ( (exp3_tmp2 + representations(r[2], ww[2][3], r[2]/2.0)) >= r[2] )

												if ( (exp3_tmp2 + representations(r[2], ww[2][3], ww[2][3]/2.0)) > r[2] ) {
													ww[2][2] = ww[2][3]/2.0;
												} else {
													sup = r[2]/2.0;
													inf = ww[2][3]/2.0;
													cur = (sup + inf)/2.0;
													while (sup - inf > precision){
														if (exp3_tmp2 + representations(r[2], ww[2][3], cur) > r[2]) {
															sup = cur;
														} else {
															inf = cur;
														}
														cur = (sup + inf)/2.0;
													}
													ww[2][2] = sup;
												}
												double lbound_r3, ubound_r3;
												if (psup.r[3] < 0.0) {
													lbound_r3 = 0.0;
													ubound_r3 = 0.0;
												} else if (pinf.r[3] > 1.0 - R - r[1] - r[2]) {
													lbound_r3 = 1.0 - R - r[1] - r[2];
													ubound_r3 = 1.0 - R - r[1] - r[2];
												} else {
													lbound_r3 = fmax(0.0, pinf.r[3]);
													ubound_r3 = fmin(1.0 - R - r[1] - r[2], psup.r[3]);
												}
												auto optim_r3 = [this, exp4_tmp3, pinf, psup](double r3){
													r[3] = r3;
													r[4] = 1.0 - R - r[1] - r[2] - r[3];

													double sup, inf, cur;
													double C_cur, C_min = 10;
													BothMay params_C_min(*this);
													double lbound_w34, ubound_w34;
													if (psup.ww[3][4] < fmax(0.0, t - p[4] - ww[1][4] - ww[2][4] - r[4])) {
														lbound_w34 = fmax(0.0, t - p[4] - ww[1][4] - ww[2][4] - r[4]);
														ubound_w34 = fmax(0.0, t - p[4] - ww[1][4] - ww[2][4] - r[4]);
													} else if (pinf.ww[3][4] > fmin(r[3]/2.0, t - p[4] - ww[1][4] - ww[2][4])) {
														lbound_w34 = fmin(r[3]/2.0, t - p[4] - ww[1][4] - ww[2][4]);
														ubound_w34 = fmin(r[3]/2.0, t - p[4] - ww[1][4] - ww[2][4]);
													} else {
														lbound_w34 = fmax(fmax(0.0, t - p[4] - ww[1][4] - ww[2][4] - r[4]), pinf.ww[3][4]);
														ubound_w34 = fmin(fmin(r[3]/2.0, t - p[4] - ww[1][4] - ww[2][4]), psup.ww[3][4]);
													}

													for (ww[3][4] = lbound_w34 ; ww[3][4] <= ubound_w34 ; ww[3][4] += step){
														ww[4][4] = t - p[4] - ww[1][4] - ww[2][4] - ww[3][4];

														// we always have ( (exp4_tmp3 + representations(r[3], ww[3][4], r[3]/2.0)) < r[3] )

														if ( (exp4_tmp3 + representations(r[3], ww[3][4], ww[3][4]/2.0)) > r[3] ) {
															ww[3][3] = ww[3][4]/2.0;
														} else {
															sup = r[3]/2.0;
															inf = ww[3][4]/2.0;
															cur = (sup + inf)/2.0;
															while (sup - inf > precision){
																if (exp4_tmp3 + representations(r[3], ww[3][4], cur) > r[3]) {
																	sup = cur;
																} else {
																	inf = cur;
																}
																cur = (sup + inf)/2.0;
															}
															ww[3][3] = sup;
														}
														C_cur = complexity();
														if (C_cur < C_min){
															C_min = C_cur;
															params_C_min = *this;
														}
													}
													*this = params_C_min;
													return C_min;
												};

												double *opt = golden_section(optim_r3, lbound_r3, ubound_r3, precision);
												if (opt[0] < C_min){
													C_min = opt[0];
													params_C_min = *this;
												}
												free(opt);
											}
										}
										*this = params_C_min;
										return C_min;
									};

									double *opt = golden_section(optim_r2, lbound_r2, ubound_r2, precision);
									if (opt[0] < C_min){
										C_min = opt[0];
										params_C_min = *this;
									}
									free(opt);
								}
							}
						}
						*this = params_C_min;
						return C_min;
					};

					double *opt = golden_section(optim_r1, lbound_r1, ubound_r1, precision);
					if (opt[0] < C_min){
						C_min = opt[0];
						params_C_min = *this;
						cerr << "-->" << C_min << endl;
					}
					free(opt);
				}
			}
		}

		*this = params_C_min;
		return C_min;
	};


	double *opt = golden_section(optim_p4, lbound_p4, ubound_p4, precision);
	double res = opt[0];
	free(opt);
	return res;
}
*/

// reduced bound and use of golden section technique for optimizing parameters p4, p3, p2, p1
double BothMay::optimization_tmp(){
	verbose = false;

	BothMay pinf(*this);
	BothMay psup(*this);

	for (int i = 0 ; i < 5 ; ++i) {
		pinf.r[i] -= radius;
		psup.r[i] += radius;
		pinf.p[i] -= radius;
		psup.p[i] += radius;
		for (int j = 1 ; j < i ; ++j) {
			pinf.ww[j][i] -= radius;
			psup.ww[j][i] += radius;
		}
	}

	double lbound_p4, ubound_p4;
	if (psup.p[4] < 0.0) {
		lbound_p4 = 0.0;
		ubound_p4 = 0.0;
	} else if (pinf.p[4] > fmin(R/2.0, t)) {
		lbound_p4 = fmin(R/2.0, t);
		ubound_p4 = fmin(R/2.0, t);
	} else {
		lbound_p4 = fmax(0.0, pinf.p[4]);
		ubound_p4 = fmin(fmin(R/2.0, t), psup.p[4]);
	}
	auto optim_p4 = [this, pinf, psup](double p4){
		p[4] = p4;

		double lbound_p3, ubound_p3;
		if (psup.p[3] < p[4]/2.0) {
			lbound_p3 = p[4]/2.0;
			ubound_p3 = p[4]/2.0;
		} else if (pinf.p[3] > p[4]) {
			lbound_p3 = p[4];
			ubound_p3 = p[4];
		} else {
			lbound_p3 = fmax(p[4]/2.0, pinf.p[3]);
			ubound_p3 = fmin(p[4], psup.p[3]);
		}

		auto optim_p3 = [this, pinf, psup](double p3){
			p[3] = p3;

			double exp4_tmp1 = representations(R, p[4], p[3]);
			double lbound_p2, ubound_p2;
			if (psup.p[2] < p[3]/2.0) {
				lbound_p2 = p[3]/2.0;
				ubound_p2 = p[3]/2.0;
			} else if (pinf.p[2] > p[3]) {
				lbound_p2 = p[3];
				ubound_p2 = p[3];
			} else {
				lbound_p2 = fmax(p[3]/2.0, pinf.p[2]);
				ubound_p2 = fmin(p[3], psup.p[2]);
			}

			auto optim_p2 = [this, exp4_tmp1, pinf, psup](double p2){
				p[2] = p2;

				double exp3_tmp1 = representations(R, p[3], p[2]);
				double lbound_p1, ubound_p1;
				if (psup.p[1] < p[2]/2.0) {
					lbound_p1 = p[2]/2.0;
					ubound_p1 = p[2]/2.0;
				} else if (pinf.p[1] > p[2]) {
					lbound_p1 = p[2];
					ubound_p1 = p[2];
				} else {
					lbound_p1 = fmax(p[2]/2.0, pinf.p[1]);
					ubound_p1 = fmin(p[2], psup.p[1]);
				}
				auto optim_p1 = [this, exp4_tmp1, exp3_tmp1, pinf, psup](double p1){
					p[1] = p1;
					p[0] = p[1]/2.0;

					double sup, inf, cur;
					double C_cur, C_min = 10;
					BothMay params_C_min(*this);

					double exp4_tmp2;
					double exp3_tmp2;
					double exp4_tmp3;
					double exp2_tmp1 = representations(R, p[2], p[1]);

					double lbound_r1, ubound_r1;
					if (psup.r[1] < 0.0) {
						lbound_r1 = 0.0;
						ubound_r1 = 0.0;
					} else if (pinf.r[1] > 1.0 - R) {
						lbound_r1 = 1.0 - R;
						ubound_r1 = 1.0 - R;
					} else {
						lbound_r1 = fmax(0.0, pinf.r[1]);
						ubound_r1 = fmin(1.0 - R, psup.r[1]);
					}
					for (r[1] = lbound_r1 ; r[1] <= ubound_r1 ; r[1] += step){

						double lbound_w14, ubound_w14;
						if (psup.ww[1][4] < 0.0) {
							lbound_w14 = 0.0;
							ubound_w14 = 0.0;
						} else if (pinf.ww[1][4] > fmin(r[1]/2.0, t - p[4])) {
							lbound_w14 = fmin(r[1]/2.0, t - p[4]);
							ubound_w14 = fmin(r[1]/2.0, t - p[4]);
						} else {
							lbound_w14 = fmax(0.0, pinf.ww[1][4]);
							ubound_w14 = fmin(fmin(r[1]/2.0, t - p[4]), psup.ww[1][4]);
						}
						for (ww[1][4] = lbound_w14 ; ww[1][4] <= ubound_w14 ; ww[1][4] += step){

							double lbound_w13, ubound_w13;
							if (psup.ww[1][3] < ww[1][4]/2.0) {
								lbound_w13 = ww[1][4]/2.0;
								ubound_w13 = ww[1][4]/2.0;
							} else if (pinf.ww[1][3] > ww[1][4]) {
								lbound_w13 = ww[1][4];
								ubound_w13 = ww[1][4];
							} else {
								lbound_w13 = fmax(ww[1][4]/2.0, pinf.ww[1][3]);
								ubound_w13 = fmin(ww[1][4], psup.ww[1][3]);
							}
							for (ww[1][3] = lbound_w13 ; ww[1][3] <= ubound_w13 ; ww[1][3] += step){
								exp4_tmp2 = exp4_tmp1 + representations(r[1], ww[1][4], ww[1][3]) - r[1];

								double lbound_w12, ubound_w12;
								if (psup.ww[1][2] < ww[1][3]/2.0) {
									lbound_w12 = ww[1][3]/2.0;
									ubound_w12 = ww[1][3]/2.0;
								} else if (pinf.ww[1][2] > ww[1][3]) {
									lbound_w12 = ww[1][3];
									ubound_w12 = ww[1][3];
								} else {
									lbound_w12 = fmax(ww[1][3]/2.0, pinf.ww[1][2]);
									ubound_w12 = fmin(ww[1][3], psup.ww[1][2]);
								}
								for (ww[1][2] = lbound_w12 ; ww[1][2] <= ubound_w12 ; ww[1][2] += step){

									exp3_tmp2 = exp3_tmp1 + representations(r[1], ww[1][3], ww[1][2]) - r[1]; //is greater than 0.0

									if ( (exp2_tmp1 + representations(r[1], ww[1][2], r[1]/2.0)) < r[1] ) {
										break;
									}
									if ( (exp2_tmp1 + representations(r[1], ww[1][2], ww[1][2]/2.0)) > r[1] ) {
										ww[1][1] = ww[1][2]/2.0;
									} else {
										sup = r[1]/2.0;
										inf = ww[1][2]/2.0;
										cur = (sup + inf)/2.0;
										while (sup - inf > precision){
											if ( (exp2_tmp1 + representations(r[1], ww[1][2], cur)) > r[1] ) {
												sup = cur;
											} else {
												inf = cur;
											}
											cur = (sup + inf)/2.0;
										}
										ww[1][1] = sup;
									}

									double lbound_r2, ubound_r2;
									if (psup.r[2] < 0.0) {
										lbound_r2 = 0.0;
										ubound_r2 = 0.0;
									} else if (pinf.r[2] > 1.0 - R - r[1]) {
										lbound_r2 = 1.0 - R - r[1];
										ubound_r2 = 1.0 - R - r[1];
									} else {
										lbound_r2 = fmax(0.0, pinf.r[2]);
										ubound_r2 = fmin(1.0 - R - r[1], psup.r[2]);
									}
									for (r[2] = lbound_r2 ; r[2] <= ubound_r2 ; r[2] += step){

										double lbound_w24, ubound_w24;
										if (psup.ww[2][4] < 0.0) {
											lbound_w24 = 0.0;
											ubound_w24 = 0.0;
										} else if (pinf.ww[2][4] > fmin(r[2]/2.0, t - p[4] - ww[1][4])) {
											lbound_w24 = fmin(r[2]/2.0, t - p[4] - ww[1][4]);
											ubound_w24 = fmin(r[2]/2.0, t - p[4] - ww[1][4]);
										} else {
											lbound_w24 = fmax(0.0, pinf.ww[2][4]);
											ubound_w24 = fmin(fmin(r[2]/2.0, t - p[4] - ww[1][4]), psup.ww[2][4]);
										}
										for (ww[2][4] = lbound_w24 ; ww[2][4] <= ubound_w24 ; ww[2][4] += step){

											double lbound_w23, ubound_w23;
											if (psup.ww[2][3] < ww[2][4]/2.0) {
												lbound_w23 = ww[2][4]/2.0;
												ubound_w23 = ww[2][4]/2.0;
											} else if (pinf.ww[2][3] > ww[2][4]) {
												lbound_w23 = ww[2][4];
												ubound_w23 = ww[2][4];
											} else {
												lbound_w23 = fmax(ww[2][4]/2.0, pinf.ww[2][3]);
												ubound_w23 = fmin(ww[2][4], psup.ww[2][3]);
											}
											for (ww[2][3] = lbound_w23 ; ww[2][3] <= ubound_w23 ; ww[2][3] += step){
												exp4_tmp3 = exp4_tmp2 + representations(r[2], ww[2][4], ww[2][3]) - r[2]; // is greater than 0.0

												if ( (exp3_tmp2 + representations(r[2], ww[2][3], r[2]/2.0)) < r[2] ) {
													break;
												}
												if ( (exp3_tmp2 + representations(r[2], ww[2][3], ww[2][3]/2.0)) > r[2] ) {
													ww[2][2] = ww[2][3]/2.0;
												} else {
													sup = r[2]/2.0;
													inf = ww[2][3]/2.0;
													cur = (sup + inf)/2.0;
													while (sup - inf > precision){
														if (exp3_tmp2 + representations(r[2], ww[2][3], cur) > r[2]) {
															sup = cur;
														} else {
															inf = cur;
														}
														cur = (sup + inf)/2.0;
													}
													ww[2][2] = sup;
												}
												double lbound_r3, ubound_r3;
												if (psup.r[3] < 0.0) {
													lbound_r3 = 0.0;
													ubound_r3 = 0.0;
												} else if (pinf.r[3] > 1.0 - R - r[1] - r[2]) {
													lbound_r3 = 1.0 - R - r[1] - r[2];
													ubound_r3 = 1.0 - R - r[1] - r[2];
												} else {
													lbound_r3 = fmax(0.0, pinf.r[3]);
													ubound_r3 = fmin(1.0 - R - r[1] - r[2], psup.r[3]);
												}

												for  (r[3] = lbound_r3 ; r[3] <= ubound_r3 ; r[3] += step){
													r[4] = 1.0 - R - r[1] - r[2] - r[3];

													double lbound_w34, ubound_w34;
													if (psup.ww[3][4] < fmax(0.0, t - p[4] - ww[1][4] - ww[2][4] - r[4])) {
														lbound_w34 = fmax(0.0, t - p[4] - ww[1][4] - ww[2][4] - r[4]);
														ubound_w34 = fmax(0.0, t - p[4] - ww[1][4] - ww[2][4] - r[4]);
													} else if (pinf.ww[3][4] > fmin(r[3]/2.0, t - p[4] - ww[1][4] - ww[2][4])) {
														lbound_w34 = fmin(r[3]/2.0, t - p[4] - ww[1][4] - ww[2][4]);
														ubound_w34 = fmin(r[3]/2.0, t - p[4] - ww[1][4] - ww[2][4]);
													} else {
														lbound_w34 = fmax(fmax(0.0, t - p[4] - ww[1][4] - ww[2][4] - r[4]), pinf.ww[3][4]);
														ubound_w34 = fmin(fmin(r[3]/2.0, t - p[4] - ww[1][4] - ww[2][4]), psup.ww[3][4]);
													}

													for (ww[3][4] = lbound_w34 ; ww[3][4] <= ubound_w34 ; ww[3][4] += step){
														ww[4][4] = t - p[4] - ww[1][4] - ww[2][4] - ww[3][4];

														if ( (exp4_tmp3 + representations(r[3], ww[3][4], r[3]/2.0)) < r[3] ) {
															break;
														}
														if ( (exp4_tmp3 + representations(r[3], ww[3][4], ww[3][4]/2.0)) > r[3] ) {
															ww[3][3] = ww[3][4]/2.0;
														} else {
															sup = r[3]/2.0;
															inf = ww[3][4]/2.0;
															cur = (sup + inf)/2.0;
															while (sup - inf > precision){
																if (exp4_tmp3 + representations(r[3], ww[3][4], cur) > r[3]) {
																	sup = cur;
																} else {
																	inf = cur;
																}
																cur = (sup + inf)/2.0;
															}
															ww[3][3] = sup;
														}
														C_cur = complexity();
														if (C_cur < C_min){
															C_min = C_cur;
															params_C_min = *this;
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
					*this = params_C_min;
					return C_min;
				};
				double *opt = golden_section(optim_p1, lbound_p1, ubound_p1, precision);
				double res = opt[0];
				free(opt);
				return res;
			};
			double *opt = golden_section(optim_p2, lbound_p2, ubound_p2, precision);
			double res = opt[0];
			free(opt);
			cerr << "--> " << res << endl;
			return res;
		};

		double *opt = golden_section(optim_p3, lbound_p3, ubound_p3, precision);
		double res = opt[0];
		free(opt);
		return res;
	};


	double *opt = golden_section(optim_p4, lbound_p4, ubound_p4, precision);
	double res = opt[0];
	free(opt);
	return res;
}

double BothMay::optimization(){

	double cost, cost_pred = M_INF;
	BothMay pred(*this);

	do {
		cerr << "radius = " << radius << "     step = " << step << "     precision = " << precision << endl;
		cost = optimization_tmp();
		if (cost >= cost_pred) {
			radius /= 2.0;
			step /= 2.0;
			*this = pred;
		} else {
			cost_pred = cost;
			pred = *this;
		}
		cerr << pred << endl;
		cerr << "Cost = " << cost_pred << endl << endl;
	} while (step >= limit);

	return cost_pred;
}

