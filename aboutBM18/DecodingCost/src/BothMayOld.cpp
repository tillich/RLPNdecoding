/*
 * BothMayOld.cpp
 *
 *  Created on: 27 août 2022
 *      Author: Kévin Carrier
 */

#include "BothMayOld.h"

BothMayOld::BothMayOld(double R, double t) : Decoder(R, t){
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

BothMayOld::BothMayOld(double R, double t, double r[5], double p[5], double ww[5][5]) : Decoder(R,t){
	this->p[0] = p[0];
	for (int i = 1 ; i < 5 ; ++i){
		this->r[i] = r[i];
		this->p[i] = p[i];
		for (int j = 1 ; j <= i ; ++j){
			this->ww[j][i] = ww[j][i];
		}
	}
}

BothMayOld::BothMayOld(const BothMayOld &bm) : Decoder(bm.R, bm.t){
	p[0] = bm.p[0];
	for (int i = 1 ; i < 5 ; ++i){
		r[i] = bm.r[i];
		p[i] = bm.p[i];
		for (int j = 1 ; j <= i ; ++j){
			ww[j][i] = bm.ww[j][i];
		}
	}
}

void BothMayOld::operator=(const BothMayOld &bm) {
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


double BothMayOld::get_r(int i) const {
	return r[i];
}

double BothMayOld::get_p(int i) const {
	return p[i];
}

double BothMayOld::get_ww(int j, int i) const {
	return ww[j][i];
}

void BothMayOld::set_params(double *r, double *p, double **ww) {
	this->p[0] = p[0];
	for (int i = 1 ; i < 5 ; ++i){
		this->r[i] = r[i];
		this->p[i] = p[i];
		for (int j = 1 ; j <= i ; ++j){
			this->ww[j][i] = ww[j][i];
		}
	}
}

void BothMayOld::set_optim_params(double radius, double step, double limit, double precision){
	this->radius = radius;
	this->step = step;
	this->limit = limit;
	this->precision = precision;
}


ostream& operator<<(ostream& os, const BothMayOld& c) {
	os.setf(ios::fixed);
	os.precision(12);
	os << "Both-May parameters (old version):" << endl;
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


double BothMayOld::correctness(int l) const {
	double res = representations(R, p[l], p[l - 1]);
	for (int i = 1 ; i < l ; ++i) {
		res += (representations(r[i], ww[i][l], ww[i][l - 1]) - r[i]);
	}
	return res;
}

double BothMayOld::complexity() const{
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

		//old version of [BM18]
		L[l] = BINOM(R, p[l]) + BINOM(r[l], ww[l][l]) -  r[l];
		for (i = 1 ; i < l ; ++i){
			L[l] += ( representations(r[i], ww[i][l], ww[i][l-1]) + BINOM(r[i], ww[i][l]) - 2.0 * BINOM(r[i], ww[i][l - 1]) );
		}
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


// very large bounds
double BothMayOld::optimization_tmp(){
	verbose = false;

	double exp4_tmp1, exp4_tmp2, exp4_tmp3;
	list<double> list_ww33;
	double exp3_tmp1, exp3_tmp2;
	list<double> list_ww22;
	double exp2_tmp1;
	list<double> list_ww11;
	double sup, inf, cur;
	double C_cur, C_min = M_INF;
	BothMayOld params_C_min(*this);

	BothMayOld pinf(*this);
	BothMayOld psup(*this);

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


double BothMayOld::optimization(){

	double cost, cost_pred = M_INF;
	BothMayOld pred(*this);

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
