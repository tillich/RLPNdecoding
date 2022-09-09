/*
 * BothMay.h
 *
 *  Created on: 10 août 2022
 *      Author: Kévin Carrier
 */

#ifndef BOTHMAY_H_
#define BOTHMAY_H_

#include "Tools.h"
#include <list>
#include "Decoder.h"

using namespace std;
using namespace tools;

class BothMay : public Decoder{
private:
	double r[5] = {0,0,0,0,0};
	double p[5] = {0,0,0,0,0};
	double ww[5][5] = {{0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0}};

	// optimization parameters
	double radius = 0.1;
	double step = 0.01;
	double limit =0.01;
	double precision = 0.001;


	double correctness(int l) const;
	double optimization_tmp();
public:
	BothMay(double R, double t);
	BothMay(double R, double t, double r[5], double p[5], double ww[5][5]);
	BothMay(const BothMay &bm);
	void operator=(const BothMay &bm);
	virtual ~BothMay(){return;};

	double complexity() const;
	double optimization();

	double get_r(int i) const;
	double get_p(int i) const;
	double get_ww(int j, int i) const;
	void set_params(double *r, double *p, double **ww);
	void set_optim_params(double radius, double step, double limit, double precision);

	friend ostream& operator<<(ostream& os, const BothMay& c);

};

#endif /* BOTHMAY_H_ */
