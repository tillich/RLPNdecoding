/*
 * MayOzerov.h
 *
 *  Created on: 27 août 2022
 *      Author: Kévin Carrier
 */

#ifndef MAYOZEROV_H_
#define MAYOZEROV_H_

#include "Tools.h"
#include <list>
#include "Decoder.h"

using namespace std;
using namespace tools;

class MayOzerov : public Decoder{
private:
	double ell[3] = {0,0,0};
	double p[4] = {0,0,0,0};

	// optimization parameters
	double precision = 0.0000001;

public:
	MayOzerov(double R, double t);
	MayOzerov(double R, double t, double *ell, double *p);
	MayOzerov(const MayOzerov &mo);
	void operator=(const MayOzerov &mo);
	virtual ~MayOzerov(){return;};

	double complexity() const;
	double optimization();

	double get_ell(int i) const;
	double get_p(int i) const;
	void set_params(double *ell, double *p);
	void set_optim_param(double precision);

	friend ostream& operator<<(ostream& os, const MayOzerov& c);

};

#endif /* MAYOZEROV_H_ */
