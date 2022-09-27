/*
 * SternMayOzerov.h
 *
 *  Created on: 10 août 2022
 *      Author: Kévin Carrier
 */

#ifndef STERNMAYOZEROV_H_
#define STERNMAYOZEROV_H_

#include "Decoder.h"
#include "Tools.h"
#define SMO_PRECISION 0.0000000001

using namespace std;
using namespace tools;

class SternMayOzerov : public Decoder {
private:
	double p = 0;
public:
	SternMayOzerov(double R, double t) : Decoder(R, t){}
	SternMayOzerov(double R, double t, double p) : Decoder(R,t), p(p){}
	virtual ~SternMayOzerov(){}
	double complexity() const;
	double optimization();

	double get_p() const;
	void set_p(double p);

	friend ostream& operator<<(ostream& os, const SternMayOzerov& c);
};

#endif /* STERNMAYOZEROV_H_ */
