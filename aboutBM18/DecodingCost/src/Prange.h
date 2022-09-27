/*
 * Prange.h
 *
 *  Created on: 10 août 2022
 *      Author: Kévin Carrier
 */

#ifndef PRANGE_H_
#define PRANGE_H_

#include "Decoder.h"
#include "Tools.h"

using namespace std;
using namespace tools;

class Prange : public Decoder {
public:
	Prange(double R, double t) : Decoder(R,t){return;}
	virtual ~Prange(){return;}
	double complexity() const;
	double optimization();

	friend ostream& operator<<(ostream& os, const Prange& c);
};

#endif /* PRANGE_H_ */
