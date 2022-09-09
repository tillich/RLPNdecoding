/*
 * Decoder.h
 *
 *  Created on: 10 août 2022
 *      Author: Kévin Carrier
 */

#ifndef DECODER_H_
#define DECODER_H_

#include <iostream>

using namespace std;

class Decoder {
protected:
	double R;
	double t;
	bool large;
	bool verbose = false;
public:
	Decoder(double R, double t);
	virtual ~Decoder(){return;}

	double get_R() const;
	void set_R(double R);
	double get_t() const;
	void set_t(double t);
	void set_verbose(bool v);

	virtual double complexity() const = 0;
	virtual double optimization() = 0;

	friend ostream& operator<<(ostream& os, const Decoder& c);
};

#endif /* DECODER_H_ */
