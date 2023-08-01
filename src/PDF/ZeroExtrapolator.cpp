#ifndef ZERO_EXTRAPOLATOR_H
#define ZERO_EXTRAPOLATOR_H

#include "LHAPDF/Extrapolator.h"

class ZeroExtrapolator : public LHAPDF::Extrapolator {
	double extrapolateXQ2(int, double, double) const override {
		return 0.0;
	}
};

#endif