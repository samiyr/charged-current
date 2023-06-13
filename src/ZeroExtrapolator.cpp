#ifndef ZERO_EXTRAPOLATOR_H
#define ZERO_EXTRAPOLATOR_H

#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/Extrapolator.h"

class ZeroExtrapolator : public LHAPDF::Extrapolator {
	double extrapolateXQ2(int, double, double) const {
		return 0.0;
	}
};

#endif