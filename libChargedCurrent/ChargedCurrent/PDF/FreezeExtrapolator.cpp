#ifndef FREEZE_EXTRAPOLATOR_H
#define FREEZE_EXTRAPOLATOR_H

#include "LHAPDF/Extrapolator.h"
#include "LHAPDF/GridPDF.h"
#include "LHAPDF/KnotArray.h"

#include <vector>
#include <algorithm>

class FreezeExtrapolator : public LHAPDF::Extrapolator {
	public:

	double extrapolateXQ2(int id, double x, double q2) const {
		if (Q2_min < 0.0 && Q2_max < 0.0) {
			Q2_min = std::pow(pdf().info().get_entry_as<double>("QMin"), 2);
			Q2_max = std::pow(pdf().info().get_entry_as<double>("QMax"), 2);
		}
		const double clamped_Q2 = std::clamp(q2, Q2_min, Q2_max);
		return pdf().interpolator().interpolateXQ2(id, x, clamped_Q2);
	}

	private:

	mutable double Q2_min = -1.0;
	mutable double Q2_max = -1.0;
};


#endif