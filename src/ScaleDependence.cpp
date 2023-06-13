#ifndef SCALE_DEPENDENCE_H
#define SCALE_DEPENDENCE_H

#include "TRFKinematics.cpp"

namespace ScaleDependence {
	static auto trivial = [](const TRFKinematics &kinematics) { return kinematics.Q2; };
}

#endif