#ifndef SIDIS_FUNCTIONS_SCALE_H
#define SIDIS_FUNCTIONS_SCALE_H

#include "Common/Constants.cpp"

#include "SIDIS/Coefficients/EvaluationParameters.cpp"

namespace SIDISFunctions::Scale {
	constexpr double quark_to_quark(const EvaluationParameters &p) {
		const double term1_delta = p.pdf * p.ff * ((2.0 * std::log(1.0 - p.z) + 1.5) * p.fragmentation_scale_log + (2.0 * std::log(1.0 - p.x) + 1.5) * p.factorization_scale_log);
		const double term2_xi = p.factorization_scale_log * ((1.0 + std::pow(p.xi, 2)) * p.pdf_hat * p.ff - 2.0 * p.pdf * p.ff) / (1.0 - p.xi);
		const double term3_xip = p.fragmentation_scale_log * ((1.0 + std::pow(p.xip, 2)) * p.pdf * p.ff_hat - 2.0 * p.pdf * p.ff) / (1.0 - p.xip);

		const double result = term1_delta / ((1.0 - p.x) * (1.0 - p.z)) + term2_xi / (1.0 - p.z) + term3_xip / (1.0 - p.x);
		return 2.0 * Constants::C_F * result / p.z;
	}

	constexpr double quark_to_gluon(const EvaluationParameters &p) {
		const double result = p.fragmentation_scale_log * p.pdf * p.ff_hat * (1.0 + std::pow(1.0 - p.xip, 2)) / p.xip;
		return 2.0 * Constants::C_F * result / p.z;
	}

	constexpr double gluon_to_quark(const EvaluationParameters &p) {
		const double result = p.factorization_scale_log * p.pdf_hat * p.ff * (std::pow(p.xi, 2) + std::pow(1.0 - p.xi, 2));
		return 2.0 * Constants::T_R * result / p.z;
	}
}

#endif