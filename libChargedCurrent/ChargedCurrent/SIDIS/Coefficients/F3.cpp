#ifndef SIDIS_FUNCTIONS_F3_H
#define SIDIS_FUNCTIONS_F3_H

#include "SIDIS/Coefficients/F2.cpp"
#include "SIDIS/Coefficients/Helper.cpp"
#include "SIDIS/Coefficients/EvaluationParameters.cpp"

namespace SIDISFunctions::F3 {
	namespace LO {
		constexpr double quark_to_quark(const EvaluationParameters &p) {
			return 2.0 * p.pdf * p.ff / (p.x * p.z);
		}

		constexpr double quark_to_gluon(const EvaluationParameters &) { return 0.0; }
		constexpr double gluon_to_quark(const EvaluationParameters &) { return 0.0; }
	}

	namespace NLO_NLP {
		constexpr double quark_to_quark(const EvaluationParameters &p) {
			return F2::NLO_NLP::quark_to_quark(p) / p.x;
		}
		constexpr double quark_to_gluon(const EvaluationParameters &) { return 0.0; }
		constexpr double gluon_to_quark(const EvaluationParameters &) { return 0.0; }
	}

	namespace NLO {
		constexpr double quark_to_quark(const EvaluationParameters &p) {
			const double F2_value = SIDISFunctions::F2::NLO::quark_to_quark(p);
			const double F3_term = Constants::C_F * p.pdf_hat * p.ff_hat * (6.0 * p.xi * p.xip + 2.0 * (1.0 - p.xi - p.xip));
			return F2_value / p.x - 2.0 * F3_term / (p.x * p.z);
		}

		constexpr double quark_to_gluon(const EvaluationParameters &p) {
			const double F2_value = SIDISFunctions::F2::NLO::quark_to_gluon(p);
			const double F3_term = Constants::C_F * p.pdf_hat * p.ff_hat * (4.0 * p.xi * (1.0 - p.xip) + 2.0 * (1.0 - p.xi) * p.xip);
			return F2_value / p.x - 2.0 * F3_term / (p.x * p.z);
		}

		constexpr double gluon_to_quark(const EvaluationParameters &p) {
			const double F2_value = SIDISFunctions::F2::NLO::gluon_to_quark(p);
			const double F3_term = Constants::T_R * p.pdf_hat * p.ff_hat * (12.0 * p.xi * (1.0 - p.xi) + 2.0 * (1.0 - 2.0 * p.xi * (1.0 - p.xi)) / p.xip - 2.0);
			return F2_value / p.x - 2.0 * F3_term / (p.x * p.z);
		}
	}

	namespace NNLO_NLP {
		constexpr double LP_CF(const EvaluationParameters &p) {
			return F2::NNLO_NLP::LP_CF(p) / p.x;
		}
		constexpr double NLP_CF(const EvaluationParameters &p) {
			return F2::NNLO_NLP::NLP_CF(p) / p.x;
		}
		constexpr double LP_CA(const EvaluationParameters &p) {
			return F2::NNLO_NLP::LP_CA(p) / p.x;
		}
		constexpr double LP_Nf(const EvaluationParameters &p) {
			return F2::NNLO_NLP::LP_Nf(p) / p.x;
		}

		constexpr double quark_to_quark(const EvaluationParameters &p) {
			return F2::NNLO_NLP::quark_to_quark(p) / p.x;
		}
		constexpr double quark_to_gluon(const EvaluationParameters &) { return 0.0; }
		constexpr double gluon_to_quark(const EvaluationParameters &) { return 0.0; }
	}
}

#endif