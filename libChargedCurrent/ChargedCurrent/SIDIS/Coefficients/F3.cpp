#ifndef SIDIS_FUNCTIONS_F3_H
#define SIDIS_FUNCTIONS_F3_H

#include "SIDIS/Coefficients/F2.cpp"
#include "SIDIS/Coefficients/Helper.cpp"
#include "SIDIS/Coefficients/EvaluationParameters.cpp"

namespace SIDISFunctions::F3 {
	namespace LO {
		// constexpr double integrand(const EvaluationParameters &p) {
		// 	return 2.0 * p.xq * p.zq / (p.x * p.z);
		// }

		constexpr double quark_to_quark(const EvaluationParameters &p) {
			return 2.0 * p.pdf * p.ff / (p.x * p.z);
		}

		constexpr double quark_to_gluon(const EvaluationParameters &p) { return 0.0; }
		constexpr double gluon_to_quark(const EvaluationParameters &p) { return 0.0; }
	}

	namespace NLO_NLP {
		constexpr double total_integrand(const EvaluationParameters &p) {
			return F2::NLO_NLP::total_integrand(p) / p.x;
		}
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

		// constexpr double delta_integrand(const EvaluationParameters &p) {
		// 	return F2::NLO::delta_integrand(p) / p.x;
		// }

		// constexpr double xi_integrand(const EvaluationParameters &p) {
		// 	return F2::NLO::xi_integrand(p) / p.x;
		// }

		// constexpr double xip_integrand(const EvaluationParameters &p) {
		// 	return F2::NLO::xip_integrand(p) / p.x;
		// }

		// constexpr double xi_xip_integrand(const EvaluationParameters &p) {
		// 	const double F2_value = F2::NLO::xi_xip_integrand(p);

		// 	const double term1 = Constants::C_F * p.xq_hat * p.zq_hat * (6.0 * p.xi * p.xip + 2.0 * (1.0 - p.xi - p.xip));
		// 	const double term2 = Constants::C_F * p.xq_hat * p.zg_hat * (4.0 * p.xi * (1.0 - p.xip) + 2.0 * (1.0 - p.xi) * p.xip);
		// 	const double term3 = Constants::T_R * p.xg_hat * p.zq_hat * (12.0 * p.xi * (1.0 - p.xi) + 2.0 * (1.0 - 2.0 * p.xi * (1.0 - p.xi)) / p.xip - 2.0);

		// 	const double value = F2_value / p.x - 2.0 * (term1 + term2 + term3) / (p.x * p.z);

		// 	return value;
		// }

		// constexpr double total_integrand(const EvaluationParameters &p) {
		// 	return Helper::make_nlo_integrand(
		// 		F3::NLO::delta_integrand, F3::NLO::xi_integrand, F3::NLO::xip_integrand, F3::NLO::xi_xip_integrand, 
		// 		p
		// 	);
		// }
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

		constexpr double total_integrand(const EvaluationParameters &p) {
			return F2::NNLO_NLP::total_integrand(p) / p.x;
		}
	}
}

#endif