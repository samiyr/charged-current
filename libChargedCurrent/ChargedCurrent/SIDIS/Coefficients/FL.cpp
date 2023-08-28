#ifndef SIDIS_FUNCTIONS_FL_H
#define SIDIS_FUNCTIONS_FL_H

#include "SIDIS/Coefficients/F2.cpp"
#include "SIDIS/Coefficients/Helper.cpp"
#include "SIDIS/Coefficients/EvaluationParameters.cpp"

namespace SIDISFunctions::FL {
	namespace LO {
		constexpr double integrand(const EvaluationParameters &p) {
			const double prefactor = p.m2 / (p.m2 + p.Q2);
			return prefactor * SIDISFunctions::F2::LO::integrand(p);
		}
	}

	namespace NLO_NLP {
		constexpr double total_integrand(const EvaluationParameters &p) {
			const double mass_prefactor = p.m2 / (p.m2 + p.Q2);
			return mass_prefactor * F2::NLO_NLP::total_integrand(p);
		}
	}

	namespace NLO {
		constexpr double delta_integrand(const EvaluationParameters &p) {
			const double prefactor = p.m2 / (p.m2 + p.Q2);
			return prefactor * SIDISFunctions::F2::NLO::delta_integrand(p);
		}

		constexpr double xi_integrand(const EvaluationParameters &p) {
			const double prefactor = p.m2 / (p.m2 + p.Q2);
			return prefactor * SIDISFunctions::F2::NLO::xi_integrand(p);
		}
		constexpr double xip_integrand(const EvaluationParameters &p) {
			const double prefactor = p.m2 / (p.m2 + p.Q2);
			return prefactor * SIDISFunctions::F2::NLO::xip_integrand(p);
		}
		constexpr double xi_xip_integrand(const EvaluationParameters &p) {
			const double term1 = p.xq_hat_zq_hat * Constants::C_F * 4.0 * p.xi * p.xip;
			const double term2 = p.xq_hat_zg_hat * Constants::C_F * 4.0 * p.xi * (1.0 - p.xip);
			const double term3 = p.xg_hat_zq_hat * Constants::T_R * 8.0 * p.xi * (1.0 - p.xi);

			const double value = term1 + term2 + term3;
			const double longitudinal_prefactor = p.Q2 / (p.m2 + p.Q2);
			const double longitudinal_contribution = 2.0 * longitudinal_prefactor * value / p.z;

			const double mass_prefactor = p.m2 / (p.m2 + p.Q2);
			const double mass_contribution = mass_prefactor * SIDISFunctions::F2::NLO::xi_xip_integrand(p);

			return longitudinal_contribution + mass_contribution;
		}

		constexpr double total_integrand(const EvaluationParameters &p) {
			return Helper::make_nlo_integrand(
				FL::NLO::delta_integrand, FL::NLO::xi_integrand, FL::NLO::xip_integrand, FL::NLO::xi_xip_integrand, 
				p
			);
		}
	}

	namespace NNLO_NLP {
		constexpr double LP_CF(const EvaluationParameters &p) {
			const double mass_prefactor = p.m2 / (p.m2 + p.Q2);
			return mass_prefactor * F2::NNLO_NLP::LP_CF(p);
		}
		constexpr double NLP_CF(const EvaluationParameters &p) {
			const double mass_prefactor = p.m2 / (p.m2 + p.Q2);
			return mass_prefactor * F2::NNLO_NLP::NLP_CF(p);
		}
		constexpr double LP_CA(const EvaluationParameters &p) {
			const double mass_prefactor = p.m2 / (p.m2 + p.Q2);
			return mass_prefactor * F2::NNLO_NLP::LP_CA(p);
		}
		constexpr double LP_Nf(const EvaluationParameters &p) {
			const double mass_prefactor = p.m2 / (p.m2 + p.Q2);
			return mass_prefactor * F2::NNLO_NLP::LP_Nf(p);
		}

		constexpr double total_integrand(const EvaluationParameters &p) {
			const double mass_prefactor = p.m2 / (p.m2 + p.Q2);
			return mass_prefactor * F2::NNLO_NLP::total_integrand(p);
		}
	}
}

#endif