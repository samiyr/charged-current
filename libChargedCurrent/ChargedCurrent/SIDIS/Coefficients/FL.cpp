#ifndef SIDIS_FUNCTIONS_FL_H
#define SIDIS_FUNCTIONS_FL_H

#include "SIDIS/Coefficients/F2.cpp"
#include "SIDIS/Coefficients/Helper.cpp"
#include "SIDIS/Coefficients/EvaluationParameters.cpp"

namespace SIDISFunctions::FL {
	namespace LO {
		constexpr double quark_to_quark(const EvaluationParameters &p) {
			const double prefactor = p.m2 / (p.m2 + p.Q2);
			return prefactor * SIDISFunctions::F2::LO::quark_to_quark(p);
		}
		constexpr double quark_to_gluon(const EvaluationParameters &p) {
			const double prefactor = p.m2 / (p.m2 + p.Q2);
			return prefactor * SIDISFunctions::F2::LO::quark_to_gluon(p);
		}
		constexpr double gluon_to_quark(const EvaluationParameters &p) {
			const double prefactor = p.m2 / (p.m2 + p.Q2);
			return prefactor * SIDISFunctions::F2::LO::gluon_to_quark(p);
		}
	}

	namespace NLO_NLP {
		constexpr double quark_to_quark(const EvaluationParameters &p) {
			const double mass_prefactor = p.m2 / (p.m2 + p.Q2);
			return mass_prefactor * F2::NLO_NLP::quark_to_quark(p);
		}
		constexpr double quark_to_gluon(const EvaluationParameters &) { return 0.0; }
		constexpr double gluon_to_quark(const EvaluationParameters &) { return 0.0; }
	}

	namespace NLO {
		constexpr double quark_to_quark(const EvaluationParameters &p) {
			const double longitudinal_prefactor = p.Q2 / (p.m2 + p.Q2);
			const double mass_prefactor = p.m2 / (p.m2 + p.Q2);

			const double longitudinal_term = p.pdf_hat * p.ff_hat * Constants::C_F * 4.0 * p.xi * p.xip;
			const double transverse_term = F2::NLO::quark_to_quark(p);

			const double longitudinal_contribution = 2.0 * longitudinal_prefactor * longitudinal_term / p.z;
			const double transverse_contribution = mass_prefactor * transverse_term;

			return longitudinal_contribution + transverse_contribution;
		}

		constexpr double quark_to_gluon(const EvaluationParameters &p) {
			const double longitudinal_prefactor = p.Q2 / (p.m2 + p.Q2);
			const double mass_prefactor = p.m2 / (p.m2 + p.Q2);

			const double longitudinal_term = p.pdf_hat * p.ff_hat * Constants::C_F * 4.0 * p.xi * (1.0 - p.xip);
			const double transverse_term = F2::NLO::quark_to_quark(p);

			const double longitudinal_contribution = 2.0 * longitudinal_prefactor * longitudinal_term / p.z;
			const double transverse_contribution = mass_prefactor * transverse_term;

			return longitudinal_contribution + transverse_contribution;
		}

		constexpr double gluon_to_quark(const EvaluationParameters &p) {
			const double longitudinal_prefactor = p.Q2 / (p.m2 + p.Q2);
			const double mass_prefactor = p.m2 / (p.m2 + p.Q2);

			const double longitudinal_term = p.pdf_hat * p.ff_hat * Constants::T_R * 8.0 * p.xi * (1.0 - p.xi);
			const double transverse_term = F2::NLO::quark_to_quark(p);

			const double longitudinal_contribution = 2.0 * longitudinal_prefactor * longitudinal_term / p.z;
			const double transverse_contribution = mass_prefactor * transverse_term;

			return longitudinal_contribution + transverse_contribution;
		}

		constexpr double total(const EvaluationParameters &p) {
			return quark_to_quark(p) + quark_to_gluon(p) + gluon_to_quark(p);
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

		constexpr double quark_to_quark(const EvaluationParameters &p) {
			const double mass_prefactor = p.m2 / (p.m2 + p.Q2);
			return mass_prefactor * F2::NNLO_NLP::quark_to_quark(p);
		}

		constexpr double quark_to_gluon(const EvaluationParameters &) { return 0.0; }
		constexpr double gluon_to_quark(const EvaluationParameters &) { return 0.0; }
	}
}

#endif