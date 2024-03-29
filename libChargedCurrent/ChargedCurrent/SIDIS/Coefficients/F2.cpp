#ifndef SIDIS_FUNCTIONS_F2_H
#define SIDIS_FUNCTIONS_F2_H

#include "SIDIS/Coefficients/Scale.cpp"
#include "SIDIS/Coefficients/Helper.cpp"
#include "SIDIS/Coefficients/NearThreshold.cpp"
#include "SIDIS/Coefficients/EvaluationParameters.cpp"

namespace SIDISFunctions::F2 {
	namespace LO {
		constexpr double quark_to_quark(const EvaluationParameters &p) {
			return 2.0 * p.pdf * p.ff / p.z;
		}

		constexpr double quark_to_gluon(const EvaluationParameters &) { return 0.0; }
		constexpr double gluon_to_quark(const EvaluationParameters &) { return 0.0; }
	}

	namespace NLO_NLP {
		constexpr double quark_to_quark(const EvaluationParameters &p) {
			using namespace NearThreshold;
			using namespace Constants;
			using namespace std::numbers;

			const auto evaluate = [&p](const double xq, const double xq_hat, const double zq, const double zq_hat) -> double {
				const double dx = d(p.x, xq);
				const double dz = d(p.z, zq);

				const auto Dx = [&](const std::size_t n) { return D(n, p.x, p.xi, xq, xq_hat); };
				const auto Dz = [&](const std::size_t n) { return D(n, p.z, p.xip, zq, zq_hat); };

				const auto lx = [&](const std::size_t n) -> double { return l(n, p.xi, xq_hat); };
				const auto lz = [&](const std::size_t n) -> double { return l(n, p.xip, zq_hat); };

				const double term1 = 2.0 * (dx * Dz(1) + dz * Dx(1) + Dx(0) * Dz(0) - 4.0 * dx * dz);

				// The zq_hat and xq_hat multiplication is because the first two terms are the only ones where there is only a x-term or a z-term,
				// which in this approach means that the PDF/FF value is missed otherwise.
				const double term2 = -2.0 * (Dx(0) * zq_hat + xq_hat * Dz(0));
				
				const double term3 = -2.0 * (dx * lz(1) + dz * lx(1));

				const double coefficient_contribution = term1 + term2 + term3;

				const double scale_1 = dx * (2.0 * Dz(0) + 1.5 * dz) + dz * (2.0 * Dx(0) + 1.5 * dx);

				const double scale_contribution = -p.factorization_scale_log * scale_1;

				return coefficient_contribution + scale_contribution;
			};

			return 2.0 * evaluate(p.pdf, p.pdf_hat, p.ff, p.ff_hat) / p.z;
		}

		constexpr double quark_to_gluon(const EvaluationParameters &) { return 0.0; }
		constexpr double gluon_to_quark(const EvaluationParameters &) { return 0.0; }
	}

	namespace NLO {
		constexpr double quark_to_quark(const EvaluationParameters &p) {
			const double term1_delta = p.pdf * p.ff * (std::pow(std::log(1.0 - p.x) + std::log(1.0 - p.z), 2) - 8.0);

			const double term2_xi = p.pdf_hat * p.ff * ((1.0 - p.xi) * (1.0 + std::log((1.0 - p.xi) / p.xi) + std::log(1.0 - p.z)) - (2.0 * p.xi * std::log(p.xi)) / (1.0 - p.xi));
			const double term3_xi = (2.0 * (std::log(1.0 - p.xi) + std::log(1.0 - p.z))) * ((p.xi * p.pdf_hat * p.ff - p.pdf * p.ff) / (1.0 - p.xi));

			const double term4_xip = p.pdf * p.ff_hat * ((1.0 - p.xip) * (1.0 + std::log(p.xip * (1.0 - p.xip)) + std::log(1.0 - p.x)) + (2.0 * p.xip * std::log(p.xip)) / (1.0 - p.xip));
			const double term5_xip = (2.0 * (std::log(1.0 - p.xip) + std::log(1.0 - p.x))) * ((p.xip * p.pdf * p.ff_hat - p.pdf * p.ff) / (1.0 - p.xip));

			const double term6_xi_xip = (1.0 - p.xi) / (1.0 - p.xip) * (p.pdf_hat * p.ff_hat - p.pdf_hat * p.ff);
			const double term7_xi_xip = (1.0 - p.xip) / (1.0 - p.xi) * (p.pdf_hat * p.ff_hat - p.pdf * p.ff_hat);
			const double term8_xi_xip = 2.0 * (p.xi * p.xip * p.pdf_hat * p.ff_hat - p.xi * p.pdf_hat * p.ff - p.xip * p.pdf * p.ff_hat + p.pdf * p.ff) / ((1.0 - p.xi) * (1.0 - p.xip));
			const double term9_xi_xip = 6.0 * p.xi * p.xip * p.pdf_hat * p.ff_hat;

			const double result = term1_delta / ((1.0 - p.x) * (1.0 - p.z)) + (term2_xi + term3_xi) / (1.0 - p.z) + (term4_xip + term5_xip) / (1.0 - p.x) + (term6_xi_xip + term7_xi_xip + term8_xi_xip + term9_xi_xip);
			const double coefficient_contribution = 2.0 * Constants::C_F * result / p.z;

			const double scale_contribution = (p.factorization_scale_log == 0 && p.fragmentation_scale_log == 0) ? 0 : Scale::quark_to_quark(p);
			
			return coefficient_contribution + scale_contribution;
		}

		constexpr double quark_to_gluon(const EvaluationParameters &p) {
			const double term1_xip = p.pdf * p.ff_hat * (p.xip + std::log(p.xip * (1.0 - p.xip)) * (1 + std::pow(1.0 - p.xip, 2)) / p.xip);
			const double term2_xip = p.pdf * p.ff_hat * (std::log(1.0 - p.x) * (p.xip + 2.0 * (1.0 - p.xip) / p.xip));

			const double term3_xi_xip = p.pdf_hat * p.ff_hat * (6.0 * p.xi * (1.0 - p.xip) + (1.0 - p.xi) / p.xip);
			const double term4_xi_xip = (p.pdf_hat * p.ff_hat * (p.xip + 2.0 * p.xi * (1.0 - p.xip) / p.xip) - p.pdf * p.ff_hat * (p.xip + 2.0 * (1.0 - p.xip) / p.xip)) / (1.0 - p.xi);

			const double result = (term1_xip + term2_xip) / (1.0 - p.x) + (term3_xi_xip + term4_xi_xip);
			const double coefficient_contribution = 2.0 * Constants::C_F * result / p.z;

			const double scale_contribution = (p.factorization_scale_log == 0 && p.fragmentation_scale_log == 0) ? 0 : Scale::quark_to_gluon(p);
			
			return coefficient_contribution + scale_contribution;
		}

		constexpr double gluon_to_quark(const EvaluationParameters &p) {
			const double term1_xi = p.pdf_hat * p.ff * (1.0 - (std::pow(p.xi, 2) + std::pow(1.0 - p.xi, 2)) * (1.0 - std::log((1.0 - p.xi) / p.xi)));
			const double term2_xi = p.pdf_hat * p.ff * (std::log(1.0 - p.z) * (1.0 - 2.0 * p.xi * (1.0 - p.xi)));

			const double term3_xi_xip = p.pdf_hat * p.ff_hat * (12.0 * p.xi * (1.0 - p.xi) + (1.0 - p.xip) / p.xip);
			const double term4_xi_xip = (p.pdf_hat * p.ff_hat * (p.xip - 2.0 * p.xi * (1.0 - p.xi) / p.xip) - p.pdf_hat * p.ff * (1.0 - 2.0 * p.xi * (1.0 - p.xi))) / (1.0 - p.xip);

			const double result = (term1_xi + term2_xi) / (1.0 - p.z) + (term3_xi_xip + term4_xi_xip);
			const double coefficient_contribution = 2.0 * Constants::T_R * result / p.z;

			const double scale_contribution = (p.factorization_scale_log == 0 && p.fragmentation_scale_log == 0) ? 0 : Scale::gluon_to_quark(p);
			
			return coefficient_contribution + scale_contribution;
		}
	}

	namespace NNLO_NLP {
		constexpr double LP_CF(const EvaluationParameters &p) {
			using namespace NearThreshold;
			using namespace Constants;
			using namespace std::numbers;

			const auto evaluate = [&p](const double xq, const double xq_hat, const double zq, const double zq_hat) -> double {
				const double dx = d(p.x, xq);
				const double dz = d(p.z, zq);

				const auto Dx = [&](const std::size_t n) { return D(n, p.x, p.xi, xq, xq_hat); };
				const auto Dz = [&](const std::size_t n) { return D(n, p.z, p.xip, zq, zq_hat); };

				const double term1 = dx * Dz(3) + dz * Dx(3);
				const double term2 = Dx(0) * Dz(2) + Dz(0) * Dx(2) + 2.0 * Dx(1) * Dz(1);
				const double term3 = Dx(0) * Dz(0) + dx * Dz(1) + dz * Dx(1);
				const double term4 = dx * Dz(0) + dz * Dx(0);
				const double term5 = dx * dz * (511.0 / 16.0 - 15.0 * zeta3 + 29.0 * std::pow(pi, 2) / 12.0 - 7.0 * std::pow(pi, 4) / 90.0);

				const double coefficient_contribution = 2.0 * term1 + 6.0 * term2 - (16.0 + 4.0 * std::pow(pi, 2) / 3.0) * term3 + 8.0 * zeta3 * term4 + term5;

				const double scale_1 = 4.0 * term3 + 6.0 * term4 + dx * dz * (9.0 / 2.0 - 2.0 * std::pow(pi, 2) / 3.0);
				const double scale_2 = -6.0 * (dx * Dz(2) + dz * Dx(2) + 2.0 * Dx(0) * Dz(1) + 2.0 * Dz(0) * Dx(1) + Dx(0) * Dz(0) + dx * Dz(1) + dz * Dx(1));
				const double scale_3 = (16.0 + 4.0 * std::pow(pi, 2) / 3.0) * term4 + dx * dz * (-20.0 * zeta3 + std::pow(pi, 2) + 93.0 / 4.0);

				const double scale_contribution = std::pow(p.factorization_scale_log, 2) * scale_1 + p.factorization_scale_log * (scale_2 + scale_3);

				return coefficient_contribution + scale_contribution;
			};

			return 2.0 * evaluate(p.pdf, p.pdf_hat, p.ff, p.ff_hat) / p.z;
		}
		constexpr double NLP_CF(const EvaluationParameters &p) {
			using namespace NearThreshold;
			using namespace Constants;
			using namespace std::numbers;

			const auto evaluate = [&p](const double xq, const double xq_hat, const double zq, const double zq_hat) -> double {
				const double dx = d(p.x, xq);
				const double dz = d(p.z, zq);

				const auto Dx = [&](const std::size_t n) -> double { return D(n, p.x, p.xi, xq, xq_hat); };
				const auto Dz = [&](const std::size_t n) -> double { return D(n, p.z, p.xip, zq, zq_hat); };

				const auto lx = [&](const std::size_t n) -> double { return l(n, p.xi, xq_hat); };
				const auto lz = [&](const std::size_t n) -> double { return l(n, p.xip, zq_hat); };

				// The zq_hat and xq_hat multiplication is because the first two terms are the only ones where there is only a x-term or a z-term,
				// which in this approach means that the PDF/FF value is missed otherwise.
				const double term1 = Dx(2) * zq_hat + Dz(2) * xq_hat + 2.0 * Dx(1) * lz(1) + 2.0 * Dz(1) * lx(1) + Dx(0) * lz(2) + Dz(0) * lx(2);
				const double term2 = dx * lz(3) + dz * lx(3);

				return -6.0 * term1 - 2.0 * term2;
			};

			return 2.0 * evaluate(p.pdf, p.pdf_hat, p.ff, p.ff_hat) / p.z;
		}
		constexpr double LP_CA(const EvaluationParameters &p) {
			using namespace NearThreshold;
			using namespace Constants;
			using namespace std::numbers;

			const auto evaluate = [&p](const double xq, const double xq_hat, const double zq, const double zq_hat) -> double {
				const double dx = d(p.x, xq);
				const double dz = d(p.z, zq);

				const auto Dx = [&](const std::size_t n) { return D(n, p.x, p.xi, xq, xq_hat); };
				const auto Dz = [&](const std::size_t n) { return D(n, p.z, p.xip, zq, zq_hat); };

				const double term1 = dx * Dz(2) + dz * Dx(2) + 2.0 * Dx(0) * Dz(1) + 2.0 * Dz(0) * Dx(1);
				const double term2 = Dx(0) * Dz(0) + dx * Dz(1) + dz * Dx(1);
				const double term3 = dx * Dz(0) + dz * Dx(0);
				const double term4 = dx * dz * (43.0 * zeta3 / 3.0 + 17.0 * std::pow(pi, 4) / 180.0 - 1535.0 / 48.0 - 269.0 * std::pow(pi, 2) / 108.0);

				const double coefficient_contribution_1 = -11.0 / 6.0 * term1 + (67.0 / 9.0 - std::pow(pi, 2) / 3.0) * term2;
				const double coefficient_contribution_2 = (7.0 * zeta3 + 11.0 * std::pow(pi, 2) / 18.0 - 101.0 / 13.5) * term3 + term4;
				const double coefficient_contribution = coefficient_contribution_1 + coefficient_contribution_2;

				const double scale_1 = term3 + 1.5 * dx * dz;
				const double scale_2 = -(67.0 / 9.0 - std::pow(pi, 2) / 3.0) * term3 + dx * dz * (6.0 * zeta3 - 11.0 * std::pow(pi, 2) / 9.0 - 17.0 / 12.0);

				const double scale_contribution = 11.0 / 6.0 * scale_1 * std::pow(p.factorization_scale_log, 2) + scale_2 * p.factorization_scale_log;

				return coefficient_contribution + scale_contribution;
			};

			return 2.0 * evaluate(p.pdf, p.pdf_hat, p.ff, p.ff_hat) / p.z;
		}
		constexpr double LP_Nf(const EvaluationParameters &p) {
			using namespace NearThreshold;
			using namespace Constants;
			using namespace std::numbers;

			const auto evaluate = [&p](const double xq, const double xq_hat, const double zq, const double zq_hat) -> double {
				const double dx = d(p.x, xq);
				const double dz = d(p.z, zq);

				const auto Dx = [&](const std::size_t n) { return D(n, p.x, p.xi, xq, xq_hat); };
				const auto Dz = [&](const std::size_t n) { return D(n, p.z, p.xip, zq, zq_hat); };

				const double term1 = dx * Dz(2) + dz * Dx(2) + 2.0 * Dx(0) * Dz(1) + 2.0 * Dz(0) * Dx(1);
				const double term2 = Dx(0) * Dz(0) + dx * Dz(1) + dz * Dx(1);
				const double term3 = dx * Dz(0) + dz * Dx(0);
				const double term4 = dx * dz * (zeta3 / 1.5 + 19.0 * std::pow(pi, 2) / 54.0 + 127.0 / 24.0);

				const double coefficient_contribution = term1 / 3.0 - 5.0 / 4.5 * term2 + term3 * (7.0 / 6.75 - std::pow(pi, 2) / 9.0) + term4;

				const double scale_1 = term3 + 1.5 * dx * dz;
				const double scale_2 = 5.0 / 4.5 * term3 + dx * dz * (1.0 / 6.0 + std::pow(pi, 2) / 4.5);

				const double scale_contribution = std::pow(p.factorization_scale_log, 2) * scale_1 / 3.0 + p.factorization_scale_log * scale_2;

				return coefficient_contribution + scale_contribution;
			};

			return 2.0 * evaluate(p.pdf, p.pdf_hat, p.ff, p.ff_hat) / p.z;
		}

		constexpr double quark_to_quark(const EvaluationParameters &p) {
			const double CF = LP_CF(p) + NLP_CF(p);
			const double CA = LP_CA(p);
			const double Nf = LP_Nf(p);

			const double nlo = F2::NLO::quark_to_quark(p) + F2::NLO::quark_to_gluon(p) + F2::NLO::gluon_to_quark(p);

			const double LP = Constants::C_F * (Constants::C_F * CF + Constants::C_A * CA + Constants::N_f * Nf);
			const double NLP = Constants::C_F * nlo * p.renormalization_scale_log * (11.0 * Constants::C_A - 2.0 * Constants::N_f) / 3.0;

			return LP + NLP;
		}

		constexpr double quark_to_gluon(const EvaluationParameters &) { return 0.0; }
		constexpr double gluon_to_quark(const EvaluationParameters &) { return 0.0; }
	}
}

#endif