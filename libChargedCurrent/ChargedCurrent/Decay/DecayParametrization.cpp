#ifndef DECAY_PARAMETRIZATION_H
#define DECAY_PARAMETRIZATION_H

#include "Utility/Math.cpp"
#include "Utility/Utility.cpp"

struct DecayParametrization {
	constexpr static DecayParametrization fit1() noexcept {
		return DecayParametrization(7.365, 1.4, 2.276, 2.04);
	}
	constexpr static DecayParametrization fit2() noexcept {
		return DecayParametrization(4.62698, 1.17383, 2.0603, 2.0565);
	}

	constexpr static std::vector<DecayParametrization> fit2_set() noexcept {
		const std::vector<std::vector<double>> parameter_set = {
			{4.13538, 1.11426, 2.01064, 2.05302}, 
			{3.52344, 1.05632, 1.86745, 2.08813}, 
			{5.49159, 1.26031, 2.17415, 2.02828}, 
			{4.65219, 1.13986, 2.16716, 2.02684}, 
			{4.81272, 1.15788, 2.23616, 2.00853}, 
			{1.81207, 0.803126, 1.33138, 2.2215}, 
			{2.48744, 0.884559, 1.66393, 2.14038}, 
			{2.16707, 0.830299, 1.60635, 2.13693}, 
			{6.0422, 1.30085, 2.15497, 2.057}, 
			{3.08592, 0.997451, 1.81307, 2.09004}, 
			{2.11793, 0.835383, 1.56551, 2.15002}, 
			{2.55209, 0.903038, 1.72966, 2.10214}, 
			{5.17884, 1.24278, 2.0241, 2.08523}, 
			{2.2788, 0.845369, 1.63586, 2.13261}, 
			{3.1376, 0.964036, 2.02487, 2.01567}, 
			{3.64627, 1.06438, 1.89462, 2.089}, 
			{2.6077, 0.926395, 1.64708, 2.14895}, 
			{5.47691, 1.23839, 2.17317, 2.03417}, 
			{12.0232, 1.59407, 2.71792, 1.95302}, 
			{4.61904, 1.1588, 2.14615, 2.0207}, 
			{1.67603, 0.705604, 1.4712, 2.15601}, 
			{2.61834, 0.941332, 1.70141, 2.10443}, 
			{4.58568, 1.12576, 2.18259, 2.02607}, 
			{1.09022, 0.532911, 1.18128, 2.24018}, 
			{3.12688, 1.00719, 1.73969, 2.13056}, 
			{3.7796, 1.0744, 1.91543, 2.08653}, 
			{6.49828, 1.31533, 2.33764, 2.0021}, 
			{6.61863, 1.33353, 2.31087, 2.0101}, 
			{4.062, 1.12687, 1.91276, 2.0932}
		};

		std::vector<DecayParametrization> parametrizations;

		for (const std::vector<double> &parameters : parameter_set) {
			parametrizations.emplace_back(parameters[0], parameters[1], parameters[2], parameters[3]);
		}

		return parametrizations;
	}

	double N;
	double alpha;
	double beta;
	double gamma;

	double beta_term_1p_alpha_beta;
	double beta_term_2p_alpha_beta;

	double gamma_prefactor_term;

	constexpr DecayParametrization() noexcept
	: N(0.0), alpha(0.0), beta(0.0), gamma(0.0), 
	beta_term_1p_alpha_beta(0.0), beta_term_2p_alpha_beta(0.0), gamma_prefactor_term(0.0) {}

	constexpr DecayParametrization(
		const double _N,
		const double _alpha,
		const double _beta,
		const double _gamma
	) noexcept
	: N(_N),
	alpha(_alpha),
	beta(_beta),
	gamma(_gamma),
	beta_term_1p_alpha_beta(Math::beta(1 + alpha, 1 + beta)),
	beta_term_2p_alpha_beta(Math::beta(2 + alpha, 1 + beta)),
	gamma_prefactor_term(std::pow(gamma, -1 - alpha))
	{ }

};


#endif