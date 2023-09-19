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

	constexpr static std::vector<DecayParametrization> fit_set_1() noexcept {
		const std::vector<std::vector<double>> parameter_set = {
			{4.62698, 1.17383, 2.0603, 2.0565},
			{2.08719, 0.826902, 1.52356, 2.16617}, 
			{3.34704, 1.01933, 1.95326, 2.04748}, 
			{4.69181, 1.19151, 1.96797, 2.09349}, 
			{5.87475, 1.29053, 2.10034, 2.07548}, 
			{6.17874, 1.29798, 2.15311, 2.06994}, 
			{8.98027, 1.44098, 2.61326, 1.95797}, 
			{6.58788, 1.34046, 2.18092, 2.06322}, 
			{7.59974, 1.40511, 2.31695, 2.02535}, 
			{20.0175, 1.7428, 3.48563, 1.82029}, 
			{2.52759, 0.913191, 1.66408, 2.13062}, 
			{4.26517, 1.15535, 1.93272, 2.09176}, 
			{13.4197, 1.59451, 3.09742, 1.86636}, 
			{3.15425, 0.997212, 1.76983, 2.11385}, 
			{5.71975, 1.24304, 2.30654, 2.00268}, 
			{1.01172, 0.515651, 1.10411, 2.26537}, 
			{3.96763, 1.10612, 1.98074, 2.0528}, 
			{10.6719, 1.53328, 2.63587, 1.97292}, 
			{7.24444, 1.36957, 2.38498, 1.99751}, 
			{3.26233, 0.993048, 1.93223, 2.05932}, 
			{1.63023, 0.69704, 1.48406, 2.14328}, 
			{1.23217, 0.608482, 1.20006, 2.24114}, 
			{2.83692, 0.93288, 1.87074, 2.06598}, 
			{4.10808, 1.11869, 2.02653, 2.04485}, 
			{3.55282, 1.0515, 1.88293, 2.08353}, 
			{2.43913, 0.895598, 1.67811, 2.11014}, 
			{14.1133, 1.63118, 3.00019, 1.89926}, 
			{3.34674, 1.04442, 1.76653, 2.13247}, 
			{2.11637, 0.829057, 1.55449, 2.14926}, 
			{2.47148, 0.895346, 1.63257, 2.13784}, 
			{10.3606, 1.51536, 2.59392, 1.98142}, 
			{2.90272, 0.966343, 1.7705, 2.10535}, 
			{4.06582, 1.11516, 1.91318, 2.10145}, 
			{2.66543, 0.913082, 1.73577, 2.10845}, 
			{12.8256, 1.60034, 2.90148, 1.90832}
		};

		std::vector<DecayParametrization> parametrizations;

		for (const std::vector<double> &parameters : parameter_set) {
			parametrizations.emplace_back(parameters[0], parameters[1], parameters[2], parameters[3]);
		}

		return parametrizations;
	}

	constexpr static std::vector<DecayParametrization> fit_set_2() noexcept {
		const std::vector<std::vector<double>> parameter_set = {
			{5.12847, 1.17875, 2.41374, 1.93202},
			{2.01662, 0.757534, 1.89883, 1.97315}, 
			{4.27705, 1.12525, 2.04765, 2.03673}, 
			{10.2546, 1.49811, 2.71334, 1.93124}, 
			{9.0773, 1.41954, 2.62996, 1.95571}, 
			{35.0458, 1.93277, 4.75843, 1.59646}, 
			{9.14676, 1.45585, 2.57359, 1.9642}, 
			{17.5179, 1.73202, 3.26352, 1.83223}, 
			{13176.3, 4.08285, 42.3206, 0.472469}, 
			{1.68427, 0.702519, 1.60323, 2.08147}, 
			{150.764, 2.47761, 8.17915, 1.27098}, 
			{0.398694, 0.0681734, 0.941842, 2.15817}, 
			{8.27331, 1.33459, 3.11974, 1.79294}, 
			{36.1691, 1.99678, 4.26504, 1.69565}, 
			{19.9958, 1.7495, 3.78158, 1.72287}, 
			{2.23154, 0.763394, 2.08861, 1.93595}, 
			{4.6589, 1.13215, 2.4543, 1.88867}, 
			{31.7077, 1.93301, 4.02469, 1.73183}, 
			{1.6939, 0.688397, 1.66044, 2.05487}, 
			{1.12861, 0.526251, 1.40021, 2.09187}, 
			{491.072, 2.95971, 10.6222, 1.15221}, 
			{2.09251, 0.809091, 1.62012, 2.1177}, 
			{2.98916, 0.943473, 1.87156, 2.06668}, 
			{113.503, 2.43144, 6.50442, 1.43207}
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