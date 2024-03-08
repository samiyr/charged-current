#ifndef DECAY_PARAMETRIZATION_H
#define DECAY_PARAMETRIZATION_H

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "Utility/Math.cpp"
#include "Utility/Utility.cpp"

struct DecayParametrization {
	constexpr static DecayParametrization fit1() noexcept {
		return DecayParametrization(3.1574, 1.00882, 1.79948, 2.1012);
	}
	constexpr static DecayParametrization fit2() noexcept {
		return DecayParametrization(7.365, 1.4, 2.276, 2.04);
	}

	constexpr static std::vector<DecayParametrization> fit_set_1() noexcept {
		const std::vector<std::vector<double>> parameter_set = {
			{3.1574, 1.00882, 1.79948, 2.1012}, 
			{1.0852, 0.522549, 1.2748, 2.16643}, 
			{2.99367, 0.990715, 1.7187, 2.1248}, 
			{2.27886, 0.858437, 1.63605, 2.12405}, 
			{2.01257, 0.808877, 1.52716, 2.15508}, 
			{3.67072, 1.0876, 1.79766, 2.12218}, 
			{2.2159, 0.850715, 1.59706, 2.1352}, 
			{1.47918, 0.695579, 1.31615, 2.1905}, 
			{1.40373, 0.614203, 1.40683, 2.16624},
			{1.75786, 0.735873, 1.46332, 2.16339}, 
			{2.42844, 0.882453, 1.64343, 2.13974}, 
			{2.2414, 0.820832, 1.65875, 2.12311}, 
			{2.80787, 0.944272, 1.76284, 2.09367}, 
			{2.20082, 0.816352, 1.72669, 2.0847}, 
			{2.23345, 0.835243, 1.63034, 2.12873}, 
			{1.91023, 0.760071, 1.59017, 2.12779}, 
			{1.40097, 0.646527, 1.36797, 2.17189}, 
			{4.03962, 1.10974, 1.99722, 2.055}, 
			{3.31209, 1.03005, 1.7951, 2.10459}, 
			{3.0577, 0.995663, 1.75082, 2.11586}, 
			{1.6773, 0.729657, 1.42702, 2.17119}, 
			{1.91028, 0.753372, 1.62664, 2.10887}, 
			{4.77146, 1.14158, 2.09896, 2.06351}, 
			{3.26849, 1.03133, 1.77903, 2.10798}, 
			{3.39915, 1.02117, 1.86321, 2.09542}, 
			{2.19529, 0.802467, 1.68457, 2.11803}, 
			{2.58173, 0.925979, 1.68291, 2.10991}, 
			{2.02873, 0.788068, 1.62485, 2.11086}, 
			{1.37372, 0.631282, 1.33972, 2.18393}, 
			{1.31275, 0.614161, 1.31491, 2.19367}, 
			{1.55543, 0.687687, 1.42167, 2.1585}
		};

		std::vector<DecayParametrization> parametrizations;

		for (const std::vector<double> &parameters : parameter_set) {
			parametrizations.emplace_back(parameters[0], parameters[1], parameters[2], parameters[3]);
		}

		return parametrizations;
	}

	constexpr static std::vector<DecayParametrization> fit_set_2() noexcept {
		const std::vector<std::vector<double>> parameter_set = {
			{3.1574, 1.00882, 1.79948, 2.1012}, 
			{3.74331, 1.09564, 1.86571, 2.09676}, 
			{1.75362, 0.726757, 1.52309, 2.14312}, 
			{1.70545, 0.738959, 1.42823, 2.16847}, 
			{1.3182, 0.620801, 1.2853, 2.19745}, 
			{3.36087, 1.03892, 1.86266, 2.08902}, 
			{3.91863, 1.10028, 1.91145, 2.09615}, 
			{1.83727, 0.761498, 1.53056, 2.14285}, 
			{2.05715, 0.810314, 1.61141, 2.12249}, 
			{2.47936, 0.881724, 1.7066, 2.11219}, 
			{2.23254, 0.863187, 1.58721, 2.14126}, 
			{1.8297, 0.77265, 1.41263, 2.19573}, 
			{2.58339, 0.936175, 1.66147, 2.12457}, 
			{2.0836, 0.863185, 1.44137, 2.18807}, 
			{2.35467, 0.866441, 1.6876, 2.10544}, 
			{0.860502, 0.424842, 1.02071, 2.29071}, 
			{3.8822, 1.10648, 1.9417, 2.07152}, 
			{3.04734, 0.996381, 1.71923, 2.12696}, 
			{4.4568, 1.19932, 1.92331, 2.10618}, 
			{1.85753, 0.787908, 1.45047, 2.16843}, 
			{2.86038, 0.938158, 1.8137, 2.09204}, 
			{6.03848, 1.32499, 2.14404, 2.06105}, 
			{1.21365, 0.537383, 1.35656, 2.17465}, 
			{2.333, 0.895201, 1.59601, 2.12666}, 
			{2.34758, 0.877585, 1.65023, 2.12143}, 
			{4.08599, 1.15672, 1.84725, 2.11355}, 
			{3.25639, 1.01466, 1.85446, 2.08897}, 
			{1.33516, 0.628823, 1.29621, 2.19018}, 
			{2.7111, 0.949754, 1.6512, 2.14458}, 
			{4.98269, 1.21467, 2.06849, 2.06784}, 
			{3.60883, 1.06823, 1.90357, 2.07327}
		};

		std::vector<DecayParametrization> parametrizations;

		for (const std::vector<double> &parameters : parameter_set) {
			parametrizations.emplace_back(parameters[0], parameters[1], parameters[2], parameters[3]);
		}

		return parametrizations;
	}

	static std::vector<DecayParametrization> fit_set_3() {
		std::ifstream file("fit_set_3.dat");
		if (!file) {
			throw std::runtime_error("Fit set file 'fit_set_3.dat' could not be opened.");
		}

		std::vector<std::vector<double>> parameter_set;

		std::size_t current_index = 0;

		std::string line;
		while (std::getline(file, line)) {
			current_index++;

			if (current_index == 1) { continue; }

			std::vector<double> parameters;

			std::istringstream stream(line);
			std::string value;

			while (std::getline(stream, value, ',')) {
				parameters.push_back(std::stod(value));
			}

			parameter_set.push_back(parameters);
		}

		std::vector<DecayParametrization> parametrizations;

		for (const std::vector<double> &parameters : parameter_set) {
			std::cout << parameters[0] << IO::endl;
			parametrizations.emplace_back(parameters[0], parameters[1], parameters[2], parameters[3]);
		}

		return parametrizations;
	}
	static std::vector<DecayParametrization> fit_set_4() {
		std::ifstream file("fit_set_4.dat");
		if (!file) {
			throw std::runtime_error("Fit set file 'fit_set_4.dat' could not be opened.");
		}

		std::vector<std::vector<double>> parameter_set;

		std::size_t current_index = 0;

		std::string line;
		while (std::getline(file, line)) {
			current_index++;

			if (current_index == 1) { continue; }

			std::vector<double> parameters;

			std::istringstream stream(line);
			std::string value;

			while (std::getline(stream, value, ',')) {
				parameters.push_back(std::stod(value));
			}

			parameter_set.push_back(parameters);
		}

		std::vector<DecayParametrization> parametrizations;

		for (const std::vector<double> &parameters : parameter_set) {
			std::cout << parameters[0] << IO::endl;
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