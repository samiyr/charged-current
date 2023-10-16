#ifndef ANALYSIS_CONSTANTS_H
#define ANALYSIS_CONSTANTS_H

#include <vector>

#include "Common/Process.cpp"

enum class AnalysisSet {
	NuTeV, NuTeV_old, CCFR, NOMAD
};

namespace AnalysisConstants {
	namespace NuTeV {
		namespace Old {
			namespace Neutrino {
				const static std::vector<double> y_bins = {0.334, 0.573, 0.790};
				const static std::vector<double> E_bins = {90.18, 174.37, 244.72};
			}
			namespace Antineutrino {
				const static std::vector<double> y_bins = {0.356, 0.586, 0.788};
				const static std::vector<double> E_bins = {78.98, 146.06, 222.14};
			}
		}
		namespace New {
			namespace Neutrino {
				const static std::vector<double> y_bins = {0.324, 0.558, 0.771};
				const static std::vector<double> E_bins = {88.29, 174.29, 247.0};
			}
			namespace Antineutrino {
				const static std::vector<double> y_bins = {0.349, 0.579, 0.776};
				const static std::vector<double> E_bins = {77.9, 143.7, 226.8};
			}
		}
	}
	namespace CCFR {
		namespace Neutrino {
			const static std::vector<double> y_bins = {0.32, 0.57, 0.795};
			const static std::vector<double> E_bins = {109.46, 209.89, 332.7};
		}
		namespace Antineutrino {
			const static std::vector<double> y_bins = {0.355, 0.596, 0.802};
			const static std::vector<double> E_bins = {87.40, 160.52, 265.76};
		}
	}
	namespace NOMAD {
		const static std::vector<double> E_bins = {
			15.91, 24.38, 28.85, 32.88, 37.31, 41.78, 46.23, 51.17, 56.73, 62.87, 69.70, 77.29, 85.78, 95.01, 105.4, 117.6, 133.0, 155.4, 205.5
		};
	}

	inline std::vector<double> get_y_bins(AnalysisSet set, Process process) {
		switch (set) {
		case AnalysisSet::NuTeV:
			return process.type == Process::Type::NeutrinoToLepton ? NuTeV::New::Neutrino::y_bins : NuTeV::New::Antineutrino::y_bins;		
		case AnalysisSet::NuTeV_old:
			return process.type == Process::Type::NeutrinoToLepton ? NuTeV::Old::Neutrino::y_bins : NuTeV::Old::Antineutrino::y_bins;		
		case AnalysisSet::CCFR:
			return process.type == Process::Type::NeutrinoToLepton ? CCFR::Neutrino::y_bins : CCFR::Antineutrino::y_bins;
		default: return {};
		}
	}
	inline std::vector<double> get_E_bins(AnalysisSet set, Process process) {
		switch (set) {
		case AnalysisSet::NuTeV:
			return process.type == Process::Type::NeutrinoToLepton ? NuTeV::New::Neutrino::E_bins : NuTeV::New::Antineutrino::E_bins;		
		case AnalysisSet::NuTeV_old:
			return process.type == Process::Type::NeutrinoToLepton ? NuTeV::Old::Neutrino::E_bins : NuTeV::Old::Antineutrino::E_bins;		
		case AnalysisSet::CCFR:
			return process.type == Process::Type::NeutrinoToLepton ? CCFR::Neutrino::E_bins : CCFR::Antineutrino::E_bins;
		case AnalysisSet::NOMAD:
			return NOMAD::E_bins;
		default: return {};
		}
	}
}

#endif