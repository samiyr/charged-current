#ifndef ANALYSIS_CONSTANTS_H
#define ANALYSIS_CONSTANTS_H

#include <vector>

namespace AnalysisConstants {
	namespace NuTeV {
		namespace Old {
			const static std::vector<double> y_bins = {0.334, 0.573, 0.790};
			const static std::vector<double> E_bins = {90.18, 174.37, 244.72};
		}
		namespace New {
			const static std::vector<double> y_bins = {0.324, 0.558, 0.771};
			const static std::vector<double> E_bins = {88.29, 174.29, 247.0};
		}
	}
	namespace CCFR {
		const static std::vector<double> y_bins = {0.32, 0.57, 0.795};
		const static std::vector<double> E_bins = {109.46, 209.89, 332.7};
	}
}

enum class AnalysisSet {
	NuTeV, NuTeV_old, CCFR
};

#endif