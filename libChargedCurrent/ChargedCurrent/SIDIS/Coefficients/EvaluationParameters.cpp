#ifndef SIDIS_FUNCTIONS_EVALUATION_PARAMETERS_H
#define SIDIS_FUNCTIONS_EVALUATION_PARAMETERS_H

namespace SIDISFunctions {
	struct EvaluationParameters {
		const double xi;
		const double xip; 
		const double x;
		const double z;

		const double pdf;
		const double pdf_hat;

		const double ff;
		const double ff_hat;

		const double renormalization_scale_log;
		const double factorization_scale_log;
		const double fragmentation_scale_log;

		const double m2;
		const double Q2;
	};
}

#endif