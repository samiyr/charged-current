#ifndef SIDIS_FUNCTIONS_EVALUATION_PARAMETERS_H
#define SIDIS_FUNCTIONS_EVALUATION_PARAMETERS_H

namespace SIDISFunctions {
	struct EvaluationParameters {
		const double xi;
		const double xip; 
		const double x;
		const double z;

		const double renormalization_scale_log;
		const double factorization_scale_log;
		const double fragmentation_scale_log;

		const double xq;
		const double xq_hat;
		const double xg_hat;

		const double zq;
		const double zq_hat;
		const double zg_hat;

		const double anti_xq;
		const double anti_xq_hat;

		const double anti_zq;
		const double anti_zq_hat;

		const double sign;

		// const double xq_zq;
		// const double xq_hat_zq; 
		// const double xq_zq_hat;
		// const double xq_hat_zq_hat;
		// const double xq_zg_hat;
		// const double xg_hat_zq;
		// const double xq_hat_zg_hat;
		// const double xg_hat_zq_hat;

		// const double log1mx;
		// const double log1mz;
		// const double logxi;
		// const double logxip;
		// const double log1mxi;
		// const double log1mxip;

		const double m2;
		const double Q2;

		const double pdf;
		const double pdf_hat;

		const double ff;
		const double ff_hat;
	};
}

#endif