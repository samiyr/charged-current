// #ifndef PDF_INTERFACE_H
// #define PDF_INTERFACE_H

// #include "Flavor.cpp"
// #include "CKM.cpp"
// #include "Process.cpp"



// class PDFInterface {
// 	public:
// 	virtual void evaluate(const double x, const double Q2) {
// 		throw std::runtime_error("Cannot call evaluate on PDFInterface directly");
// 	}
// 	virtual double xf_evaluate(const FlavorType flavor, const double x, const double Q2) {
// 		throw std::runtime_error("Cannot call xf_evaluate on PDFInterface directly");
// 	}
// 	virtual double xf(const FlavorType flavor) const {
// 		throw std::runtime_error("Cannot call xf on PDFInterface directly");
// 	}
// 	virtual double alpha_s(const double Q2) const {
// 		throw std::runtime_error("Cannot call alpha_s on PDFInterface directly");
// 	}

// 	double xf_sum(const FlavorVector &flavors, const FlavorVector &reflected) {
// 		double sum = 0;
// 		for (const auto flavor : flavors) {
// 			const double ckm_contribution = CKM::sum_of_squares(reflected, flavor);
// 			sum += ckm_contribution * xf(flavor);
// 		}
// 		return sum;
// 	}
// 	double xg() {
// 		return xf(Flavor::Gluon);
// 	}
// 	double xg_sum(const FlavorVector &flavors, const FlavorVector &reflected) {
// 		double ckm_sum = 0;
// 		for (const auto flavor : flavors) {
// 			const double ckm_contribution = CKM::sum_of_squares(reflected, flavor);
// 			ckm_sum += ckm_contribution;
// 		}
// 		return 2 * ckm_sum * xg();
// 	}
// 	double xq_sum(const FlavorVector &upper_flavors, 
// 	const FlavorVector &lower_flavors, 
// 	const FlavorVector &upper_antiflavors, 
// 	const FlavorVector &lower_antiflavors, 
// 	const bool quark_minus, 
// 	const Process process) {
// 		const bool positive_W = process.positive_W();
// 		const double xq_sum = positive_W ? xf_sum(lower_flavors, upper_flavors) : xf_sum(upper_flavors, lower_flavors);
// 		const double anti_xq_sum = positive_W ? xf_sum(upper_antiflavors, lower_antiflavors) : xf_sum(lower_antiflavors, upper_antiflavors);
// 		const double xq = quark_minus ? xq_sum - anti_xq_sum : xq_sum + anti_xq_sum;
// 		return xq;
// 	}
// };

// #endif