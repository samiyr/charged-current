#ifndef INTEGRATION_PARAMETERS_H
#define INTEGRATION_PARAMETERS_H

#include <iostream>

struct IntegrationParameters {
	struct Cuba {
		struct Vegas {
			int starting_iterations = 1000;
			int continuing_iterations = 1000;
			int batch_size = 1000;
			int grid_number = 0;

			friend std::ostream &operator <<(std::ostream &os, const Vegas &vegas) {
				os << "#cuba.vegas.starting_iterations = " 				<< vegas.starting_iterations 			<< IO::endl;
				os << "#cuba.vegas.continuing_iterations = " 			<< vegas.continuing_iterations 			<< IO::endl;
				os << "#cuba.vegas.stabatch_sizerting_iterations = " 	<< vegas.batch_size 					<< IO::endl;
				os << "#cuba.vegas.grid_number = " 						<< vegas.grid_number 					<< IO::endl;
				return os;
			}
		};

		struct Suave {
			int new_subdivision_evaluations = 500;
			int minimum_samples = 20;
			double flatness = 100.0;

			friend std::ostream &operator <<(std::ostream &os, const Suave &suave) {
				os << "#cuba.suave.new_subdivision_evaluations = " 		<< suave.new_subdivision_evaluations 	<< IO::endl;
				os << "#cuba.suave.minimum_samples = " 					<< suave.minimum_samples 				<< IO::endl;
				os << "#cuba.suave.flatness = " 						<< suave.flatness 						<< IO::endl;
				return os;
			}
		};

		struct Cuhre {
			int cubature_rule = 0;

			friend std::ostream &operator <<(std::ostream &os, const Cuhre &cuhre) {
				os << "#cuba.cuhre.cubature_rule = " 					<< cuhre.cubature_rule 					<< IO::endl;
				return os;
			}
		};

		Vegas vegas;
		Suave suave;
		Cuhre cuhre;

		int seed = 0;
		int minimum_evaluations = 0;
		int maximum_evaluations = 1'000'000;

		double maximum_relative_error = 1e-2;
		double maximum_absolute_error = 1e-12;

		bool use_all_samples = true;
		bool apply_smoothing = true;
		bool delete_state_file_after_integration = true;
		bool reset_state = false;
		int random_number_generator_level = 0;
	};

	struct GSL {
		std::size_t points = 200'000;
		double max_chi_squared_deviation = 0.2;
		double max_relative_error = 1e-3;
		unsigned int iter_max = 5;

		bool grid_warmup = true;

		friend std::ostream &operator <<(std::ostream &os, const GSL &gsl) {
			os << "#gsl.points = " 						<< gsl.points 						<< IO::endl;
			os << "#gsl.max_chi_squared_deviation = " 	<< gsl.max_chi_squared_deviation 	<< IO::endl;
			os << "#gsl.max_relative_error = " 			<< gsl.max_relative_error 			<< IO::endl;
			os << "#gsl.iter_max = " 					<< gsl.iter_max 					<< IO::endl;
			return os;
		}
	};

	Cuba cuba;
	GSL gsl;

	friend std::ostream &operator <<(std::ostream &os, const IntegrationParameters &p) {
		os << "#cuba.seed = "									<< p.cuba.seed 										<< IO::endl;
		os << "#cuba.minimum_evaluations = "					<< p.cuba.minimum_evaluations 						<< IO::endl;
		os << "#cuba.maximum_evaluations = "					<< p.cuba.maximum_evaluations 						<< IO::endl;

		os << "#cuba.maximum_relative_error = "					<< p.cuba.maximum_relative_error 					<< IO::endl;
		os << "#cuba.maximum_absolute_error = "					<< p.cuba.maximum_absolute_error 					<< IO::endl;
		
		os << "#cuba.use_all_samples = "						<< p.cuba.use_all_samples 							<< IO::endl;
		os << "#cuba.apply_smoothing = "						<< p.cuba.apply_smoothing 							<< IO::endl;
		os << "#cuba.delete_state_file_after_integration = "	<< p.cuba.delete_state_file_after_integration 		<< IO::endl;
		os << "#cuba.reset_state = "							<< p.cuba.reset_state 								<< IO::endl;
		os << "#cuba.random_number_generator_level = "			<< p.cuba.random_number_generator_level 			<< IO::endl;

		os << p.cuba.vegas;
		os << p.cuba.suave;
		os << p.cuba.cuhre;
		os << p.gsl;

		return os;
	}
};

#endif