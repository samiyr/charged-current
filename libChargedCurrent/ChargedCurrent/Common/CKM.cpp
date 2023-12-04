#ifndef CKM_H
#define CKM_H

#include <array>

#include "Common/Flavor.cpp"
#include "Common/Constants.cpp"

// Convenient access to the squared CKM matrix elements.
struct CKM {
	private:
	// The 2D CKM matrix in a 1D form in such a way that invalid elements, such as V_dd or V_ds, which don't belong in the CKM matrix, are zero. 
	// The elements are all squared.
	constexpr static std::array<double, 36> _ckm_full_matrix = {
		/* d: 				(d) 					u 						(s) 					c 						(b) 					t 				*/
							0,				Constants::V_ud,				0,				Constants::V_cd,				0,				Constants::V_td	,
		/* u: 				d 						(u)						s 						(c)						b 						(t) 			*/
					Constants::V_ud,				0,				Constants::V_us,				0,				Constants::V_ub,				0		,
		/* s:				(d) 					u 						(s) 					c 						(b) 					t 				*/
							0,				Constants::V_us,				0,				Constants::V_cs,				0,				Constants::V_ts	,
		/* c: 				d 						(u)						s 						(c)						b 						(t) 			*/
					Constants::V_cd,				0,				Constants::V_cs,				0,				Constants::V_cb,				0		,
		/* b: 				(d) 					u 						(s) 					c 						(b) 					t 				*/
							0,				Constants::V_ub,				0,				Constants::V_cb,				0,				Constants::V_tb	,
		/* t: 				d 						(u)						s 						(c)						b 						(t) 			*/
					Constants::V_td,				0,				Constants::V_ts,				0,				Constants::V_tb,				0		
	};
	public:
	// Returns the squared CKM matrix element corresponding to the two flavors. The order of the flavors doesn't matter, and they can also be antiflavors
	// (i.e. with a negative sign). If an invalid flavor combination is given, such as dd or ds, returns zero. If flavors out of the range [-6, -1] + [1, 6]
	// are provided, the result is undefined behaviour.
	constexpr static double squared(const FlavorType flavor1, FlavorType flavor2) {
		const std::size_t f1 = static_cast<std::size_t>(std::abs(flavor1));
		const std::size_t f2 = static_cast<std::size_t>(std::abs(flavor2));
		const std::size_t row = f1 - 1;
		const std::size_t column = f2 - 1;
		const std::size_t column_count = 6; 
		return _ckm_full_matrix[row * column_count + column];
	}
};


#endif