#ifndef CKM_H
#define CKM_H

#include "Flavor.cpp"
#include "Constants.cpp"

struct CKM {
	private:
	constexpr static double _ckm_full_matrix[36] = {
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
	constexpr static double squared(const FlavorType flavor1, FlavorType flavor2) {
		const size_t f1 = static_cast<size_t>(std::abs(flavor1));
		const size_t f2 = static_cast<size_t>(std::abs(flavor2));
		const size_t row = f1 - 1;
		const size_t column = f2 - 1;
		const size_t column_count = 6; 
		return _ckm_full_matrix[row * column_count + column];
	}
};


#endif