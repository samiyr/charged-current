#ifndef LHASET_INTERFACE_H
#define LHASET_INTERFACE_H

#include <string>

#include "LHAPDF/LHAPDF.h"

#include "PDF/Interfaces/LHAInterface.cpp"
#include "PDF/PDFConcept.cpp"

template <typename explicit_isospin = std::false_type, std::derived_from<LHAPDF::Extrapolator> Extrapolator = ZeroExtrapolator>
class LHASetInterface {
	public:
	using size_type = std::size_t;

	const std::string set_name;

	double Z = 1.0;
	double A = 1.0;

	private:
	LHASetInterface(
		const std::string set_name,		
		const bool _use_multipliers, 
		const std::vector<double> _multipliers,
		const bool _use_global_multiplier = false,
		const double _global_multiplier = 1.0) noexcept 
		: set_name(set_name), 
		use_multipliers(_use_multipliers), 
		multipliers(_multipliers), 
		use_global_multiplier(_use_global_multiplier), 
		global_multiplier(_global_multiplier) {
			
		std::unique_ptr<LHAPDF::PDFInfo> info(LHAPDF::mkPDFInfo(set_name, 0));
		member_count = info->get_entry_as<size_type>("NumMembers");
	}

	public:
	LHASetInterface(std::string _set_name, const std::vector<double> _multipliers) noexcept
	: LHASetInterface(_set_name, _multipliers.size() == TOTAL_FLAVORS, _multipliers) { }

	LHASetInterface(std::string _set_name) noexcept
	: LHASetInterface(_set_name, false, {}) { }

	LHAInterface<explicit_isospin, Extrapolator> operator[](const size_type member) const {
		LHAInterface<explicit_isospin, Extrapolator> pdf_member(set_name, member, use_multipliers, multipliers, use_global_multiplier, global_multiplier);
		pdf_member.Z = Z;
		pdf_member.A = A;

		return pdf_member;
	}

	LHAInterface<explicit_isospin, Extrapolator> central() const {
		return this->operator[](0);
	}

	double quark_mass(const FlavorType flavor) const {
		return central().quark_mass(flavor);
	}

	double Q2_min() const {
		return central().Q2_min();
	}
	
	double Q2_max() const {
		return central().Q2_max();
	}

	class iterator {
		public:
		iterator(const int member, const LHASetInterface<explicit_isospin, Extrapolator> *set) : member(member), set(set) {}
		iterator operator++() { member++; return *this; }
		bool operator!=(const iterator &other) { return member != other.member; }
		const LHAInterface<explicit_isospin, Extrapolator> operator*() const { return set->operator[](member); }

		private:
		int member;
		const LHASetInterface<explicit_isospin, Extrapolator> *set;
	};

	size_type size() const {
		return member_count;
	}

	iterator begin() const {
		return iterator(0, this);
	}

	iterator end() const {
		return iterator(member_count - 1, this);
	}

	private:
	size_type member_count;

	const bool use_multipliers;
	const std::vector<double> multipliers;

	bool use_global_multiplier = false;
	double global_multiplier = 1.0;
};

#endif