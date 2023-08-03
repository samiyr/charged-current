#ifndef PYTHIA_SIDIS_H
#define PYTHIA_SIDIS_H

#include "Pythia8/Pythia.h"
#include "Pythia8/PythiaParallel.h"
#include <vector>
#include <boost/format.hpp>
#include <boost/histogram.hpp>
#include <iostream>
#include <string>


namespace PythiaSIDIS {
	void cross_section(const std::size_t events, const std::string filename) {
		using namespace Pythia8;
		using namespace boost::histogram;
		PythiaParallel pythia;

		pythia.readString("Beams:eCM = 318.");
		pythia.readString("Beams:idA = 2212");
		pythia.readString("Beams:idB = 13");
		pythia.readString("WeakBosonExchange:ff2ff(t:W) = on");
		pythia.settings.parm("PhaseSpace:Q2Min", 10.0);
		pythia.readString("SpaceShower:dipoleRecoil = on");
		pythia.readString("SpaceShower:pTmaxMatch = 2");
		pythia.readString("PDF:lepton = off");
		pythia.readString("TimeShower:QEDshowerByL = off");
		pythia.readString("Next:numberCount = 100000");

		pythia.init();

		auto hist = make_histogram(axis::regular<>(40, 0.01, 1.0, "x"), axis::regular<>(40, 0.01, 1.0, "z"), axis::regular<>(9, 10.0, 100.0, "Q2"));

		pythia.run(events, [&](Pythia *ptr) {
			Event &event = ptr->event;

			std::vector<Particle> D_mesons;

			for (std::size_t i_particle = 1; i_particle < event.size(); i_particle++) {
				const auto particle = event[i_particle];

				if (particle.idAbs() == 421) {
					D_mesons.push_back(particle);
				}
			}

			if (D_mesons.empty()) { return; }

			const Vec4 hadron_momentum = event[1].p();
			const Vec4 lepton_in_momentum = event[2].p();
			const Vec4 lepton_out_momentum = event[6].p();
			const Vec4 q = lepton_in_momentum - lepton_out_momentum;

			const double Q2 = -q.m2Calc();
			const double x = Q2 / (2.0 * hadron_momentum * q);

			for (auto const &particle : D_mesons) {
				const Vec4 D_momentum = particle.p();
				const double z = (hadron_momentum * D_momentum) / (hadron_momentum * q);

				// std::cout << x << ", " << z << ", " << Q2 << IO::endl;

				hist(x, z, Q2);
			}
		});

		const double factor = pythia.sigmaGen() / pythia.weightSum();

		std::ofstream file(filename);

		for (auto &&x : indexed(hist)) {
			const double value = *x * factor / (x.bin(0).width() * x.bin(1).width() * x.bin(2).width());
			file << x.bin(0).center() << ", " << x.bin(1).center() << ", " << x.bin(2).center() << ", " << value << IO::endl;
		}

		file.close();
	}
}

#endif