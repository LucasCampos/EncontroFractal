#ifndef QUANTITIES_HPP
#define QUANTITIES_HPP

#include <vector>
#include "Particle.hpp"

template <class Particle>
double CalculateKinEnergy(const std::vector<Particle>& p) {
	double Ek = 0;
	for (auto a: p)
		Ek += a.SumVels2();
	Ek = Ek*0.5;
	return Ek;
}

template <class Particle, class Vector>
double CalculatePotEnergy(const std::vector<Particle>& p, double box) {
	const int n = p.size();
	double Ep =0;
	for (int i=0; i<(n-1); i++) {
		for (int j=i+1; j<n; j++) {
			Vector d = p[i].r-p[j].r;
			d = CorrectBox(box, d);
			const double r2 = d.MagnitudeSquared();
			const double r6_ = 1.0/(r2*r2*r2);
			const double e = 4.0*(r6_*r6_ - r6_);
			Ep += e;
		}
	}
	return Ep;
}

template <class Particle, class Vector>
double CalculateEnergy(const std::vector<Particle>& p, double box) {
	double Ek = CalculateKinEnergy(p);
	double Ep = CalculatePotEnergy<Particle,Vector>(p,box);
	return Ep + Ek;
}
#endif
