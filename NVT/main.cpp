#include <iostream>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <tuple>
#include <limits>
#include "RNG.hpp"
#include "Particle.hpp"
#include "ExtVar.hpp"
#include "Vector3D.hpp" 
#include "Quantities.hpp"

using namespace std;

typedef ParticleGeneral<Vector3D> Particle;

double CalculateTemperature(const vector<Particle>& p) {

	const double K = CalculateKinEnergy(p);
	return 2.0*K/(3.0*p.size());

}

vector<Vector3D> CalculateInteraction(const vector<Particle>& p, double box) {

	vector<Vector3D> forc(p.size());

	const int n = p.size();
	for (int i=0; i<(n-1); i++) {
		for (int j=i+1; j<n; j++) {
			Vector3D d = p[i].r-p[j].r;
			d = CorrectBox(box, d);
			const double r2 = d.MagnitudeSquared();
			const double r6_ = 1.0/(r2*r2*r2);
			const double f = 12*(2.0*r6_*r6_ - r6_)/r2;
			d *= f;

			forc[i] += d;			
			forc[j] -= d;			
		}
	}
	return forc;
}

void InitializeSystem(vector<Particle>& p, ExtVar& s, double box, double temperature, int nSideX, int nSideY, int nSideZ) {

	int N=nSideX*nSideY*nSideZ;
	p.resize(N);
	RNG rng;

	s.r = 0.0;
	s.v = 0.0;
	s.m = 1.0;

	const Vector3D centerer = Vector3D(-0.5,-0.5, -0.5)*box;
	const double stepX = box/nSideX;
	const double stepY = box/nSideY;
	const double stepZ = box/nSideZ;
	int idx = 0;
	for (int i=0; i<nSideX; i++) {
		for (int j=0; j<nSideY; j++) {
			for (int k=0; k<nSideZ; k++) {
				p[idx].r = Vector3D(i*stepX,j*stepY,k*stepZ) + centerer;
				p[idx].v = Vector3D(rng.UniCentered(),rng.UniCentered(),rng.UniCentered())*sqrt(3*temperature);
				/* p[idx].v = Vector3D(1.0, 1.0,1.0)*temperature; */
				idx++;
			}
		}
	}
}

void UpdateSystem(vector<Particle>& p, ExtVar& s, double dt, double box, double temperature) {

	const int n = p.size();

	double p1 = -3*n*temperature;
	double p2 = 0;
	for (auto a:p) {
		p2 += a.SumVels2();
	}
	s.v = (p1 + p2)/s.m;
	s.r += s.v*dt;

	auto f = CalculateInteraction(p, box);
	for (int i=0; i<n; i++) {
		p[i].v += (f[i] - p[i].v*s.r)*dt;
		p[i].r = CorrectBox(box, p[i].r + p[i].v*dt);
	}

}

int main() {

	int nSideX = 5;
	int nSideY = 5;
	int nSideZ = 5;
	vector<Particle> p;
	ExtVar s;
	double box = 10.0;
	double dt = 0.0001;
	double temperature = 1.00;

	InitializeSystem(p,s,box, temperature, nSideX, nSideY, nSideZ);

	ofstream initial("initial.dat");
	for (auto a:p) {
		initial << a.r << endl;
	}	

	ofstream trans("trans.dat");
	ofstream temp("temp.dat");
	ofstream frict("friction.dat");
	ofstream energy("energy.dat");

	for (int step =0;step <5e6; step++) {
		UpdateSystem(p, s, dt, box, temperature); 
		if (step %2000 == 0) {
			cout << step << endl;
			double v2 =0;
			for (auto a:p) {
				trans << a.r << endl;
				v2+=a.SumVels2();
			}
			temp << step << " " << CalculateTemperature(p) << endl;
			energy << step << " " << CalculateEnergy<Particle, Vector3D>(p,box) << endl;
			frict << step << " " << s.r << endl;
		}
	}

	ofstream fin("final.dat");
	for (auto a:p) {
		fin << a.r << a.v << endl;
	}	
}
