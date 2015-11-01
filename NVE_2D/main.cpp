#include <iostream>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <tuple>
#include <limits>
#include "RNG.hpp"
#include "Particle.hpp"
#include "Vector2D.hpp" 
#include "Quantities.hpp"

using namespace std;

typedef ParticleGeneral<Vector2D> Particle;

double CalculateTemperature(const vector<Particle>& p) {

	const double K = CalculateKinEnergy(p);
	return 2.0*K/(3.0*p.size());

}

vector<Vector2D> CalculateInteraction(const vector<Particle>& p, double box) {

	vector<Vector2D> forc(p.size());

	const int n = p.size();
	for (int i=0; i<(n-1); i++) {
		for (int j=i+1; j<n; j++) {
			Vector2D d = p[i].r-p[j].r;
			d = CorrectBox(box, d);
			const double r2 = d.MagnitudeSquared();
			const double r6_ = 1.0/(r2*r2*r2);
			const double f = 24*(2.0*r6_*r6_ - r6_)/r2;
			d *= f;

			forc[i] += d;			
			forc[j] -= d;			
		}
	}
	return forc;
}

void InitializeSystem(vector<Particle>& p, double box, int nSideX, int nSideY) {

	int N=nSideX*nSideY;
	p.resize(N);
	RNG rng;

	const Vector2D centerer = Vector2D(-0.5,-0.5)*box;
	const double stepX = box/nSideX;
	const double stepY = box/nSideY;
	int idx = 0;
	for (int i=0; i<nSideX; i++) {
		for (int j=0; j<nSideY; j++) {
			p[idx].r = Vector2D(i*stepX,j*stepY) + centerer;
			p[idx].v = Vector2D(rng.UniCentered(), rng.UniCentered());
			idx++;
		}
	}
}

void UpdateSystem(vector<Particle>& p, double dt, double box) {

	const int n = p.size();

	auto f = CalculateInteraction(p, box);
	for (int i=0; i<n; i++) {
		p[i].v += f[i]*dt;
		p[i].r = CorrectBox(box, p[i].r + p[i].v*dt);
	}

}

int main() {

	int nSideX = 3;
	int nSideY = 3;
	vector<Particle> p;
	double box = 10.0;
	double dt = 1e-5;

	InitializeSystem(p,box, nSideX, nSideY);

	ofstream initial("initial.dat");
	for (auto a:p) {
		initial << a.r << endl;
	}	

	ofstream trans("trans.dat");
	ofstream temp("temp.dat");
	ofstream energy("energy.dat");

	for (int step =0;step <5e6; step++) {
		UpdateSystem(p, dt, box); 
		if (step %2000 == 0) {
			cout << step << endl;
			double v2 =0;
			for (auto a:p) {
				trans << a.r << endl;
				v2+=a.SumVels2();
			}
			temp << step << " " << CalculateTemperature(p) << endl;
			energy << step << " " << CalculateEnergy<Particle, Vector2D>(p,box) << endl;
		}
	}

	ofstream fin("final.dat");
	for (auto a:p) {
		fin << a.r << a.v << endl;
	}	
}
