#include <iostream>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <tuple>
#include <limits>
#include "RNG.hpp"
#include "Particle.hpp"
#include "ExtVar.hpp"
#include "Vector2D.hpp" 

using namespace std;

typedef ParticleGeneral<Vector2D> Particle;

vector<Vector2D> CalculateInteraction(const vector<Particle>& p, double box) {

	vector<Vector2D> forc(p.size());

	const int n = p.size();
	for (int i=0; i<(n-1); i++) {
		for (int j=i+1; j<n; j++) {
			Vector2D d = p[i].r-p[j].r;
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

void InitializeSystem(vector<Particle>& p, ExtVar& s, double box, double temperature, int nSideX, int nSideY) {

	int N=nSideX*nSideY;
	p.resize(N);
	RNG rng;

	s.r = 0.0;
	s.v = 0.0;
	s.m = 1.0;

	const Vector2D centerer = Vector2D(-0.5,-0.5)*box;
	const double stepX = box/nSideX;
	const double stepY = box/nSideY;
	int idx = 0;
	for (int i=0; i<nSideX; i++) {
		for (int j=0; j<nSideY; j++) {
			p[idx].r = Vector2D(i*stepX,j*stepY) + centerer;
			p[idx].v = Vector2D(rng.UniCentered(),rng.UniCentered())*sqrt(2*temperature);
			/* p[idx].v = Vector2D(1.0, 1.0,1.0)*temperature; */
			idx++;
		}
	}
}

void UpdateSystem(vector<Particle>& p, ExtVar& s, double dt, double box, double temperature) {

	const int n = p.size();

	double p1 = -2*n*temperature;
	double p2 = 0;
	for (auto a:p) {
		p2 += a.SumVels2();
	}
	s.v = (p1 + p2)/s.m;
	s.r += s.v*dt;
	/* cout << "ds/dt = " << p1 << " + " << p2 << endl; */
	/* cout << s.r << endl; */
	auto f = CalculateInteraction(p, box);
	for (int i=0; i<n; i++) {
		/* cout << p[i].v << endl; */
		p[i].v += (f[i] - p[i].v*s.r)*dt;
		p[i].r = CorrectBox(box, p[i].r + p[i].v*dt);

	}
	/* cout << s.r << endl; */
	/* exit(0); */

}

int main() {

	int nSideX = 14;
	int nSideY = 14;
	vector<Particle> p;
	ExtVar s;
	double box = 20.0;
	double dt = 0.0001;
	double temperature = 1.00;

	InitializeSystem(p,s,box, temperature, nSideX, nSideY);

	ofstream initial("initial.dat");
	for (auto a:p) {
		initial << a.r << endl;
	}	

	ofstream trans("trans.dat");
	ofstream temp("temp.dat");
	ofstream frict("friction.dat");

	for (int step =0;step <5e6; step++) {
		UpdateSystem(p, s, dt, box, temperature); 
		if (step %2000 == 0) {
			cout << step << endl;
			double v2 =0;
			for (auto a:p) {
				trans << a.r << endl;
				v2+=a.SumVels2();
			}
			temp << step << " " << v2/(2*p.size()) << endl;
			frict << step << " " << s.r << endl;
		}
	}

	ofstream fin("final.dat");
	fin << p.size() << endl;
	fin << s.r << " " << s.v  << " " << s.m << endl;
	fin << temperature << endl;
	fin << box << endl;
	fin << endl;
	for (auto a:p) {
		fin << a.r << a.v << endl;
	}	
}
