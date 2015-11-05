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

Vector3D CalculateInteractionSingle(Vector3D d) {
	const double r2 = d.MagnitudeSquared();
	const double r6_ = 1.0/(r2*r2*r2);
	const double f = 24*(2.0*r6_*r6_ - r6_)/r2;
	d *= f;
	return d;
}

vector<Vector3D> CalculateInteraction(const vector<Particle>& p, double box) {

	vector<Vector3D> forc(p.size());

	const int n = p.size();
	for (int i=0; i<(n-1); i++) {
		for (int j=i+1; j<n; j++) {
			Vector3D d = p[i].r-p[j].r;
			d = CorrectBox(box, d);
			d = CalculateInteractionSingle(d);

			forc[i] += d;			
			forc[j] -= d;			
		}
	}
	return forc;
}

double CalculateTemperature(const vector<Particle>& p) {

	const double K = CalculateKinEnergy(p);
	return 2.0*K/(3.0*p.size());

}

double CalculateVirial(const vector<Particle>& p, double box) {

	const int n = p.size();
	double a1 = 0;
	for (int i=0; i<(n-1); i++) {
		for (int j=i+1; j<n; j++) {
			Vector3D d = p[i].r-p[j].r;
			d = CorrectBox(box, d);
			Vector3D f = CalculateInteractionSingle(d);
			a1 += f*d;
		}
	}

	return a1;
}

double CalculatePressure(const vector<Particle>& p, const ExtVar& V, double temperature, double box) {

	const int n = p.size();
	const double virial = CalculateVirial(p,box);

	//cout <<"Pressure: "<< (n*temperature + virial/(3*n))/V.r << endl;
	//cout << n*temperature << endl;
	//cout <<virial/(3*n) << endl;
	//exit(0);
	
	return (n*temperature + virial/(3.0))/V.r;

}


void InitializeSystem(vector<Particle>& p, ExtVar& V, double box, int nSideX, int nSideY, int nSideZ) {

	int N=nSideX*nSideY*nSideZ;
	p.resize(N);
	RNG rng;

	V.r = box*box*box;
	V.v = 0;
	V.a = 0;
	V.m = 10.0;

	const Vector3D centerer = Vector3D(-0.5,-0.5, -0.5)*box;
	const double stepX = box/nSideX;
	const double stepY = box/nSideY;
	const double stepZ = box/nSideZ;
	int idx = 0;
	for (int i=0; i<nSideX; i++) {
		for (int j=0; j<nSideY; j++) {
			for (int k=0; k<nSideZ; k++) {
				p[idx].r = Vector3D(i*stepX,j*stepY,k*stepZ) + centerer;
				p[idx].v = Vector3D(rng.UniCentered(), rng.UniCentered(), rng.UniCentered());
				idx++;
			}
		}
	}
}

double AccelVolume(const vector<Particle>& p, const ExtVar& V, double pressure, double box) {

	const int n = p.size();

	const double virial = CalculateVirial(p,box);
	const double temperature = CalculateTemperature(p);

	return -pressure + (n*temperature + virial/3.0)/V.r;

}

void UpdateSystem(vector<Particle>& p, ExtVar& V, double pressure, double dt, double& box) {

	const int n = p.size();
	auto f = CalculateInteraction(p, box);

	V.a = AccelVolume(p,V,pressure, box)/V.m;

	V.v += V.a*dt;
	V.r += V.v*dt;
	
	const double dlnV = V.v/(3.0*V.r);

	//cout << V.r << " " << V.v << " " << V.a << endl;
	//cout << dlnV << endl;
	box = pow(V.r,1.0/3.0);
	
	for (int i=0; i<n; i++) {
		p[i].v += (f[i]   - p[i].v*dlnV)*dt;
		p[i].r += (p[i].v + p[i].r*dlnV)*dt;
		p[i].r = CorrectBox(box, p[i].r);
	}


}

int main() {

	int nSideX = 5;
	int nSideY = 5;
	int nSideZ = 3;
	vector<Particle> p;
	ExtVar V;
	double box = 4.0;
	double dt = 0.0001;
	double pressure = 200;

	InitializeSystem(p, V, box, nSideX, nSideY, nSideZ);

	ofstream initial("initial.dat");
	for (auto a:p) {
		initial << a.r << endl;
	}	

	ofstream trans("trans.dat");
	ofstream temp("temp.dat");
	ofstream energy("energy.dat");
	ofstream press("pressure.dat");
	ofstream volume("volume.dat");

	for (int step =0;step <5e6; step++) {
		UpdateSystem(p, V, pressure, dt, box); 
		if (step %2000 == 0) {
			cout << step << endl;
			double v2 =0;
			for (auto a:p) {
				trans << a.r << endl;
				v2+=a.SumVels2();
			}
			double temperature = CalculateTemperature(p);
			temp << step << " " << temperature << endl;
			energy << step << " " << CalculateEnergy<Particle, Vector3D>(p,box) << endl;
			press << step << " " << CalculatePressure(p,V,temperature, box) << endl;
			volume << step << " " << V.r << endl;
		}
	}

	ofstream fin("final.dat");
	for (auto a:p) {
		fin << a.r << a.v << endl;
	}	
}
