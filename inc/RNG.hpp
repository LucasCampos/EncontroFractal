#ifndef RANDOM_NUMBERS_HPP
#define RANDOM_NUMBERS_HPP

#include <cmath>
#include <random>
#include <chrono>

class RNG {
	
	std::default_random_engine rng;
	std::uniform_real_distribution<double> uni;
	std::uniform_real_distribution<double> uniD;
	std::normal_distribution<double> gau;

	double stddev;

	void InitializeAll(unsigned seed, double dev) {
		SetSeed(seed);
		SetDev(dev);
	}

	public:

	RNG(){
		InitializeAll(std::chrono::system_clock::now().time_since_epoch().count(), 1.0);
	}

	RNG(double dev) {
		InitializeAll(std::chrono::system_clock::now().time_since_epoch().count(), dev);
	}

	RNG(double dev, int seed) {
		InitializeAll(seed, dev);
	}

	void SetSeed(unsigned seed) {
		rng = std::default_random_engine(seed);
		uni = std::uniform_real_distribution<double>(0.0, 1.0);
	}

	void SetDev(double dev) {
		uniD = std::uniform_real_distribution<double>(0.0, dev);
		gau = std::normal_distribution<double>(0.0,dev);
		stddev = dev;
	}

	double inline GetDev() const {
		return gau.stddev();
	}

	double inline Uniform() {
		return uniD(rng);
	}

	double inline UniCentered() {
		return 2.0*(uniD(rng)-0.5);
	}

	double inline Gaussian() 
	{
		return gau(rng);
	}

	//Algorithm proposed at DOI 10.1007/s10955-011-0364-y, by 
	//N.M.A. Mohamed
	double MaxwellE() {

		double y;
		double rat;
		do {
			const double r1 = uni(rng);
			const double r2 = uni(rng);

			rat = r2/r1;
			y = -2*log(r1);
		} while (y*2.712609 < rat*rat);
		return y*stddev;

	}

	double MaxwellV() {
		const double E = MaxwellE();
		return sqrt(2*E);
	}
};


#endif

