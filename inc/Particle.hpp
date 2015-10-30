#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "Vector2D.hpp"
#include "Vector3D.hpp"

template <class Vector>
struct ParticleGeneral {
	Vector r,v,a;
	double SumVels2() const {
		return v.MagnitudeSquared();
	}
};

double CorrectBox(double box, double p) {
	if (p > box/2)
		p-=box;
	else if (p <-box/2)
		p+=box;
	return p;
}

Vector3D CorrectBox(double box, Vector3D r) {
	r.x = CorrectBox(box,r.x);
	r.y = CorrectBox(box,r.y);
	r.z = CorrectBox(box,r.z);
	return r;
}

Vector2D CorrectBox(double box, Vector2D r) {
	r.x = CorrectBox(box,r.x);
	r.y = CorrectBox(box,r.y);
	return r;
}

#endif
