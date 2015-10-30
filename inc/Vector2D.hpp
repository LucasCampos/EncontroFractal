#ifndef VECTORND_H
#define VECTORND_H

#include <iostream>
#include <cmath>

class Vector2D {
	public:
		double x,y;

		Vector2D(): x(0),y(0){};
		Vector2D(double _x, double _y): x(_x), y(_y) {};
		
		void SetElem(int i, double a) {
			switch(i) {
				case 0: x=a; break;
				case 1: y=a; break;
			}
		}

		Vector2D Sum(const Vector2D& v) const {
			return Vector2D(x+v.x, y+v.y);
		}

		Vector2D Subtract(const Vector2D& v) const {
			return Vector2D(x-v.x, y-v.y);
		}
		
		double InternalProduct(const Vector2D& v) const{
			return x*v.x+y*v.y;
		}

		double MagnitudeSquared() const{
			return x*x + y*y;
		}

		double Magnitude() const{
			return sqrt(MagnitudeSquared());
		}

		Vector2D Normalize() const{
			double mag = Magnitude();
			return Vector2D(x/mag, y/mag);
		}

		Vector2D MultScalar(double a) const{
			return Vector2D(a*x, a*y);
		}

		Vector2D operator+ (const Vector2D& v) const{
			return Sum(v);
		}

		Vector2D operator- (const Vector2D& v) const{
			return Subtract(v);
		}

		double operator* (const Vector2D& v) const{
			return InternalProduct(v);
		}

		Vector2D operator* (double a) const{
			return MultScalar(a);
		}

		Vector2D operator/ (double a) const{
			return MultScalar(1.0/a);
		}

		Vector2D& operator+= (const Vector2D& v) {
			x+=v.x;
			y+=v.y;
			return (*this);
		}

		Vector2D& operator-= (const Vector2D& v) {
			x-=v.x;
			y-=v.y;
			return (*this);
		}

		Vector2D& operator*= (double a) {
			x*=a;
			y*=a;
			return (*this);
		}


		friend std::ostream& operator<<(std::ostream& stream,const Vector2D& v){

			stream << v.x << "   " << v.y;
			return stream;
		}
		
};


#endif
