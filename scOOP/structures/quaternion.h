/** @file quaternion.h*/

#ifndef QUATERNION_H
#define QUATERNION_H

#include <iostream>
#include <math.h>

using namespace std;

class Quat{
public:
    double w,x,y,z;

    Quat(double w,double x, double y, double z): w(w), x(x), y(y), z(z) {}

    Quat(): w(0.0), x(0.0), y(0.0), z(0.0){}

    set( double axis_x, double axis_y, double axis_z, double angle ){
        double vc, vs;
        vc = cos(angle);
        angle > 0.0 ? vs = sqrt(1.0 - vc*vc) : vs = -sqrt(1.0 - vc*vc);
        w = vc;
        x = axis_x * vs;
        y = axis_y * vs;
        z = axis_z * vs;
    }

    void print() {
        cout <<"( "<<w<<", "<<x<<", "<<y<<", "<<z<<" )"<< endl;
    }
};

#endif // QUATERNION_H
