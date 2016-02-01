/** @file quaternion.h*/

#ifndef QUATERNION_H
#define QUATERNION_H

#include <iostream>

#include "Vector.h"

using namespace std;

class Vector;

class Quat{
public:
    double w,x,y,z;

    Quat(double w,double x, double y, double z) :
            w( w ),
            x( x ),
            y( y ),
            z( z )
    {}

    Quat( Vector &vec, double vc, double vs ) :
        w( vc ),
        x( vec.x * vs ),
        y( vec.y * vs ),
        z( vec.z * vs )
    {}

//    Quat( Vector &vec, double angle ){
//        double  w = cos( angle ),
//                vs = sqrt( 1 - w * w );
//        x = vec.x * vs;
//        y = vec.y * vs;
//        z = vec.z * vs;
//    }

//    void chaneRotationSign(){
//        x *= -1.0;
//        y *= -1.0;
//        z *= -1.0;
//    }

    void print() {
        cout <<"( "<<w<<", "<<x<<", "<<y<<", "<<z<<" )"<< endl;
    }
};

#endif // QUATERNION_H
