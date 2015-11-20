/** @file Vector.h*/

#ifndef VECTOR_H
#define VECTOR_H

#include <math.h>
#include "quaternion.h"
#include <sstream>
#include "../mc/randomGenerator.h"
#include "macros.h"



#ifdef RAN2
    extern Ran2 ran2;
#else
  #ifdef DSFMT
    extern Dsfmt ran2;
  #else
    extern MersenneTwister ran2;
  #endif
#endif

class Quat;

class Vector {
public:
    double x, y, z;
    Vector(){}
    Vector(double x, double y, double z) : x(x), y(y), z(z) {}

    std::string info() {
        std::ostringstream o;
        o << "(" <<x << ", " << y << ", " << z <<")";
        return o.str();
    }

    std::string infoRaw() {
        std::ostringstream o;
        o <<x << " " << y << " " << z;
        return o.str();
    }

    /**
     * @brief Return a norm of Vector
     * @return
     */
    inline double size() {
        return sqrt( pow(this->x,2) + pow(this->y,2) + pow(this->z,2));
    }

    /**
     * @brief normalise Normalise a vector to have unit length.  For speed during heavy use, it is
       not checked that the supplied vector has non-zero length.
     */
    inline void normalise() {
        double tot = size();
        if (tot !=0.0) {
            tot = 1.0 / tot;
            x *= tot;
            y *= tot;
            z *= tot;
        }
    }

    /**
     * @brief dotProduct == Scallar product
     * @param other
     * @return
     */
    inline double dot(Vector& other)  {
        return x*other.x + y*other.y + z*other.z;
    }

    inline void scale(double scale) {
        x=x*scale; y=y*scale, z=z*scale;
    }

    inline Vector operator- (const Vector& o) const {
        return Vector(x-o.x, y-o.y,z-o.z);
    }

    inline void operator-= (const Vector& o) {
        x-=o.x;
        y-=o.y;
        z-=o.z;
    }

    inline void operator+= (const Vector& o) {
        x+=o.x;
        y+=o.y;
        z+=o.z;
    }

    inline bool operator==(Vector& o) {
        if(x == o.x && y == o.y && z == o.z) return true;
        else return false;
    }

    inline bool operator!= (Vector& o) {
        if(o.x != x || o.y != y || o.z != z ) return true;
        return false;
    }

    inline Vector operator* (double scale) {
        return Vector(this->x*scale, this->y*scale, this->z*scale);
    }

    inline void operator*= (double scale) {
        this->x*=scale, this->y*=scale, this->z*=scale;
    }

    friend inline Vector operator* (double,Vector&);

    inline Vector cross(Vector& B) {
        return Vector(this->y*B.z - this->z*B.y, -this->x*B.z + this->z*B.x, this->x*B.y - this->y*B.x);
    }

    inline void ortogonalise(Vector& B) {
        double dp(dot(B));    this->x -= dp*B.x;    this->y -= dp*B.y;    this->z -= dp*B.z;
    }

    // using Vector&, double, double instead of Quat -> circular include
    inline void rotate(Vector& vec, double vc, double vs) {
        double t2,t3,t4,t5,t6,t7,t8,t9,t10,newx,newy,newz, qw(vc), qx(vec.x * vs), qy(vec.y * vs), qz(vec.z * vs);

        /*    t1 = quat.w * quat.w; */
        t2 =  qw * qx;
        t3 =  qw * qy;
        t4 =  qw * qz;
        t5 = -qx * qx;
        t6 =  qx * qy;
        t7 =  qx * qz;
        t8 = -qy * qy;
        t9 =  qy * qz;
        t10 = -qz * qz;

        newx = 2.0 * ( (t8+t10) * x + (t6-t4)  * y + (t3+t7) * z ) + x;
        newy = 2.0 * ( (t4+t6)  * x + (t5+t10) * y + (t9-t2) * z ) + y;
        newz = 2.0 * ( (t7-t3)  * x + (t2+t9)  * y + (t5+t8) * z ) + z;

        x = newx;
        y = newy;
        z = newz;
    }

    inline void rotate(Quat& q) {
        double t2,t3,t4,t5,t6,t7,t8,t9,t10,newx,newy,newz;

        /*    t1 = quat.w * quat.w; */
        t2 =  q.w * q.x;
        t3 =  q.w * q.y;
        t4 =  q.w * q.z;
        t5 = -q.x * q.x;
        t6 =  q.x * q.y;
        t7 =  q.x * q.z;
        t8 = -q.y * q.y;
        t9 =  q.y * q.z;
        t10 = -q.z * q.z;

        newx = 2.0 * ( (t8+t10) * x + (t6-t4)  * y + (t3+t7) * z ) + x;
        newy = 2.0 * ( (t4+t6)  * x + (t5+t10) * y + (t9-t2) * z ) + y;
        newz = 2.0 * ( (t7-t3)  * x + (t2+t9)  * y + (t5+t8) * z ) + z;

        x = newx;
        y = newy;
        z = newz;
    }

    /**
     * @brief ranvec    Returns an evenly distributed random unit vector2 of unit length.
                        See Allen & Tildesley p349 or Frenkel & Smit p410.
     * @return    RANDOM vector2 ON UNIT SPHERE
     */
    inline void randomUnitSphere() {
        double a, b, xi1, xi2;

        do {
            xi1 = 1.0 - 2.0*ran2();
            xi2 = 1.0 - 2.0*ran2();

            a = xi1*xi1 + xi2*xi2;
        } while (a > 1.0);

        b = 2.0 * sqrt(1.0 - a);

        x = xi1 * b;
        y = xi2 * b;
        z = 1.0 - 2.0*a;
    }

    /**
     * @brief ranUnitSphereFaunus for comparison of grandcanonical with faunus
     */
    inline void ranUnitSphereFaunus() {
      double r2;
      do {
        this->x = 2*ran2()-1;
        this->y = 2*ran2()-1;
        this->z = 2*ran2()-1;
        r2 = this->dot(*this);
      } while (r2>1);
      this->scale(1/sqrt(r2));
    }

    inline void randomUnitCube() {
        x = ran2();
        y = ran2();
        z = ran2();
    }

    static inline Vector getRandomUnitCube() {
        return Vector(ran2(), ran2(), ran2());
    }

    static inline Vector getRandomUnitSphere() {
        Vector vec;
        vec.randomUnitSphere();
        return vec;
    }

    static inline Vector getOrthogonalVector(Vector input){
        if (input.x == 0.0 && input.y == 0.0)                   // Theoreticaly if even input.z == 0 then input is zeroth vector ... form definition of orthogonalyty any vector is perpendicular to Vector(0, 0, 0)....
            return Vector(0.0, 1.0, 0.0);
        return Vector(-input.y, input.x, 0.0);
    }

    static inline Vector getRandomOrthogonalVector(Vector input){
        Vector  vec = getOrthogonalVector(input);               // get orthogonal vector to input vector

        double  cosRanfomAngle = cos(ran2()*2*PI),              // get random rotation around input vector
                sinRandomAngle = sqrt(1 - cosRanfomAngle*cosRanfomAngle);

        vec.rotate(input, cosRanfomAngle, sinRandomAngle);      // rotate vector around input vector

        return vec;
    }

    static inline Vector getRandomUnitCone(Vector axis, const double maxangle){
        Vector  vec     = axis,                                 //returned vector in cone
                axis2   = getRandomUnitSphere();                //orthogonal vector to vec

        double  cosAngle = cos(maxangle*ran2()),                // get cosinus of random number in interval [0, maxangle]
                sinAngle = sqrt(1 - cosAngle*cosAngle);         // get sinus of same angle

        axis2.ortogonalise(axis);                               // now make vector orthogonal to cone axis
        axis2.normalise();                                      // just normalise vector

        vec.rotate(axis2, cosAngle, sinAngle);                  // roatate vector pointing in axis of cone by random angle in interval

        return vec;
    }

    static inline Vector getRandomUnitConeUniform(Vector axis, const double maxangle){
        Vector  vec     = axis,                                 // Vector that will be returned, at first initialized as unit vector in direction of cone axis
                axis2   = getRandomOrthogonalVector(axis);      // Randomly distributed vector in plane perpendicular to cone axis

<<<<<<< HEAD
        double  cosAngle = cos(maxangle),                       // Calculate cos() of maxangle defined by angle of cone
                multiplaier = ran2()*(1-cosAngle)+cosAngle;     // Get coordinate of output vector in direction of cone axis
=======
        double  multiplaier = ran2()*sin(maxangle);                // multiplaier = (r*sin(angle))*rand[0, 1] ==> [r*sin(angle), 0] ==> get multiplaier to get vectors to project on angle arch
>>>>>>> 7addede0107eacb3fe9cf53689749deda80c6ebf

        vec.normalise();                                        // Before we start working we want to normalise both axis
        axis2.normalise();

        vec *= multiplaier;
        vec += (sqrt(1-multiplaier*multiplaier))*axis2;

        return vec;
    }

    static inline Vector getRandomUnitSphereFaunus() {
        Vector vec;
        vec.ranUnitSphereFaunus();
        return vec;
    }
};

inline Vector operator* (double scale, Vector& vec) {
    return Vector(vec.x*scale, vec.y*scale, vec.z*scale);
}




#endif // VECTOR_H
