#ifndef TOOLS_H
#define TOOLS_H

#include <cmath>
#include <vector>
#include <iostream>



const double PI = 3.14159265358979;



class Coord2D {
    public:
    double x, y;

    Coord2D(double ix, double iy) {
        x = ix;
        y = iy;
    }

    double length() const {
        return sqrt(x * x + y * y);
    }

    Coord2D  operator-() const {
        return Coord2D(-x,  -y);
    }

    Coord2D& operator+=(const Coord2D& p) {
        x += p.x;
        y += p.y;
        return *this;
    }

    Coord2D& operator-=(const Coord2D& p) {
        x -= p.x;
        y -= p.y;
        return *this;
    }

    Coord2D& operator*=(const double a) {
        x *= a;
        y *= a;
        return *this;
    }

    Coord2D& operator/=(const double a) {
        x /= a;
        y /= a;
        return *this;
    }

    Coord2D rotateRight(void)      const { return Coord2D(y,-x); }
    Coord2D rotateLeft(void)      const { return Coord2D(-y,x); }
};

inline Coord2D operator+ (const Coord2D& p, const Coord2D& q)
                { return Coord2D(p.x+q.x, p.y+q.y); }
inline Coord2D operator- (const Coord2D& p, const Coord2D& q)
                { return Coord2D(p.x-q.x, p.y-q.y); }
inline Coord2D operator* (const Coord2D& p, const double a)
                { return Coord2D(a*p.x, a*p.y); }
inline Coord2D operator* (const double a, const Coord2D& p)
                { return Coord2D(a*p.x, a*p.y); }
inline Coord2D operator/ (const Coord2D& p, const double a)
                { return Coord2D(p.x/a, p.y/a); }

// dot product
inline double dotprod( const Coord2D& a, const Coord2D& b ) { return a.x*b.x + a.y*b.y; }
// cross product
inline double crossprod( const Coord2D& a, const Coord2D& b ) { return a.x*b.y - a.y*b.x; }
// solve for [x y]*r = b
inline Coord2D CalcCombination( const Coord2D& b, const Coord2D& x, const Coord2D& y )
{
    Coord2D res;
    res.x = b.x * y.y - b.y * y.x;
    res.y = b.y * x.x - b.x * x.y;
    res /= crossprod( x, y );
    return res;
}

// output
inline ostream& operator << ( ostream& s, Coord2D x )
{
    s << '(' << x.x << ',' << x.y << ')';
    return s;
}

/* ======================================================================= */
/*                                 Coord3D                                 */
/* ======================================================================= */

class Coord3D {
    public:
    double x, y, z;

    Coord3D(){}
    Coord3D(double ix, double iy, double iz) {x=ix; y=iy; z=iz;}
    double len() const { return sqrt(x*x + y*y); }
    Coord3D  operator- ()           const { return Coord3D(-x,-y,-z); }
    Coord3D& operator+=(const Coord3D& p) { x+=p.x; y+=p.y; z+=p.z; return *this; }
    Coord3D& operator-=(const Coord3D& p) { x-=p.x; y-=p.y; z-=p.z; return *this; }
    Coord3D& operator*=(const double a)   { x*=a;   y*=a;   z*=a;   return *this; }
    Coord3D& operator/=(const double a)   { x/=a;   y/=a;   z/=a;   return *this; }
    Coord3D  PermuteRight(void)     const { return Coord3D(z,y,x); }
    Coord3D  PermuteLeft (void)     const { return Coord3D(y,z,x); }
};

inline Coord3D operator+ (const Coord3D& p, const Coord3D& q)
                { return Coord3D(p.x+q.x, p.y+q.y, p.z+q.z); }
inline Coord3D operator- (const Coord3D& p, const Coord3D& q)
                { return Coord3D(p.x-q.x, p.y-q.y, p.z-q.z); }
inline Coord3D operator* (const Coord3D& p, const double a)
                { return Coord3D(a*p.x, a*p.y, a*p.z); }
inline Coord3D operator* (const double a, const Coord3D& p)
                { return Coord3D(a*p.x, a*p.y, a*p.z); }
inline Coord3D operator/ (const Coord3D& p, const double a)
                { return Coord3D(p.x/a, p.y/a, p.z/a); }

// dot product
inline double dotprod( const Coord3D& a, const Coord3D& b ) { return a.x*b.x + a.y*b.y + a.z*b.z; }
// cross product
inline Coord3D crossprod( const Coord3D& a, const Coord3D& b )
{
    return Coord3D( a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x );
}

// output
inline ostream& operator << ( ostream& s, Coord3D x )
{
    s << '(' << x.x << ',' << x.y << ',' << x.z << ')';
    return s;
}

/* ======================================================================= */
/*                               Use Coord2D                               */
/* ======================================================================= */

typedef Coord2D Coord;

/* ======================================================================= */

int search( vector<int> v, int x, int start=-1, int end=-1 );

#endif
