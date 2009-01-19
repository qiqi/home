// Automatic type conversion functions between C++ and Python
//
// Copyright (C) 2007, 2007 Qiqi Wang (qiqi.wang+stanford@gmail.com)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.



#ifndef WANGQ_UTIL_VECTOR_UTIL_H
#define WANGQ_UTIL_VECTOR_UTIL_H

#include <cmath>



namespace wangq
{
    class Position {
        public:
        // Position() : x(0.0), y(0.0), z(0.0) {}
        Position(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
        double x, y, z;
    };
    
    class Vector {
        public:
        // Vector() : x(0.0), y(0.0), z(0.0) {}
        Vector(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
        double x, y, z;
    };


    
    inline Vector operator + (const Vector& v1, const Vector& v2)
    {
        return Vector(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
    }
    
    inline Vector operator - (const Vector& v1, const Vector& v2)
    {
        return Vector(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
    }
    
    inline Vector operator - (const Vector& v)
    {
        return Vector(- v.x, - v.y, - v.z);
    }
    
    inline Vector operator - (const Position& p1, const Position& p2)
    {
        return Vector(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
    }
    
    inline Vector operator * (const Vector& v, double a)
    {
        return Vector(v.x * a, v.y * a, v.z * a);
    }
    
    inline Vector operator * (double a, const Vector& v)
    {
        return Vector(v.x * a, v.y * a, v.z * a);
    }
    
    inline Vector operator / (const Vector& v, double a)
    {
        return Vector(v.x / a, v.y / a, v.z / a);
    }


    
    inline Position operator + (const Position& p, const Vector& dp)
    {
        return Position(p.x + dp.x, p.y + dp.y, p.z + dp.z);
    }
    
    inline Position operator + (const Vector& dp, const Position& p)
    {
        return Position(p.x + dp.x, p.y + dp.y, p.z + dp.z);
    }


    
    inline double norm (const Vector& v)
    {
        return sqrt (v.x * v.x + v.y * v.y + v.z * v.z);
    }
    
    inline Vector normalize (const Vector& v)
    {
        return v / norm(v);
    }
    
    inline double dot (const Vector& u, const Vector& v)
    {
        return u.x * v.x + u.y * v.y + u.z * v.z;
    }
    
}   // namespace wangq

#endif
