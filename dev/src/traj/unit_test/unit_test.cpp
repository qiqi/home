// Boost Python interface for cpadj module.
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



#include <iostream>

#include "util/VectorUtil.h"
#include "traj/Trajectory.h"



using namespace std;
using namespace wangq;
using namespace wangq::traj;

void test_perigee_gradient1()
{
    cout << "\ntest_perigee_gradient1..." << endl;
    double h = EARTH_RADIUS;
    double vr = 2.0E+3;
    double vt = 8.0E+3;
    double p0 = perigee (h, vr, vt);
    double p1 = perigee (h + 1000, vr, vt);
    double p2 = perigee (h, vr + 100, vt);
    double p3 = perigee (h, vr, vt + 100);
    double grad_h, grad_vr, grad_vt;
    perigee_gradient (grad_h, grad_vr, grad_vt, h, vr, vt);
    double dp1 = grad_h * 1000;
    double dp2 = grad_vr * 100;
    double dp3 = grad_vt * 100;
    cout << p1 - p0 << ", " << dp1 << endl;
    cout << p2 - p0 << ", " << dp2 << endl;
    cout << p3 - p0 << ", " << dp3 << endl;
}

void test_perigee_gradient2()
{
    cout << "\ntest_perigee_gradient2..." << endl;
    double x0 = EARTH_RADIUS, y0 = EARTH_RADIUS * 2, z0 = -0.5 * EARTH_RADIUS;
    double vx0 = 3.0E+3, vy0 = -2.0E+3, vz0 = 5.0E+3;
    double p0 = perigee (Position(x0, y0, z0), Vector(vx0, vy0, vz0));
    double p1 = perigee (Position(x0 + 1000, y0, z0), Vector(vx0, vy0, vz0));
    double p2 = perigee (Position(x0, y0 + 1000, z0), Vector(vx0, vy0, vz0));
    double p3 = perigee (Position(x0, y0, z0 + 1000), Vector(vx0, vy0, vz0));
    double p4 = perigee (Position(x0, y0, z0), Vector(vx0 + 100, vy0, vz0));
    double p5 = perigee (Position(x0, y0, z0), Vector(vx0, vy0 + 100, vz0));
    double p6 = perigee (Position(x0, y0, z0), Vector(vx0, vy0, vz0 + 100));
    Vector grad_x(0,0,0), grad_v(0,0,0);
    perigee_gradient (grad_x, grad_v, Position(x0, y0, z0), 
                                      Vector(vx0, vy0, vz0));
    double dp1 = grad_x.x * 1000;
    double dp2 = grad_x.y * 1000;
    double dp3 = grad_x.z * 1000;
    double dp4 = grad_v.x * 100;
    double dp5 = grad_v.y * 100;
    double dp6 = grad_v.z * 100;
    cout << p1 - p0 << ", " << dp1 << endl;
    cout << p2 - p0 << ", " << dp2 << endl;
    cout << p3 - p0 << ", " << dp3 << endl;
    cout << p4 - p0 << ", " << dp4 << endl;
    cout << p5 - p0 << ", " << dp5 << endl;
    cout << p6 - p0 << ", " << dp6 << endl;
}

void test_d_gravity()
{
    cout << "\ntest_d_gravity..." << endl;
    Position x(EARTH_RADIUS, 2 * EARTH_RADIUS, 0.0);
    Vector dx(1000.0, 500.0, -200.0);
    Vector dg0 = gravity(x + dx) - gravity(x);
    Vector dg1 = d_gravity(x, dx);
    cout << dg0.x << " = " << dg1.x << endl;
    cout << dg0.y << " = " << dg1.y << endl;
    cout << dg0.z << " = " << dg1.z << endl;
}

void test_adjoint()
{
    cout << "\ntest_adjoint..." << endl;
    Position x0(EARTH_RADIUS, 0.0, 0.0);
    Vector v0(2.0E+3, 8.0E+3, 0.0);
    double m0 = 1.0E+5;
    double isp = 500;
    double dt = 1000.0;
    vector<double> fx, fy, fz;
    for (int i = 0; i < 2; ++ i) {
        fx.push_back(0.0E+6);
        fy.push_back(0.0E+6);
        fz.push_back(0.0E+6);
    }
    Trajectory t0(x0, v0, m0, isp, dt, fx, fy, fz);
    double p0 = perigee(t0.trajectory.back(), t0.velocity.back());

    Position x1(EARTH_RADIUS + 100, 0.0, 0.0);
    Trajectory t1(x1, v0, m0, isp, dt, fx, fy, fz);
    double p1 = perigee(t1.trajectory.back(), t1.velocity.back());

    double dp0 = dot(t0.grad_traj[0], t1.trajectory[0] - t0.trajectory[0]);
    double dp1 = dot(t0.grad_traj[1], t1.trajectory[1] - t0.trajectory[1]) +
                 dot(t0.grad_vel[1], t1.velocity[1] - t0.velocity[1]);
    double dp2 = dot(t0.grad_traj[2], t1.trajectory[2] - t0.trajectory[2]) +
                 dot(t0.grad_vel[2], t1.velocity[2] - t0.velocity[2]);
    cout << p1 - p0 << endl;
    cout << "0: " << dp0 << endl;
    cout << "1: " << dp1 << endl;
    cout << "2: " << dp2 << endl;
}

int main()
{
    test_perigee_gradient1();
    test_perigee_gradient2();
    test_d_gravity();
    test_adjoint();
    return 0;
}

