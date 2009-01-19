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



#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "util/VectorUtil.h"
#include "traj/Trajectory.h"

using namespace std;
using namespace wangq;
using namespace wangq::traj;



double wangq::traj::perigee (double h, double vr, double vt, double gm)
{
    double potential = - gm / h;
    double kinetic = 0.5 * (vr * vr + vt * vt);
    double e = potential + kinetic;  // total energy
    double w = h * vt;               // angular momentum
    // solve for: v^2/2 - gm/h = e
    //                     v*h = w
    // eleimnate v we get: h^2 + (gm/e) h - w^2/(2e) = 0
    double sqrt_delta = sqrt( (gm/e)*(gm/e) + 2 * w*w/e );
    return 0.5 * (- gm / e - sqrt_delta);
}

void wangq::traj::perigee_gradient (double& grad_h, double& grad_vr,
                                    double& grad_vt, double h, double vr,
                                    double vt, double gm)
{
    double potential = - gm / h;
    double kinetic = 0.5 * (vr * vr + vt * vt);
    double e = potential + kinetic;  // total energy
    double w = h * vt;               // angular momentum
    // calculate derivative of the solution of: h^2 + (gm/e) h - w^2/(2e) = 0
    double sqrt_delta = sqrt( (gm/e)*(gm/e) + 2 * w*w/e );
    double grad_e = 0.5/(e*e) * (gm + (gm*gm/e + w*w) / sqrt_delta);
    double grad_w = - w/e / sqrt_delta;
    // calculate derivative to h and vt, vr
    grad_h = grad_w * vt + grad_e * gm / (h*h);
    grad_vt = grad_w * h + grad_e * vt;
    grad_vr = grad_e * vr;
}

double wangq::traj::perigee (const Position& x, const Vector& v, double gm,
                             const Position& center)
{
    double h = norm(x - center);
    double vr = dot(x - center, v) / h;
    double speed = norm(v);
    double vt = sqrt(speed * speed - vr * vr);
    return perigee (h, vr, vt);
}

void wangq::traj::perigee_gradient (Vector& grad_x, Vector& grad_v,
                                    const Position& x, const Vector& v,
                                    double gm, const Position& center)
{
    Vector d = x - center;
    double h = norm(d);
    double vr = dot(d, v) / h;
    double speed = norm(v);
    double vt = sqrt(speed * speed - vr * vr);
    // calculate gradient with respect to the scalar quantities
    double grad_h, grad_vr, grad_vt;
    perigee_gradient (grad_h, grad_vr, grad_vt, h, vr, vt, gm);
    // calculate the vector gradient
    double dvt_dvx = dot(v, Vector(d.y*d.y + d.z*d.z, -d.x*d.y, -d.x*d.z));
    double dvt_dvy = dot(v, Vector(-d.x*d.y, d.x*d.x + d.z*d.z, -d.y*d.z));
    double dvt_dvz = dot(v, Vector(-d.x*d.z, -d.y*d.z, d.x*d.x + d.y*d.y));
    grad_x.x = grad_h * d.x/h + (grad_vr - grad_vt * vr/vt) * dvt_dvx / (h*h*h);
    grad_x.y = grad_h * d.y/h + (grad_vr - grad_vt * vr/vt) * dvt_dvy / (h*h*h);
    grad_x.z = grad_h * d.z/h + (grad_vr - grad_vt * vr/vt) * dvt_dvz / (h*h*h);
    grad_v.x = grad_vr * d.x / h + grad_vt * dvt_dvx / (vt * h*h);
    grad_v.y = grad_vr * d.y / h + grad_vt * dvt_dvy / (vt * h*h);
    grad_v.z = grad_vr * d.z / h + grad_vt * dvt_dvz / (vt * h*h);
}



Vector wangq::traj::gravity (const Position x, double gm,
                             const Position& center)
{
    double r = norm(center - x);
    return gm * (center - x) / (r * r * r);
}

Vector wangq::traj::gravity (const Position& x1, const Position& x2,
                             double gm, const Position& center)
{
    return 0.5 * (gravity(x1, gm, center) + gravity(x2, gm, center));
}

Vector wangq::traj::d_gravity (const Position& x, const Vector& dx,
                               double gm, const Position& center)
{
    Vector d = center - x;
    double r = norm(d);
    double xr = d.x / r, yr = d.y / r, zr = d.z / r;
    double dgx = -dot(dx, Vector(1 - 3*xr*xr, -3*xr*yr, -3*xr*zr));
    double dgy = -dot(dx, Vector(-3*xr*yr, 1 - 3*yr*yr, -3*yr*zr));
    double dgz = -dot(dx, Vector(-3*xr*zr, -3*yr*zr, 1 - 3*zr*zr));
    return gm * Vector(dgx, dgy, dgz) / (r*r*r);
}

Vector wangq::traj::d_gravity (const Position& x1, const Position& x2,
                               const Vector& dx1, const Vector& dx2,
                               double gm, const Position& center)
{
    return 0.5 * (d_gravity(x1, dx1, gm, center) +
                  d_gravity(x2, dx2, gm, center));
}



int advance_trajectory (Position& x, Vector& v, double& m, double isp,
                double dt, const Vector& f)
{
    const double TOLERANCE = 1.0E-12;
    Position x0 = x;
    Vector v0 = v;
    double m0 = m;
    m -= norm(f) * dt / isp / G0;
    if (m <= 0.0) {  // too much fuel burned, cannot happen.
        m = m0;
        return -1;
    }
    Vector af = f / sqrt(m * m0);
    // inner iteration
    Vector v1 = v;
    Position x1 = x;
    int n_iter = 0;
    do {
        ++ n_iter;
        v1 = v, x1 = x;
        Vector g = gravity(x0 ,x);
        v = v0 + dt * (g + af);
        x = x0 + dt * (v + v0) / 2.0;
        if (++ n_iter > 1000) {
            throw string("Inner iteration diverge, reduce dt.");
        }
    } while(norm(v1 - v) / norm(v) > TOLERANCE or
            norm(x1 - x) / norm(x - x0) > TOLERANCE);
    return 0;
}

Vector retrograde_trajectory (Vector& grad_x, Vector& grad_v, double& grad_m,
                              const Position& x, const Position& x0,
                              double m0, double isp, double dt, const Vector& f)
{
    const double TOLERANCE = 1.0E-12;
    Vector grad_x0 = grad_x;
    Vector grad_v0 = grad_v;
    // inner iteration
    Vector grad_x1 = grad_x;
    Vector grad_v1 = grad_v;
    int n_iter = 0;
    do {
        grad_v1 = grad_v, grad_x1 = grad_x;
        grad_x = grad_x0 + d_gravity(x, x0, grad_v, Vector(0,0,0)) * dt;
        grad_v = grad_v0 + 0.5 * grad_x * dt;
        if (++ n_iter > 1000) {
            throw string("Inner iteration diverge, reduce dt.");
        }
    } while(norm(grad_v1 - grad_v) / norm(grad_v) > TOLERANCE or
            norm(grad_x1 - grad_x) / norm(grad_x) > TOLERANCE);
    // explicit part
    grad_x0 = grad_x, grad_v0 = grad_v;
    grad_x = grad_x + d_gravity(x, x0, Vector(0,0,0), grad_v0) * dt;
    grad_v = grad_v + 0.5 * grad_x0 * dt;
    Vector grad_a = grad_v0 * dt;
    // mass part
    double m = m0 - norm(f) * dt / isp / G0;
    double ma = sqrt(m * m0);
    grad_m = grad_m - dot(grad_a, f) * 0.5 / (m * ma);
    Vector grad_f = grad_a / ma;
    if (norm(f) > 0.0) grad_f = grad_f - grad_m * f / norm(f) * dt / isp / G0;
    else grad_f = grad_f - grad_m * grad_f / norm(grad_f) * dt / isp / G0;
    grad_m = grad_m - dot(grad_a, f) * 0.5 / (m0 * ma);
    return grad_f;
}



template<class T>
void reverse (vector<T>& v)
{
    for (int i = 0; i < v.size() / 2; ++ i) {
        int j = v.size() - 1 - i;
        T temp = v[i];
        v[i] = v[j];
        v[j] = temp;
    }
}





Trajectory::Trajectory (Position x, Vector v, double m, double isp,
                        double dt, vector<double> fx, vector<double> fy,
                        vector<double> fz)
{
    assert (fx.size() == fy.size() and fx.size() == fz.size());
    // forward sweep
    cout << "Calculating trajectory." << endl;
    trajectory.push_back(x);
    velocity.push_back(v);
    mass.push_back(m);
    time.push_back(0.0);
    for (int i = 0; i < fx.size(); ++ i) {
        Vector f(fx[i], fy[i], fz[i]);
        int res = advance_trajectory (x, v, m, isp, dt, f);
        assert (res >= 0);
        trajectory.push_back(x);
        velocity.push_back(v);
        mass.push_back(m);
        time.push_back(time.back() + dt);
        force.push_back(f);
    }
    // adjoint sweep
    cout << "Calculating adjoint." << endl;
    Vector grad_x(0,0,0), grad_v(0,0,0);
    double grad_m = 0.0;
    perigee_gradient (grad_x, grad_v, trajectory.back(), velocity.back());
    grad_traj.push_back (grad_x);
    grad_vel.push_back (grad_v);
    grad_mass.push_back (grad_m);
    for (int i = fx.size() - 1; i >= 0; -- i) {
        Vector grad_f = retrograde_trajectory (grad_x, grad_v, grad_m,
                trajectory[i+1], trajectory[i], mass[i], isp, dt, force[i]);
        grad_force.push_back(grad_f);
        grad_traj.push_back (grad_x);
        grad_vel.push_back (grad_v);
        grad_mass.push_back (grad_m);
    }
    // reverse gradients
    reverse(grad_force);
    reverse(grad_traj);
    reverse(grad_vel);
    reverse(grad_mass);
}



vector<double> Trajectory::trajectory_x()
{
    vector<double> result;
    for (int i = 0; i < trajectory.size(); ++ i)
        result.push_back(trajectory[i].x);
    return result;
}

vector<double> Trajectory::trajectory_y()
{
    vector<double> result;
    for (int i = 0; i < trajectory.size(); ++ i)
        result.push_back(trajectory[i].y);
    return result;
}

vector<double> Trajectory::trajectory_z()
{
    vector<double> result;
    for (int i = 0; i < trajectory.size(); ++ i)
        result.push_back(trajectory[i].z);
    return result;
}



vector<double> Trajectory::velocity_x()
{
    vector<double> result;
    for (int i = 0; i < velocity.size(); ++ i)
        result.push_back(velocity[i].x);
    return result;
}

vector<double> Trajectory::velocity_y()
{
    vector<double> result;
    for (int i = 0; i < velocity.size(); ++ i)
        result.push_back(velocity[i].y);
    return result;
}

vector<double> Trajectory::velocity_z()
{
    vector<double> result;
    for (int i = 0; i < velocity.size(); ++ i)
        result.push_back(velocity[i].z);
    return result;
}



vector<double> Trajectory::get_mass()
{ return mass; }

vector<double> Trajectory::get_time()
{ return time; }

vector<double> Trajectory::get_grad_m()
{ return grad_mass; }



vector<double> Trajectory::grad_fx()
{
    vector<double> result;
    for (int i = 0; i < grad_force.size(); ++ i)
        result.push_back(grad_force[i].x);
    return result;
}

vector<double> Trajectory::grad_fy()
{
    vector<double> result;
    for (int i = 0; i < grad_force.size(); ++ i)
        result.push_back(grad_force[i].y);
    return result;
}

vector<double> Trajectory::grad_fz()
{
    vector<double> result;
    for (int i = 0; i < grad_force.size(); ++ i)
        result.push_back(grad_force[i].z);
    return result;
}



vector<double> Trajectory::grad_tx()
{
    vector<double> result;
    for (int i = 0; i < grad_traj.size(); ++ i)
        result.push_back(grad_traj[i].x);
    return result;
}

vector<double> Trajectory::grad_ty()
{
    vector<double> result;
    for (int i = 0; i < grad_traj.size(); ++ i)
        result.push_back(grad_traj[i].y);
    return result;
}

vector<double> Trajectory::grad_tz()
{
    vector<double> result;
    for (int i = 0; i < grad_traj.size(); ++ i)
        result.push_back(grad_traj[i].z);
    return result;
}



vector<double> Trajectory::grad_vx()
{
    vector<double> result;
    for (int i = 0; i < grad_vel.size(); ++ i)
        result.push_back(grad_vel[i].x);
    return result;
}

vector<double> Trajectory::grad_vy()
{
    vector<double> result;
    for (int i = 0; i < grad_vel.size(); ++ i)
        result.push_back(grad_vel[i].y);
    return result;
}

vector<double> Trajectory::grad_vz()
{
    vector<double> result;
    for (int i = 0; i < grad_vel.size(); ++ i)
        result.push_back(grad_vel[i].z);
    return result;
}

