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



#ifndef WANGQ_TRAJ_TRAJECTORY_H
#define WANGQ_TRAJ_TRAJECTORY_H

#include <vector>

#include <util/VectorUtil.h>



namespace wangq {
namespace traj {
    const double GRAVITY = 6.67428E-11;  // gravitational constant in m^3/kg/s^2
    const double EARTH_MASS = 5.9736E+24;  // mass of earth in kg
    const double EARTH_GM = EARTH_MASS * GRAVITY;
    const double EARTH_RADIUS = 6372.8E+3;
    const double G0 = EARTH_GM / (EARTH_RADIUS * EARTH_RADIUS);
    
    
    
    class Trajectory {
        public:
        std::vector<wangq::Position> trajectory;
        std::vector<wangq::Vector> velocity;
        std::vector<double> mass, time;
        std::vector<wangq::Vector> force;
        std::vector<wangq::Vector> grad_traj, grad_vel, grad_force;
        std::vector<double> grad_mass;
    
        Trajectory (wangq::Position x0, wangq::Vector v0, double m0, double isp,
                    double dt, std::vector<double> fx, std::vector<double> fy,
                    std::vector<double> fz);

        std::vector<double> get_mass();
        std::vector<double> get_time();
        std::vector<double> trajectory_x();
        std::vector<double> trajectory_y();
        std::vector<double> trajectory_z();
        std::vector<double> velocity_x();
        std::vector<double> velocity_y();
        std::vector<double> velocity_z();
        std::vector<double> grad_fx();
        std::vector<double> grad_fy();
        std::vector<double> grad_fz();
        std::vector<double> grad_tx();
        std::vector<double> grad_ty();
        std::vector<double> grad_tz();
        std::vector<double> grad_vx();
        std::vector<double> grad_vy();
        std::vector<double> grad_vz();
        std::vector<double> get_grad_m();
    };
    
    
    
    double perigee (double h, double vr, double vt, double gm = EARTH_GM);
    
    double perigee (const wangq::Position& x, const wangq::Vector& v,
                    double gm = EARTH_GM, const wangq::Position& center =
                                          wangq::Position(0,0,0));
    
    void perigee_gradient (double& grad_h, double& grad_vr, double& grad_vt,
                           double h, double vr, double vt,
                           double gm = EARTH_GM);
    
    void perigee_gradient (wangq::Vector& grad_x, wangq::Vector& grad_v,
                           const wangq::Position& x, const wangq::Vector& v,
                           double gm = EARTH_GM, const wangq::Position& center
                                             = wangq::Position(0,0,0));
    
    
    
    wangq::Vector gravity (const wangq::Position x, double gm = EARTH_GM,
                           const wangq::Position& center =
                                        wangq::Position(0,0,0));
    
    wangq::Vector gravity (const wangq::Position& x1, const wangq::Position& x2,
                    double gm = EARTH_GM, const wangq::Position& center =
                                          wangq::Position(0,0,0));
    
    wangq::Vector d_gravity (const wangq::Position& x, const wangq::Vector& dx,
                             double gm = EARTH_GM,
                             const wangq::Position& center =
                                   wangq::Position(0,0,0));
    
    wangq::Vector d_gravity (const wangq::Position& x1,
                             const wangq::Position& x2,
                             const wangq::Vector& dx1, const wangq::Vector& dx2,
                             double gm = EARTH_GM,
                             const wangq::Position& center =
                                   wangq::Position(0,0,0));

} } // namespace wangq::traj

#endif
