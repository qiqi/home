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



#include <boost/python.hpp>  // this must be put first

#include <vector>

#include "util/BoostPythonConversion.h"
#include "traj/Trajectory.h"

using namespace std;
using namespace boost::python;
using namespace wangq;
using namespace wangq::traj;



double (*perigee_simple)
    (const Position&, const Vector&, double, const Position&) = &perigee;
BOOST_PYTHON_FUNCTION_OVERLOADS(perigee_overloads, perigee, 2, 4)

BOOST_PYTHON_MODULE(trajectory)
{
    def("perigee", perigee_simple, perigee_overloads());

    class_<Position>("Position", init<double, double, double>())
        .def_readonly("x", &Position::x)
        .def_readonly("y", &Position::y)
        .def_readonly("z", &Position::z);

    class_<Vector>("Vector", init<double, double, double>())
        .def_readonly("x", &Vector::x)
        .def_readonly("y", &Vector::y)
        .def_readonly("z", &Vector::z);
        
    class_<Trajectory>("Trajectory", init<Position, Vector, double, double,
                          double, std::vector<double>, std::vector<double>,
                          std::vector<double> >())
        .add_property("m", &Trajectory::get_mass)
        .add_property("t", &Trajectory::get_time)
        .add_property("x", &Trajectory::trajectory_x)
        .add_property("y", &Trajectory::trajectory_y)
        .add_property("z", &Trajectory::trajectory_z)
        .add_property("vx", &Trajectory::velocity_x)
        .add_property("vy", &Trajectory::velocity_y)
        .add_property("vz", &Trajectory::velocity_z)
        .add_property("grad_fx", &Trajectory::grad_fx)
        .add_property("grad_fy", &Trajectory::grad_fy)
        .add_property("grad_fz", &Trajectory::grad_fz)
        .add_property("grad_vx", &Trajectory::grad_vx)
        .add_property("grad_vy", &Trajectory::grad_vy)
        .add_property("grad_vz", &Trajectory::grad_vz)
        .add_property("grad_tx", &Trajectory::grad_tx)
        .add_property("grad_ty", &Trajectory::grad_ty)
        .add_property("grad_tz", &Trajectory::grad_tz)
        .add_property("grad_m", &Trajectory::get_grad_m);
}

