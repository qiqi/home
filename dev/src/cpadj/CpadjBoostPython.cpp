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

#include <util/BoostPythonConversion.h>
#include <cpadj/CpadjPde.h>
#include <cpadj/CpadjBurgers.h>

using namespace std;
using namespace boost::python;
using namespace wangq;
using namespace wangq::cpadj;



// default argument for CpadjPde::retrograde
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(CpadjPde_overloads, retrograde, 1, 2)

// avoid overloaded member function 'advance'
void (CpadjPde::*advance_simple)(double, double) = &CpadjPde::advance;
void (CpadjPde::*advance_complex)(const std::vector<double>& cfl_numbers)
    = &CpadjPde::advance;

BOOST_PYTHON_MODULE(cpadj)
{
    class_<CpadjPde>("Pde", init<const vector<double>&, int>())
        .def(init<const vector<double>&>())
        .def("advance", advance_simple)
        .def("advance", advance_complex)
        .def("retrograde", &CpadjPde::retrograde, CpadjPde_overloads())
        .add_property("checkPoints", &CpadjPde::checkPoints)
        .add_property("checkPointLevels", &CpadjPde::checkPointLevels)
        .add_property("solution", &CpadjPde::solution)
        .add_property("adjoint", &CpadjPde::adjoint)
        .add_property("memorySize", &CpadjPde::memorySize);

    class_<CpadjBurgers, bases<CpadjPde> >("Burgers",
        init<const vector<double>&, int>())
        .def(init<const vector<double>&>());

    // Automatic conversion of std::vector<double> and std::vector<int>
    // to and from python list;
    // included from util/BoostPythonConversion.h.
//    to_python_converter<vector<double>, stl_vector_to_python_list<double> >();
//    to_python_converter<vector<int>, stl_vector_to_python_list<int> >();
//    stl_vector_from_python_list<double>();
//    stl_vector_from_python_list<int>();
}

