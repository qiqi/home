// Provides automatic type conversions between C++ and Python
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



#include <boost/python.hpp>

#include <vector>

#include <util/BoostPythonConversion.h>

using namespace std;
using namespace boost::python;
using namespace wangq;



BOOST_PYTHON_MODULE(conversion)
{
    to_python_converter<vector<double>, stl_vector_to_python_list<double> >();
    to_python_converter<vector<int>, stl_vector_to_python_list<int> >();
    stl_vector_from_python_list<double>();
    stl_vector_from_python_list<int>();
}

