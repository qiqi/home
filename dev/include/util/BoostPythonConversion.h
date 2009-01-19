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



#ifndef WANGQ_UTIL_BOOSTPYTHONCONVERSION_H
#define WANGQ_UTIL_BOOSTPYTHONCONVERSION_H

#include <boost/python.hpp>

#include <vector>



namespace wangq
{
    template<class T>
    struct stl_vector_to_python_list
    {
        static PyObject* convert(std::vector<T> const& v)
        {
            boost::python::list result;
            for (int i = 0; i < v.size(); ++i) {
                result.append(v[i]);
            }
            return boost::python::incref(result.ptr());
        }
    };
    
    
    
    template<class T>
    struct stl_vector_from_python_list
    {
        stl_vector_from_python_list()
        {
            boost::python::converter::registry::push_back(
                &convertible,
                &construct,
                boost::python::type_id<std::vector<T> >());
        }
    
        static void* convertible (PyObject* obj_ptr)
        {
            return obj_ptr;
        }
    
        static void construct (PyObject* obj_ptr,
            boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            boost::python::object obj(boost::python::borrowed(obj_ptr));
            void* storage = (
              (boost::python::converter::rvalue_from_python_storage<
                  std::vector<T> >*)data)->storage.bytes;
            new (storage) std::vector<T>;
            ((std::vector<T>*)storage)->resize (PyObject_Length(obj.ptr()));
            for (int i = 0; i < PyObject_Length(obj.ptr()); ++i) {
                (*(std::vector<T>*)storage)[i] = static_cast<T>(
                    boost::python::extract<double>(obj[i]));
            }
            data->convertible = storage;
        }
    };

}   // namespace wangq

#endif

