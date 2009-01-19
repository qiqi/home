// Checkpointing adjoint solver for Burgers' equations.
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



#ifndef WANGQ_CPADJ_CPADJBURGERS_H
#define WANGQ_CPADJ_CPADJBURGERS_H

#include <list>
#include <vector>

#include <cpadj/CpadjPde.h>



namespace wangq
{
    namespace cpadj
    {
        class CpadjBurgers : public CpadjPde {
            public:
                CpadjBurgers(const std::vector<double>& solution_initial,
                             int memory_size=5);
            private:
                virtual double advanceSolution(double cfl);
                virtual double retrogradeAdjoint(double cfl);
                double deltaT(double cfl);
                std::vector<double> fluxes_, adjoint_fluxes_;
        };

    }   // namespace cpadj
}   // namespace wangq



#endif

