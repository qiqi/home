// Monte Carlo adjoint solver for Burgers' equation.
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



#ifndef WANGQ_MCADJ_MCADJBURGERS_H
#define WANGQ_MCADJ_MCADJBURGERS_H

#include <mcadj/McadjPde.h>



namespace wangq
{
    namespace mcadj
    {
        class McadjBurgers : public McadjPde {
            public:
                McadjBurgers(const std::vector<double>& solution,
                             const std::vector<int>& num_walkers,
                             bool is_fast=false);
            private:
                virtual double advanceSolution(double cfl);
                virtual void assignWeights(double cfl);
                std::vector<double> fluxes_, lastSolution_;
                double deltaT_;
        };

    }   // namespace mcadj
}   // namespace wangq



#endif

