// Monte Carlo adjoint solver for 1 dimensional partial differential equations.
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



#ifndef WANGQ_MCADJ_MCADJPDE_H
#define WANGQ_MCADJ_MCADJPDE_H

#include <mcadj/Mcadj1d.h>



namespace wangq
{
    namespace mcadj
    {
        class McadjPde {
            public:
                McadjPde(const std::vector<double>& solution,
                         const std::vector<int>& num_walkers,
                         bool is_fast=false);
                virtual ~McadjPde();
                std::vector<double> solution();
                std::vector<double>
                    concludeRandomWalk(const std::vector<double>& input);
                void advance(double cfl, int num_steps);
            protected:
                std::vector<double> solution_;
                std::vector<double> leftWeights_, centerWeights_, rightWeights_;
            private:
                virtual double advanceSolution(double cfl);
                virtual void assignWeights(double cfl);
                Mcadj1d* mcadj1d_;
        };

    }   // namespace mcadj
}   // namespace wangq



#endif

