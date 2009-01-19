// Monte Carlo adjoint solver in 1 dimensional grid.
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



#ifndef WANGQ_MCADJ_MCADJ1D_H
#define WANGQ_MCADJ_MCADJ1D_H

#include <vector>

#include <mcadj/Mcadj.h>



namespace wangq
{
    namespace mcadj
    {
        class Mcadj1d {
            public:
                enum Side {SIDE_LEFT, SIDE_CENTER, SIDE_RIGHT};
                Mcadj1d(int num_grids);
                void assignWeights(const std::vector<double>& weights_left,
                                   const std::vector<double>& weights_center,
                                   const std::vector<double>& weights_right);
                void initializeRandomWalk(const std::vector<int>& num_walkers,
                                          bool is_fast=false);
                void advanceRandomWalk();
                void concludeRandomWalk(const std::vector<double>& input,
                                        std::vector<double>& output);
            private:
                std::vector<McadjNode> nodes_;
                std::vector<McadjWalker> walkers_;
        };

    }   // namespace mcadj
}   // namespace wangq



#endif

