// Checkpointing adjoint solver for partial differential equations.
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



#ifndef WANGQ_CPADJ_CPADJPDE_H
#define WANGQ_CPADJ_CPADJPDE_H

#include <list>
#include <vector>

#include <cpadj/CpadjDynamicCheckPoint.h>



namespace wangq
{
    namespace cpadj
    {
        class CpadjPde {
            public:
                CpadjPde(const std::vector<double>& solution_initial,
                         int memory_size=5);
                virtual ~CpadjPde() {;}

                void advance(double cfl, double time);
                void advance(const std::vector<double>& cfl_numbers);
                void retrograde(const std::vector<double>& adjoint_initial,
                                int num_steps = 0);

                std::vector<int> checkPoints();
                std::vector<int> checkPointLevels();

                std::vector<double> solution() {
                    return solution_;
                }
                std::vector<double> adjoint() {
                    return adjoint_;
                }
                double memorySize() {
                    return memory_.memorySize();
                }
            protected:
                std::vector<double> solution_, adjoint_;
            private:
                virtual double advanceSolution(double cfl);
                virtual double retrogradeAdjoint(double cfl);
                class Memory {
                    public:
                        Memory() {}
                        Memory(const std::vector<double>& solution_, double cfl)
                            : solution(solution_),
                              cflNumbers(std::vector<double>(1, cfl)) {}
                        void mergeWith(const Memory& next);
                        std::vector<double> solution;
                        std::vector<double> cflNumbers;
                };
                DynamicCheckPoints<Memory> memory_;
        };

    }   // namespace cpadj
}   // namespace wangq



#endif

