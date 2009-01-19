// Dynamic checkpointing adjoint solver for partial differential equations.
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



#include <cassert>
#include <iostream>
#include <list>
#include <vector>

#include <cpadj/CpadjPde.h>

using namespace std;
using namespace wangq::cpadj;



CpadjPde::CpadjPde(const vector<double>& solution_initial, int memory_size)
    : solution_(solution_initial), memory_(memory_size)
{;}



vector<int> CpadjPde::checkPointLevels()
{
    vector<int> levels;
    for (int i = 0; i <= memory_.numCheckPoints(); ++i)
        levels.push_back (memory_.level(i));
    return levels;
}



vector<int> CpadjPde::checkPoints()
{
    vector<int> check_points;
    // int position = 0;
    // for (int i = 0; i < memory_.numCheckPoints(); ++i) {
    //     check_points.push_back(position);
    //     position += memory_[i].cflNumbers.size();
    // }
    for (int i = 0; i <= memory_.numCheckPoints(); ++i)
        check_points.push_back (memory_.checkPoint(i));
    return check_points;
}



void CpadjPde::advance(double cfl, double time)
{
    double t = 0.0;
    while (t < time) {
        memory_.push(Memory(solution_, cfl));
        double dt = advanceSolution(cfl);
        t += dt;
        if (t + dt > time) {
            cfl = cfl * (time - t) / dt + 1.0E-8;
        }
    }
}



void CpadjPde::advance(const vector<double>& cfl_numbers)
{
    for (int i = 0; i < cfl_numbers.size(); ++i) {
        memory_.push(Memory(solution_, cfl_numbers[i]));
        advanceSolution(cfl_numbers[i]);
    }
}



void CpadjPde::retrograde(const vector<double>& adjoint_initial, int num_steps)
{
    if (num_steps == 0)
        num_steps = memory_.length();
    assert (memory_.length() >= num_steps);
    adjoint_ = adjoint_initial;

    for (int i = 0; i < num_steps; ++i) {
        Memory latest_memory = memory_.pop();

        solution_ = latest_memory.solution;
        double cfl = latest_memory.cflNumbers.back();
        latest_memory.cflNumbers.pop_back();
        advance(latest_memory.cflNumbers);

        retrogradeAdjoint(cfl);
    }
}



void CpadjPde::Memory::mergeWith(const CpadjPde::Memory& next)
{
    for (int i = 0; i < next.cflNumbers.size(); ++i) {
        cflNumbers.push_back (next.cflNumbers[i]);
    }
}



double CpadjPde::advanceSolution(double cfl)
{
    return cfl;
}



double CpadjPde::retrogradeAdjoint(double cfl)
{
    return cfl;
}

