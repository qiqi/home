// Dynamic checkpointing adjoint solver for Burgers' equation.
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
#include <cmath>
#include <iostream>
#include <vector>

#include <cpadj/CpadjBurgers.h>

using namespace std;
using namespace wangq::cpadj;


CpadjBurgers::CpadjBurgers(const vector<double>& solution_initial,
                           int memory_size)
    : CpadjPde(solution_initial, memory_size)
{
    int grid_size = solution_initial.size();
    fluxes_.resize(grid_size + 1, 0.0);
    adjoint_fluxes_.resize(grid_size + 1, 0.0);
}



double CpadjBurgers::advanceSolution(double cfl)
{
    double delta_t = deltaT(cfl);

    // compute fluxes
    if (solution_.front() >= 0.0)
        fluxes_.front() = 0.0;
    else
        fluxes_.front() = solution_.front() * solution_.front();

    if (solution_.back() <= 0.0)
        fluxes_.back() = 0.0;
    else
        fluxes_.back() = solution_.back() * solution_.back();

    for (int i = 1; i < solution_.size(); ++i) {
        if (solution_[i - 1] + solution_[i] >= 0.0)
            fluxes_[i] = solution_[i - 1] * solution_[i - 1];
        else
            fluxes_[i] = solution_[i] * solution_[i];
    }

    // advance solution
    for (int i = 0; i < solution_.size(); ++i) {
        solution_[i] += delta_t * (fluxes_[i] - fluxes_[i + 1]);
    }

    return delta_t / solution_.size();
}



double CpadjBurgers::retrogradeAdjoint(double cfl)
{
    assert (solution_.size() == adjoint_.size());
    double delta_t = deltaT(cfl);

    // compute adjoint fluexes
    adjoint_fluxes_.assign(solution_.size() + 1, 0.0);
    for (int i = 0; i < solution_.size(); ++i) {
        adjoint_fluxes_[i] += adjoint_[i] * delta_t;
        adjoint_fluxes_[i + 1] -= adjoint_[i] * delta_t;
    }

    // advance adjoint
    if (solution_.front() >= 0.0)
        ;
    else
        adjoint_.front() += 2 * solution_.front() * adjoint_fluxes_.front();

    if (solution_.back() <= 0.0)
        ;
    else
        adjoint_.back() += 2 * solution_.back() * adjoint_fluxes_.back();

    for (int i = 1; i < solution_.size(); ++i) {
        if (solution_[i - 1] + solution_[i] >= 0.0)
            adjoint_[i - 1] += 2 * solution_[i - 1] * adjoint_fluxes_[i];
        else
            adjoint_[i] += 2 * solution_[i] * adjoint_fluxes_[i];
    }

    return delta_t / solution_.size();
}



double CpadjBurgers::deltaT(double cfl)
{
    assert (cfl > 0.0 && cfl <= 1.0);

    // compute maximum wave speed and time step
    double max_speed = 0.0;
    for (int i = 0; i < solution_.size(); ++i) {
        if (2 * fabs(solution_[i]) > max_speed) {
            max_speed = 2 * fabs(solution_[i]);
        }
    }
    return cfl / max(max_speed, 1.0E-12);
}

