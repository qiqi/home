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



#include <algorithm>
#include <cassert>
#include <cmath>

#include <mcadj/McadjBurgers.h>

using namespace std;
using namespace wangq::mcadj;



McadjBurgers::McadjBurgers(const vector<double>& solution,
                           const vector<int>& num_walkers,
                           bool is_fast)
    : McadjPde(solution, num_walkers, is_fast)
{
    int grid_size = solution.size();
    fluxes_.resize(grid_size + 1, 0.0);
    lastSolution_.resize(grid_size, 0.0);
}



double McadjBurgers::advanceSolution(double cfl)
{
    assert (cfl > 0.0 && cfl <= 1.0);

    // compute maximum wave speed and time step
    double max_speed = 0.0;
    for (int i = 0; i < solution_.size(); ++i) {
        if (2 * fabs(solution_[i]) > max_speed) {
            max_speed = 2 * fabs(solution_[i]);
        }
    }
    deltaT_ = cfl / max(max_speed, 1.0E-12);

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
    lastSolution_ = solution_;
    for (int i = 0; i < solution_.size(); ++i) {
        solution_[i] += deltaT_ * (fluxes_[i] - fluxes_[i + 1]);
    }
}



void McadjBurgers::assignWeights(double cfl)
{
    for (int i = 0; i < solution_.size(); ++i) {
        leftWeights_[i] = 0.0;
        centerWeights_[i] = 1.0;
        rightWeights_[i] = 0.0;
    }

    if (lastSolution_.front() < 0.0)
        centerWeights_.front() += 2 * deltaT_ * lastSolution_.front();

    if (lastSolution_.back() > 0.0)
        centerWeights_.back() -= 2 * deltaT_ * lastSolution_.back();

    for (int i = 1; i < lastSolution_.size(); ++i) {
        if (lastSolution_[i - 1] + lastSolution_[i] >= 0.0) {
            centerWeights_[i - 1] -= 2 * deltaT_ * lastSolution_[i - 1];
            rightWeights_[i - 1] += 2 * deltaT_ * lastSolution_[i - 1];
        }
        else {
            leftWeights_[i] -= 2 * deltaT_ * lastSolution_[i];
            centerWeights_[i] += 2 * deltaT_ * lastSolution_[i];
        }
    }

    leftWeights_.front() = 1.0 - centerWeights_.front() - rightWeights_.front();
    rightWeights_.back() = 1.0 - centerWeights_.back() - leftWeights_.back();
}

