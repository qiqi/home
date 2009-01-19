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



#include <cassert>

#include <mcadj/Mcadj1d.h>
#include <mcadj/McadjPde.h>

using namespace std;
using namespace wangq::mcadj;



McadjPde::McadjPde(const vector<double>& solution,
                   const vector<int>& num_walkers, bool is_fast)
{
    int grid_size = solution.size();
    assert (grid_size == num_walkers.size());

    solution_.resize (grid_size, 0.0);
    solution_ = solution;
    
    mcadj1d_ = new Mcadj1d(grid_size);
    mcadj1d_->initializeRandomWalk(num_walkers, is_fast);

    leftWeights_.resize (grid_size, 0.0);
    centerWeights_.resize (grid_size, 1.0);
    rightWeights_.resize (grid_size, 0.0);
}



McadjPde::~McadjPde()
{
    delete mcadj1d_;
}



vector<double> McadjPde::solution()
{
    return solution_;
}



vector<double> McadjPde::concludeRandomWalk(const vector<double>& input)
{
    vector<double> output(input.size(), 0.0);
    mcadj1d_->concludeRandomWalk(input, output);
    return output;
}



void McadjPde::advance(double cfl, int num_steps)
{
    for (int i = 0; i < num_steps; ++i) {
        advanceSolution(cfl);
        assignWeights(cfl);
        mcadj1d_->assignWeights (leftWeights_, centerWeights_, rightWeights_);
        mcadj1d_->advanceRandomWalk();
    }
}



double McadjPde::advanceSolution(double cfl)
{
    assert (cfl > 0.0 && cfl <= 1.0);
}



void McadjPde::assignWeights(double cfl)
{
    for (int i = 0; i < solution_.size(); ++i) {
        leftWeights_[i] = 0.0;
        centerWeights_[i] = 1.0;
        rightWeights_[i] = 0.0;
    }
}

