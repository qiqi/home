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



#include <cassert>

#include <mcadj/Mcadj.h>
#include <mcadj/Mcadj1d.h>

using namespace std;
using namespace wangq::mcadj;



Mcadj1d::Mcadj1d(int num_grids)
{
    nodes_.resize(num_grids);

    for (int i = 0; i < num_grids; ++i) {
        nodes_[i].id = i;
        nodes_[i].outNodes.resize(3);
        nodes_[i].outWeights.resize(3);

        if (i != 0) {
            nodes_[i].outNodes[SIDE_LEFT] = & nodes_[i - 1];
        }
        if (i != num_grids - 1) {
            nodes_[i].outNodes[SIDE_RIGHT] = & nodes_[i + 1];
        }
        nodes_[i].outNodes[SIDE_CENTER] = & nodes_[i];
    }
}



void Mcadj1d::assignWeights(const vector<double>& weights_left,
                             const vector<double>& weights_center,
                             const vector<double>& weights_right)
{
    assert (weights_left.size() == nodes_.size());
    assert (weights_center.size() == nodes_.size());
    assert (weights_right.size() == nodes_.size());

    for (int i = 0; i < nodes_.size(); ++i) {
        nodes_[i].outWeights[SIDE_LEFT] = weights_left[i];
        nodes_[i].outWeights[SIDE_CENTER] = weights_center[i];
        nodes_[i].outWeights[SIDE_RIGHT] = weights_right[i];
    }
}



void Mcadj1d::initializeRandomWalk(const vector<int>& num_walkers, bool is_fast)
{
    assert (num_walkers.size() == nodes_.size());

    int total_walkers = 0;
    for (int i = 0; i < nodes_.size(); ++i) {
        walkers_.resize(total_walkers + num_walkers[i]);
        for (int j = total_walkers; j < total_walkers + num_walkers[i]; ++j) {
            walkers_[j].currentState = & nodes_[i];
            walkers_[j].initialState = & nodes_[i];
            walkers_[j].lastState = 0;
            walkers_[j].weight = 1.0;
            walkers_[j].isFast = is_fast;
        }
        total_walkers += num_walkers[i];
    }

    assert (total_walkers == walkers_.size());
}



void Mcadj1d::advanceRandomWalk()
{
    for (int i = 0; i < walkers_.size(); ++i) {
        walkers_[i].randomWalk();
    }
}



void Mcadj1d::concludeRandomWalk(const vector<double>& input,
                                  vector<double>& output)
{
    assert (input.size() == nodes_.size());
    assert (output.size() == nodes_.size());

    vector<int> num_walkers;
    vector<double> sums;
    num_walkers.resize (nodes_.size(), 0);
    sums.resize (nodes_.size(), 0);

    for (int i = 0; i < walkers_.size(); ++i) {
        int init_id = walkers_[i].initialState->id;

        int final_id;
        if (walkers_[i].currentState != TERMINAL_STATE) {
            final_id = walkers_[i].currentState->id;
            sums[init_id] += input[final_id] * walkers_[i].weight;
        }

        num_walkers[init_id] ++;
    }

    for (int i = 0; i < nodes_.size(); ++i) {
        output[i] = sums[i] / double(num_walkers[i]);
    }
}

