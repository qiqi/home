// Basic Monte Carlo adjoiont solver functionalities.
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



#include <cstdlib>
#include <cmath>

#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <mcadj/Mcadj.h>

using namespace wangq::mcadj;




int McadjNode::randomOutIndex(double& probability)
{
    // Initialize uniform(0, 1) random number generator uni.
    static boost::mt19937 generator;
    static boost::uniform_real<double> uni_dist(0, 1);
    static boost::variate_generator
           <boost::mt19937&, boost::uniform_real<double> >
           uni(generator, uni_dist);

    double weight_sum = 0.0;
    for (int i = 0; i < outNodes.size(); ++i) {
        weight_sum += fabs(outWeights[i]);
    }

    double select = uni() * weight_sum;
    for (int i = 0; i < outNodes.size(); ++i) {
        select -= fabs(outWeights[i]);
        if (select < 0) {
            probability = fabs(outWeights[i]) / weight_sum;
            return i;
        }
    }
}



int McadjNode::randomOutIndex_fast(double& probability)
{
    // Fast version of randomOutIndex, based on the assumption that
    // the weight sum is 1, and all weights are non-negative.

    double select = rand() / double(RAND_MAX);
    for (int i = 0; i < outNodes.size(); ++i) {
        select -= outWeights[i];
        if (select < 0) {
            probability = outWeights[i];
            return i;
        }
    }
}



void McadjWalker::randomWalk()
{
    double probability;
    if (currentState != TERMINAL_STATE) {
        int i;
        if (isFast)
            i = currentState->randomOutIndex_fast(probability);
        else
            i = currentState->randomOutIndex(probability);

        weight *= currentState->outWeights[i] / probability;

        lastState = currentState;
        currentState = currentState->outNodes[i];
    }
}

