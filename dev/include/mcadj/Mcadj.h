// Basic Monte Carlo adjoint solver functionality.
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



#ifndef WANGQ_MCADJ_MCADJ_H
#define WANGQ_MCADJ_MCADJ_H

#include <vector>

#define TERMINAL_STATE 0



namespace wangq
{
    namespace mcadj
    {
        class McadjNode;
        
        class McadjNode{
            public:
                int id;
                std::vector<McadjNode*> outNodes;
                std::vector<double> outWeights;
                int randomOutIndex(double& probability);
                int randomOutIndex_fast(double& probability);
        };
        
        class McadjWalker{
            public:
                McadjWalker() {
                    isFast = false;
                    weight = 1.0;
                }
                McadjNode* currentState;
                McadjNode* lastState;
                McadjNode* initialState;
                double weight;
                bool isFast;
                void randomWalk();
        };

    }   // namespace mcadj
}   // namespace wangq



#endif
