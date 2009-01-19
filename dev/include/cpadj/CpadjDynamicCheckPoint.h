// DynamicCheckPoints template that contains the optimal checkpointing
// algorithm.
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



#ifndef WANGQ_CPADJ_CPADJDYNAMICCHECKPOINT_H
#define WANGQ_CPADJ_CPADJDYNAMICCHECKPOINT_H

#include <vector>
#include <iostream>



namespace wangq
{
    namespace cpadj
    {
        template<class C>
        class DynamicCheckPoints {
            public:
                DynamicCheckPoints(int memory_size=5);

                const C& operator[](int index);
                void push(const C& c);
                const C& pop();

                int length() {
                    return length_;
                }
                int numCheckPoints() {
                    assert (cpCheckPoints_.size() == cpLevels_.size());
                    assert (cpIndices_.size() + 1 == cpLevels_.size());
                    return cpCheckPoints_.size() - 1;
                }
                int memorySize() {
                    return memorySize_;
                }
                int checkPoint(int i) {
                    assert (i >= 0 and i <= numCheckPoints());
                    return cpCheckPoints_[i];
                }
                int level(int i) {
                    assert (i >= 0 and i <= numCheckPoints());
                    return cpLevels_[i];
                }

            private:
                std::vector<C> cpMemory_;
                std::vector<int> cpIndices_;
                std::vector<int> cpCheckPoints_;
                std::vector<int> cpLevels_;

                int memorySize_;
                int length_;
        };



        // Implimentation
        
        template<class C>
        DynamicCheckPoints<C>::DynamicCheckPoints(int memory_size)
        {
            assert (memory_size >= 1);
            memorySize_ = memory_size;
            length_ = 0;

            cpCheckPoints_.push_back(0);
            cpLevels_.push_back(0);
        }
        
        
        
        template<class C>
        void DynamicCheckPoints<C>::push(const C& c)
        {
            assert (cpIndices_.size() + 1 == cpCheckPoints_.size());
            assert (cpIndices_.size() + 1 == cpLevels_.size());

            std::vector<bool> dispensable(cpLevels_.size() - 1, false);
            int max_level = cpLevels_.back();
            for (int i = cpLevels_.size() - 2; i > 0; -- i) {
                dispensable[i] = (max_level > cpLevels_[i]);
                max_level = std::max (max_level, cpLevels_[i]);
            }

            ++length_;
            cpCheckPoints_.push_back(length_);

            // Enough memory for more checkpoints.
            if (cpCheckPoints_.size() <= memorySize_ + 1) {
                cpLevels_.push_back (0);
                if (cpMemory_.size() <= memorySize_) {
                    cpIndices_.push_back (cpMemory_.size());
                    cpMemory_.push_back (c);
                }
                else {
                    // Find an empty memory spot.
                    std::vector<bool> empty(memorySize_, true);
                    for (int i = 0; i < cpIndices_.size(); ++ i)
                        empty [cpIndices_[i]] = false;
                    int i_empty = 0;
                    while (not empty[i_empty]) ++ i_empty;

                    assert (i_empty < memorySize_);
                    cpIndices_.push_back(i_empty);
                    cpMemory_[i_empty] = c;
                }
            }
            else {
                // Search for a dispensable checkpoint.
                int i_recycle = 0;
                for (int i = 1; i < cpLevels_.size(); ++ i) {
                    if (dispensable[i])
                        i_recycle = i;
                }

                if (i_recycle > 0) {
                    cpLevels_.push_back(0);
                }
                else {
                    // no dispensable checkpoint, add higher level checkpoint.
                    int lvl = cpLevels_.back();
                    i_recycle = cpLevels_.size() - 1;
                    // new checkpoint should be of one level higher.
                    cpLevels_.push_back(lvl + 1);
                }

                // Recycle the old checkpoint.
                cpCheckPoints_.erase (cpCheckPoints_.begin() + i_recycle);
                cpLevels_.erase (cpLevels_.begin() + i_recycle);

                if (i_recycle < memorySize_) {
                    int i_deleted = cpIndices_[i_recycle];
                    int i_before_deleted = cpIndices_[i_recycle - 1];
                    cpMemory_[i_before_deleted].mergeWith(cpMemory_[i_deleted]);
                    cpMemory_[i_deleted] = c;

                    cpIndices_.erase (cpIndices_.begin() + i_recycle);
                    cpIndices_.push_back (i_deleted);
                }
                else {  // Second last checkpoint deleted.
                    cpMemory_[cpIndices_.back()].mergeWith(c);
                }

                // Update the level of the initial checkpoint.
                if (cpLevels_.back() > cpLevels_[0])
                    cpLevels_[0] = cpLevels_.back();
            }

            // std::cout << length_ << "    ";
            // for (int i = 0; i < cpCheckPoints_.size(); ++ i) {
            //     std::cout << cpCheckPoints_[i] << ' ';
            // }
            // std::cout << '\n';
        }
        
        
        
        template<class C>
        const C& DynamicCheckPoints<C>::pop()
        {
            assert (not cpIndices_.empty());

            assert (cpIndices_.size() + 1 == cpCheckPoints_.size());
            assert (cpIndices_.size() + 1 == cpLevels_.size());

            static C c;
            c = cpMemory_[cpIndices_.back()];

            cpIndices_.pop_back();
            cpCheckPoints_.pop_back();
            cpLevels_.pop_back();

            length_ = cpCheckPoints_.back();

            // std::cout << length_ << " -- ";
            // for (int i = 0; i < cpCheckPoints_.size(); ++ i) {
            //     std::cout << cpCheckPoints_[i] << ' ';
            // }
            // std::cout << '\n';

            return c;
        }



        template<class C>
        const C& DynamicCheckPoints<C>::operator[](int index)
        {
            assert (index < cpIndices_.size());
            return cpMemory_[cpIndices_[index]];
        }

    }   // namespace cpadj
}   // namespace wangq



#endif

