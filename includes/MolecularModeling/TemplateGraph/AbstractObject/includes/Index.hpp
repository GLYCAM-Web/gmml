#ifndef ABSTRACTOBJECT_INCLUDES_INDEX_HPP
#define ABSTRACTOBJECT_INCLUDES_INDEX_HPP

#include <string>
#include <vector>

namespace abstrab
{
    class Index
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        Index()
        {
            this->setIndex(this->generateIndex());
        }

        Index(unsigned int index)
        {
            this->setIndex(index);
        }

        // copy constructor
        inline Index(const Index& rhs) : index_m(rhs.index_m)
        {}

        // move constructor
        inline Index(Index&& rhs) : index_m(rhs.index_m)
        {}

        inline Index& operator=(const Index& rhs)
        {
            return *this = Index(rhs);
        }

        // move assignment
        inline Index& operator=(Index&& rhs)
        {
            this->index_m = rhs.index_m;
            return *this;
        }

        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        inline unsigned int getIndex() const
        {
            return index_m;
        }

        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        inline void setIndex(unsigned int index)
        {
            index_m = index;
        }

      private:
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        inline unsigned int generateIndex()
        {
            static unsigned int s_NodeIndex =
                0; // static keyword means it is created only once and persists beyond scope of code block.
            return s_NodeIndex++; // makes copy of index, increments the real index, then returns the value in the copy
        }

        //////////////////////////////////////////////////////////
        //                       ATTRIBUTES                     //
        //////////////////////////////////////////////////////////
        unsigned int index_m;
    };

} // namespace abstrab
#endif // INDEX_HPP
