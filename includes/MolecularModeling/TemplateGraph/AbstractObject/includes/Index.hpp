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
	Index(unsigned long long index)
	{
		this->setIndex(index);
	}

	//copy constructor
	inline Index(const Index &rhs) :
			index_m(rhs.index_m)
	{
	}

	//move constructor
	inline Index(Index &&rhs) :
			index_m(rhs.index_m)
	{
	}

	//copy assignment
	inline Index& operator=(const Index &rhs)
	{
		this->index_m = rhs.index_m;
		return *this;
	}

	//move assignment
	inline Index& operator=(Index &&rhs)
	{
		this->index_m = rhs.index_m;
		return *this;
	}

	//////////////////////////////////////////////////////////
	//                       ACCESSOR                       //
	//////////////////////////////////////////////////////////
	inline unsigned long long getIndex() const
	{
		return index_m;
	}
	//////////////////////////////////////////////////////////
	//                       MUTATOR                        //
	//////////////////////////////////////////////////////////
	inline void setIndex(unsigned long long index)
	{
		index_m = index;
	}
private:
	//////////////////////////////////////////////////////////
	//                       FUNCTIONS                      //
	//////////////////////////////////////////////////////////
	inline unsigned long long generateIndex()
	{
		static unsigned long long s_NodeIndex = 0; // static keyword means it is created only once and persists beyond scope of code block.
		return s_NodeIndex++; // makes copy of index, increments the real index, then returns the value in the copy
	}
	//////////////////////////////////////////////////////////
	//                       ATTRIBUTES                     //
	//////////////////////////////////////////////////////////
	unsigned long long index_m;
};

}
#endif // INDEX_HPP
