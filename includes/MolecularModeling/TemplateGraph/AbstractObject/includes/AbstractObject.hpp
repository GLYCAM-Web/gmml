#ifndef ABSTRACTOBJECT_INCLUDES_GENERIC_OBJECT_HPP
#define ABSTRACTOBJECT_INCLUDES_GENERIC_OBJECT_HPP

#include "Labels.hpp"
#include "Index.hpp"

namespace abstrab
{
/* TODO: Figure out a better name and what data we want to include. Do
 * 			we want/need index counting? This could be a pain when
 * 			deletions and insertions come into play. All my algos (except
 * 			subgraph matching) ignore these signifiers.
 * 			We could name it "generic identifiers" or something of the
 * 			sort.
 */
class AbstractObject : public Labels, public Index
{
public:
	/************************************************
	 *  CONSTRUCTORS/DESTRUCTORS
	 ***********************************************/

	inline AbstractObject(std::string name) :
		Labels(name), Index()
	{
	}

	inline AbstractObject(std::string name, std::string label) :
		Labels(name, label), Index()
	{
	}

	inline AbstractObject(std::string name, unsigned long long index) :
			Labels(name), Index(index)
	{
	}

	inline AbstractObject(std::string name, std::vector<std::string> labels) :
		Labels(name, labels), Index()
	{
	}

	//copy constructor
	inline AbstractObject(const AbstractObject &rhs) :
		Labels(rhs.getName(), rhs.getLabels()), Index(rhs.getIndex())
	{
	}

	//move constructor
	inline AbstractObject(AbstractObject &&rhs) :
			Labels(rhs.getName(), rhs.getLabels()), Index(rhs.getIndex())
	{
	}

	//copy assignment
	inline AbstractObject& operator=(const AbstractObject &rhs)
	{
		return *this;
	}

	//move assignment
	inline AbstractObject& operator=(AbstractObject &&rhs)
	{
		return *this;
	}


	virtual ~AbstractObject();


	/************************************************
	 *  FUNCTIONS
	 ***********************************************/
	inline bool compareLabels(const std::vector<std::string> otherLabels)
	{
		return this->Labels::compareLabels(otherLabels, 1);
	}

private:
	/************************************************
	 *  ATTRIBUTES
	 ***********************************************/


};

inline AbstractObject::~AbstractObject()
{
}

}
#endif // ABSTRACT_OBJECT_HPP
