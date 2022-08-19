#ifndef INCLUDES_CODEUTILS_TEMPLATEDSELECTIONS_HPP_
#define INCLUDES_CODEUTILS_TEMPLATEDSELECTIONS_HPP_

#include <vector>

namespace codeUtils
{
template <class RandomAccessIterator, class T>
bool isElementPresent( RandomAccessIterator first, RandomAccessIterator last, const T& value )
{
	while( first != last )
	{
		if( *first == value )
		{
			return true;
		}
		++first;
	}
	return false;
}

// Must pass in T in order for the return type to be deduced by the compiler. There may be a better way but I can't find it. Tried passing in T typed iterators, but no deal as the containers are templated too
// This works fine, but the extra T thing makes it real ugly and you get warnings about it not being used.
// My compromise is to use the function below, but it loses the ability to pass in iterators to any container type, now it's just vector.
// Perhaps if I can figure out how to pass in a lamda?(nope, member functions need binding I think) indicating the member function to use e.g. getName vs getNumber, it could
// then deduce the type of T from that and we could go back to this one, and futher reduce code duplication.
//template <class RandomAccessIterator, class T>
//std::vector<T*> getElementsWithName( RandomAccessIterator first, RandomAccessIterator last, std::vector<std::string> queryNames, const T& value)

template <class T>
std::vector<T*> getElementsWithNames(const std::vector<T*>& inputVector, const std::vector<std::string>& queryNames)
{
	std::vector<T*> results;
	for (auto & element : inputVector)
	{
		if (std::find(queryNames.begin(), queryNames.end(), element->getName()) != queryNames.end())
		{
			results.push_back(element);
		}
	}
	return results;
}

template <class T>
const T* findElementWithName(const std::vector<const T*>& inputVector, const std::string& queryName)
{
	for (const auto & element : inputVector)
	{
		if ( element->getName() == queryName )
		{
			return element;
		}
	}
	return nullptr;
}

template <class T>
T* findElementWithName(const std::vector<T*>& inputVector, const std::string& queryName)
{
	for (auto & element : inputVector)
	{
		if ( element->getName() == queryName )
		{
			return element;
		}
	}
	return nullptr;
}

template <class T>
const T* findElementWithNumber(const std::vector<const T*>& inputVector, const int& queryNumber)
{
	for (const auto & element : inputVector)
	{
		if ( element->getNumber() == queryNumber )
		{
			return element;
		}
	}
	return nullptr;
}

} // namespace

#endif /* INCLUDES_CODEUTILS_TEMPLATEDSELECTIONS_HPP_ */
