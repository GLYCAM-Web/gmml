#ifndef INCLUDES_CODEUTILS_TEMPLATEDSELECTIONS_HPP_
#define INCLUDES_CODEUTILS_TEMPLATEDSELECTIONS_HPP_

#include <vector>
#include <algorithm> // find
//ToDo move this to CentralDataStructure/Selections and rename the namespace.
// Must pass in T in order for the return type to be deduced by the compiler. There may be a better way but I can't find it. Tried passing in T typed iterators, but no deal as the containers are templated too
// This works fine, but the extra T thing makes it real ugly and you get warnings about it not being used.
// My compromise is to use the function below, but it loses the ability to pass in iterators to any container type, now it's just vector.
// Perhaps if I can figure out how to pass in a lamda?(nope, member functions need binding I think) indicating the member function to use e.g. getName vs getNumber, it could
// then deduce the type of T from that and we could go back to this one, and further reduce code duplication.
//template <class RandomAccessIterator, class T>
//std::vector<T*> getElementsWithName( RandomAccessIterator first, RandomAccessIterator last, std::vector<std::string> queryNames, const T& value)
namespace codeUtils
{
//template <class T> // Here type deduction is funky, because T can be both const and non-const in inputs, so it can't deduce properly. Couldn't find a reasonable solution.
//bool isElementPresent( const std::vector<T>& inputVector, const T query )
//{
//    for (auto & element : inputVector)
//    {
//		if( element == query )
//		{
//			return true;
//		}
//	}
//	return false;
//}

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

// ToDo: This const version may be unnecessary. Can pass in const types to the next one.
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
T* findElementWithNumber(const std::vector<T*>& inputVector, const int& queryNumber)
{
	for (auto & element : inputVector)
	{
		if ( element->getNumber() == queryNumber )
		{
			return element;
		}
	}
	return nullptr;
}

template <class T>
std::vector<T> findElementsNotInVector(const std::vector<T>& inputVector, const std::vector<T>& excludeElements)
{
    std::vector<T> elementsInInputVectorButNotInQueryElements;
    for (auto element : inputVector)
    { // if element is not in the exclude list
        if (std::find(excludeElements.begin(), excludeElements.end(), element) == excludeElements.end())
        {
            elementsInInputVectorButNotInQueryElements.push_back(element);
        }
    }
    return elementsInInputVectorButNotInQueryElements;
}

} // namespace

#endif /* INCLUDES_CODEUTILS_TEMPLATEDSELECTIONS_HPP_ */
