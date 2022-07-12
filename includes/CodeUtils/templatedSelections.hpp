#ifndef INCLUDES_CODEUTILS_TEMPLATEDSELECTIONS_HPP_
#define INCLUDES_CODEUTILS_TEMPLATEDSELECTIONS_HPP_

#include <vector>

namespace codeUtils
{
template <class RandomAccessIterator, class T>
bool isThingPresentInContainer( RandomAccessIterator first, RandomAccessIterator last, const T& value )
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

//template <class RandomAccessIterator, class T>
//std::vector<class *T> getThingsWithNameInContainer( RandomAccessIterator first, RandomAccessIterator last, const std::string& queryName )
//{
//	std::vector<class *T> foundThingsWithName;
//	while( first != last )
//	{
//		if( *first->getName() == queryName )
//		{
//			foundThingsWithName.push_back(*first);
//		}
//		++first;
//	}
//	return foundThingsWithName;
//}


//template <class RandomAccessIterator, class T>
//class T* findThingInContainer( RandomAccessIterator first, RandomAccessIterator last, const T& value )
//{
//    while( first != last ) {
//        if( *first == value ) {
//            return *first;
//        }
//        ++first;
//    }
//    return nullptr;
//}

} // namespace

#endif /* INCLUDES_CODEUTILS_TEMPLATEDSELECTIONS_HPP_ */
