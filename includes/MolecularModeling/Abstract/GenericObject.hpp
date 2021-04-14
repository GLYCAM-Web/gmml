#ifndef GENERIC_OBJECT_HPP
#define GENERIC_OBJECT_HPP

#include "Labels.hpp"
#include "Index.hpp"
#include "Visitors.hpp"

namespace Abstract
{ // An inheritable composite that gives Labels, Index, Visitors attributes etc.
	class GenericObject 
	{
	public:
		//////////////////////////////////////////////////////////
		//                       CONSTRUCTOR                    //
		//////////////////////////////////////////////////////////
		GenericObject() {};
		GenericObject(std::string label) {this->AddLabel(label);}
		//////////////////////////////////////////////////////////
		//                       ACCESSOR                       //
		//////////////////////////////////////////////////////////
		// Visitors
		inline bool GetIsVisitedBy(std::string visitor) {return visitors_.GetIsVisitedBy(visitor);}
		// Labels
		inline std::vector<std::string> GetLabels() {return labels_.GetLabels();}
		inline std::string GetLabel() {return labels_.GetLabel();}
		// Index
		inline unsigned long long GetIndex() {return index_.GetIndex();}
		//////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        // Visitors
		inline void AddVisitor(std::string visitor = "") {visitors_.AddVisitor(visitor);}
		inline void RemoveVisitor(std::string visitor) {visitors_.RemoveVisitor(visitor);}
		// Labels
		inline void SetLabels(std::vector<std::string> labels) {labels_.SetLabels(labels);}
		inline void AddLabel(std::string label) {labels_.AddLabel(label);}
 		// Index
		inline void SetIndex(unsigned long long index) {index_.SetIndex(index);}
		//////////////////////////////////////////////////////////
        //                      FUNCTIONS                       //
        //////////////////////////////////////////////////////////
		// Labels
		inline bool CompareLabels(const std::vector<std::string> otherLabels) {return labels_.CompareLabels(otherLabels);}
		inline std::string FindLabelContaining(const std::string query) {return labels_.FindLabelContaining(query);}

		inline std::string Print() {return this->GetLabel();}
		//////////////////////////////////////////////////////////
        //                  OPERATOR OVERLOADING                //
        //////////////////////////////////////////////////////////
        bool operator== ( GenericObject& rhs)  { return (this->GetIndex() == rhs.GetIndex());}
        bool operator!= ( GenericObject& rhs)  { return (this->GetIndex() != rhs.GetIndex());}
	private:
		//////////////////////////////////////////////////////////
        //                       ATTRIBUTES                     //
        //////////////////////////////////////////////////////////
		Index index_;
		Labels labels_;
		Visitors visitors_;
		
	};
}
#endif // GENERIC_OBJECT_HPP