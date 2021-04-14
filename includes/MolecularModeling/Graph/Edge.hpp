#ifndef T_EDGE_HPP
#define T_EDGE_HPP

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include "includes/MolecularModeling/Abstract/GenericObject.hpp"

namespace TemplateGraph
{
	template <class T> class Node;
	template <class T> class Edge; // forward declare 
	template <class T>
	class Edge : public Abstract::GenericObject
	{
	public:
		//////////////////////////////////////////////////////////
		//                       CONSTRUCTOR                    //
		//////////////////////////////////////////////////////////
		Edge() : GenericObject {} {};
		Edge(Node<T>* source, Node<T>* target, std::string label = ""); // Only Node should call this. How to enforce?
		~Edge() { std::cout << "Edge labeled " << this->GetLabel() << ", with index " << this->GetIndex() << " destroyed\n";}
		//////////////////////////////////////////////////////////
		//                       ACCESSOR                       //
		//////////////////////////////////////////////////////////	
		inline Node<T>* GetSource() {return source_;}
		inline Node<T>* GetTarget() {return target_;}
		//////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
		bool CompareEdgeAndNodeLabels(Edge<T>* otherEdge);
	private:
		//////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
		inline void SetSource(Node<T>* source) {source_ = source;}
		inline void SetTarget(Node<T>* target) {target_ = target;}       		
		//////////////////////////////////////////////////////////
        //                       ATTRIBUTES                     //
        //////////////////////////////////////////////////////////
        Node<T>* source_;
		Node<T>* target_;
		//////////////////////////////////////////////////////////
        //                       FRIENDS                        //
        //////////////////////////////////////////////////////////
		friend class Node<T>; // Allows Node to access private stuff like GetSource().
	};
	//////////////////////////////////////////////////////////
    //                       DEFINITIONS                    //
    //////////////////////////////////////////////////////////
template <typename T> 
	Edge<T>::Edge(Node<T>* source, Node<T>* target, std::string label)
	{
		this->SetSource(source);
		this->SetTarget(target);
		this->AddLabel(label);
		//std::cout << "Edge labeled " << this->GetLabel() << ", with index " << this->GetIndex() << " constructed\n";
	}

template <typename T>  
	bool Edge<T>::CompareEdgeAndNodeLabels(Edge<T>* otherEdge)
	{ // If any label here matches any in other label, return true
		if(this->CompareLabels(otherEdge->GetLabels()) 
			&& this->GetSource()->CompareLabels(otherEdge->GetSource()->GetLabels())
			&& this->GetTarget()->CompareLabels(otherEdge->GetTarget()->GetLabels()) )
				return true;
		return false;
	}

}// TemplateGraph namespace
#endif // T_EDGE_HPP