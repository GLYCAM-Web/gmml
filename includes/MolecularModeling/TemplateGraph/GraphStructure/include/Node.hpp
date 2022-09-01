#ifndef TEMPLATEGRAPH_GRAPHSTRUCTURE_INCLUDE_NODE_HPP
#define TEMPLATEGRAPH_GRAPHSTRUCTURE_INCLUDE_NODE_HPP

//#include "../../LazyPrints/LazyPrinters.hpp"
#include "./Edge.hpp"
#include "./GenericGraphObject.hpp"

#include <memory>
#include <unordered_set>

namespace glygraph
{
  template<class T>
  class Node : public GenericGraphObject, public std::enable_shared_from_this<Node<T>>
  {
  public:
    /************************************************
     *  CONSTRUCTORS/DESTRUCTORS
     ***********************************************/
    Node();
    Node(std::string name_t);
    Node(std::string name_t, std::string label_t);

    // copy constructor
    Node(const Node<T> &rhs);

    // move constructor
    Node(Node<T> &&rhs);

    // copy assignment
    Node<T> &operator=(const Node<T> &rhs);

    // move assignment
    Node<T> &operator=(Node<T> &&rhs);

    virtual ~Node();

    /************************************************
     *  GETTER/SETTER
     ***********************************************/

    inline T *getDeriviedClass()
    {
      auto derived = static_cast<T*>(this);
      return derived;
    }

    std::vector<T*> getNeighbors();

    std::vector<Edge<T>*> getEdges() const;
    std::vector<Edge<T>*> getOutEdges() const;
    std::vector<Edge<T>*> getInEdges() const;

    /************************************************
     *  MUTATORS
     ***********************************************/
    void addEdge(Edge<T> *edgeToAdd_t);
    /* TODO: Finalize how we would like to add nodes to one another. The addNeighbor will just be a
     * 			wrapper for our addChild. Please note how I avoided having an "addEdge" because from
     * 			my understanding and what I am pretty sure our use will be every single edge is
     * 			guaranteed to have both a source and sink node.
     */
    void addNeighbor(std::string edgeName_t, T* const &newNeighbor_t);

    void addChild(std::string edgeName_t, T* const &childNode_t);

    void addParent(std::string edgeName_t, T* const &parentNode_t);
    /* NOTE: We MUST remove the edge from out "child" node BEFORE deleting the edge by deleting the
     * 			unique_ptr that owns it. This is handled in edge's destructor, all we have to worry
     * 			about is deleting the unique ptr that owns our edge.
     */
    void removeEdgeBetween(T* const &otherNode_t);

    void removeInEdge(Edge<T> *edgeToRemove_t);
    void removeOutEdge(Edge<T> *edgeToRemove_t);

    /************************************************
     *  FUNCTIONS
     ***********************************************/
    bool     isNeighbor(T* const &otherNode_t);
    Edge<T> *getConnectingEdge(T* const &otherNode_t);

	/************************************************
	 *  LAMBDAS
	 ***********************************************/
	// Implemented for sorting the vector of incoming edges by the source objects < operator. Uses a lambda function.
	inline void sortInEdgesBySourceTObjectComparator()
	{
		std::sort(inEdges_m.begin(), inEdges_m.end(),
				[](Edge<T>* e1, Edge<T>* e2)
				{ // Lambda function for doing the sort.
							return ( *(e1->getSourceNode()->getDeriviedClass()) > *(e2->getSourceNode()->getDeriviedClass()) );
				});
		return;
	}

    std::vector<T* > getChildren();
    std::vector<T*> getParents();

  private:
    /************************************************
     *  CONSTRUCTORS/DESTRUCTORS
     ***********************************************/
// ToDo: Make all of these private and have T be friend. // Safeguards mistake like class prepAtom : public Node<pdbAtom>;
    //    Node (): GenericGraphObject("INVALID NODE") {};
//    friend T;
    /************************************************
     *  ATTRIBUTES
     ***********************************************/
    std::vector<std::unique_ptr<Edge<T>>> outEdges_m;
    std::vector<Edge<T> *>                inEdges_m;

    /************************************************
     *  FUNCTIONS
     ***********************************************/
    void edgeConnectionUpdate();

    bool isChildOf(T* const &possibleParent_t);
    bool isParentOf(T* const &possibleChild_t);

    friend class Edge<T>;
  };

  template<class T>
  inline Node<T>::Node() : GenericGraphObject("INVALID NODE") {};

  template<class T>
  inline Node<T>::Node(std::string name_t) : GenericGraphObject(name_t)
  {
    // lazyInfo(__LINE__, __func__,
    //		"Created node with name <" + this->getName() + ">");
  }

  template<class T>
  inline Node<T>::Node(std::string name_t, std::string label_t) : GenericGraphObject(name_t, label_t)
  {
    // lazyInfo(__LINE__, __func__,
    //		"Created node with name <" + this->getName()
    //				+ ">\n\tAnd with label <" + this->getLabel() + ">");
  }

  template<class T>
  inline Node<T>::~Node()
  {
    std::vector<Edge<T> *> tempInEdge = this->inEdges_m;
    //lazyInfo(__LINE__, __func__, "Destroying Node: " + this->getName());
    // go through and hit all our parents, i.e. the ones that own the incoming edge and delete them
    // TODO: Do this but not lazy
    this->outEdges_m.clear();
    for (Edge<T> *currInEdge : tempInEdge)
      {
        currInEdge->getSourceNode()->removeOutEdge(currInEdge);
      }
    tempInEdge.clear();
    // gmml::log(__LINE__, __FILE__, gmml::INF, ("Node labelled " + this->getLabel() + " destroyed\n"));
  }

  // Copy constructor
  template<class T>
  inline Node<T>::Node(const Node<T> &rhs)
      : GenericGraphObject(rhs.getName(), rhs.getLabels(), rhs.getConnectivityTypeIdentifier())
  {
    // std::cout << "\n\tGiven object ptr: " << rhs.objectPtr << "\n\n";
    // lazyInfo(__LINE__, __func__, "Calling copy constructor");
    for (Edge<T> const *currInEdge : rhs.inEdges_m)
      {
    	//ToDo Check if the "new" is best. Shouldn't Edge class copy constructor not handle some of this?
        std::unique_ptr<Edge<T>> tempIn(new Edge<T>(*currInEdge));

        tempIn.get()->setTargetNode(this->getDeriviedClass());

        this->inEdges_m.push_back(tempIn.get());
        tempIn.get()->getSourceNode()->outEdges_m.push_back(std::move(tempIn));
      }
    for (std::unique_ptr<Edge<T>> const &currOutEdge : rhs.outEdges_m)
      {
        std::unique_ptr<Edge<T>> tempOut(new Edge<T>(*currOutEdge.get()));

        tempOut.get()->setSourceNode(this->getDeriviedClass());

        tempOut.get()->getTargetNode()->inEdges_m.push_back(tempOut.get());
        this->outEdges_m.push_back(std::move(tempOut));
      }
  }

  // move constructor
  template<class T>
  inline Node<T>::Node(Node<T> &&rhs)
      : GenericGraphObject(rhs.getName(), rhs.getLabels(), rhs.getConnectivityTypeIdentifier()),
        outEdges_m(std::move(rhs.outEdges_m)), inEdges_m(std::move(rhs.inEdges_m))
  {
    // lazyInfo(__LINE__, __func__, "Calling node move constructor");
    this->edgeConnectionUpdate();
  }

  // copy assignment
  template<class T>
  inline Node<T> &Node<T>::operator=(const Node<T> &rhs)
  {
    return *this = Node<T>(rhs);
  }

  // move assignment
  template<class T>
  inline Node<T> &Node<T>::operator=(Node<T> &&rhs)
  {
    // lazyInfo(__LINE__, __func__, "Calling node move assignment");
    this->setName(rhs.getName());
    this->setLabels(rhs.getLabels());
    this->setConnectivityTypeIdentifier(rhs.getConnectivityTypeIdentifier());

    this->inEdges_m  = std::move(rhs.inEdges_m);
    this->outEdges_m = std::move(rhs.outEdges_m);
    this->edgeConnectionUpdate();

    return *this;
  }

  template<class T>
  inline std::vector<T* > Node<T>::getNeighbors()
  {
    std::vector<T* > childrenVec = this->getChildren();
    std::vector<T* > parentsVec  = this->getParents();
    parentsVec.insert(parentsVec.end(), childrenVec.begin(), childrenVec.end());

    return parentsVec;
  }

  template<class T>
  inline std::vector<Edge<T> *> Node<T>::getEdges() const
  {
    std::vector<Edge<T> *> outEdgesVec = this->getOutEdges();
    std::vector<Edge<T> *> inEdgesVec  = this->getInEdges();
    outEdgesVec.insert(outEdgesVec.end(), inEdgesVec.begin(), inEdgesVec.end());
    return outEdgesVec;
  }

  template<class T>
  inline std::vector<Edge<T> *> Node<T>::getOutEdges() const
  {
    std::vector<Edge<T> *> outEdgeVecToReturn;
    for (std::unique_ptr<Edge<T>> const &currOutEdge : this->outEdges_m)
      {
        outEdgeVecToReturn.push_back(currOutEdge.get());
      }
    return outEdgeVecToReturn;
  }

  template<class T>
  inline std::vector<Edge<T> *> Node<T>::getInEdges() const
  {
    return this->inEdges_m;
  }

  template<class T>
  inline void Node<T>::addNeighbor(std::string edgeName_t, T* const &newNeighbor_t)
  {
    this->addChild(edgeName_t, newNeighbor_t);
  }

  template<class T>
  inline void Node<T>::addChild(std::string edgeName_t, T *const &childNode_t)
  {
    if (this->isNeighbor(childNode_t))
      {
        //badBehavior(__LINE__, __func__, "Trying to make create an edge between two nodes that are already neighbors");
      }
    else if (this == childNode_t)
      {
        //badBehavior(__LINE__, __func__, "Trying to add self as child, stop that!");
      }
    else
      {
        std::unique_ptr<Edge<T>> tempEdge(new Edge<T>(edgeName_t, this->getDeriviedClass(), childNode_t));

        childNode_t->inEdges_m.push_back(tempEdge.get());

        this->outEdges_m.push_back(std::move(tempEdge));
      }
  }

  template<class T>
  inline void Node<T>::addParent(std::string edgeName_t, T* const &parentNode_t)
  {
    parentNode_t->addChild(edgeName_t, this->getDeriviedClass());
  }

  template<class T>
  inline void Node<T>::removeEdgeBetween(T* const &otherNode_t)
  {
    if (this->isNeighbor(otherNode_t))
      {
        if (this->isChildOf(otherNode_t))
          {
            Edge<T> *edgeToRemove = this->getConnectingEdge(otherNode_t);
            otherNode_t->removeOutEdge(edgeToRemove);
          }
        else if (this->isParentOf(otherNode_t))
          {
            otherNode_t->removeEdgeBetween(this->getDeriviedClass());
          }
        else
          {
            //badBehavior(__LINE__, __func__, "Tried to remove an edge for some reason it isnt a parent or child");
          }
      }
    else
      {
        //badBehavior(__LINE__, __func__,
        //            "Tried to remove an edge between node <" + this->getName() + "> and node <" +
        //                otherNode_t->getName() + ">");
      }
  }

  template<class T>
  inline bool Node<T>::isNeighbor(T* const &otherNode_t)
  {
    return this->isChildOf(otherNode_t) || this->isParentOf(otherNode_t);
  }

  template<class T>
  bool Node<T>::isChildOf(T* const &possibleParent_t)
  {
    for (Edge<T> *currInEdge : this->inEdges_m)
      {
        if (currInEdge->getSourceNode() == possibleParent_t)
          {
            return true;
          }
      }
    return false;
  }

  template<class T>
  inline bool Node<T>::isParentOf(T* const &possibleChild_t)
  {
    for (std::unique_ptr<Edge<T>> const &currOutEdge : this->outEdges_m)
      {
        if (currOutEdge.get()->getTargetNode() == possibleChild_t)
          {
            return true;
          }
      }
    return false;
  }

  template<class T>
  inline std::vector<T* > Node<T>::getChildren()
  {
    std::vector<T* > childrenVecToReturn;
    for (std::unique_ptr<Edge<T>> const &currOutEdge : this->outEdges_m)
      {
        childrenVecToReturn.push_back(currOutEdge.get()->getTargetNode());
        // lazyInfo(__LINE__, __func__,
        //		currOutEdge.get()->getTargetNode()->getName());
        // childrenVecToReturn.push_back(currOutEdge.get()->getTargetNode().shared_from_this());
      }
    return childrenVecToReturn;
  }

  template<class T>
  inline void Node<T>::edgeConnectionUpdate()
  {
    // So even tho we moved all the references to an edge, we have
    //	not updated said edges with their new vertex appropriately
    for (Edge<T> *currInEdge : this->inEdges_m)
      {
        currInEdge->setTargetNode(this->getDeriviedClass());
      }
    for (std::unique_ptr<Edge<T>> &currOutEdge : this->outEdges_m)
      {
        currOutEdge.get()->setSourceNode(this->getDeriviedClass());
      }
  }

  template<class T>
  inline std::vector<T*> Node<T>::getParents()
  {
    std::vector<T*> parentsVecToReturn;
    for (Edge<T> *currInEdge : this->inEdges_m)
      {
        parentsVecToReturn.push_back(currInEdge->getSourceNode());
      }
    return parentsVecToReturn;
  }

  template<class T>
  inline void Node<T>::removeInEdge(Edge<T> *edgeToRemove_t)
  {
    // lazyInfo(__LINE__, __func__,
    //		"removing in edge <" + edgeToRemove->getName()
    //				+ "> from node named <" + this->getName() + ">");
    this->inEdges_m.erase(std::remove(this->inEdges_m.begin(), this->inEdges_m.end(), edgeToRemove_t),
                          this->inEdges_m.end());
  }

  template<class T>
  inline void Node<T>::removeOutEdge(Edge<T> *edgeToRemove_t)
  {
    // lazyInfo(__LINE__, __func__,
    //			"removing out edge <" + edgeToRemove->getName()
    //	 				+ "> from node named <" + this->getName() + ">");
    for (unsigned int outIndex = 0; outIndex <= this->outEdges_m.size(); outIndex++)
      {
        if (this->outEdges_m[outIndex].get() == edgeToRemove_t)
          {
            this->outEdges_m.erase(this->outEdges_m.begin() + outIndex);
          }
      }
  }

  template<class T>
  inline Edge<T> *glygraph::Node<T>::getConnectingEdge(T* const &otherNode_t)
  {
    if (this->isNeighbor(otherNode_t))
      {
        if (this->isChildOf(otherNode_t))
          {
            for (Edge<T> *currInEdge : this->inEdges_m)
              {
                if (currInEdge->getSourceNode() == otherNode_t)
                  {
                    return currInEdge;
                  }
              }
          }
        else if (this->isParentOf(otherNode_t))
          {
            for (std::unique_ptr<Edge<T>> const &currOutEdge : this->outEdges_m)
              {
                if (currOutEdge.get()->getTargetNode() == otherNode_t)
                  {
                    return currOutEdge.get();
                  }
              }
          }

        // NOTE: Should never hit here
        //badBehavior(__LINE__, __func__, "Found that two nodes were neihgbors but could not find the connecting edge");
      }
    else
      {
        //badBehavior(__LINE__, __func__, "Tried to get a connecting edge whene the two nodes are not connected");
        return NULL;
      }
    return NULL;
  }

} // namespace temp_graph

#endif // end TEMPLATEGRAPH_GRAPHSTRUCTURE_INCLUDE_NODE_HPP
