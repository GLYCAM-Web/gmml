#ifndef TEMPLATEGRAPH_GRAPHSTRUCTURE_INCLUDE_EDGE_HPP
#define TEMPLATEGRAPH_GRAPHSTRUCTURE_INCLUDE_EDGE_HPP

//#include "../../LazyPrints/LazyPrinters.hpp"
#include "./GenericGraphObject.hpp"

#include <memory>

namespace glygraph
{
  template<class T>
  class Edge : public GenericGraphObject
  {
  public:
    /************************************************
     *  CONSTRUCTORS/DESTRUCTORS
     ***********************************************/
    Edge();
    Edge(std::string name_t, T* const &sourceNode_t, T* const &targetNode_t);
    Edge(std::string name_t, std::vector<std::string> labels_t, T* const &sourceNode_t,
         T* const &targetNode_t);

    // copy constructor
    Edge(const Edge<T> &rhs);

    // move constructor
    Edge(Edge<T> &&rhs);

    // copy assignment
    Edge<T> &operator=(const Edge<T> &rhs);

    // move assignment
    Edge<T> &operator=(Edge<T> &&rhs);

    virtual ~Edge();

    T* getTargetNode();
    T* getSourceNode();
    const T* getTargetNode() const;
    const T* getSourceNode() const;

    /************************************************
     *  GETTER/SETTER
     ***********************************************/
    // NOTE: Using shared pointer to get our source and sink in order to ensure
    // 			source and sink are good and alive.
    void setSourceNode(T *source_t);
    void setTargetNode(T *target_t);

  private:
    /************************************************
     *  ATTRIBUTES
     ***********************************************/
    // NOTE: Source node = the node that has a unique_ptr to this edge
    T* sourceNode_m;
    // NOTE: Sink node = the node that has a raw pointer to this edge
    T* targetNode_m;
  };

  template<class T>
  inline Edge<T>::Edge()
  {
    //badBehavior(__LINE__, __func__, "Warning calling default constructor");
    this->targetNode_m = nullptr;
    this->sourceNode_m = nullptr;
  }

  template<class T>
  inline Edge<T>::Edge(std::string name_t, T *const &sourceNode_t, T *const &targetNode_t)
      : GenericGraphObject(name_t)
  {
    this->targetNode_m = targetNode_t;
    this->sourceNode_m = sourceNode_t;
  }

  template<class T>
  inline Edge<T>::Edge(std::string name_t, std::vector<std::string> labels_t, T *const &sourceNode_t,
                       T *const &targetNode_t)
      : GenericGraphObject(name_t, labels_t)
  {
    this->targetNode_m = targetNode_t;
    this->sourceNode_m = sourceNode_t;
  }

  template<class T>
  inline Edge<T>::~Edge()
  {
    // have our edge destructor remove itself from our inList then let die
    this->targetNode_m->removeInEdge(this);
    //to help match outputs
    // std::stringstream logss;
    // logss << "Edge labeled " << this->getLabel() << ", with index " << this->getIndex() << " destroyed\n";
    // gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
    // lazyInfo(__LINE__, __func__,
    //			"Edge with name <" + this->getName() + "> deleted");
  }

  // copy constructor
  template<class T>
  inline Edge<T>::Edge(const Edge<T> &rhs)
      : GenericGraphObject(rhs.getName(), rhs.getLabels(), rhs.getConnectivityTypeIdentifier()),
        sourceNode_m(rhs.getSourceNode()), targetNode_m(rhs.getTargetNode())
  {
    // lazyInfo(__LINE__, __func__,
    //		"Calling copy constructor on " + this->getName());
  }

  // move constructor
  template<class T>
  inline Edge<T>::Edge(Edge<T> &&rhs)
      : GenericGraphObject(rhs.getName(), rhs.getLabels(), rhs.getConnectivityTypeIdentifier()),
        sourceNode_m(rhs.getSourceNode()), targetNode_m(rhs.getTargetNode())
  {
    // wanted data has been yoinked so we go ahead and delete this edge that we dont care about
    //	anymore. As stated in move assignment we dont care what state we leave our rhs in after a move
    // lazyInfo(__LINE__, __func__,
    //			"Calling move constructor on " + this->getName());

    rhs.getSourceNode()->removeOutEdge(rhs);
  }

  // copy assignment
  template<class T>
  inline Edge<T> &Edge<T>::operator=(const Edge<T> &rhs)
  {
    return *this = Edge<T>(rhs);
  }

  // move assignment
  template<class T>
  inline Edge<T> &Edge<T>::operator=(Edge<T> &&rhs)
  {
    // lazyInfo(__LINE__,__func__, "Edge move assignment");
    // Please note that in order to help prevent some bad behavior due to moving an edge
    //	causing bad connectivity to arise (i.e. multigraph creation, etc.) I am using
    //	the delete on move paradigm in order to help prevent this. Keep in mind move
    //	implies that we dont care about what happens to our rhs.
    this->sourceNode_m = rhs.sourceNode_m;
    this->targetNode_m = rhs.targetNode_m;
    this->setName(rhs.getName());
    this->setLabels(rhs.getLabels());

    // after we yoink data wanted from our rhs we go ahead and delete it
    rhs.getSourceNode()->removeOutEdge(rhs);

    return *this;
  }

  template<class T>
  inline void Edge<T>::setSourceNode(T *source_t)
  {
    this->sourceNode_m = source_t;
  }

  template<class T>
  inline void Edge<T>::setTargetNode(T *target_t)
  {
    this->targetNode_m = target_t;
  }

  template<class T>
  inline T* Edge<T>::getTargetNode()
  {
	  return this->targetNode_m;
  }

  template<class T>
  inline T* Edge<T>::getSourceNode()
  {
	  return this->sourceNode_m;
  }

  template<class T>
  inline const T* Edge<T>::getTargetNode() const
  {
    return this->targetNode_m;
  }

  template<class T>
  inline const T* Edge<T>::getSourceNode() const
  {
    return this->sourceNode_m;
  }

} // namespace temp_graph
#endif // end TEMPLATEGRAPH_GRAPHSTRUCTURE_INCLUDE_EDGE_HPP
