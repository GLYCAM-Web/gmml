#ifndef TEMPLATEGRAPH_GRAPHSTRUCTURE_INCLUDE_HALF_ADJACENCY_MATRIX_HPP
#define TEMPLATEGRAPH_GRAPHSTRUCTURE_INCLUDE_HALF_ADJACENCY_MATRIX_HPP

#include <vector>

//#include "../../LazyPrints/LazyPrinters.hpp"

namespace glygraph
{
  // TODO:Remove forward declare, why throwing errors?
  template<class T>
  class Node;

  template<class T>
  class HalfAdjacencyMatrix
  {
  public:
    /************************************************
     *  CONSTRUCTORS/DESTRUCTORS
     ***********************************************/
    HalfAdjacencyMatrix();
    HalfAdjacencyMatrix(std::vector<Node<T> *> const &nodeList_t);
    // copy constructor
    HalfAdjacencyMatrix(const HalfAdjacencyMatrix<T> &rhs);

    ~HalfAdjacencyMatrix() {};

    /************************************************
     *  GETTER/SETTER
     ***********************************************/
    unsigned int getNumEdges();
    unsigned int getNumNodes();

    /************************************************
     *  MUTATORS
     ***********************************************/
    void connect(unsigned int aNodeIndex_t, unsigned int bNodeIndex_t);
    void disconnect(unsigned int aNodeIndex_t, unsigned int bNodeIndex_t);

    /************************************************
     *  FUNCTIONS
     ***********************************************/
    bool isConnected(unsigned int aNodeIndex_t, unsigned int bNodeIndex_t);

    // TODO: Actually fix our constructors so we do not have to use this initialize workaround stuff
    void initializeWorkaround(const HalfAdjacencyMatrix<T> &rhs);
    void initializeWorkaround(std::vector<Node<T> *> const &nodeList_t);
    void emptyInitializeWorkaround(const HalfAdjacencyMatrix<T> &rhs);

    /************************************************
     *  OPERATOR OVERLOADS
     ***********************************************/
    inline bool operator()(unsigned int nodeAIndex_t, unsigned int nodeBIndex_t)
    {
      return isConnected(nodeAIndex_t, nodeBIndex_t);
    }

    // our XOR overload
    inline HalfAdjacencyMatrix<T> operator^(const HalfAdjacencyMatrix<T> &rhs) const
    {
      // TODO: Figure this out, issue is not passing in our
      //			type. Tries to pass in the int.
      // HalfAdjacencyMatrix<T> result(this->numNodes);
      HalfAdjacencyMatrix<T> result(rhs);
      result.emptyInitializeWorkaround(rhs);
      if (this->numNodes_m != rhs.numNodes_m)
        {
          //badBehavior(__LINE__, __func__, "WARNING NUMBER OF NODES NOT EQUAL XOR");
          return result;
        }
      else
        {
          for (unsigned int bitListIndex = 0; bitListIndex < this->bitList_m.size(); bitListIndex++)
            {
              if ((this->bitList_m[bitListIndex] || rhs.bitList_m[bitListIndex]) &&
                  (this->bitList_m[bitListIndex] != rhs.bitList_m[bitListIndex]))
                {
                  result.bitList_m[bitListIndex] = 1;
                  ++result.numEdges_m;
                }
            }
          return result;
        }
    }

    // our XOR-EQUALS overload
    inline HalfAdjacencyMatrix<T> &operator^=(const HalfAdjacencyMatrix<T> &rhs)
    {
      if (this->numNodes_m != rhs.numNodes_m)
        {
          //badBehavior(__LINE__, __func__, "WARNING NUMBER OF NODES NOT EQUAL XOREQUALS");
        }
      numEdges_m = 0;
      for (unsigned int bitListIndex = 0; bitListIndex < this->bitList_m.size(); bitListIndex++)
        {
          if ((this->bitList_m[bitListIndex] || rhs.bitList_m[bitListIndex]) &&
              (this->bitList_m[bitListIndex] != rhs.bitList_m[bitListIndex]))
            {
              this->bitList_m[bitListIndex] = 1;
              numEdges_m++;
            }
          else
            this->bitList_m[bitListIndex] = 0;
        }
      return *this;
    }

    // our equals operator, well assignment equals
    inline HalfAdjacencyMatrix<T> &operator=(const HalfAdjacencyMatrix<T> &rhs)
    {
      if (this->numNodes_m != rhs.numNodes_m)
        {
          //badBehavior(__LINE__, __func__, "WARNING NUMBER OF NODES NOT EQUAL EQUAL");
        }
      this->bitList_m  = rhs.bitList_m;
      this->numEdges_m = rhs.numEdges_m;
      return *this;
    }

  private:
    // Allows for proper index lookup. We use a half adjacency matrix because
    // 	our data is mirrored on the diagonal and we would like to save space
    unsigned int index(const unsigned int aNodeIndex_t, const unsigned int bNodeIndex_t);
    // Our actual connections/Adj matrix
    std::vector<bool> bitList_m;
    unsigned int      numEdges_m;
    unsigned int      numNodes_m;
    long long         indexFactor_m;
  };

  template<class T>
  HalfAdjacencyMatrix<T>::HalfAdjacencyMatrix()
  {
    //badBehavior(__LINE__, __func__, "Warning default constructor called");
    this->numEdges_m    = 0;
    this->numNodes_m    = 0;
    this->indexFactor_m = 0;
  }

  template<class T>
  HalfAdjacencyMatrix<T>::HalfAdjacencyMatrix(std::vector<Node<T> *> const &nodeList_t)
  {
    unsigned int numNodes = nodeList_t.size();
    this->bitList_m.assign(((numNodes * (numNodes - 1)) / 2), 0);
    this->numNodes_m    = numNodes;
    this->numEdges_m    = 0;
    this->indexFactor_m = (1 + 2 * (numNodes - 2));
  }

  template<class T>
  HalfAdjacencyMatrix<T>::HalfAdjacencyMatrix(const HalfAdjacencyMatrix<T> &rhs)
  {
    this->bitList_m     = rhs.bitList_m;
    this->numEdges_m    = rhs.numEdges_m;
    this->numNodes_m    = rhs.numNodes_m;
    this->indexFactor_m = rhs.indexFactor_m;
  }

  template<class T>
  void HalfAdjacencyMatrix<T>::initializeWorkaround(const HalfAdjacencyMatrix<T> &rhs)
  {
    this->bitList_m     = rhs.bitList_m;
    this->numEdges_m    = rhs.numEdges_m;
    this->numNodes_m    = rhs.numNodes_m;
    this->indexFactor_m = rhs.indexFactor_m;
  }

  template<class T>
  void HalfAdjacencyMatrix<T>::initializeWorkaround(std::vector<Node<T> *> const &nodeList_t)
  {
    unsigned int numNodes = nodeList_t.size();
    this->bitList_m.assign(((numNodes * (numNodes - 1)) / 2), 0);
    this->numNodes_m    = numNodes;
    this->numEdges_m    = 0;
    this->indexFactor_m = (1 + 2 * (numNodes - 2));
  }

  template<class T>
  unsigned int HalfAdjacencyMatrix<T>::getNumEdges()
  {
    return this->numEdges_m;
  }

  template<class T>
  void HalfAdjacencyMatrix<T>::connect(unsigned int aNodeIndex_t, unsigned int bNodeIndex_t)
  {
    if (this->bitList_m[index(aNodeIndex_t, bNodeIndex_t)])
      {
        //badBehavior(__LINE__, __func__, "TRYING TO ADD A CONNECTION THAT WAS ALREADY THERE");
      }
    else
      {
        this->bitList_m[index(aNodeIndex_t, bNodeIndex_t)] = true;
        this->numEdges_m++;
      }
  }

  template<class T>
  void HalfAdjacencyMatrix<T>::disconnect(unsigned int aNodeIndex_t, unsigned int bNodeIndex_t)
  {
    if (!(this->bitList_m[index(aNodeIndex_t, bNodeIndex_t)]))
      {
        //badBehavior(__LINE__, __func__, "TRYING TO REMOVE A CONNECTION THAT WAS NOT THERE");
      }
    else
      {
        this->bitList_m[index(aNodeIndex_t, bNodeIndex_t)] = false;
        this->numEdges_m--;
      }
  }

  template<class T>
  bool HalfAdjacencyMatrix<T>::isConnected(unsigned int aNodeIndex_t, unsigned int bNodeIndex_t)
  {
    if (aNodeIndex_t == bNodeIndex_t)
      return false;
    return this->bitList_m[index(aNodeIndex_t, bNodeIndex_t)];
  }

  template<class T>
  unsigned int HalfAdjacencyMatrix<T>::getNumNodes()
  {
    return this->numNodes_m;
  }

  template<class T>
  void HalfAdjacencyMatrix<T>::emptyInitializeWorkaround(const HalfAdjacencyMatrix<T> &rhs)
  {
    unsigned int numNodes = rhs.numNodes_m;
    this->bitList_m.assign(((numNodes * (numNodes - 1)) / 2), 0);
    this->numNodes_m    = numNodes;
    this->numEdges_m    = 0;
    this->indexFactor_m = (1 + 2 * (numNodes - 2));
  }

  template<class T>
  unsigned int HalfAdjacencyMatrix<T>::index(const unsigned int aNodeIndex_t, const unsigned int bNodeIndex_t)
  {
    if ((aNodeIndex_t < this->numNodes_m) && (bNodeIndex_t < this->numNodes_m) && (aNodeIndex_t != bNodeIndex_t))
      {
        long long aLongIndex = aNodeIndex_t;
        long long bLongIndex = bNodeIndex_t;
        if (aNodeIndex_t < bNodeIndex_t)
          {
            return (bLongIndex - aLongIndex * (aLongIndex - indexFactor_m) / 2) - 1;
          }
        else
          {
            return (aLongIndex - bLongIndex * (bLongIndex - indexFactor_m) / 2) - 1;
          }
      }
    else
      {
        //badBehavior(__LINE__, __func__, "ERROR DUE TO OUR INDEX PICKS");
        return -1;
      }
  }

} // namespace temp_graph
#endif // end TEMPLATEGRAPH_GRAPHSTRUCTURE_INCLUDE_HALF_ADJACENCY_MATRIX_HPP
