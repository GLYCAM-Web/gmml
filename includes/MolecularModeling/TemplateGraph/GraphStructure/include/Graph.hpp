#ifndef TEMPLATEGRAPH_GRAPHSTRUCTURE_INCLUDE_GRAPH_HPP
#define TEMPLATEGRAPH_GRAPHSTRUCTURE_INCLUDE_GRAPH_HPP

//#include "../../LazyPrints/LazyPrinters.hpp"
#include "HalfAdjacencyMatrix.hpp"
#include "Node.hpp"

#include <algorithm>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>

namespace glygraph
{

  template<class T>
  class Graph
  {
  public:
    /************************************************
     *  CONSTRUCTORS/DESTRUCTORS
     ***********************************************/
    Graph();
    // TODO: Ensure we would like this functionality, current idea is pass root node then get all traversable from this
    // node and store in our set
    Graph(Node<T> *const &initialNode_t);
    Graph(std::vector<Node<T> *> const &nodeList_t);

    // copy constructor
    Graph(const Graph<T> &rhs);

    // move constructor
    Graph(Graph<T> &&rhs);

    // copy assignment
    Graph<T> &operator=(const Graph<T> &rhs);

    // move assignment
    Graph<T> &operator=(Graph<T> &&rhs);

    virtual ~Graph();

    /************************************************
     *  GETTER/SETTER
     ***********************************************/
    // TODO: Finalize what we would like to pass. Using weak_ptr is nice because we can easily check
    // 			if our node is even useful/alive still but this would get annoying to constantly pass
    // 			weak_ptr<Node<T>>. As of now, for our getNodes I will be using a raw ptr due to the fact
    // 			that there should be no deletions etc. when we run our algos. I will do a check prior to
    // 			returning the vector to ensure that all nodes in our node-list are still alive.
    //
    std::vector<Node<T> *> getNodes();
    HalfAdjacencyMatrix<T> getAdjMatrix() const;

    unsigned int getIndexFromNode(Node<T> *const &queryNode_t);
    Node<T> *    getNodeFromIndex(unsigned int const &queryIndex_t);
    /************************************************
     *  MUTATORS
     ***********************************************/

    /************************************************
     *  FUNCTIONS
     ***********************************************/
    std::string getGraphvizLink();
    //rewalk our graph and update our allNodes_m to what the current walk is. Have both dfs and bfs
    // available. Please note that startNode_t must be within allNodes_m before we start a walk
    void dfsWalk(Node<T> *const&startNode_t);
    void bfsWalk(Node<T> *const&startNode_t);

  private:
    /************************************************
     *  ATTRIBUTES
     ***********************************************/
    HalfAdjacencyMatrix<T> adjMatrix_m;
    std::vector<Node<T> *> allNodes_m;
    // TODO: Ensure the correct hashing function is being used. Must be 100% sure, am only somewhat sure.
    std::unordered_map<unsigned int, Node<T> *> nodeLookup_m;
    std::unordered_map<Node<T> *, unsigned int> indexLookup_m;

    /************************************************
     *  HELPER FUNCTIONS
     ***********************************************/
    void populateAdjacencyMatrix();
    void populateLookups();
    void lazyExpiredFixer();

	//Used to assist in recursion when running bfs or dfs
    void dfsHelper(Node<T> *const &currentNode_t, std::unordered_set<Node<T> *> &visitedNodeSet_t,
            				std::vector<Node<T> *> &reachableNodes_t);
    void bfsHelper(Node<T> *const &currentNode_t, std::unordered_set<Node<T> *> &visitedNodeSet_t,
    						std::vector<Node<T> *> &reachableNodes_t);

    std::vector<Node<T> *> getReachableNodes(Node<T> *const &startingNode_t);
    // NOTE: To be used when we are passed solely a root node.
    void getReachableHelper(Node<T> *const &currentNode_t, std::unordered_set<Node<T> *> &visistedNodeSet_t,
                            std::vector<Node<T> *> &reachableNodes_t);
  }; // end graph class

  template<class T>
  inline Graph<T>::Graph()
  {
    //badBehavior(__LINE__, __func__, "Warning calling default graph constructor");
  }

  template<class T>
  inline Graph<T>::Graph(Node<T> *const &initialNode_t)
  {
	this->allNodes_m.push_back(initialNode_t);
	this->dfsWalk(initialNode_t);

	// populate our lookups
    this->populateLookups();
    this->populateAdjacencyMatrix();
  }

  template<class T>
  inline Graph<T>::Graph(std::vector<Node<T> *> const &nodeList_t)
  {
    if (nodeList_t.size() > 0)
      {
        // Lazy way to prevent dupes, again need to come up with
        // a more efficient way to actually prevent our dupes
        std::unordered_set<Node<T> *> tempNodeSet(nodeList_t.begin(), nodeList_t.end());
        for (Node<T> *currNode : tempNodeSet)
          {
            this->allNodes_m.push_back(currNode);
          }

        this->populateLookups();
        this->populateAdjacencyMatrix();
      }
    else
      {
        //badBehavior(__LINE__, __func__, "Was passed a nodelist of size 0");
      }
  }

  template<class T>
  inline Graph<T>::~Graph()
  {
    //lazyInfo(__LINE__, __func__, "Graph deleted");
  }

  template<class T>
  inline std::vector<Node<T> *> Graph<T>::getNodes()
  {
    return this->allNodes_m;
  }

  template<class T>
  inline HalfAdjacencyMatrix<T> Graph<T>::getAdjMatrix() const
  {
    return this->adjMatrix_m;
  }

  template<class T>
  inline void Graph<T>::populateAdjacencyMatrix()
  {
    if ((this->allNodes_m.size() > 0) && (this->indexLookup_m.size()))
      {
        this->adjMatrix_m.initializeWorkaround(this->getNodes());
        for (Node<T> *currNode : this->allNodes_m)
          {
            for (Node<T> *currNodeNeighbor : currNode->getNeighbors())
              {
                if (!(this->adjMatrix_m.isConnected(this->indexLookup_m[currNode],
                                                    this->indexLookup_m[currNodeNeighbor])))
                  {
                    this->adjMatrix_m.connect(this->indexLookup_m[currNode], this->indexLookup_m[currNodeNeighbor]);
                  }
              }
          }
      }
    else
      {
        //badBehavior(__LINE__, __func__, "no nodes present! cannot populate adj matrix");
      }
  }

  // NOTE: Must call before we try to create our adj matrix
  template<class T>
  inline void Graph<T>::populateLookups()
  {
    if (this->allNodes_m.size() > 0)
      {
        this->nodeLookup_m.clear();
        this->indexLookup_m.clear();

        int currIndex = 0;
        for (Node<T> *currNode : this->allNodes_m)
          {

            this->nodeLookup_m.insert({currIndex, currNode});

            this->indexLookup_m.insert({currNode, currIndex});

            currIndex++;
          }
      }
    else
      {
        //badBehavior(__LINE__, __func__, "Warning tried to populate lookups with no node list");
      }
  }

  template<class T>
  inline std::vector<Node<T> *> Graph<T>::getReachableNodes(Node<T> *const &startingNode_t)
  {
    std::unordered_set<Node<T> *> visitedNodes;
    // TODO: Please note that this current method does increase the size of our call stack a good bit due to the use of
    // recursion. 			Depending on how large of graphs we are dealing with this could become an issue and it may be
    // a better 			call to use a different method.

    std::vector<Node<T> *> reachableVecToReturn;

    this->getReachableHelper(startingNode_t, visitedNodes, reachableVecToReturn);
    return reachableVecToReturn;
  }

  // Should be correct. Passing pointer by reference
  template<class T>
  inline unsigned int Graph<T>::getIndexFromNode(Node<T> *const &queryNode_t)
  {
    return this->indexLookup_m[queryNode_t];
  }

  template<class T>
  inline Node<T> *Graph<T>::getNodeFromIndex(unsigned int const &queryIndex_t)
  {
    return this->nodeLookup_m[queryIndex_t];
  }

  // TODO: Find a better way to remove all expired ptrs that will always work
  // 			I am worried that just iterating through index will not work well
  // 			in all cases. Could just do currIndex-- once we find an expired
  // 			and remove it at the end of our current loop iteration tho.
  //
  template<class T>
  inline void Graph<T>::lazyExpiredFixer()
  {
    // Possibly a good way, need to run through some tests.
    unsigned int ogSize = this->allNodes_m.size();
    for (unsigned int currIndex = 0; currIndex < this->allNodes_m.size(); currIndex++)
      {
        if (this->allNodes_m[currIndex].expired())
          {
            this->allNodes_m.erase(this->allNodes_m.begin() + currIndex);
            currIndex--;
          }
      }

    //	Extremely heavy way but def works
    //
    // unsigned int ogSize = this->allNodes.size();
    // std::vector<std::weak_ptr<Node<T>>> dustyList = this->allNodes;
    // this->allNodes.clear();
    // for (std::weak_ptr<Node<T>> currDusty : dustyList)
    // {
    // if (!(currDusty.expired()))
    // {
    // this->m_allNodes.push_back(currDusty);
    // }
    // }
    if (ogSize != this->allNodes_m.size())
      {
        this->populateLookups();
        this->populateAdjacencyMatrix();
      }
  }

  template<class T>
  inline std::string Graph<T>::getGraphvizLink()
  {
    std::string connectionArrow = "%20-%3E%20";
    std::string newLine         = "%0A%09";
    std::string baseURL         = "https://dreampuf.github.io/GraphvizOnline/#digraph%20G%20%7B%0A%09";
    std::string endBracket      = "%0A%7D";

    // first make a collection of all of our node connections, just use set cause less chars lazy
    std::map<Node<T> *, std::set<Node<T> *>> nodeNeighs;
    for (Node<T> *currNode : this->getNodes())
      {
        for (Node<T> *currNeigh : currNode->getNeighbors())
          {
            // ensure that the current connection is not already present
            //	We know we do not have a specific connection if either
            // 		A) We do NOT have our neighbor as a key in the node neighs
            // 						OR
            // 		B) If we DO have our neighbor as a key in the node neighs, then we do NOT have
            // 				the currNode as a member
            if ((nodeNeighs.count(currNeigh) == 0) || (nodeNeighs[currNeigh].count(currNode) == 0))
              {
                if (nodeNeighs[currNode].count(currNeigh) == 0)
                  {
                    nodeNeighs[currNode].insert(currNeigh);
                  }
                else
                  {
                    //badBehavior(__LINE__, __func__, "ARGH MATEY");
                  }
              }
          } // end 4
      }     // end 44

    std::string endURL;

    for (std::pair<Node<T> *, std::set<Node<T> *>> currPair : nodeNeighs)
      {
        for (Node<T> *currNeigh : currPair.second)
          {
            endURL += currPair.first->getName() + "->" + currNeigh->getName() + newLine;
          }
      }
    endURL += endBracket;
    baseURL += endURL;
    return baseURL;
  }
  /*
   adjMatrix;
   std::vector<Node<T>*> allNodes;
   // TODO: Ensure the correct hashing function is being used. Must be 100% sure, am only somewhat sure.
   std::unordered_map<unsigned int, Node<T>*> nodeLookup;
   std::unordered_map<Node<T>*, unsigned int> m_indexLookup;
   */

  // copy constructor
  template<class T>
  inline Graph<T>::Graph(const Graph<T> &rhs)
      : adjMatrix_m(rhs.adjMatrix_m), nodeLookup_m(rhs.nodeLookup_m), indexLookup_m(rhs.indexLookup_m),
        allNodes_m(rhs.allNodes_m)
  {
  }

  // move constructor
  template<class T>
  inline Graph<T>::Graph(Graph<T> &&rhs)
      : adjMatrix_m(rhs.adjMatrix_m), nodeLookup_m(rhs.nodeLookup_m), indexLookup_m(rhs.indexLookup_m),
        allNodes_m(rhs.allNodes_m)
  {
  }

  // copy assignment
  template<class T>
  inline Graph<T> &Graph<T>::operator=(const Graph<T> &rhs)
  {
    this->adjMatrix_m   = rhs.adjMatrix_m;
    this->allNodes_m    = rhs.allNodes_m;
    this->indexLookup_m = rhs.indexLookup_m;
    this->nodeLookup_m  = rhs.nodeLookup_m;
    return *this;
  }

  // move assignment
  template<class T>
  inline Graph<T> &Graph<T>::operator=(Graph<T> &&rhs)
  {
    this->adjMatrix_m   = rhs.adjMatrix_m;
    this->allNodes_m    = rhs.allNodes_m;
    this->indexLookup_m = rhs.indexLookup_m;
    this->nodeLookup_m  = rhs.nodeLookup_m;
    return *this;
  }

template<class T>
inline void Graph<T>::dfsWalk(Node<T> *const&startNode_t)
{
	if (std::find(this->allNodes_m.begin(), this->allNodes_m.end(), startNode_t) != this->allNodes_m.end())
	{
		std::unordered_set<Node<T> *> visitedNodes;
		// TODO: Please note that this current method does increase the size of our call stack a good bit due to the use of
		// recursion. 			Depending on how large of graphs we are dealing with this could become an issue and it may be
		// a better 			call to use a different method.

		std::vector<Node<T> *> reachableVecToReturn;

		//NOTE: reachable vector to return redundant due to visited nodes
		this->dfsHelper(startNode_t, visitedNodes, reachableVecToReturn);
		this->allNodes_m.clear();
		this->allNodes_m = reachableVecToReturn;
	}
}

template<class T>
inline void Graph<T>::dfsHelper(Node<T> *const&currentNode_t,
		std::unordered_set<Node<T>*> &visitedNodeSet_t,
		std::vector<Node<T>*> &reachableNodes_t)
{
    reachableNodes_t.push_back(currentNode_t);

    visitedNodeSet_t.insert(currentNode_t);
    for (Node<T> *currNeighbor : currentNode_t->getNeighbors())
      {
        if (!(visitedNodeSet_t.count(currNeighbor)))
          {
            this->dfsHelper(currNeighbor, visitedNodeSet_t, reachableNodes_t);
          }
      }
}

template<class T>
inline void Graph<T>::bfsWalk(Node<T> *const&startNode_t)
{
}

template<class T>
inline void Graph<T>::bfsHelper(Node<T> *const&currentNode_t,
		std::unordered_set<Node<T>*> &visitedNodeSet_t,
		std::vector<Node<T>*> &reachableNodes_t)
{
}

  template<class T>
  inline void Graph<T>::getReachableHelper(Node<T> *const &currentNode_t, std::unordered_set<Node<T> *> &visitedNodeSet_t,
                                           std::vector<Node<T> *> &reachableNodes_t)
  {

    reachableNodes_t.push_back(currentNode_t);

    visitedNodeSet_t.insert(currentNode_t);
    for (Node<T> *currNeighbor : currentNode_t->getNeighbors())
      {
        if (!(visitedNodeSet_t.count(currNeighbor)))
          {
            this->getReachableHelper(currNeighbor, visitedNodeSet_t, reachableNodes_t);
          }
      }
  }

} // namespace temp_graph
#endif // TEMPLATEGRAPH_GRAPHSTRUCTURE_INCLUDE_GRAPH_HPP
