#ifndef TEMPLATEGRAPH_ALGORITHMS_INCLUDE_CONNECTIVITYIDENTIFIER_HPP
#define TEMPLATEGRAPH_ALGORITHMS_INCLUDE_CONNECTIVITYIDENTIFIER_HPP

#include <memory>
#include <vector>

#include "../../GraphStructure/include/Graph.hpp"
#include "../../GraphStructure/include/Node.hpp"

//#include "../../LazyPrints/LazyPrinters.hpp"

namespace
{

  //	Please note that this algo only labels our bridge EDGES, not our bridge nodes. We will
  //		do a little extra lifting to label all of our bridge nodes in our main function and
  // 		not this one that is hidden away in an anonymous namespace.

  template<class T>
  void bridgeDetectHelperDFS(glygraph::Node<T> *currNode_t, glygraph::Node<T> *nextNode_t,
                             glygraph::Graph<T> &runningGraph_t, std::vector<int> &preTime_t,
                             std::vector<int> &lowestTime_t, unsigned int &counter_t)
  {
    unsigned int nextIndex = runningGraph_t.getIndexFromNode(nextNode_t);

    preTime_t[nextIndex]    = counter_t++;
    lowestTime_t[nextIndex] = preTime_t[nextIndex];

    // now we run our dfs recursion on our neighbors.
    for (glygraph::Node<T> *currNeigh : nextNode_t->getNeighbors())
      {
        unsigned int currNeighIndex = runningGraph_t.getIndexFromNode(currNeigh);

        //	We know we have not run our algo on said node yet so we need to
        // 		go ahead and use our friend recurion to run through the nodes.
        if (preTime_t[currNeighIndex] == -1)
          {
            bridgeDetectHelperDFS(nextNode_t, currNeigh, runningGraph_t, preTime_t, lowestTime_t, counter_t);
            lowestTime_t[nextIndex] = std::min(lowestTime_t[nextIndex], lowestTime_t[currNeighIndex]);
            if (lowestTime_t[currNeighIndex] == preTime_t[currNeighIndex])
              {
                //	Now we know all edges between our two nodes are bridge edges so
                //		we must label them as such. Since we are not dealing with a multigraph
                //  	we know we only have 1 edge to label for now.
                glygraph::Edge<T> *connectingEdge = currNeigh->getConnectingEdge(nextNode_t);
                connectingEdge->setConnectivityTypeIdentifier(glygraph::ConnectivityType::BRIDGE);
              } // end if where we lable our bridge edges

          } // end if
        else if (currNeighIndex != runningGraph_t.getIndexFromNode(currNode_t))
          {
            lowestTime_t[nextIndex] = std::min(lowestTime_t[nextIndex], preTime_t[currNeighIndex]);
          }
      } // end bridgeDetectHelperDFS
  }
} // namespace

namespace id_connectivity
{

  // Now this is a little more odd. We could just remove an edge, run dfs/bfs and
  // 	make sure all our nodes are still connected but that gives us a time complexity
  // 	of O(E*(N+E)) which is not very good especially considering that we want to be fast
  // 	and stay zoomin. We can do something better which is based off the idea of articulation
  // 	points.
  //
  // 	A node is an articulation point if and only if removing it and the connected edges disconnects
  // 	our graph. Please note that this is slightly different than the definition of a bridge due
  // 	to the fact that this is a little bit "looser" of a definition in comparison to a bridge.
  // 	We run DFS and in our call tree an edge is considered a bridge if there are no other
  // 	alternatives to reach our "parent" node (node we are looking at) or an ancestor of said parent
  // 	node from our "child" node (next node we are looking at).
  template<class T>
  void bridgeDetect(glygraph::Graph<T> &graphToBridgeDetect_t)
  {
    unsigned int numOfNodes = graphToBridgeDetect_t.getNodes().size();

    // our lowest time for a certain node to be visited
    std::vector<int> lowestTime(numOfNodes, -1);
    // keep track of the time to visit our preorder
    std::vector<int> preTime(numOfNodes, -1);

    unsigned int counter = 0;

    for (glygraph::Node<T> *currNode : graphToBridgeDetect_t.getNodes())
      {
        // If our "time" for our previous time is -1 we know that we have not checked said
        // 	node yet thus we must run our algo.
        if (preTime[graphToBridgeDetect_t.getIndexFromNode(currNode)] == -1)
          {
            bridgeDetectHelperDFS(currNode, currNode, graphToBridgeDetect_t, preTime, lowestTime, counter);
          }
      } // end labeling all of our bridge EDGES

    // Now we hit all of our nodes, we only label a node as a bridge if all of its
    // 	edges are bridge edges.

    for (glygraph::Node<T> *currNode : graphToBridgeDetect_t.getNodes())
      {
        if (!(currNode->getConnectivityTypeIdentifier() == glygraph::ConnectivityType::LEAF))
          {
            bool possibleBridgeNode = false;

            for (glygraph::Edge<T> *currEdge : currNode->getEdges())
              {
                // If we ever hit an edge that is NOT a bridge edge we know we dont want to insert
                if (!(currEdge->getConnectivityTypeIdentifier() == glygraph::ConnectivityType::BRIDGE))
                  {
                    possibleBridgeNode = false;
                    break;
                  }
                else
                  {
                    possibleBridgeNode = true;
                  }
              }
            if (possibleBridgeNode)
              {
                currNode->setConnectivityTypeIdentifier(glygraph::ConnectivityType::BRIDGE);
              }
          }
      }

  } // end our bridge detect algo

  // This is easy to do, we know that a node is a leaf if it has 1 neighbor. Thus we label our node
  // 	that is being checked as a leaf and the edge associated with it as a bridge edge.
  template<class T>
  void leafDetect(glygraph::Graph<T> &graphToLeafDetect_t)
  {
    for (glygraph::Node<T> *currNode : graphToLeafDetect_t.getNodes())
      {
        // we are at leaf
        if (currNode->getNeighbors().size() == 1)
          {
            // Set node as leaf
            currNode->setConnectivityTypeIdentifier(glygraph::ConnectivityType::LEAF);
            int errorEdgeCounter = 0;
            for (glygraph::Edge<T> *currEdge : currNode->getEdges())
              {
                if (errorEdgeCounter > 1)
                  {
                    //badBehavior(__LINE__, __func__,
                    //            "WE FOUND A NODE THAT HAS MORE THAN ONE EDGE BUT ONLY HAS ONE NEIGHBOR");
                    break;
                  }
                else
                  {
                    // set edge as bridge
                    currEdge->setConnectivityTypeIdentifier(glygraph::ConnectivityType::BRIDGE);
                  }
                errorEdgeCounter++;
              } // end currEdge for
          }
      } // end all nodes for
  }     // end leaf detect

  //	The only issue is I dislike this naming, also I would like to have this only
  // 		label the nodes within a cycle instead of having to run all of our algos.
  // 		We know if a node/edge is neither a bridge or a node then it MUST be in a
  // 		cycle. I checked out a good bit of algos and even thew one together myself
  // 		but it seems like the best approach is done by using bridge detection.
  //
  // 		My self-developed algo was lazy and would just run dfs through and check
  // 		for a leading edge to the node we called then we could label said node as
  // 		incycle then also label the first edge we hit leaving our node and the one
  // 		hit to go back to the node as in cycle but this must be run on every single
  // 		node. Another way is we just label all nodes as if they are in a cycle then
  // 		change each edge between 2 nodes that are in cycle to a cycle edge. Finally
  // 		we could keep a stack of what edges & nodes we visited and then if we hit a
  // 		cycle then all in the stack are considered incycle but this is basically the
  // 		same as the bridge algo so we would be stupid to not just call the leaf &
  // 		bridge detection then just label all unknown as incycle. Let me know if we
  // 		would prefer to do something else for cycle signifier.
  template<class T>
  void cycleDetect(glygraph::Graph<T> &graphToCycleDetect_t)
  {
    leafDetect(graphToCycleDetect_t);
    bridgeDetect(graphToCycleDetect_t);
    for (glygraph::Node<T> *currNode : graphToCycleDetect_t.getNodes())
      {
        if (currNode->getConnectivityTypeIdentifier() == glygraph::ConnectivityType::UNKNOWN)
          {
            currNode->setConnectivityTypeIdentifier(glygraph::ConnectivityType::INCYCLE);
            //	now we need to hit our connecting edges to the
            // 		newly labeled node and do the same logic.
            // 		I dont like this it is slow and gross.
            for (glygraph::Edge<T> *currEdge : currNode->getEdges())
              {
                if (currEdge->getConnectivityTypeIdentifier() == glygraph::ConnectivityType::UNKNOWN)
                  {
                    currEdge->setConnectivityTypeIdentifier(glygraph::ConnectivityType::INCYCLE);
                  }
              } // end our edge labeling for loop
          }
      } // end node for loop

    int badNodeCounter = 0;
    int badEdgeCounter = 0;

    // Lazy post check. We should have all labeled as of now.
    for (glygraph::Node<T> *currNode : graphToCycleDetect_t.getNodes())
      {
        if (currNode->getConnectivityTypeIdentifier() == glygraph::ConnectivityType::UNKNOWN)
          {
            badNodeCounter++;
          }
        for (glygraph::Edge<T> *currEdge : currNode->getEdges())
          {
            if (currEdge->getConnectivityTypeIdentifier() == glygraph::ConnectivityType::UNKNOWN)
              {
                badEdgeCounter++;
              }
          }
      }

    if ((badNodeCounter > 0) || (badEdgeCounter > 0))
      {
        //badBehavior(__LINE__, __func__,
        //            "Warning after running our total lableing there are " + std::to_string(badNodeCounter) +
        //                " UNKNOWN nodes and " + std::to_string(badEdgeCounter) + " UNKNOWN edges");
      }
  }

  //	Until I finalize our algos, this is just a lazy wrapper to make
  // 		our other coders happy and to show that this does total
  // 		connectivity type labeling.
  //
  //
  template<class T>
  void identifyConnectivity(glygraph::Graph<T> &graphToConnectivityDetect_t)
  {
    cycleDetect(graphToConnectivityDetect_t);
  }
} // namespace connect_id

#endif // end TEMPLATEGRAPH_ALGORITHMS_INCLUDE_CONNECTIVITYIDENTIFIER_HPP
