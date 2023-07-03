#ifndef TEMPLATEGRAPH_ALGORITHMS_INCLUDE_SUBGRAPHMATCHING_HPP
#define TEMPLATEGRAPH_ALGORITHMS_INCLUDE_SUBGRAPHMATCHING_HPP

#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../../GraphStructure/include/Graph.hpp"
#include "../../GraphStructure/include/HalfAdjacencyMatrix.hpp"
#include "../../GraphStructure/include/Node.hpp"

// #include "../../LazyPrints/LazyPrinters.hpp"

namespace
{
    //	This is used to extract the patterns of our query graph. A pattern is defined by a key and a set of node
    // signifiers.
    // This is used to "step through" what our query graph is and see if the steps relate to our actual graph. As of
    // now, our key and node signifiers are the names of the nodes. This can be changed but please note it will require
    // us to cascade changes throughout our algo.
    template<class T>
    std::map<std::string, std::vector<std::string>> patternExtractor(glygraph::Graph<T>& patternGraph_t)
    {
        std::map<std::string, std::vector<std::string>> foundPatterns;

        for (unsigned int indexA = 0; indexA < patternGraph_t.getNodes().size(); indexA++)
        {
            // our "A-node"'s name will be our key for patterns
            std::pair<std::string, std::vector<std::string>> currentPattern;
            // we now place current nodes name as key for our current pattern
            currentPattern.first = patternGraph_t.getNodeFromIndex(indexA)->getName();
            // now we hit all other nodes and if the 2 nodes are connected we throw the
            //	secondary node signifier in our pattern
            for (unsigned int indexB = 0; indexB < patternGraph_t.getNodes().size(); indexB++)
            {
                // now we must check if the two nodes are connected. We can do this faster.
                // Fix up later
                if (patternGraph_t.getAdjMatrix().isConnected(indexA, indexB))
                {
                    currentPattern.second.push_back(patternGraph_t.getNodeFromIndex(indexB)->getName());
                }
            }
            // once we populate all our patterns to the current key, we throw the pair
            // into our total pattern set, then go to next node.
            foundPatterns.insert(currentPattern);
        }
        return foundPatterns;
    } // end found patterns

    // forward declare because dumb
    template<class T>
    void searchMatches(std::vector<glygraph::Node<T>*> matches_t,
                       std::map<std::string, std::vector<std::string>> patterns_t,
                       std::pair<std::vector<glygraph::Node<T>*>, std::vector<glygraph::Edge<T>*>>& resultsPair_t,
                       std::unordered_set<glygraph::Node<T>*>& visitedKeys_t, glygraph::Node<T>* previousNode_t,
                       glygraph::Graph<T>& graphSearch_t);

    template<class T>
    int searchForPatterns(unsigned int currNodeIndex_t, std::map<std::string, std::vector<std::string>> patterns_t,
                          std::pair<std::vector<glygraph::Node<T>*>, std::vector<glygraph::Edge<T>*>>& resultsPair_t,
                          std::unordered_set<glygraph::Node<T>*>& visitedKeys_t, glygraph::Node<T>* previousNode_t,
                          glygraph::Graph<T>& graphSearch_t)
    {
        glygraph::Node<T>* currNode = graphSearch_t.getNodeFromIndex(currNodeIndex_t);

        // we check if we have a corresponding pattern to this node
        if (patterns_t.count(currNode->getName()) && !visitedKeys_t.count(currNode))
        {
            visitedKeys_t.insert(currNode);
            // When we are searching and we hit a match in our graphSearch (i.e. a parital
            // 	match onto our query graph) we want to put them in our matched nodes, then
            // 	we want to check search our patterns for each match.
            std::vector<glygraph::Node<T>*> foundMatches;

            for (glygraph::Node<T>* interestingNode : graphSearch_t.getNodes())
            {
                // for now we just want to check if we have our current node
                // we are checking out in our results vector.
                bool isInterestingInResults = (std::find(resultsPair_t.first.begin(), resultsPair_t.first.end(),
                                                         interestingNode) != resultsPair_t.first.end());
                // if we dont have current node in our results we want to check it out
                if (!(visitedKeys_t.count(interestingNode) && !isInterestingInResults))
                {
                    // now ensure the two nodes are actually connected
                    int interestingNodeIndex = graphSearch_t.getIndexFromNode(interestingNode);

                    if (graphSearch_t.getAdjMatrix().isConnected(currNodeIndex_t, interestingNodeIndex))
                    {
                        // after confirmed connect we check if our interesting node
                        // is within the pattern requirements of our "current" node
                        // recall our keys for our pattern matching is the name of the node
                        std::vector<std::string> tempPatternReqs = patterns_t[currNode->getName()];

                        bool isInterestingInReqs = (std::find(tempPatternReqs.begin(), tempPatternReqs.end(),
                                                              interestingNode->getName()) != tempPatternReqs.end());

                        if (isInterestingInReqs)
                        {
                            // We now know our interesting node is not in our results,
                            // 		is connected to the current node we are checking out,
                            //		and that it is a pattern requirement for the current
                            // 		node we are checking out. Thus we put it on our matches
                            // 		that we are going to have to search through in this same
                            // 		manner.
                            foundMatches.push_back(interestingNode);
                        }
                    } // end the if for if teh 2 nodes are connected
                }
            }
            // end our for loop, this finds all matches we want to check and possibly throw in our
            // pairResults

            unsigned int currNodeReqsLength = patterns_t[currNode->getName()].size();
            // now we want to make sure that we have AT LEAST the same amount of matches
            // 	as we do for the current node requirements. Keep in mind we need to hit
            // 	all of each nodes pattern requirements to continue with a valid/matching
            // 	subgraph. Also we need to ensure our matches size is larger than 0 because
            // 	if we have 0 matches and 0 req length we are at a leaf
            if ((foundMatches.size() >= currNodeReqsLength) && (foundMatches.size() > 0))
            {
                // we now prune our matches in order to make the rest of our resursive journey
                // easier and faster. Once we match a pattern who cares about it anymore for this current run
                patterns_t.erase(currNode->getName());

                //	Now we get the edge between the previous node and this node to insert.
                //
                //  NOTE: For some odd reason when I take out this size check I end up
                //			trying to get a connecting edge between 2 nodes (many times) when
                //			they do NOT have an edge between them. Was due to never updating the
                //			"firstNode" value which would stay the same value as the node returned
                //			by checking key "0" from our graph.
                // if (resultsPair.first.size() > 0)
                //{
                // Since I have a small brain and need to check whats getting algo mad I wanna
                //	see if my idea is correct. It was, we were trying to get an edge between
                //	our node and itself (a "loop" edge).
                // if (!(currNode == resultsPair.first.back()))
                if (!(previousNode_t == currNode))
                {
                    resultsPair_t.second.push_back(currNode->getConnectingEdge(previousNode_t));
                }
                //}				//end our bit that inserts edge.

                // add the current node we are checking out to our results since we are good so far
                resultsPair_t.first.push_back(currNode);
                previousNode_t = currNode;
                searchMatches(foundMatches, patterns_t, resultsPair_t, visitedKeys_t, previousNode_t, graphSearch_t);

                // return 2 to designate we are continuing our traversal
                return 2;
            }
            else if ((foundMatches.size() == 0) && (currNodeReqsLength == 0))
            {
                // we have hit a leaf, update our patterns
                patterns_t.erase(currNode->getName());

                if (resultsPair_t.first.size() > 0)
                {
                    // Since I have a small brain and need to check whats getting algo mad I wanna
                    //	see if my idea is correct. It was, we were trying to get an edge between
                    //	our node and itself (a "loop" edge).
                    if (!(currNode == resultsPair_t.first.back()))
                    {
                        resultsPair_t.second.push_back(currNode->getConnectingEdge(resultsPair_t.first.back()));
                    }
                } // end our bit that inserts edge.

                resultsPair_t.first.push_back(currNode);

                // return 1 for our leaf case
                return 1;
            }
        }
        // return 0 since we dont hit a match
        return 0;
    } // end search patterns

    // This is used to search all the matches we previously found. Works by calling our
    // 	search patterns function for each matched node, then we find THOSE matches
    // 	then hit this function again.
    //
    template<class T>
    void searchMatches(std::vector<glygraph::Node<T>*> matches_t,
                       std::map<std::string, std::vector<std::string>> patterns_t,
                       std::pair<std::vector<glygraph::Node<T>*>, std::vector<glygraph::Edge<T>*>>& resultsPair_t,
                       std::unordered_set<glygraph::Node<T>*>& visitedKeys_t, glygraph::Node<T>* previousNode_t,
                       glygraph::Graph<T>& graphSearch_t)
    {
        for (glygraph::Node<T>* currMatch : matches_t)
        {
            // as before we need to check if our node we are checking out is
            // already in our results i.e. already been run through
            bool isMatchInResults = (std::find(resultsPair_t.first.begin(), resultsPair_t.first.end(), currMatch) !=
                                     resultsPair_t.first.end());
            // We want to make sure we havent used the current node as an "entry point"
            // 		for our search pattern function AND we want to make sure she doesnt
            // 		exist in our results yet.
            //
            if (!isMatchInResults && !(visitedKeys_t.count(currMatch)))
            {
                int searchState = searchForPatterns(graphSearch_t.getIndexFromNode(currMatch), patterns_t,
                                                    resultsPair_t, visitedKeys_t, previousNode_t, graphSearch_t);
                // if we hit a leaf we get out
                if (searchState == 1)
                {
                    break;
                }
            }
        }
    } // end search matches

} // namespace

// end anon namespace

//	As with the issue brought up within the cycle decomposition algo, the same issue of
//		us accidently get an induced structure arises here. Granted it may be more
//		difficult to accidently get an induced subgraph but it is totally possible.
//		This would be extremely annoying because we would have to completely hunt down
//		that it is this algo itself and the way we traverse the data is causing us to
//		return bad data.
//
//	TODO: Check if this works for cycles and we correctly return the needed data,
//			for instance how the cycle decomp returns all cycles and all edges. I am worried
//			this will end up missing an edge.
namespace subgraph_match
{
    template<class T>
    std::unordered_map<glygraph::Node<T>*, std::pair<std::vector<glygraph::Node<T>*>, std::vector<glygraph::Edge<T>*>>>
    findSubgraphs(glygraph::Graph<T>& mainGraph_t, glygraph::Graph<T>& queryGraph_t)
    {
        // This will be what is returned, we be slowly built using what is right below it.
        std::unordered_map<glygraph::Node<T>*,
                           std::pair<std::vector<glygraph::Node<T>*>, std::vector<glygraph::Edge<T>*>>>
            subgraphEdgeNodeResults;

        std::pair<std::vector<glygraph::Node<T>*>, std::vector<glygraph::Edge<T>*>> pairedResult;

        // first we want to grab all our patterns that we will use to match on our main graph
        std::map<std::string, std::vector<std::string>> patternsToMatch = patternExtractor(queryGraph_t);

        std::vector<glygraph::Node<T>*> currResults;

        std::unordered_set<glygraph::Node<T>*> keyVisitTracker;

        std::unordered_map<glygraph::Node<T>*, std::vector<glygraph::Node<T>*>> finalResults;

        // now we run a search for each node
        for (unsigned int searchStartNodeIndex = 0; searchStartNodeIndex < mainGraph_t.getNodes().size();
             searchStartNodeIndex++)
        {
            glygraph::Node<T>* firstNode = mainGraph_t.getNodeFromIndex(searchStartNodeIndex);
            searchForPatterns(searchStartNodeIndex, patternsToMatch, pairedResult, keyVisitTracker, firstNode,
                              mainGraph_t);

            if (pairedResult.first.size() != 0)
            {
                subgraphEdgeNodeResults.insert({mainGraph_t.getNodeFromIndex(searchStartNodeIndex), pairedResult});
                // after we recursively hit all patterns we want to go ahead and clear out
                //	our results so we hit all again. From my understanding this is, in our case,
                //	the best way to clear our pairedResult pair.
                pairedResult = {};
            }
            keyVisitTracker.clear();
        }
        return subgraphEdgeNodeResults;
    }

} // namespace subgraph_match

#endif // END TEMPLATEGRAPH_ALGORITHMS_INCLUDE_SUBGRAPHMATCHING_HPP
