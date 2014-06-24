#ifndef EDGE_HPP
#define EDGE_HPP

#include <string>
#include <iostream>

namespace Geometry
{
    template <class T>
    class Graph;
    class EdgeAttribute;
    template <class T>
    class Edge
    {
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            Edge();

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the dst
              * @return dst_ attribute of the current object of this class
              */
            Graph<T>* GetDst();
            /*! \fn
              * An accessor function in order to access to the properties
              * @return properties_ attribute of the current object of this class
              */
            EdgeAttribute* GetProperties();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the dst of the current object
              * Set the dst_ attribute of the current graph
              * @param dst The dst attribute of the current object
              */
            void SetDst(Graph<T>* dst);
            /*! \fn
              * A mutator function in order to set the properties of the current object
              * Set the properties_ attribute of the current graph
              * @param properties The properties attribute of the current object
              */
            void SetProperties(EdgeAttribute* properties);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the edge contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            Graph<T>* dst_;
            EdgeAttribute* properties_;
    };
}

#endif // EDGE_HPP
