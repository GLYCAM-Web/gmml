#ifndef EDGEATTRIBUTE_HPP
#define EDGEATTRIBUTE_HPP

#include <string>
#include <iostream>

#include "../../common.hpp"

namespace Geometry
{
    class EdgeAttribute
    {
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            EdgeAttribute();

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the force
              * @return force_ attribute of the current object of this class
              */
            double GetForce();
            /*! \fn
              * An accessor function in order to access to the weight
              * @return weight_ attribute of the current object of this class
              */
            double GetWeight();
            /*! \fn
              * An accessor function in order to access to the length
              * @return length_ attribute of the current object of this class
              */
            double GetLength();
            /*! \fn
              * An accessor function in order to access to the edge type
              * @return edge_type_ attribute of the current object of this class
              */
            gmml::EdgeTypeDescriptor GetEdgeType();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the force of the current object
              * Set the force_ attribute of the current edge attribute
              * @param force The force attribute of the current object
              */
            void SetForce(double force);
            /*! \fn
              * A mutator function in order to set the weight of the current object
              * Set the weight_ attribute of the current edge attribute
              * @param weight The weight attribute of the current object
              */
            void SetWeight(double weight);
            /*! \fn
              * A mutator function in order to set the length of the current object
              * Set the length_ attribute of the current edge attribute
              * @param length The weight attribute of the current object
              */
            void SetLength(double length);
            /*! \fn
              * A mutator function in order to set the structure of the current object
              * Set the edge_type_ attribute of the current edge attribute
              * @param edge_type The edge type attribute of the current object
              */
            void SetEdgeType(gmml::EdgeTypeDescriptor edge_type);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the edge attribute contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            double force_;
            double weight_;
            double length_;
            gmml::EdgeTypeDescriptor edge_type_;
    };
}

#endif // EDGEATTRIBUTE_HPP
