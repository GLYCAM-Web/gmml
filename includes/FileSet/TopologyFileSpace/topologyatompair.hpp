#ifndef TOPOLOGYATOMPAIR_HPP
#define TOPOLOGYATOMPAIR_HPP

#include <string>
#include <iostream>
#include <map>

namespace TopologyFileSpace
{
    class TopologyAtomPair
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            typedef std::map<std::string, double> TopologyCoefficientMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            TopologyAtomPair();

            TopologyAtomPair(std::string atom_type, TopologyCoefficientMap coefficient_a_map, TopologyCoefficientMap coefficient_b_map);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////           
            /*! \fn
              * An accessor function in order to access to the atom type
              * @return atom_type_ attribute of the current object of this class
              */
            std::string GetAtomType();
            /*! \fn
              * An accessor function in order to access to the coefficient a
              * @return coefficient_a_ attribute of the current object of this class
              */
            TopologyCoefficientMap GetCoefficientA();
            /*! \fn
              * An accessor function in order to access to the coefficient b
              * @return coefficient_b_ attribute of the current object of this class
              */
            TopologyCoefficientMap GetCoefficientB();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the atom type of the current object
              * Set the atom_type_ attribute of the current topology atom type
              * @param atom_type The atom type attribute of the current object
              */
            void SetAtomType(std::string atom_type);
            /*! \fn
              * A mutator function in order to set the coefficient a of the current object
              * Set the coefficient_a_ attribute of the current topology atom type
              * @param coefficient_a The coefficient a attribute of the current object
              */
            void SetCoefficientA(TopologyCoefficientMap coefficient_a);
            /*! \fn
              * A mutator function in order to set the coefficient b of the current object
              * Set the coefficient_b_ attribute of the current topology atom type
              * @param coefficient_b The coefficient b attribute of the current object
              */
            void SetCoefficientB(TopologyCoefficientMap coefficient_b);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the atom type contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string atom_type_;
            TopologyCoefficientMap coefficient_a_;
            TopologyCoefficientMap coefficient_b_;

    };
}

#endif // TOPOLOGYATOMPAIR_HPP
