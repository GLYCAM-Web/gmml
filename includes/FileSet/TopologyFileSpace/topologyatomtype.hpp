#ifndef TOPOLOGYATOMTYPE_HPP
#define TOPOLOGYATOMTYPE_HPP

#include <string>
#include <iostream>
#include <map>

namespace TopologyFileSpace
{
    class TopologyAtomType
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            typedef std::map<int, double> TopologyCoefficientAMap;
            typedef std::map<int, double> TopologyCoefficientBMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            TopologyAtomType();

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the atom type index
              * @return atom_type_index_ attribute of the current object of this class
              */
            int GetAtomTypeIndex();
            /*! \fn
              * An accessor function in order to access to the index
              * @return index_ attribute of the current object of this class
              */
            int GetIndex();
            /*! \fn
              * An accessor function in order to access to the coefficient a
              * @return coefficient_a_ attribute of the current object of this class
              */
            TopologyCoefficientAMap GetCoefficientA();
            /*! \fn
              * An accessor function in order to access to the coefficient b
              * @return coefficient_b_ attribute of the current object of this class
              */
            TopologyCoefficientBMap GetCoefficientB();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the atom type index of the current object
              * Set the atom_type_index_ attribute of the current topology atom type
              * @param atom_type_index The atom type attribute of the current object
              */
            void SetAtomTypeIndex(int atom_type_index);
            /*! \fn
              * A mutator function in order to set the index of the current object
              * Set the index_ attribute of the current topology atom type
              * @param index The index attribute of the current object
              */
            void SetIndex(int index);
            /*! \fn
              * A mutator function in order to set the coefficient a of the current object
              * Set the coefficient_a_ attribute of the current topology atom type
              * @param coefficient_a The coefficient a attribute of the current object
              */
            void SetCoefficientA(TopologyCoefficientAMap coefficient_a);
            /*! \fn
              * A mutator function in order to set the coefficient b of the current object
              * Set the coefficient_b_ attribute of the current topology atom type
              * @param coefficient_b The coefficient b attribute of the current object
              */
            void SetCoefficientB(TopologyCoefficientBMap coefficient_b);

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
            int atom_type_index_;
            int index_;
            TopologyCoefficientAMap coefficient_a_;
            TopologyCoefficientBMap coefficient_b_;

    };
}

#endif // TOPOLOGYATOMTYPE_HPP
