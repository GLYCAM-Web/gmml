#ifndef TOPOLOGYATOMTYPE_HPP
#define TOPOLOGYATOMTYPE_HPP

#include <string>
#include <iostream>

namespace TopologyFileSpace
{
    class TopologyAtomType
    {
        public:
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
              * An accessor function in order to access to the atom type
              * @return atom_type_ attribute of the current object of this class
              */
            std::string GetAtomType();
            /*! \fn
              * An accessor function in order to access to the atom type index
              * @return index_ attribute of the current object of this class
              */
            int GetIndex();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the atom type of the current object
              * Set the atom_type_ attribute of the current topology atom type
              * @param atom_type The atom type attribute of the current object
              */
            void SetAtomType(const std::string atom_type);
            /*! \fn
              * A mutator function in order to set the index of the current object
              * Set the index_ attribute of the current topology atom type
              * @param index The index attribute of the current object
              */
            void SetIndex(int index);

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
            int index_;

    };
}

#endif // TOPOLOGYATOMTYPE_HPP
