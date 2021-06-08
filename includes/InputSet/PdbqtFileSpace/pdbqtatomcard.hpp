#ifndef PDBQTATOMCARD_HPP
#define PDBQTATOMCARD_HPP

#include <string>
#include <sstream>
#include <map>
#include <iostream>

namespace PdbqtFileSpace
{
    class PdbqtAtom;
    class PdbqtAtomCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * A mapping between atom serial number and its belonging factors and information in a model in a pdbqt file
              */
            typedef std::map<int, PdbqtAtom*> PdbqtAtomMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbqtAtomCard();
            /*! \fn
              * A constructor that get a stream block of atom card and parse the whole block to fill the related fields
              * @param stream_block A whole block of atoms belonging to a model in a pdbqt file
              */
            PdbqtAtomCard(std::ifstream& stream_block);
            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the record name in a atom card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the atoms in a atom card
              * @return atoms_ attribute of the current object of this class
              */
            PdbqtAtomMap GetAtoms();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current qt atom card
              * @param record_name The record name attribute of the current object
              */
            void SetRecordName(const std::string record_name);
            /*! \fn
              * A mutator function in order to set the atoms of the current object
              * Set the atoms_ attribute of the current qt atom card
              * @param atoms Atom attribute of the current object
              */
            void SetAtoms(PdbqtAtomMap atoms);
	    void AddAtom(PdbqtAtom* atom);  //Yao added 11/21/2019

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the qt atom card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;           /*!< Record name of qt atom card: "ATOM" */
            PdbqtAtomMap atoms_;                /*!< Map of all atoms informatin that belong to a specific model in a qt pdb file by their serial numbers */

    };
}

#endif // PDBQTATOMCARD_HPP
