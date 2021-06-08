// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia, Dave Montgomery

#ifndef PDBDISULFIDEBONDSECTION_HPP
#define PDBDISULFIDEBONDSECTION_HPP

#include <string>
#include <map>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbDisulfideResidueBond;
    class PdbDisulfideBondSection
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * A mapping between disulfide residue serial number and its belonging records in a disulfide residue bond card in a pdb file
              */
            typedef std::map<int, PdbDisulfideResidueBond*> DisulfideResidueBondMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbDisulfideBondSection();
            /*! \fn
              * Constructor with required parameters
              * @param record_name
              */
            PdbDisulfideBondSection(const std::string& record_name);
            /*! \fn
              * A constructor that get a stream block of disulfide bond card and parse the whole block to fill the related fields
              * @param stream_block A whole block of disulfide bond card in a pdb file
              */
            PdbDisulfideBondSection(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the record name in a disulfide bond card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the disulfide residue bonds in a disulfide bond card
              * @return residue_name_ attribute of the current object of this class
              */
            DisulfideResidueBondMap GetDisulfideResidueBonds();
/** @}*/
            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current disulfide bond card
              * @param record_name The record name of the current object
              */
            void SetRecordName(const std::string record_name);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the disulfide residue bond card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string record_name_;                           /*!< Record name of disulfide bond in a pdb file: "SSBOND" */
            DisulfideResidueBondMap disulfide_residue_bonds_;   /*!< Map of disulfide bonds in a disulfide bond card in a pdb file by their serial numbers */
    };
}

#endif // PDBDISULFIDEBONDSECTION_HPP
