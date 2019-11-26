// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia, Dave Montgomery

#ifndef PDBCOMPOUNDSECTION_HPP
#define PDBCOMPOUNDSECTION_HPP

#include <string>
#include <map>
#include <iostream>

namespace PdbFileSpace
{
    class PdbCompoundSpecification;
    class PdbCompoundSection
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * A mapping between a molecule id and its assigned compound specifications in a compound card of a pdb file
              */
            typedef std::map<std::string, PdbCompoundSpecification*> PdbCompoundSpecificationMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbCompoundSection();
            /*! \fn
              * Constructor with required parameters
              * @param record_name
              */
            PdbCompoundSection(const std::string& record_name);
            /*! \fn
              * A constructor that get a stream block of compound card and parse the whole block to fill the related fields
              * @param stream_block A whole block of compound card in a pdb file
              */
            PdbCompoundSection(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the record name in a compound card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to map of compound specifications of current object
              * @return compound_specifications_ attribute of the current object of this class
              */
            PdbCompoundSpecificationMap GetCompoundSpecifications();
/** @*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current compund card
              * @param record_name The record name of the current object
              */
            void SetRecordName(const std::string record_name);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the compound card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;                               /*!< Record name of compound card in a pdb file: "COMPND" */
            PdbCompoundSpecificationMap compound_specifications_;   /*!< Map of all compound specifications of a compound record by its molecule id */

    };
}

#endif // PDBCOMPOUNDSECTION_HPP
