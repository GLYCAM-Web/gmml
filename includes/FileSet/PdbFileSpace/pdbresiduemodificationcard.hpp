// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia// Author: Alireza Khatamian

#ifndef PDBRESIDUEMODIFICATIONCARD_HPP
#define PDBRESIDUEMODIFICATIONCARD_HPP

#include <map>
#include <string>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbResidueModification;
    class PdbResidueModificationCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::map<std::string, PdbResidueModification*> ResidueModificationMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbResidueModificationCard();
            /*! \fn
              * Constructor with required parameters
              * @param record_name
              */
            PdbResidueModificationCard(const std::string& record_name);
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbResidueModificationCard(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the record name in a residue modification card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the residue modifications in a residue modification card
              * @return residue_modifications_ attribute of the current object of this class
              */
            ResidueModificationMap GetResidueModifications();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current residue modification card
              * @param record_name The record name of the current object
              */
            void SetRecordName(const std::string record_name);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb residue modification card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string record_name_;
            ResidueModificationMap residue_modifications_;
    };
}

#endif // PDBRESIDUEMODIFICATIONCARD_HPP
