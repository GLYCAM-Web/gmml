// Created by: Dave Montgomery

#ifndef PDBOBSOLETESECTION_HPP
#define PDBOBSOLETESECTION_HPP

#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbObsoleteCard;
    class PdbObsoleteSection
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * Vector of PDB Obsolete Cards
              */
            typedef std::vector<PdbObsoleteCard*> ObsoleteCardVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbObsoleteSection();
            /*! \fn
              * Constructor with required parameters
              * @param record_name
              */
            PdbObsoleteSection(const std::string& record_name,
                                const std::string& continuation,
                                const std::string& replacement_date,
                                const std::vector<std::string>& identifier_codes);
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbObsoleteSection(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the record name in a obsolete card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the continuation in a obsolete card
              * @return continuation_ attribute of the current object of this class
              */
            std::string GetContinuation();
            /*! \fn
              * An accessor function in order to access to the replacement date in a obsolete card
              * @return replacement_date_ attribute of the current object of this class
              */
            std::string GetReplacementDate();
            /*! \fn
              * An accessor function in order to access to the identifier codes in a obsolete card
              * @return identifier_codes_ attribute of the current object of this class
              */
            std::vector<std::string> GetIdentifierCodes();
            /*! \fn
              * An accessor function in order to access to the obsoletes in a obsolete card
              * @return obsolete_cards_ attribute of the current object of this class
              */
            ObsoleteCardVector GetObsoleteCards();
/** @}*/
            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{pdbobsolete
              */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current obsolete card
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
              * A function to print out the obsolete section contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string record_name_;             /*!< Record name of obsolete card in a pdb file >*/
            std::string continuation_;            /*!< Continuation obsolete card in a pdb file >*/
            std::string replacement_date_;        /*!< Replacement Date of obsolete card in a pdb file >*/
            std::vector<std::string> identifier_codes_;        /*!< ID Codes of replacement PDB entries of obsolete card in a pdb file >*/
            ObsoleteCardVector obsolete_cards_;   /*!< Vector of obsolete cards >*/
    };
}
#endif // PDBOBSOLETESECTION_HPP
