// Created by: Dave Montgomery

#ifndef PDBREMARKSECTION_HPP
#define PDBREMARKSECTION_HPP

#include <string>
#include <map>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    // class PdbRemark;

    class PdbRemarkSection
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * Mapping between model serial number and the model
              */
            // typedef std::map<int, PdbRemark*> PdbRemarkMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbRemarkSection();
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbRemarkSection(std::stringstream& stream_block);
            /*! \fn
              * Constructor with required parameters
              * @param record_name Name for a title card record which appears in the first column of each line in a pdb file
              * @param title Title of a pdb file
              */
            PdbRemarkSection(const std::string& record_name, const std::string& remark_cards);


            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the record name in a Remark card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the Remarks in a Remark card
              * @return remark_cards_ attribute of the current object of this class
              */
            std::string GetRemarks();
            /*! \fn
              * An accessor function in order to access to the Resolution in a Remark card
              * @return resolution_ attribute of the current object of this class
              */
            float GetResolution();
            /*! \fn
              * An accessor function in order to access to the B Factor in a Remark card
              * @return b_factor_ attribute of the current object of this class
              */
            float GetBFactor();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current Remark card
              * @param record_name The record name attribute of the current object
              */
            void SetRecordName(const std::string record_name);
            /*! \fn
              * A mutator function in order to set the Remarks of the current object
              * Set the remark_cards_ attribute of the current remark card
              * @param remark_cards The remark attribute of the current object
              */
            void SetRemarks(const std::string remark_cards);
            /*! \fn
              * A mutator function in order to set the Resolution of the current object
              * Set the resolution_ attribute of the current remark card
              * @param resolution The Resolution attribute of the current object
              */
            void SetResolution(const float resolution);
            /*! \fn
              * A mutator function in order to set the B Factor of the current object
              * Set the b_factor_ attribute of the current remark card
              * @param b_factor The B Factor attribute of the current object
              */
            void SetBFactor(const float b_factor);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the remark card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;           /*!< Name of remark card record which is in the first column of each line of a pdb file >*/
            std::string remark_cards_;          /*!< Remarks that are in remark card of a pdb file >*/
            float resolution_;            /*!< Resolution of PDB >*/
            float b_factor_;              /*!< B Factor of PDB >*/

    };
}

#endif // PDBREMARKSECTION_HPP
