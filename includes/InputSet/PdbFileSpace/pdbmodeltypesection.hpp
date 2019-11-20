// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia, Dave Montgomery

#ifndef PDBMODELTYPESECTION_HPP
#define PDBMODELTYPESECTION_HPP

#include <string>
#include <vector>
#include <string>
#include <iostream>

namespace PdbFileSpace
{
    class PdbModelTypeSection
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbModelTypeSection();
            /*! \fn
              * Constructor with required parameters
              * @param record_name Name of the model type card record
              * @param comments List of comments that are in the model type card
              */
            PdbModelTypeSection(const std::string& record_name, const std::vector<std::string>& comments);
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbModelTypeSection(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the record name in model type card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the comments in model type card
              * @return comments_ attribute of the current object of this class
              */
            std::vector<std::string> GetComments();
/** @}*/
            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current model type card
              * @param record_name The record name of the current object
              */
            void SetRecordName(const std::string record_name);
            /*! \fn
              * A mutator function in order to set the list of comments of the current object
              * Set the comments_ of the current model type card
              * @param comments The comments of the current object
              */
            void SetComments(const std::vector<std::string> comments);
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the model type card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string record_name_;               /*!< Record name of model type card which is the first column of each line in a pdb file >*/
            std::vector<std::string> comments_;     /*!< List of comments for a model type card >*/
    };
}

#endif // PDBMODELTYPESECTION_HPP
