// Created by: Dave Montgomery

#ifndef PDBSOURCECARD_HPP
#define PDBSOURCECARD_HPP

#include <string>
#include <iostream>


namespace PdbFileSpace
{
    class PdbSourceCard
    {
        public:

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbSourceCard();
            /*! \fn
              * Constructor with required parameters
              * @param line
              */
            PdbSourceCard(std::string& line);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the record name in a source
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the Token attribute in a source
              * @return token_ attribute of the current object of this class
              */
            std::string GetToken();
            /*! \fn
              * An accessor function in order to access to the Value attribute in a source
              * @return value_ attribute of the current object of this class
              */
            std::string GetValue();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current source
              * @param record_name The record name attribute of the current object
              */
            void SetRecordName(const std::string record_name);
            /*! \fn
              * A mutator function in order to set the Token attribute of the current object
              * Set the token_ attribute of the current source
              * @param token The Token attribute of the current object
              */
            void SetToken(std::string token_);
            /*! \fn
              * A mutator function in order to set the value attribute of the current object
              * Set the value_ attribute of the current source
              * @param value The Value attribute of the current object
              */
            void SetValue(std::string value);
/** @}*/
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the source contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;                       /*!< Record name of source card which is the first column of each line of the card >*/
            std::string token_;                             /*!< Token of source card >*/
            std::string value_;                             /*!< Value of the source card>*/
    };
}

#endif // PDBSOURCECARD_HPP
