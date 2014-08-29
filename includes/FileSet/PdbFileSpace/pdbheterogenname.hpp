// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBHETEROGENNAME_HPP
#define PDBHETEROGENNAME_HPP

#include <string>
#include <iostream>

namespace PdbFileSpace
{
    class PdbHeterogenName
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbHeterogenName();
            /*! \fn
              * Constructor with required parameters
              * @param heterogen_identifier
              * @param heterogen_name
              */
            PdbHeterogenName(const std::string& heterogen_identifier, const std::string& heterogen_name);
            /*! \fn
              * Constructor with required parameters
              * @param specification_block
              */
            PdbHeterogenName(std::stringstream& specification_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the heterogen identifier in a heterogen name class
              * @return heterogen_identifier_ attribute of the current object of this class
              */
            std::string GetHeterogenIdentifier();
            /*! \fn
              * An accessor function in order to access to the heterogen name in a heterogen name class
              * @return heterogen_name_ attribute of the current object of this class
              */
            std::string GetHeterogenName();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the heterogen identifier of the current object
              * Set the heterogen_identifier_ attribute of the current heterogen name
              * @param heterogen_identifier The heterogen identifier of the current object
              */
            void SetHeterogenIdentifier(const std::string heterogen_identifier);
            /*! \fn
              * A mutator function in order to set the heterogen name of the current object
              * Set the heterogen_name_ attribute of the current heterogen name
              * @param heterogen_name The heterogen name of the current object
              */
            void SetHeterogenName(const std::string heterogen_name);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the heterogen name contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string heterogen_identifier_;
            std::string heterogen_name_;
    };
}
#endif // PDBHETEROGENNAME_HPP
