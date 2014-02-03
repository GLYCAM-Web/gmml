// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBHETEROGENNAMECARD_HPP
#define PDBHETEROGENNAMECARD_HPP

#include <string>
#include <map>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbHeterogenName;
    class PdbHeterogenNameCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::map<std::string, PdbHeterogenName*> HeterogenNameMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbHeterogenNameCard();
            /*! \fn
              * Constructor with required parameters
              * @param record_name
              */
            PdbHeterogenNameCard(const std::string& record_name);
            PdbHeterogenNameCard(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the record name in a heterogen name card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the heterogen names in a heterogen name card
              * @return heterogen_names_ attribute of the current object of this class
              */
            HeterogenNameMap GetHeterogenNames();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current heterogen name card
              * @param record_name The record name of the current object
              */
            void SetRecordName(const std::string record_name);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string record_name_;
            HeterogenNameMap heterogen_names_;
    };
}

#endif // PDBHETEROGENNAMECARD_HPP
