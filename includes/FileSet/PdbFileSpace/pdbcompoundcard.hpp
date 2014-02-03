// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBCOMPOUNDCARD_HPP
#define PDBCOMPOUNDCARD_HPP

#include <string>
#include <map>
#include <iostream>

namespace PdbFileSpace
{
    class PdbCompoundSpecification;

    class PdbCompoundCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            typedef std::map<std::string, PdbCompoundSpecification*> PdbCompoundSpecificationMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbCompoundCard();
            /*! \fn
              * Constructor with required parameters
              * @param record_name
              */
            PdbCompoundCard(const std::string& record_name);
            PdbCompoundCard(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
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

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current compund card
              * @param record_name The record name of the current object
              */
            void SetRecordName(const std::string record_name);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;
            PdbCompoundSpecificationMap compound_specifications_;

    };
}

#endif // PDBCOMPOUNDCARD_HPP
