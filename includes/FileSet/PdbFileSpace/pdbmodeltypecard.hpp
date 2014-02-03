// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBMODELTYPECARD_HPP
#define PDBMODELTYPECARD_HPP

#include <string>
#include <vector>
#include <string>
#include <iostream>

namespace PdbFileSpace
{
    class PdbModelTypeCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbModelTypeCard();
            /*! \fn
              * Constructor with required parameters
              * @param record_name
              * @param comments
              */
            PdbModelTypeCard(const std::string& record_name, const std::vector<std::string>& comments);
            PdbModelTypeCard(std::stringstream& stream_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
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

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
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
            std::vector<std::string> comments_;
    };
}

#endif // PDBMODELTYPECARD_HPP
