#ifndef PDBQTREMARKCARD_HPP
#define PDBQTREMARKCARD_HPP

#include <string>
#include <iostream>

namespace PdbqtFileSpace
{
    class PdbqtRemarkCard
    {
        public:

            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbqtRemarkCard();
            /*! \fn
              * Constructor with required parameters
              * @param line
              */
            PdbqtRemarkCard(std::string line);


            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the record name in a pdbqt remark card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the value in a pdbqt remark card
              * @return value_ attribute of the current object of this class
              */
            std::string GetValue();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current pdbqt remark card
              * @param record_name The record name attribute of the current object
              */
            void SetRecordName(const std::string record_name);
            /*! \fn
              * A mutator function in order to set the value of the current object
              * Set the value_ attribute of the current pdbqt remark card
              * @param value The value attribute of the current object
              */
            void SetValue(const std::string value);


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
            std::string record_name_;
            std::string value_;
    };
}

#endif // PDBQTREMARKCARD_HPP
