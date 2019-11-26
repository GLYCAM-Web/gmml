#ifndef PDBQTTORSIONALDOFCARD_HPP
#define PDBQTTORSIONALDOFCARD_HPP

#include <string>
#include <iostream>

namespace PdbqtFileSpace
{
    class PdbqtTorsionalDoFCard
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
            PdbqtTorsionalDoFCard();
            /*! \fn
              * Constructor with required parameters
              * @param line
              */
            PdbqtTorsionalDoFCard(std::string line);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the record name in a pdbqt torsional degree of freedom card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the value in a pdbqt torsional degree of freedom
              * @return number_of_torsional_dof_ attribute of the current object of this class
              */
            int GetNumberofTorsionalDoF();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current pdbqt number of torsional degree of freedom card
              * @param record_name The record name attribute of the current object
              */
            void SetRecordName(const std::string record_name);
            /*! \fn
              * A mutator function in order to set the value of the current object
              * Set the value_ attribute of the current pdbqt number of torsional degree of freedom card
              * @param number_of_torsional_dof The torsional degree of freedom attribute of the current object
              */
            void SetNumberOfTorsionalDoF(int number_of_torsional_dof);


            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the number of torsional degree of freedom card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;
            int number_of_tosional_dof_;


    };
}

#endif // PDBQTTORSIONALDOFCARD_HPP
