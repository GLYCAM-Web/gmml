#ifndef PDBQTMODELCARD_HPP
#define PDBQTMODELCARD_HPP

#include <string>
#include <map>
#include <sstream>
#include <iostream>

namespace PdbqtFileSpace
{
    class PdbqtModel;
    class PdbqtModelCard
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * Mapping between model serial number and the model
              */
            typedef std::map<int, PdbqtModel*> PdbqtModelMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbqtModelCard();
            /*! \fn
              * Constructor with required parameters
              * @param stream_block
              */
            PdbqtModelCard(std::ifstream& stream_block);


            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the record name in a pdbqt model card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the models in a pdbqt model card
              * @return models_ attribute of the current object of this class
              */
            PdbqtModelMap GetModels();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current pdbqt model card
              * @param record_name The record name attribute of the current object
              */
            void SetRecordName(const std::string record_name);
            /*! \fn
              * A mutator function in order to set the models of the current object
              * Set the models_ attribute of the current pdbqt model card
              * @param models The model attribute of the current object
              */
            void SetModels(PdbqtModelMap models);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdbqt model card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;           /*!< Name of model card record which is in the first column of each line of a pdbqt file >*/
            PdbqtModelMap models_;                /*!< Models that are in model card of a pdbqt file >*/

    };
}

#endif // PDBQTMODELCARD_HPP
