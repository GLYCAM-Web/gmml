#ifndef PDBMODEL_HPP
#define PDBMODEL_HPP

#include <string>
#include <sstream>

namespace PdbFileSpace
{
    class PdbModelResidueSet;

    class PdbModel
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////


            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            PdbModel();
            PdbModel(std::stringstream& model_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            int GetModelSerialNumber();
            PdbModelResidueSet* GetModelResidueSet();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetModelSerialNumber(int model_serial_number);
            void SetModelResidueSet(PdbModelResidueSet* model_residue_set);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            int model_serial_number_;
            PdbModelResidueSet* model_residue_set_;

    };
}

#endif // PDBMODEL_HPP
