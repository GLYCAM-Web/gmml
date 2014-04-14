#ifndef PDBPREPROCESSORDISULFIDEBOND_HPP
#define PDBPREPROCESSORDISULFIDEBOND_HPP

#include <string>
#include <iostream>

namespace PdbPreprocessorSpace
{
    class PdbPreprocessorDisulfideBond
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbPreprocessorDisulfideBond();

            PdbPreprocessorDisulfideBond(char chain_id_1, char chain_id_2, int residue_sequence_number_1, int residue_sequence_number_2, double distance, bool is_bonded);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the residue chain Id 1
              * @return residue_chain_id_1_ attribute of the current object of this class
              */
            char GetResidueChainId1();
            /*! \fn
              * An accessor function in order to access to the residue chain Id 2
              * @return residue_chain_id_2_ attribute of the current object of this class
              */
            char GetResidueChainId2();
            /*! \fn
              * An accessor function in order to access to the residue sequence number 1
              * @return residue_sequence_Number_1_ attribute of the current object of this class
              */
            int GetResidueSequenceNumber1();
            /*! \fn
              * An accessor function in order to access to the residue sequence number 2
              * @return residue_sequence_Number_2_ attribute of the current object of this class
              */
            int GetResidueSequenceNumber2();
            /*! \fn
              * An accessor function in order to access to the distance
              * @return distance_ attribute of the current object of this class
              */
            double GetDistance();
            /*! \fn
              * An accessor function in order to access to the is bonded attribute
              * @return is_bonded_ attribute of the current object of this class
              */
            bool GetIsBonded();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the residue chain id 1 of the current object
              * Set the residue_chain_id_1_ attribute of the current pdb preprocessor disulfide bond
              * @param residue_chain_id_1 The residue chain id 1 attribute of the current object
              */
            void SetResidueChainId1(char residue_chain_id_1);
            /*! \fn
              * A mutator function in order to set the residue chain id 2 of the current object
              * Set the residue_chain_id_2_ attribute of the current pdb preprocessor disulfide bond
              * @param residue_chain_id_2 The residue chain id 2 attribute of the current object
              */
            void SetResidueChainId2(char residue_chain_id_2);
            /*! \fn
              * A mutator function in order to set the residue sequence number 1 of the current object
              * Set the residue_sequence_number_1_ attribute of the current pdb preprocessor disulfide bond
              * @param residue_sequence_number_1 The residue sequence number 1 attribute of the current object
              */
            void SetResidueSequenceNumber1(int residue_sequence_number_1);
            /*! \fn
              * A mutator function in order to set the residue sequence number 2 of the current object
              * Set the residue_sequence_number_2_ attribute of the current pdb preprocessor disulfide bond
              * @param residue_sequence_number_2 The residue sequence number 2 attribute of the current object
              */
            void SetResidueSequenceNumber2(int residue_sequence_number_2);
            /*! \fn
              * A mutator function in order to set the distance of the current object
              * Set the distance_ attribute of the current pdb preprocessor disulfide bond
              * @param distance The distance attribute of the current object
              */
            void SetDistance(double distance);
            /*! \fn
              * A mutator function in order to set the is bonded of the current object
              * Set the is_bonded_ attribute of the current pdb preprocessor disulfide bond
              * @param is_bonded The is bonded attribute of the current object
              */
            void SetIsBonded(bool is_bonded);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb preprocessor disulfide bond contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            char residue_chain_id_1_;
            char residue_chain_id_2_;
            int residue_sequence_number_1_;
            int residue_sequence_number_2_;
            double distance_;
            bool is_bonded_;

    };
}

#endif // PDBPREPROCESSORDISULFIDEBOND_HPP
