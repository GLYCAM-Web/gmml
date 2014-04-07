#ifndef PDBPREPROCESSORUNRECOGNIZEDRESIDUE_HPP
#define PDBPREPROCESSORUNRECOGNIZEDRESIDUE_HPP

#include <string>
#include <iostream>

namespace PdbPreprocessorSpace
{
    class PdbPreprocessorUnrecognizedResidue
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbPreprocessorUnrecognizedResidue();

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the residue chain Id
              * @return residue_chain_id_ attribute of the current object of this class
              */
            char GetResidueChainId();
            /*! \fn
              * An accessor function in order to access to the residue index
              * @return residue_index_ attribute of the current object of this class
              */
            int GetResidueIndex();
            /*! \fn
              * An accessor function in order to access to the residue name
              * @return residue_name_ attribute of the current object of this class
              */
            std::string GetResidueName();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the residue chain id of the current object
              * Set the residue_chain_id_ attribute of the current pdb preprocessor unrecognized residue
              * @param residue_chain_id The residue chain id attribute of the current object
              */
            void SetResidueChainId(char residue_chain_id);
            /*! \fn
              * A mutator function in order to set the residue index of the current object
              * Set the residue_index_ attribute of the current pdb preprocessor unrecognized residue
              * @param residue_index The residue index attribute of the current object
              */
            void SetResidueIndex(int residue_index);
            /*! \fn
              * A mutator function in order to set the residue name of the current object
              * Set the residue_name_ attribute of the current pdb preprocessor unrecognized residue
              * @param residue_name The residue name attribute of the current object
              */
            void SetResidueName(std::string residue_name);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb preprocessor unrecognized residue contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            char residue_chain_id_;
            int residue_index_;
            std::string residue_name_;

    };
}

#endif // PDBPREPROCESSORUNRECOGNIZEDRESIDUE_HPP
