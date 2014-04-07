#ifndef PDBPREPROCESSORREPLACEDHYDROGEN_HPP
#define PDBPREPROCESSORREPLACEDHYDROGEN_HPP

#include <string>
#include <iostream>

namespace PdbPreprocessorSpace
{
    class PdbPreprocessorReplacedHydrogen
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbPreprocessorReplacedHydrogen();

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the atom index
              * @return atom_index_ attribute of the current object of this class
              */
            int GetAtomIndex();
            /*! \fn
              * An accessor function in order to access to the atom name
              * @return atom_name_ attribute of the current object of this class
              */
            std::string GetAtomName();
            /*! \fn
              * An accessor function in order to access to the residue name
              * @return residue_name_ attribute of the current object of this class
              */
            std::string GetResidueName();
            /*! \fn
              * An accessor function in order to access to the residue number
              * @return residue_number_ attribute of the current object of this class
              */
            int GetResidueNumber();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the atom index of the current object
              * Set the atom_index_ attribute of the current pdb preprocessor replaced hydrogen
              * @param atom_index The residue index attribute of the current object
              */
            void SetAtomIndex(int atom_index);
            /*! \fn
              * A mutator function in order to set the atom name of the current object
              * Set the atom_name_ attribute of the current pdb preprocessor replaced hydrogen
              * @param atom_name The atom name attribute of the current object
              */
            void SetAtomName(std::string atom_name);
            /*! \fn
              * A mutator function in order to set the residue name of the current object
              * Set the residue_name_ attribute of the current pdb preprocessor replaced hydrogen
              * @param residue_name The residue name attribute of the current object
              */
            void SetResidueName(std::string residue_name);
            /*! \fn
              * A mutator function in order to set the residue number of the current object
              * Set the residue_number_ attribute of the current pdb preprocessor replaced hydrogen
              * @param residue_number The residue number attribute of the current object
              */
            void SetResidueNumber(int residue_number);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the pdb preprocessor replaced hydrogen contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            int atom_index_;
            std::string atom_name_;
            std::string residue_name_;
            int residue_number_;

    };
}


#endif // PDBPREPROCESSORREPLACEDHYDROGEN_HPP
