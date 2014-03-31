// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBLINKRESIDUE_HPP
#define PDBLINKRESIDUE_HPP

#include <string>
#include <iostream>

namespace PdbFileSpace
{
    class PdbLinkResidue
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbLinkResidue();
            /*! \fn
              * Constructor with required parameters
              * @param atom_name
              * @param alternate_location_indicator
              * @param residue_name
              * @param residue_chain_identifier
              * @param residue_sequence_number
              * @param residue_insertion_code
              * @param symmetry_operator
              */
            PdbLinkResidue(const std::string& atom_name, char alternate_location_indicator, const std::string& residue_name, char residue_chain_identifier,
                           int residue_sequence_number, char residue_insertion_code, int symmetry_operator);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the atom name in a link residue
              * @return atom_name_ attribute of the current object of this class
              */
            std::string GetAtomName();
            /*! \fn
              * An accessor function in order to access to the alternate location indicator in a link residue
              * @return alternate_location_indicator_ attribute of the current object of this class
              */
            char GetAlternateLocationIndicator();
            /*! \fn
              * An accessor function in order to access to the residue name in a link residue
              * @return residue_name_ attribute of the current object of this class
              */
            std::string GetResidueName();
            /*! \fn
              * An accessor function in order to access to the residue chain identifier in a link residue
              * @return residue_chain_identifier_ attribute of the current object of this class
              */
            char GetResidueChainId();
            /*! \fn
              * An accessor function in order to access to the residue sequence number in a link residue
              * @return residue_sequence_number_ attribute of the current object of this class
              */
            int GetResidueSequenceNumber();
            /*! \fn
              * An accessor function in order to access to the residue insertion code in a link residue
              * @return residue_insertion_code_ attribute of the current object of this class
              */
            char GetResidueInsertionCode();
            /*! \fn
              * An accessor function in order to access to the symmetry operator in a link residue
              * @return symmetry_operator_ attribute of the current object of this class
              */
            int GetSymmetryOperator();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the atom name of the current object
              * Set the atom_name_ attribute of the current link residue
              * @param atom_name The atom name of the current object
              */
            void SetAtomName(const std::string atom_name);
            /*! \fn
              * A mutator function in order to set the alternate location indicator of the current object
              * Set the alternate_location_indicator_ attribute of the current link residue
              * @param alternate_location_indicator The alternate location indicator of the current object
              */
            void SetAlternateLocationIndicator(char alternate_location_indicator);
            /*! \fn
              * A mutator function in order to set the residue name of the current object
              * Set the residue_name_ attribute of the current link residue
              * @param residue_name The residue name of the current object
              */
            void SetResidueName(const std::string residue_name);
            /*! \fn
              * A mutator function in order to set the residue chain identifier of the current object
              * Set the residue_chain_identifier_ attribute of the current link residue
              * @param residue_chain_identifier The residue chain identifier of the current object
              */
            void SetResidueChainId(char residue_chain_identifier);
            /*! \fn
              * A mutator function in order to set the residue sequence number of the current object
              * Set the residue_sequence_number_ attribute of the current link residue
              * @param residue_sequence_number The residue sequence number of the current object
              */
            void SetResidueSequenceNumber(int residue_sequence_number);
            /*! \fn
              * A mutator function in order to set the residue insertion code of the current object
              * Set the residue_insertion_code_ attribute of the current link residue
              * @param residue_insertion_code The residue insertion code of the current object
              */
            void SetResidueInsertionCode(char residue_insertion_code);
            /*! \fn
              * A mutator function in order to set the symmetry operator of the current object
              * Set the symmetry_operator_ attribute of the current link residue
              * @param symmetry_operator The symmetry operator of the current object
              */
            void SetSymmetryOperator(int symmetry_operator);

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
            std::string atom_name_;
            char alternate_location_indicator_;
            std::string residue_name_;
            char residue_chain_identifier_;
            int residue_sequence_number_;
            char residue_insertion_code_;
            int symmetry_operator_;
    };
}

#endif // PDBLINKRESIDUE_HPP
