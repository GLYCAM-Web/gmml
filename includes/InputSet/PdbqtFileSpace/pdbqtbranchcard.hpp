#ifndef PDBQTBRANCHCARD_HPP
#define PDBQTBRANCHCARD_HPP

#include <string>
#include <iostream>
#include <vector>

namespace PdbqtFileSpace
{
    class PdbqtAtomCard;
    class PdbqtBranchCard
    {
        public:

            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * A vector of pdbqt branch cards
              */
            typedef std::vector<PdbqtBranchCard*> BranchCardVector;

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbqtBranchCard();
            /*! \fn
              * A constructor that get a stream block of branch card and parse the whole block to fill the related fields
              * @param stream_block A whole block of branches belonging to a model in a pdbqt file
              */
            PdbqtBranchCard(std::ifstream& stream_block, std::vector<PdbqtAtomCard*>& ACV);
            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the record name in a pdbqt branch card
              * @return record_name_ attribute of the current object of this class
              */
            std::string GetRecordName();
            /*! \fn
              * An accessor function in order to access to the solid atom serial number in a pdbqt branch card
              * @return solid_atom_serial_number_ attribute of the current object of this class
              */
            int GetSolidAtomSerialNumber();
            /*! \fn
              * An accessor function in order to access to the rotatble atom serial number in a pdbqt branch card
              * @return rotatable_atom_serial_number_ attribute of the current object of this class
              */
            int GetRotatbleAtomSerialNumber();
            /*! \fn
              * An accessor function in order to access to the pdbqt rotatble atom set
              * @return rotatable_atom_set_ attribute of the current object of this class
              */
            PdbqtAtomCard* GetRotatableAtomSet();
            /*! \fn
              * An accessor function in order to access to the pdbqt branch cards
              * @return sub_branches_ attribute of the current object of this class
              */
            BranchCardVector GetSubBranches();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the record name of the current object
              * Set the record_name_ attribute of the current pdbqt branch card
              * @param record_name The record name attribute of the current object
              */
            void SetRecordName(const std::string record_name);
            /*! \fn
              * A mutator function in order to set the solid atom serial number of the current object
              * Set the solid_atom_serial_number_ attribute of the current pdbqt branch card
              * @param solid_atom_serial_number The solid atom serial number attribute of the current object
              */
            void SetSolidAtomSerialNumber(int solid_atom_serial_number);
            /*! \fn
              * A mutator function in order to set the rotatable atom serial number of the current object
              * Set the rotatable_atom_serial_number_ attribute of the current pdbqt branch card
              * @param rotatable_atom_serial_number The rotatable atom serial number attribute of the current object
              */
            void SetRotatableAtomSerialNumber(int rotatable_atom_serial_number);
            /*! \fn
              * A mutator function in order to set rotatble atom set of the current object
              * Set the rotatable_atom_set_ attribute of the current pdbqt branch card
              * @param rotatable_atom_set The rotatable atom set attribute of the current object
              */
            void SetRotatableAtomSet(PdbqtAtomCard* rotatable_atom_set);
            /*! \fn
              * A mutator function in order to set sub branches of the current object
              * Set the sub_branches_ attribute of the current pdbqt branch card
              * @param sub_branches The sub branches attribute of the current object
              */
            void SetSubBranches(BranchCardVector sub_branches);
            /*! \fn
              * A function in order to add the sub branch to the current object
              * Set the sub_branches_ attribute of the current pdbqt branch card
              * @param sub_branches The sub branch of the current object
              */
            void AddSubBranch(PdbqtBranchCard* sub_branch);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the branch card contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;
            int solid_atom_serial_number_;
            int rotatable_atom_serial_number_;
            PdbqtAtomCard* rotatable_atom_set_;
            BranchCardVector sub_branches_;
    };
}

#endif // PDBQTBRANCHCARD_HPP
