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
            /*! \typedef
              * A vector of pdbqt atom cards
              */
            typedef std::vector<PdbqtAtomCard*> AtomCardVector;

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            PdbqtBranchCard();

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
            AtomCardVector GetRotatableAtomSet();
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
            void SetRotatableAtomSet(AtomCardVector rotatable_atom_set);
            /*! \fn
              * A function in order to add the rotatable atom to the current object
              * Set the rotatable_atom_ attribute of the current pdbqt branch card
              * @param rotatable_atom The rotatable atom of the current object
              */
            void AddRotatableAtom(PdbqtAtomCard* rotatable_atom);
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
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string record_name_;
            int solid_atom_serial_number_;
            int rotatable_atom_serial_number_;
            AtomCardVector rotatable_atom_set_;
            BranchCardVector sub_branches_;
    };
}

#endif // PDBQTBRANCHCARD_HPP
