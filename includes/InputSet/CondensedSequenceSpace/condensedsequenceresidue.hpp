#ifndef CONDENSEDSEQUENCERESIDUE_HPP
#define CONDENSEDSEQUENCERESIDUE_HPP

#include <string>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>

namespace CondensedSequenceSpace
{
    class CondensedSequence;
    class CondensedSequenceResidue
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::map<int, std::string> DerivativeMap;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            CondensedSequenceResidue();
            CondensedSequenceResidue(std::string residue_string, CondensedSequenceSpace::CondensedSequence* condensed_sequence);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            bool GetIsTerminalSugar();
            bool GetIsTerminalAglycone();
            std::string GetIsomer();
            std::string GetConfiguration();
            std::string GetName();
            int GetAnomericCarbon();
            int GetOxygenPosition();
            DerivativeMap GetDerivatives();
            int GetParentId();
	    std::vector<int> GetChildIds();
	    int GetBondId(); //Added by Yao 08/03/2018. Bond Id is to label the index of the bond of a residue to its parent, starting from reducing end.Numbering starts from 0.
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            void SetIsTerminalSugar(bool is_terminal_sugar);
            void SetIsTerminalAglycone(bool is_terminal_aglycone);
            void SetIsomer(std::string isomer);
            void SetConfiguration(std::string configuration);
            void SetName(std::string name);
            void SetAnomericCarbon(int anomeric_carbon);
            void SetOxygenPosition(int oxygen_position);
            void SetDerivatives(DerivativeMap derivatives);
            void SetParentId(int parent_id);
	    void AddChildId(int child_id);
	    void RemoveChildId(int child_id);
	    void SetBondId(int bond_id); //Added by Yao 08/03/2018. Bond Id is to label the index of the bond of a residue to its parent, starting from reducing end.Numbering starts from 0.
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the condensed sequence residue contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            bool is_terminal_aglycone_ = false;
            bool is_terminal_sugar_ = false;
            std::string isomer_;
            std::string configuration_;
            std::string name_;
            int anomeric_carbon_;
            int oxygen_position_;
            DerivativeMap derivatives_;
            int parent_id_;
	    std::vector<int> child_ids_;
	    int bond_id_;  //Added by Yao 08/03/2018. Bond Id is to label the index of the bond of a residue to its parent, starting from reducing end.Numbering starts from 0.
    };
}

#endif // CONDENSEDSEQUENCERESIDUE_HPP
