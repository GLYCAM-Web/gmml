#ifndef CONDENSEDSEQUENCEGLYCAM06RESIDUE_HPP
#define CONDENSEDSEQUENCEGLYCAM06RESIDUE_HPP
#include <string>
#include <vector>
#include <iostream>

namespace CondensedSequenceSpace
{
    class CondensedSequenceGlycam06Residue    {
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
            CondensedSequenceGlycam06Residue();
            CondensedSequenceGlycam06Residue(std::string name);
            CondensedSequenceGlycam06Residue(std::string name, std::string anomeric_carbon, std::string parent_oxygen, bool is_derivative = false);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            std::string GetName();
            std::string GetAnomericCarbon();
            std::string GetParentOxygen();
            bool GetIsDerivative();
            int GetParentId();
	    int GetBondId();  //Added by Yao Xiao 08/03/2018
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            void SetName(std::string name);
            void SetAnomericCarbon(std::string anomeric_carbon);
            void SetParentOxygen(std::string parent_oxygen);
            void SetIsDerivative(bool is_derivative);
            void SetParentId(int parent_id);
	    void SetBondId(int bond_id);  //Added by Yao Xiao 08/03/2018
/** @}*/
            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the condensed sequence tree residue contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string name_;
            std::string anomeric_carbon_;
            std::string parent_oxygen_;
            bool is_derivative_;
            int parent_id_;
	    int bond_id_;  //Added by Yao Xiao 08/03/2018
    };
}

#endif // CONDENSEDSEQUENCEGLYCAM06RESIDUE_HPP
