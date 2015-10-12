#ifndef CONDENSEDSEQUENCEAMBERPREPRESIDUE_HPP
#define CONDENSEDSEQUENCEAMBERPREPRESIDUE_HPP

#include <string>
#include <vector>
#include <iostream>

namespace CondensedSequenceSpace
{
    class CondensedSequenceAmberPrepResidue
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
            CondensedSequenceAmberPrepResidue();
            CondensedSequenceAmberPrepResidue(std::string name);
            CondensedSequenceAmberPrepResidue(std::string name, std::string anomeric_carbon, std::string parent_oxygen);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            std::string GetName();
            std::string GetAnomericCarbon();
            std::string GetParentOxygen();
//            bool GetHasParent();
            int GetParentId();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetName(std::string name);
            void SetAnomericCarbon(std::string anomeric_carbon);
            void SetParentOxygen(std::string parent_oxygen);
//            void SetHasParent(bool has_parent);
            void SetParentId(int parent_id);

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
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string name_;
            std::string anomeric_carbon_;
            std::string parent_oxygen_;
//            bool has_parent_;
            int parent_id_;
    };
}

#endif // CONDENSEDSEQUENCEAMBERPREPRESIDUE_HPP
