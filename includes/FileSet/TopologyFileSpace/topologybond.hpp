#ifndef TOPOLOGYBOND_HPP
#define TOPOLOGYBOND_HPP

#include <string>
#include <iostream>
#include <vector>

namespace TopologyFileSpace
{
    class TopologyBondType;

    class TopologyBond
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            TopologyBond();

            TopologyBond(std::vector<std::string> bonds, std::vector<std::string> residue_names );

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the bonds
              * @return bonds_ attribute of the current object of this class
              */
            std::vector<std::string> GetBonds();
            /*! \fn
              * An accessor function in order to access to the bond type
              * @return bond_type_ attribute of the current object of this class
              */
            TopologyBondType* GetBondType();
            /*! \fn
              * An accessor function in order to access to boolean value that indicates if bond includes hydrogen
              * @return including_hydrogen_ attribute of the current object of this class
              */
            bool GetIncludingHydrogen();
            /*! \fn
              * An accessor function in order to access to the residue names
              * @return residue_names_ attribute of the current object of this class
              */
            std::vector<std::string> GetResidueNames();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A mutator function in order to set the bonds of the current object
              * Set the bonds_ attribute of the current topology bond
              * @param bonds The bonds attribute of the current object
              */
            void SetBonds(std::vector<std::string> bonds);
            /*! \fn
              * A mutator function in order to set the bond type of the current object
              * Set the bond_type_ attribute of the current topology bond type
              * @param bond_type The bond type attribute of the current object
              */
            void SetBondType(TopologyBondType* bond_type);
            /*! \fn
              * A mutator function in order to set the boolean including hydrogen attribute of the current object
              * Set the including_hydrogen_ attribute of the current topology bond type
              * @param including_hydrogen The including hydrogen attribute of the current object
              */
            void SetIncludingHydrogen(bool including_hydrogen);
            /*! \fn
              * A mutator function in order to set residue names of the current object
              * Set the residue_names_ attribute of the current topology bond type
              * @param residue_names The residue names of the current object
              */
            void SetResidueNames(std::vector<std::string> residue_names);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the topology bond contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::vector<std::string> bonds_;
            TopologyBondType* bond_type_;
            bool including_hydrogen_;
            std::vector<std::string> residue_names_;

    };
}

#endif // TOPOLOGYBOND_HPP
