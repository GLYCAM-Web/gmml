#ifndef TOPOLOGYDIHEDRAL_HPP
#define TOPOLOGYDIHEDRAL_HPP

#include <string>
#include <iostream>
#include <vector>

namespace TopologyFileSpace
{
    class TopologyDihedralType;
    class TopologyDihedral
    {
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            TopologyDihedral();
            /*! \fn
              * Constructor with required parameters
              * @param dihedral_atoms Atoms invovled in a dihedral (four bonded atoms)
              * @param residue_names Residue names of the atoms involved in a dihedral (four residue names)
              */
            TopologyDihedral(std::vector<std::string> dihedral_atoms, std::vector<std::string> residue_names);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to the dihedrals
              * @return dihedrals_ attribute of the current object of this class
              */
            std::vector<std::string> GetDihedrals();
            /*! \fn
              * An accessor function in order to access to the dihedral type
              * @return dihedral_type_ attribute of the current object of this class
              */
            TopologyDihedralType* GetDihedralType();
            /*! \fn
              * An accessor function in order to access to boolean value that indicates if dihedral is improper
              * @return is_improper_ attribute of the current object of this class
              */
            bool GetIsImproper();
            /*! \fn
              * An accessor function in order to access to boolean value that indicates if dihedral is in ignored group interaction
              * @return ignored_group_interaction_ attribute of the current object of this class
              */
            bool GetIgnoredGroupInteraction();
            /*! \fn
              * An accessor function in order to access to boolean value that indicates if dihedral includes hydrogen
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
              * A mutator function in order to set the dihedrals of the current object
              * Set the dihedrals_ attribute of the current topology dihedral
              * @param dihedrals The dihedrals attribute of the current object
              */
            void SetDihedrals(std::vector<std::string> dihedrals);
            /*! \fn
              * A mutator function in order to set the dihedral type of the current object
              * Set the dihedral_type_ attribute of the current topology dihedral
              * @param dihedral_type The dihedral type attribute of the current object
              */
            void SetDihedralType(TopologyDihedralType* dihedral_type);
            /*! \fn
              * A mutator function in order to set the boolean improper attribute of the current object
              * Set the is_improper_ attribute of the current topology dihedral
              * @param is_improper The is improper attribute of the current object
              */
            void SetIsImproper(bool is_improper);
            /*! \fn
              * A mutator function in order to set the boolean ignored group interaction attribute of the current object
              * Set the ignored_group_interaction_ attribute of the current topology dihedral
              * @param ignored_group_interaction The ignored group interactionr attribute of the current object
              */
            void SetIgnoredGroupInteraction(bool ignored_group_interaction);
            /*! \fn
              * A mutator function in order to set the boolean including hydrogen attribute of the current object
              * Set the including_hydrogen_ attribute of the current topology dihedral
              * @param including_hydrogen The ignored group interactionr attribute of the current object
              */
            void SetIncludingHydrogen(bool including_hydrogen);
            /*! \fn
              * A mutator function in order to set residue names of the current object
              * Set the residue_names_ attribute of the current topology dihedral
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
              * A function to print out the topology dihedral contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::vector<std::string> dihedrals_;            /*!< Four bonded atom names involved in a dihedral >*/
            TopologyDihedralType* dihedral_type_;           /*!< Type of the dihedral regarding to topology dihedral types >*/
            bool is_improper_;                              /*!< Indicates that the dihedral is proper or not >*/
            bool ignored_group_interaction_;                /*!< Indicates that the dihedral has to be ignored in group interaction or not >*/
            bool including_hydrogen_;                       /*!< Indicates that the dihedral has a hydrogen atom or not >*/
            std::vector<std::string> residue_names_;        /*!< Names of residues of the atoms involved in a dihedral >*/
    };
}

#endif // TOPOLOGYDIHEDRAL_HPP
