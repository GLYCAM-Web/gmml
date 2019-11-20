#ifndef RESIDUE_HPP
#define RESIDUE_HPP


#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "residueproperties.hpp"
#include "../GeometryTopology/coordinate.hpp"
#include "../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../includes/ParameterSet/PrepFileSpace/prepfileprocessingexception.hpp"



namespace MolecularModeling
{
    class Assembly;
    class Atom;
    class ResidueNode;
    //class PrepFileResidue; //This is not in the MolecularModeling namespace
    class Residue; // Forward declare for the vector typedef
    typedef std::vector<MolecularModeling::Residue*> ResidueVector;
    class Residue : public ResidueProperties
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
           // typedef std::vector<Atom*> AtomVector;
            typedef std::vector<std::string> StringVector;
        //typedef PrepFileSpace::PrepFileResidue PrepFileResidue; // What was this doing? Ah ok, it was bandaging a deeper issue.
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            Residue();
            Residue(Assembly* assembly, std::string name);
            Residue(Residue* residue);
            Residue(Residue& residue);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the assembly
              * @return assembly_ attribute of the current object of this class
              */
            Assembly* GetAssembly();
            /*! \fn
              * An accessor function in order to access to the name
              * @return name_ attribute of the current object of this class
              */
            std::string GetName();
            /*! \fn
              * An accessor function in order to access to the number
              * @return number split from id_ attribute of the current object of this class
              */
            std::string GetNumber();
            /*! \fn
              * An accessor function in order to access to the atoms
              * @return atoms_ attribute of the current object of this class
              */
            AtomVector GetAtoms();
            /*! \fn
              * An accessor function in order to access to the head atoms
              * @return head_atoms_ attribute of the current object of this class
              */
            AtomVector GetHeadAtoms();
            /*! \fn
              * An accessor function in order to access to the tail atoms
              * @return tail_atoms_ attribute of the current object of this class
              */
            AtomVector GetTailAtoms();
            /*! \fn
              * An accessor function in order to access to the chemical type
              * @return chemical_residue_ attribute of the current object of this class
              */
            std::string GetChemicalType();
            /*! \fn
              * An accessor function in order to access to the description
              * @return description_ attribute of the current object of this class
              */
            std::string GetDescription();
            /*! \fn
              * An accessor function in order to access to the id
              * @return id_ attribute of the current object of this class
              */
            std::string GetId();
            /*! \fn                                                                              //Added by ayush on 11/20/17 for residuenodes in assembly
              * An accessor function in order to access to the node
              * @return node_ attribute of the current object of this class
              */
            ResidueNode* GetNode();

            bool GetIsSugarDerivative();

            bool GetIsAglycon();
            /*! \fn
            * An accessor function in order to access to the index
            * @return index_ attribute of the current object of this class
            */
            unsigned long long GetIndex() const;
            
            bool GetIsSugar();

/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A mutator function in order to set the assembly of the current object
              * Set the assembly_ attribute of the current residue
              * @param assembly The assembly attribute of the current object
              */
            void SetAssembly(Assembly* assembly);
            /*! \fn
              * A mutator function in order to set the name of the current object
              * Set the name_ attribute of the current residue
              * @param name The name attribute of the current object
              */
            void SetName(std::string name);
            /*! \fn
              * A mutator function in order to set the atoms of the current object
              * Set the atoms_ attribute of the current residue
              * @param atoms The atoms attribute of the current object
              */
            void SetAtoms(AtomVector atoms);
            /*! \fn
              * A function in order to add the atom to the current object
              * Set the atoms_ attribute of the current residue
              * @param atom The atom of the current object
              */
            void AddAtom(Atom* atom);
            /*! \fn
              * A function in order to remove the atom from the current object
              * Set the atoms_ attribute of the current residue
              * @param atom The atom of the current object
              */
            void RemoveAtom(Atom* atom, bool remove_bonds = true);
            /*! \fn
              * A mutator function in order to set the head atoms of the current object
              * Set the head_atoms_ attribute of the current residue
              * @param head_atoms The head atoms attribute of the current object
              */
            void SetHeadAtoms(AtomVector head_atoms);
            /*! \fn
              * A function in order to add the head atom to the current object
              * Set the head_atom_ attribute of the current residue
              * @param head_atom The head atom of the current object
              */
            void AddHeadAtom(Atom* head_atom);
            /*! \fn
              * A mutator function in order to set the tail atoms of the current object
              * Set the tail_atoms_ attribute of the current residue
              * @param tail_atoms The head atoms attribute of the current object
              */
            void SetTailAtoms(AtomVector tail_atoms);
            /*! \fn
              * A function in order to add the tail atom to the current object
              * Set the tail_atoms_ attribute of the current residue
              * @param tail_atom The tail atom of the current object
              */
            void AddTailAtom(Atom* tail_atom);
            /*! \fn
              * A mutator function in order to set the chemical type of the current object
              * Set the chemical_type_ attribute of the current residue
              * @param chemical_type The chemical type attribute of the current object
              */
            void SetChemicalType(std::string chemical_type);
            /*! \fn
              * A mutator function in order to set the description of the current object
              * Set the description_ attribute of the current residue
              * @param description The description attribute of the current object
              */
            void SetDescription(std::string description);
            /*! \fn
              * A mutator function in order to set the id of the current object
              * Set the id attribute of the current residue
              * @param id The identification attribute of the current object
              */
            void SetId(std::string id);

            /*! \fn                                                                                          //Added by ayush on 11/20/17 for residuenode in assembly
              * A mutator function in order to set the node of the current object
              * Set the node_ attribute of the current residue
              * @param node The node attribute of the current object
              */
            void SetNode(ResidueNode* node);

            /*! \fn
              * A mutator function that replaces the coordinates of the atoms of the current object
              * Replace the coordinate attribute for atoms of the current residue
              * @param atoms The atom attribute of the current object
              */
            void ReplaceAtomCoordinates(AtomVector *newAtoms);

            void SetIsSugarDerivative(bool is_derivative);

            void SetIsAglycon(bool is_aglycon);
            
            
            void SetIsSugar(bool is_sugar);
/** @}*/
            //////////////////////////////////////////////////////////
            //                       FUNCTIONS                      //
            //////////////////////////////////////////////////////////
	    
            void BuildResidueFromPrepFileResidue(PrepFileSpace::PrepFileResidue *prep_residue);
            /// Check if all atoms in the residue have their element symbols --> Label directly (1st priority)
            bool CheckSymbolBasedElementLabeling();
            /// Check if all atoms in the residue have their atom type --> Element symbols come from parameter file (2nd priority)
            bool CheckParameterBasedElementLabeling();
            bool GraphElementLabeling();
            bool GraphSymbolBasedElementLabeling();
            bool GraphParameterBasedElementLabeling();
            bool GraphPredictionBasedElementLabeling();
            AtomVector GetAtomsWithLowestIntraDegree();
            double CalculateAtomicOverlaps(Assembly *assemblyB);
            double CalculateAtomicOverlaps(AtomVector assemblyBAtoms);
            bool CheckIfProtein();
            bool CheckIfWater();
            //GeometryTopology::Coordinate GetRingCenter(); disabled by OG. GetIsRing returns true for all atoms even when IsRing wasn't set.
            GeometryTopology::Coordinate GetGeometricCenter();
            Atom* GetAtom(std::string query_name);
            Atom* GetAtom(unsigned long long query_index);
            Atom* GetAtomWithId(std::string query_id);

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the md_atom contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

            void PrettyPrintHet(std::ostream& out = std::cerr);
            void PrintHetResidues(std::ostream& out = std::cerr);
            void PrintHetAtoms(std::ostream& out = std::cerr);

            void WriteHetResidues(std::ofstream& out);
            void WriteHetAtoms(std::ofstream& out);

            //////////////////////////////////////////////////////////
            //                   OVERLOADED OPERATORS               //
            //////////////////////////////////////////////////////////
            bool operator== (const Residue &otherResidue);
            bool operator!= (const Residue &otherResidue);

        private:

            unsigned long long generateIndex();

            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            Assembly* assembly_;                /*!< Pointer back to the assembly that the current residue belongs to >*/
            std::string name_;                  /*!< Name of residue >*/
            AtomVector atoms_;                  /*!< List of atoms building the residue >*/
            AtomVector head_atoms_;             /*!< List of head atoms in the residue >*/
            AtomVector tail_atoms_;             /*!< List of tail atoms in the residue >*/
            std::string chemical_type_;         /*!< A descriptor in order to describe chemical type of the residue >*/
            std::string description_;           /*!< A short description of the residue >*/
            std::string id_;                    /*!< An identifier for a residue that is generated based on the type of the given file from which the structure has to be built >*/
            unsigned long long index_;         /*!< A unqiue index for each residue in an assembly >*/

	    //Added by Yao 06/13/2018
	    bool is_sugar_derivative_ = false;
	    bool is_aglycon_ = false;
            ResidueNode* node_;                 /*!< A Pointer to a node of the graph structure that indicates this residue >*/              //Added by ayush on 11/20/17 for residuenode in assembly
            
      //Added by Dave 2/1/19
      bool is_sugar_ = false;
    };
 }


#endif // RESIDUE_HPP
