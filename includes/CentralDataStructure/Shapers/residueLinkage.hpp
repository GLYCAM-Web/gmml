#ifndef GMML_INCLUDES_CENTRAL_DATA_STRUCTURE_SHAPERS_RESIDUE_LINKAGE_HPP
#define GMML_INCLUDES_CENTRAL_DATA_STRUCTURE_SHAPERS_RESIDUE_LINKAGE_HPP
/*
 * This class figures out the rotatable bonds between two residues
 * Starts/ends at the CA atoms in proteins. Looks for cycles (as they aren't rotatable).
 * Stores each rotatable bond as a RotatableDihedral object.
 */
#include "includes/CentralDataStructure/residue.hpp"
#include "includes/CentralDataStructure/Shapers/rotatableDihedral.hpp"
#include "includes/CodeUtils/logging.hpp"
#include <vector>

namespace cds
{
    class ResidueLinkage
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        ResidueLinkage() {};
        ResidueLinkage(cds::Residue* nonReducingResidue1, cds::Residue* reducingResidue2,
                       bool reverseAtomsThatMove = true);
        ResidueLinkage(cds::Residue* nonReducingResidue1, cds::Residue* reducingResidue2,
                       std::vector<cds::Atom*> alsoMovingAtoms, bool reverseAtomsThatMove = true);
        //~ResidueLinkage() {std::cout << "Linkage dtor for " << this->GetFromThisResidue1()->getId() << " -Link- "  <<
        // this->GetToThisResidue2()->getId() << "\n";}
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        std::vector<RotatableDihedral> GetRotatableDihedrals() const;
        std::vector<RotatableDihedral> GetRotatableDihedralsWithMultipleRotamers() const;
        unsigned long int GetNumberOfRotatableDihedrals() const;
        int GetNumberOfShapes(const bool likelyShapesOnly = false) const;

        inline cds::Residue* GetFromThisResidue1() const
        {
            return from_this_residue1_;
        }

        inline cds::Residue* GetToThisResidue2() const
        {
            return to_this_residue2_;
        }

        bool CheckIfConformer() const;

        inline bool GetIfExtraAtoms() const
        {
            return isExtraAtoms_;
        }

        inline std::vector<cds::Atom*> GetExtraAtoms()
        {
            return extraAtomsThatMove_;
        }

        inline unsigned long long GetIndex() const
        {
            return index_;
        }

        std::string GetName() const;

        std::vector<cds::Residue*>& GetMovingResidues();
        std::vector<cds::Residue*>& GetFixedResidues();

        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        inline void SetRotatableDihedrals(std::vector<RotatableDihedral> rotatableDihedrals)
        {
            rotatableDihedrals_ = rotatableDihedrals;
        }

        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        void SetDefaultShapeUsingMetadata();
        void SetRandomShapeUsingMetadata(bool useRanges = true);
        void SetSpecificShapeUsingMetadata(int shapeNumber, bool useRanges = false);
        void SetSpecificShape(std::string dihedralName, std::string selectedRotamer);
        void SetCustomDihedralAngles(std::vector<double> dihedral_angles);
        void SetShapeToPrevious();
        void SetRandomDihedralAngles();
        void DetermineAtomsThatMove();
        void SimpleWiggle(std::vector<cds::Atom*>& overlapAtomSet1, std::vector<cds::Atom*>& overlapAtomSet2,
                          const int angleIncrement = 5);
        void SimpleWiggleCurrentRotamers(std::vector<cds::Atom*>& overlapAtomSet1,
                                         std::vector<cds::Atom*>& overlapAtomSet2, const int angleIncrement = 5);
        void SimpleWiggleCurrentRotamers(std::vector<cds::Residue*>& overlapSet1,
                                         std::vector<cds::Residue*>& overlapSet2, const int angleIncrement = 5);

        inline void SetIndex(unsigned long long index)
        {
            index_ = index;
        }

        //////////////////////////////////////////////////////////
        //                       DISPLAY FUNCTION               //
        //////////////////////////////////////////////////////////
        std::string Print() const;

        //////////////////////////////////////////////////////////
        //                  OPERATOR OVERLOADING                //
        //////////////////////////////////////////////////////////
        bool operator==(const ResidueLinkage& rhs) const
        {
            return (this->GetIndex() == rhs.GetIndex());
        }

        bool operator!=(const ResidueLinkage& rhs) const
        {
            return (this->GetIndex() != rhs.GetIndex());
        }

      private:
        //////////////////////////////////////////////////////////
        //                    PRIVATE FUNCTIONS                 //
        //////////////////////////////////////////////////////////
        std::vector<cds::Residue*> GetResidues() const;

        inline bool GetIfReversedAtomsThatMove() const
        {
            return reverseAtomsThatMove_;
        }

        inline cds::Atom* GetFromThisConnectionAtom1() const
        {
            return from_this_connection_atom1_;
        }

        inline cds::Atom* GetToThisConnectionAtom2() const
        {
            return to_this_connection_atom2_;
        }

        inline void SetIfReversedAtomsThatMove(bool b)
        {
            reverseAtomsThatMove_ = b;
        }

        void AddExtraAtomsThatMove(std::vector<cds::Atom*> extraAtoms);
        void InitializeClass(cds::Residue* from_this_residue1, cds::Residue* to_this_residue2,
                             bool reverseAtomsThatMove);
        bool CheckIfViableLinkage() const;
        std::vector<RotatableDihedral> FindRotatableDihedralsConnectingResidues(cds::Atom* from_this_connection_atom1,
                                                                                cds::Atom* to_this_connection_atom2);
        // std::vector<cds::Atom*> DealWithBranchesFromLinkages(std::vector<cds::Atom*> linearLinkageAtoms, Atom
        // *cycle_point1, Atom *cycle_point2);
        //  Previous function generates a list of linearly connected atoms that define the rotatable bonds
        //  This function splits that list into groups of 4 and creates RotatableDihedral objects
        std::vector<RotatableDihedral> SplitAtomVectorIntoRotatableDihedrals(std::vector<cds::Atom*> atoms);
        gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector FindMetadata() const;
        void AddMetadataToRotatableDihedrals(gmml::MolecularMetadata::GLYCAM::DihedralAngleDataVector metadata);
        void SetResidues(cds::Residue* residue1, cds::Residue* residue2);
        void SetConnectionAtoms(cds::Residue* residue1, cds::Residue* residue2);
        void SetConformerUsingMetadata(bool useRanges = false, int conformerNumber = 0);
        unsigned long long GenerateIndex();
        std::string DetermineLinkageNameFromResidueNames() const;
        void DetermineMovingResidues();

        inline void SetName(std::string name)
        {
            name_ = name;
        }

        //////////////////////////////////////////////////////////
        //                       ATTRIBUTES                     //
        //////////////////////////////////////////////////////////
        cds::Residue* from_this_residue1_      = nullptr;
        cds::Residue* to_this_residue2_        = nullptr;
        cds::Atom* from_this_connection_atom1_ = nullptr;
        cds::Atom* to_this_connection_atom2_   = nullptr;
        std::vector<RotatableDihedral> rotatableDihedrals_;
        bool reverseAtomsThatMove_ = true;
        std::vector<cds::Atom*> extraAtomsThatMove_;
        bool isExtraAtoms_        = true;
        unsigned long long index_ = 0;
        std::string name_         = "";             // e.g. "DGalpb1-6DGlcpNAc". It being empty works with GetName();
        std::vector<cds::Residue*> movingResidues_; // overlap speedups
        std::vector<cds::Residue*> fixedResidues_;  // overlap speedups
    };
} // namespace cds
#endif // GMML_INCLUDES_GEOMETRYTOPOLOGY_RESIDUELINKAGES_RESIDUE_LINKAGE_HPP
