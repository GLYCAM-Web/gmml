#ifndef INCLUDES_CENTRALDATASTRUCTURE_RESIDUE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_RESIDUE_HPP

#include "includes/CentralDataStructure/atom.hpp"
#include "includes/CentralDataStructure/coordinate.hpp"
#include "includes/CentralDataStructure/Readers/Pdb/pdbResidueId.hpp" // getId
#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"
#include "includes/CodeUtils/constants.hpp" // iNotSet

#include <vector>
#include <memory>    // unique_ptr
#include <algorithm> //swap

using cds::Coordinate;

namespace cds
{
    enum ResidueType
    {
        Protein,
        Sugar,
        Aglycone,
        Derivative,
        Solvent,
        Deoxy,
        Undefined
    };

    class Residue : public glygraph::Node<Residue>
    {
      public:
        //////////////////////////////////////////////////////////
        //                    CONSTRUCTOR                       //
        //////////////////////////////////////////////////////////
        Residue()
        {} //{std::cout << "Residue default ctor\n";}

        Residue(const std::string& residueName, const Residue* referenceResidue);
        Residue(Residue&& other) noexcept; // Move Ctor
        Residue(const Residue& other);     // Copy Ctor
        Residue& operator=(Residue other); // Move & Copy assignment operator.

        virtual ~Residue()
        {} // std::cout << "Residue dtor for " << this->getName() << std::endl;}// Dtor, virtual so that derived classes
           // dtors will get triggered if possible.

        //////////////////////////////////////////////////////////
        //                    ACCESSOR                          //
        //////////////////////////////////////////////////////////
        inline virtual const std::string& getName() const
        {
            return name_;
        }

        const std::string GetParmName() const;
        std::vector<Atom*> getAtoms() const;
        std::vector<std::string> getAtomNames() const;
        std::string getStringId(std::string moleculeNumber = constants::sNotSet) const;
        std::vector<Coordinate*> getCoordinates() const;
        const Coordinate* getGeometricCenter();

        inline ResidueType GetType() const
        {
            return type_;
        }

        inline unsigned int getNumber() const
        {
            return number_;
        }

        virtual pdb::ResidueId getId() const;

        //////////////////////////////////////////////////////////
        //                    MUTATOR                           //
        //////////////////////////////////////////////////////////
        inline void setName(const std::string& s)
        {
            name_ = s;
        }

        inline void setAtoms(std::vector<std::unique_ptr<Atom>> v)
        {
            atoms_ = std::move(v);
        }

        Atom* addAtom(std::unique_ptr<Atom> myAtom);
        // void addAtom(Atom* myAtom);
        bool deleteAtom(const Atom* atom);

        std::vector<std::unique_ptr<Atom>> extractAtoms()
        {
            return std::move(atoms_);
        }

        inline void SetType(ResidueType type)
        {
            type_ = type;
        }

        inline void setNumber(unsigned int i)
        {
            number_ = i;
        }

        //////////////////////////////////////////////////////////
        //                    FUNCTIONS                         //
        //////////////////////////////////////////////////////////
        typename std::vector<std::unique_ptr<Atom>>::iterator FindPositionOfAtom(const Atom* queryAtom);
        Atom* FindAtom(const std::string queryName) const;
        Atom* FindAtom(const int& queryNumber) const;
        bool contains(const Atom* queryAtom) const;
        std::vector<const Atom*> getAtomsConnectedToOtherResidues() const;
        void getAtomPairsConnectedToOtherResidues(std::vector<std::pair<const Atom*, const Atom*>>& foundAtoms) const;
        //        std::vector<std::pair<const Atom*, const Atom*>> getAtomPairsConnectedToOtherResidues() const;
        void MakeDeoxy(std::string oxygenNumber);
        void RemoveHydroxyHydrogen(std::string hydrogenNumber);
        const Coordinate* calculateGeometricCenter();
        ResidueType determineType(const std::string& residueName);
        //////////////////////////////////////////////////////////
        //                    DISPLAY                           //
        //////////////////////////////////////////////////////////
        virtual void Print(std::ostream& out) const;

        //////////////////////////////////////////////////////////
        //                  OPERATOR OVERLOADING                //
        //////////////////////////////////////////////////////////
        virtual bool operator==(const Residue& rhs) const
        {
            return (this->getIndex() == rhs.getIndex());
        }

        virtual bool operator!=(const Residue& rhs) const
        {
            return (this->getIndex() != rhs.getIndex());
        }

        virtual bool operator>(const Residue& rhs) const
        {
            return (this->getNumber() > rhs.getNumber());
        }

        virtual bool operator<(const Residue& rhs) const
        {
            return (this->getNumber() < rhs.getNumber());
        }

        //////////////////////////////////////////////////////////
        //               Copy-Swap                    //
        //////////////////////////////////////////////////////////
        friend void swap(Residue& lhs, Residue& rhs)
        {
            using std::swap;
            swap(static_cast<glygraph::Node<cds::Residue>&>(lhs), static_cast<glygraph::Node<cds::Residue>&>(rhs));
            swap(lhs.atoms_, rhs.atoms_);
            swap(lhs.name_, rhs.name_);
            swap(lhs.geometricCenter_, rhs.geometricCenter_);
            swap(lhs.type_, rhs.type_);
            swap(lhs.number_, rhs.number_);
        }

      private:
        //////////////////////////////////////////////////////////
        //                    ATTRIBUTES                        //
        //////////////////////////////////////////////////////////
        std::vector<std::unique_ptr<Atom>> atoms_;
        std::string name_ = "   ";
        Coordinate geometricCenter_;
        ResidueType type_ = Undefined; // enum Type. See enum above.
        unsigned int number_ =
            1; // constants::iNotSet; ToDo: For prep residues a default 1 value is good. Is there a reason not to?
    };
} // namespace cds
#endif // INCLUDES_CENTRALDATASTRUCTURE_RESIDUE_HPP
