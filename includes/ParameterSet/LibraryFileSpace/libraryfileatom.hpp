#ifndef LIBRARYFILEATOM_HPP
#define LIBRARYFILEATOM_HPP

#include <string>
#include <vector>
#include "../../Geometry/coordinate.hpp"
using namespace Geometry;

namespace LibraryFileSpace
{
    class LibraryFileAtom
    {
        public:
            /////////////////////////////// CONSTRUCTOR ////////////////////////////////
            LibraryFileAtom();
            LibraryFileAtom(const std::string& type, const std::string& name, int residue_index, int atom_index, int atomic_number, double charge);
            LibraryFileAtom(const std::string& type, const std::string& name, int residue_index, int atom_index, int atomic_number, double charge,
                            const Coordinate& coordinate, const std::vector<int> bonded_atoms_indices, int atom_order);

            ///////////////////////////////// ACCESSOR /////////////////////////////////
            std::string GetType();
            std::string GetName();
            int GetResidueIndex();
            int GetAtomIndex();
            int GetAtomicNumber();
            double GetCharge();
            Coordinate GetCoordinate();
            std::vector<int> GetBondedAtomsIndices();
            int GetAtomOrder();

            ////////////////////////////////// MUTATOR /////////////////////////////////
            void SetType(const std::string type);
            void SetName(const std::string name);
            void SetResidueIndex(int residue_index);
            void SetAtomIndex(int atom_index);
            void SetAtomicNumber(int atomic_number);
            void SetCoordinate(const Coordinate& coordinate);
            void SetBondedAtomsIndices(const std::vector<int> bonded_atoms_indices);
            void AddBondedAtomIndex(int index);
            void SetAtomOrder(int atom_order);

            //////////////////////////// DISPLAY FUNCTION //////////////////////////////
            void Print(std::ostream& out);


        private:
            /////////////////////////////// ATTRIBUTES /////////////////////////////////
            std::string type_;              // Atom type; Set by the first column of the atom section of a library file
            std::string name_;              // Atom name: Set by the 2nd column of the atom section of a library file
            int residue_index_;             // Residue index that the atom belongs to; Set by the 4th column of the atom section of a library file
                                            // (this number should indicates the index of the residue in the index section of the library file but it could be a wrong number)
            int atom_index_;                // Index of the atom in the belonging residue; Set by the 6th column of the atom section of a library file
            int atomic_number_;             // Atomic number of the atom; Set by the 7th column of the atom section of a library file
            double charge_;                 // Charge of the atom; Set by the 8th column of the atom section of a library file
            Coordinate coordinate_;         // Coordinate of the atom; Set by the corresponding line in the positions section of a library file
            std::vector<int> bonded_atoms_indices_; // List of atom indices that are bonded to the atom; Set by the corresponding atom index in the connectivity section of a library file
            int atom_order_;                // Order of atoms of the residue in the atom section; Set by a line counter iterates on lines of the atom section of a library file
    };
}

#endif // LIBRARYFILEATOM_HPP
