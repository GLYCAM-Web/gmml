#ifndef LIBRARYFILERESIDUE_HPP
#define LIBRARYFILERESIDUE_HPP

#include <string>
#include <vector>
#include <map>

namespace LibraryFileSpace
{
    class LibraryFileAtom;
    class LibraryFileResidue
    {
        public:
            ///////////////////////////////// TYPE DEFINITION ///////////////////////////////////////
            typedef std::map<int, LibraryFileAtom*> AtomMap;

            /////////////////////////////// CONSTRUCTOR ////////////////////////////////
            LibraryFileResidue();
            LibraryFileResidue(std::string& name);
            LibraryFileResidue(std::string& name, std::vector<LibraryFileAtom*>& atoms, int head_atom_index, int tail_atom_index, double box_angle = 0,
                               double box_length = 0, double box_width = 0, double box_height = 0);

            //////////////////////////////// ACCESSOR //////////////////////////////////
            std::string GetName();
            AtomMap GetAtoms();
            LibraryFileAtom* GetAtomByIndex(int index);
            LibraryFileAtom* GetAtomByOrder(int order);
            double GetBoxAngle();
            double GetBoxLength();
            double GetBoxWidth();
            double GetBoxHeight();
            int GetHeadAtomIndex();
            int GetTailAtomIndex();

            ///////////////////////////////// MUTATOR //////////////////////////////////
            void SetName(std::string& name);
            void SetAtoms(AtomMap& atoms);
            void AddAtom(LibraryFileAtom* atom);
            void SetBoxAngle(double box_angle);
            void SetBoxLength(double box_length);
            void SetBoxWidth(double box_width);
            void SetBoxHeight(double box_height);
            void SetHeadAtomIndex(int head_atom_index);
            void SetTailAtomIndex(int tail_atom_index);

            //////////////////////////// DISPLAY FUNCTION //////////////////////////////
            void Print(std::ostream& out);

        private:
            /////////////////////////////// ATTRIBUTES /////////////////////////////////
            std::string name_;          // Residue nam; Set by each line of the index section of a library file
            AtomMap atoms_;             // List of atoms belonging to the residue mapped to their order numbers; set by various sections (atom, connectivity, positions) of a library file
            double box_angle_;          // Box angle; set by the 2nd line of the boundbox section of a library file
            double box_length_;         // Box length; set by the 3rd line of the boundbox section of a library file
            double box_width_;          // Box width; set by the 4th line of the boundbox section of a library file
            double box_height_;         // Box height; set by the 5th line of the boundbox section of a library file
            int head_atom_index_;       // Head atom index; set by the 1st line of the connect section of a library file
            int tail_atom_index_;       // Tail atom index; set by the 2nd line of the connect section of a library file
    };
}

#endif // LIBRARYFILERESIDUE_HPP
