#ifndef PDBMODELRESIDUESET_HPP
#define PDBMODELRESIDUESET_HPP

#include <string>
#include <vector>
#include <sstream>
#include <iostream>

namespace PdbFileSpace
{
    class PdbAtomCard;
    class PdbHeterogenAtomCard;

    class PdbModelResidueSet
    {
        public:
            //////////////////////////////////////////////////////////
            //                       TYPE DEFINITION                //
            //////////////////////////////////////////////////////////
            typedef std::vector<PdbAtomCard*> AtomCardVector;
            typedef std::vector<PdbHeterogenAtomCard*> HeterogenAtomCardVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            PdbModelResidueSet();
            PdbModelResidueSet(std::stringstream& residue_set_block);

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            AtomCardVector GetAtoms();
            HeterogenAtomCardVector GetHeterogenAtoms();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetAtoms(const  AtomCardVector atoms);
            void AddAtom(PdbAtomCard* atom);
            void SetHeterogenAtoms(const HeterogenAtomCardVector heterogen_atoms);
            void AddHeterogenAtom(PdbHeterogenAtomCard* heterogen_atom);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            AtomCardVector atoms_;
            HeterogenAtomCardVector heterogen_atoms_;

    };
}

#endif // PDBMODELRESIDUESET_HPP
