#ifndef PDBMODELRESIDUESET_HPP
#define PDBMODELRESIDUESET_HPP

#include <string>
#include <vector>

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

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            AtomCardVector GetAtoms();
            HeterogenAtomCardVector GetHeterogenAtoms();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetAtoms(const  AtomCardVector atoms);
            void SetHeterogenAtoms(const HeterogenAtomCardVector heterogen_atoms);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            AtomCardVector atoms_;
            HeterogenAtomCardVector heterogen_atoms_;

    };
}

#endif // PDBMODELRESIDUESET_HPP
