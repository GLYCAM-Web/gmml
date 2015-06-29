#ifndef MONOSACCHARIDE_HPP
#define MONOSACCHARIDE_HPP

#include <string>
#include <vector>
#include <map>

#include "../MolecularModeling/atom.hpp"
#include "chemicalcode.hpp"
#include "sugarname.hpp"

namespace Glycam
{
    struct Monosaccharide {
            int mono_id;
            std::vector<std::vector<MolecularModeling::Atom*> > side_atoms_;
            std::vector<MolecularModeling::Atom*> cycle_atoms_;
            ChemicalCode* chemical_code_;
            SugarName sugar_name_;
            std::map<std::string, std::string> derivatives_map_;
            std::string cycle_atoms_str_;
            std::string anomeric_status_;

            Monosaccharide() {}
    } ;
}

#endif // MONOSACCHARIDE_HPP
