#ifndef GLYCOSIDICLINKAGE_HPP
#define GLYCOSIDICLINKAGE_HPP

#include <vector>
#include <string>
#include <sstream>

#include "../MolecularModeling/assembly.hpp"
#include "../utils.hpp"
#include "../MolecularModeling/residue.hpp"

namespace Glycan
{
  class Monosaccharide;
  class Oligosaccharide;
  class GlycosidicLinkage
  {
  public:
    Monosaccharide* reducing_mono_ = NULL;
    Monosaccharide* non_reducing_mono_ = NULL;
    Monosaccharide* non_reducing_mono_2_ = NULL; //In case of anomeric-anomeric linkage
    bool anomeric_anomeric_linkage_ = false;
    std::string linkage_type_; //ie "1-4"
    std::string inverse_linkage_type_; //ie "4-1" needed for writing as you travel in reverse
    double phi_angle_;
    double psi_angle_;
    double omega_angle_ = -9999;
    MolecularModeling::Atom* reducing_mono_carbon_ = NULL;
    MolecularModeling::Atom* non_reducing_mono_carbon_ = NULL;
    MolecularModeling::Atom* non_reducing_mono_2_carbon_ = NULL;
    MolecularModeling::Atom* glycosidic_oxygen_ = NULL;
    
    GlycosidicLinkage(Monosaccharide* sourceMono, Monosaccharide* targetMono, std::string source_carbon_ID, std::string target_carbon_ID);
    //FUNCTIONS
    double CalculatePhiAngle();
    double CalculatePsiAngle();
    double CalculateOmegaAngle();
    
  };
}

#endif // GLYCOSIDICLINKAGE_HPP
