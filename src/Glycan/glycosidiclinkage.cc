//Created 12/11/18
//Dave Montgomery

#include <vector>
#include <string>
#include <sstream>

#include "../../includes/MolecularModeling/assembly.hpp"
#include "../../includes/utils.hpp"
#include "../../includes/MolecularModeling/residue.hpp"
#include "../../includes/Glycan/oligosaccharide.hpp"
#include "../../includes/Glycan/monosaccharide.hpp"
#include "../../includes/Glycan/glycosidiclinkage.hpp"

using Glycan::Monosaccharide;
using Glycan::GlycosidicLinkage;
//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
GlycosidicLinkage::GlycosidicLinkage(Monosaccharide* sourceMono, Monosaccharide* targetMono, std::string source_carbon_ID, std::string target_carbon_ID)
{
  int local_debug = -1;
  if(local_debug > 0)
  {
    gmml::log(__LINE__, __FILE__, gmml::INF, sourceMono->anomeric_carbon_pointer_->GetId());
    gmml::log(__LINE__, __FILE__, gmml::INF, source_carbon_ID);
    gmml::log(__LINE__, __FILE__, gmml::INF, targetMono->anomeric_carbon_pointer_->GetId());
    gmml::log(__LINE__, __FILE__, gmml::INF, target_carbon_ID);
  }
  reducing_mono_ = NULL;
  non_reducing_mono_ = NULL;
  non_reducing_mono_2_ = NULL;
  anomeric_anomeric_linkage_ = false;


// This didn't work for all sugars, adding a check for carbon number; 1 is anomeric unless there is no 1, then 2 is.
  for(std::vector<MolecularModeling::Atom*>::iterator it = sourceMono->cycle_atoms_.begin(); it != sourceMono->cycle_atoms_.end(); it++)
  {
    MolecularModeling::Atom* thisAtom = *it;
    if(thisAtom->GetId() == source_carbon_ID)
    {
      if(thisAtom == sourceMono->anomeric_carbon_pointer_)
      {
        non_reducing_mono_ = sourceMono;
        non_reducing_mono_carbon_ = thisAtom;
      }
    }
  }
  for(std::vector<MolecularModeling::Atom*>::iterator it = targetMono->cycle_atoms_.begin(); it != targetMono->cycle_atoms_.end(); it++)
  {
    MolecularModeling::Atom* thisAtom = *it;
    if(thisAtom->GetId() == target_carbon_ID)
    {
      if(thisAtom == targetMono->anomeric_carbon_pointer_)
      {
        if(non_reducing_mono_ == NULL)
        {
          non_reducing_mono_ = targetMono;
          non_reducing_mono_carbon_ = thisAtom;
        }
        else
        {
          non_reducing_mono_2_ = targetMono;
          anomeric_anomeric_linkage_ = true;
          non_reducing_mono_2_carbon_ = thisAtom;
        }
      }
    }
  }

  if(!anomeric_anomeric_linkage_)
  {
    if(non_reducing_mono_ == sourceMono)
    {
      reducing_mono_ = targetMono;
    }
    else if(non_reducing_mono_ == targetMono)
    {
      reducing_mono_ = sourceMono;
    }
  }
  if(non_reducing_mono_ == NULL)
  {
    gmml::log(__LINE__, __FILE__, gmml::ERR, "Unable to find non_reducing_mono_");
  }
  if((reducing_mono_ == NULL) && (non_reducing_mono_2_ == NULL))
  {
    gmml::log(__LINE__, __FILE__, gmml::ERR, "Unable to find reducing_mono_ or 2nd non_reducing_mono_");
  }
  std::stringstream ss;
  ss << source_carbon_ID[1]  << "-" << target_carbon_ID[1];
  linkage_type_ = ss.str();
  inverse_linkage_type_ = linkage_type_;
  std::reverse(inverse_linkage_type_.begin(), inverse_linkage_type_.end());
  if(local_debug > 0)
  {
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
  }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
