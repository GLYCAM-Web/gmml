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
  int local_debug = 0;
  if(local_debug > 0)
  {
    if(sourceMono->anomeric_carbon_pointer_ != NULL)
      gmml::log(__LINE__, __FILE__, gmml::INF, "Source Mono Anomeric Carbon: " + sourceMono->anomeric_carbon_pointer_->GetId());
    gmml::log(__LINE__, __FILE__, gmml::INF, "Source Carbon ID " + source_carbon_ID);
    if(targetMono->anomeric_carbon_pointer_ != NULL)
      gmml::log(__LINE__, __FILE__, gmml::INF, "Target Mono Anomeric Carbon: " + targetMono->anomeric_carbon_pointer_->GetId());
    gmml::log(__LINE__, __FILE__, gmml::INF,  "Target Carbon ID " +target_carbon_ID);
  }
  reducing_mono_ = NULL;
  non_reducing_mono_ = NULL;
  non_reducing_mono_2_ = NULL;
  anomeric_anomeric_linkage_ = false;


// This didn't work for all sugars, adding a check for carbon number; 1 is anomeric unless there is no 1, then 2 is.
  for(std::vector<MolecularModeling::Atom*>::iterator it = sourceMono->cycle_atoms_.begin(); it != sourceMono->cycle_atoms_.end(); it++)
  {
    MolecularModeling::Atom* thisAtom = *it;
    if(thisAtom->GetId() == source_carbon_ID && sourceMono->anomeric_carbon_pointer_ != NULL)
    {
      if(thisAtom->GetId() == sourceMono->anomeric_carbon_pointer_->GetId())
      {
          if(local_debug > 0)
          {
            gmml::log(__LINE__, __FILE__, gmml::INF, "Source Mono Anomeric Carbon in Linkage");
          }
        non_reducing_mono_ = sourceMono;
        non_reducing_mono_carbon_ = thisAtom;
      }
    }
  }
  for(std::vector<MolecularModeling::Atom*>::iterator it = targetMono->cycle_atoms_.begin(); it != targetMono->cycle_atoms_.end(); it++)
  {
    MolecularModeling::Atom* thisAtom = *it;
    if(thisAtom->GetId() == target_carbon_ID && targetMono->anomeric_carbon_pointer_ != NULL)
    {
      if(thisAtom->GetId() == targetMono->anomeric_carbon_pointer_->GetId())
      {
        if(local_debug > 0)
        {
          gmml::log(__LINE__, __FILE__, gmml::INF, "Target Mono Anomeric Carbon in Linkage");
        }
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
  
  phi_angle_ = CalculatePhiAngle();
  psi_angle_ = CalculatePsiAngle();
  if(linkage_type_ == "1-6" || inverse_linkage_type_ == "1-6")
  {
    omega_angle_ = CalculateOmegaAngle();
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

//TODO make these functions work

double Glycan::GlycosidicLinkage::CalculatePhiAngle()
{
  int local_debug = -1;
  if(local_debug > 0)
  {
    gmml::log(__LINE__, __FILE__, gmml::INF, "CalculatingPhiAngle");
  }
  //(O5-C1-O-Cx') {Ring oxygen of child_oligo}-{child_atom_id}-{glycosidic_atom_id}-{parent_atom_id}

  MolecularModeling::Atom* O5 = NULL;
  MolecularModeling::Atom* C1 = non_reducing_mono_carbon_;
  MolecularModeling::Atom* glycosidicO = NULL;
  MolecularModeling::Atom* Cx = NULL;

  //Get Anomeric Carbon and ring oxygen
  if(non_reducing_mono_ != NULL)
  {
    for(std::vector<MolecularModeling::Atom*>::iterator it = non_reducing_mono_->cycle_atoms_.begin();
        it != non_reducing_mono_->cycle_atoms_.end(); ++it)
    {
      MolecularModeling::Atom* this_cycle_atom = (*it);
      std::string this_element = this_cycle_atom->GetElementSymbol();
      std::string this_atom_id = this_cycle_atom->GetId();

      if (this_atom_id.find("O") != std::string::npos)
      {
        O5 = this_cycle_atom;
      }
    }
  }
  if((O5 == NULL)||(C1 == NULL))
  {
    gmml::log(__LINE__, __FILE__, gmml::ERR, "O5 or C1 is null in calculate Phi Angle");
    return -9999;
  }
  if(local_debug > 0)
  {
    gmml::log(__LINE__, __FILE__, gmml::INF, "O5");
    gmml::log(__LINE__, __FILE__, gmml::INF, O5->GetId());
    gmml::log(__LINE__, __FILE__, gmml::INF, "C1");
    gmml::log(__LINE__, __FILE__, gmml::INF, C1->GetId());
  }
  //Find glycosidicO
  std::vector<MolecularModeling::Atom*> C1_neighbors = C1->GetNode()->GetNodeNeighbors();
  for (std::vector<MolecularModeling::Atom*>::iterator it = C1_neighbors.begin(); it != C1_neighbors.end(); ++it)
  {
    MolecularModeling::Atom* this_neighbor = (*it);
    std::string this_atom_id = this_neighbor->GetId();
    if (this_atom_id.find("O") != std::string::npos && !this_neighbor->GetIsCycle())
    {
      glycosidicO = this_neighbor;
    }
  }
  if(glycosidicO == NULL)
  {
    // gmml::log(__LINE__, __FILE__, gmml::INF, glycosidicO->GetId());
    gmml::log(__LINE__, __FILE__, gmml::ERR, "glycosidicO is null in calculate Phi Angle");
    return -9999;
  }
  else
  {
    glycosidic_oxygen_ = glycosidicO;
  }

  //Get Cx
  std::vector<MolecularModeling::Atom*> glycosidicO_neighbors = glycosidicO->GetNode()->GetNodeNeighbors();
  for (std::vector<MolecularModeling::Atom*>::iterator it = glycosidicO_neighbors.begin(); it != glycosidicO_neighbors.end(); ++it)
  {
    MolecularModeling::Atom* this_neighbor = (*it);
    if (this_neighbor != C1)
    {
      Cx = this_neighbor;
    }
  }
  if(Cx == NULL)
  {
    // gmml::log(__LINE__, __FILE__, gmml::INF, Cx->GetId());
    gmml::log(__LINE__, __FILE__, gmml::ERR, "Cx is null in calculate Phi Angle");
    return -9999;
  }

  if(local_debug > 0)
  {
    gmml::log(__LINE__, __FILE__, gmml::INF, "glycosidicO");
    gmml::log(__LINE__, __FILE__, gmml::INF, glycosidicO->GetId());
    gmml::log(__LINE__, __FILE__, gmml::INF, "Cx");
    gmml::log(__LINE__, __FILE__, gmml::INF, Cx->GetId());
  }
  
  double Phi = non_reducing_mono_->assembly_->CalculateTorsionAngleByAtoms(O5, C1, glycosidicO, Cx);
  if (local_debug > 0)
  {
    gmml::log(__LINE__, __FILE__, gmml::INF, std::to_string(Phi));
  }
  Phi = gmml::ConvertRadian2Degree(Phi);
  return Phi;
}

double Glycan::GlycosidicLinkage::CalculatePsiAngle()
{
  //(C1-O-Cx'-C[x-1]') {child_atom_id}-{glycosidic_atom_id}-{parent_atom_id}-{parent_atom_id - 1}

  MolecularModeling::Atom* C1 = non_reducing_mono_carbon_;
  MolecularModeling::Atom* glycosidicO = NULL;
  MolecularModeling::Atom* Cx = NULL;
  MolecularModeling::Atom* Cx_1 = NULL;

  //Get Anomeric Carbon
  if(C1 == NULL)
  {
    gmml::log(__LINE__, __FILE__, gmml::ERR, "C1 is null in calculate Psi Angle");
    return -9999;
  }

  //Find glycosidicO
  std::vector<MolecularModeling::Atom*> C1_neighbors = C1->GetNode()->GetNodeNeighbors();
  for (std::vector<MolecularModeling::Atom*>::iterator it = C1_neighbors.begin(); it != C1_neighbors.end(); ++it)
  {
    MolecularModeling::Atom* this_neighbor = (*it);
    std::string this_atom_id = this_neighbor->GetId();
    if (this_atom_id.find("O") != std::string::npos && !this_neighbor->GetIsCycle())
    {
      glycosidicO = this_neighbor;
    }
  }
  if(glycosidicO == NULL)
  {
    gmml::log(__LINE__, __FILE__, gmml::ERR, "glycosidicO is null in calculate Psi Angle");
    return -9999;
  }

  //Get Cx
  std::vector<MolecularModeling::Atom*> glycosidicO_neighbors = glycosidicO->GetNode()->GetNodeNeighbors();
  for (std::vector<MolecularModeling::Atom*>::iterator it = glycosidicO_neighbors.begin(); it != glycosidicO_neighbors.end(); ++it)
  {
    MolecularModeling::Atom* this_neighbor = (*it);
    if (this_neighbor != C1)
    {
      Cx = this_neighbor;
    }
  }
  if(Cx == NULL)
  {
    gmml::log(__LINE__, __FILE__, gmml::ERR, "Cx is null in calculate Psi Angle");
    return -9999;
  }

  //Get Cx-1
  std::vector<MolecularModeling::Atom*> Cx_neighbors = Cx->GetNode()->GetNodeNeighbors();
  for (std::vector<MolecularModeling::Atom*>::iterator it = Cx_neighbors.begin(); it != Cx_neighbors.end(); ++it)
  {
    MolecularModeling::Atom* this_neighbor = (*it);
    std::string this_atom_id = this_neighbor->GetId();
    int Cx_atom_number = Cx->GetId().at(1) - '0';
    int this_atom_number = this_atom_id.at(1) - '0';
    if ((Cx_atom_number != 1) && ((Cx_atom_number - 1) == this_atom_number))
    {
      Cx_1 = this_neighbor;
    }
    else if((Cx_atom_number == 1) && ((Cx_atom_number + 1) == this_atom_number))
    {
      Cx_1 = this_neighbor;
    }
  }
  if(Cx_1 == NULL)
  {
    // gmml::log(__LINE__, __FILE__, gmml::INF, Cx->GetId());
    gmml::log(__LINE__, __FILE__, gmml::ERR, "Cx_1 is null in calculate Psi Angle");
    return -9999;
  }

  double Psi = non_reducing_mono_->assembly_->CalculateTorsionAngleByAtoms(C1, glycosidicO, Cx, Cx_1);
  Psi = gmml::ConvertRadian2Degree(Psi);
  return Psi;
}

double Glycan::GlycosidicLinkage::CalculateOmegaAngle()
{
  int local_debug = -1;
  //(O-C6'-C5'-O5') {glycosidic_atom_id}-{parent_atom_id}-{Carbon 5 in parent oligo}-{Ring oxygen in parent_oligo}

  MolecularModeling::Atom* glycosidicO = NULL;
  MolecularModeling::Atom* C6prime = NULL;
  MolecularModeling::Atom* C5prime = NULL;
  MolecularModeling::Atom* O5prime = NULL;

  //Get C5' and O5' first
  if (local_debug > 0)
  {
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Get C5' and O5' first");
  }
  for(std::vector<MolecularModeling::Atom*>::iterator it = reducing_mono_->cycle_atoms_.begin();
      it != reducing_mono_->cycle_atoms_.end(); ++it)
  {
    MolecularModeling::Atom* this_cycle_atom = (*it);
    std::string this_atom_id = this_cycle_atom->GetId();
    int this_atom_number = this_atom_id.at(1) - '0';
    std::string this_element = this_cycle_atom->GetElementSymbol();
    if (this_atom_number == 5)
    {
      if (this_atom_id.find("C") != std::string::npos)
      {
        C5prime = this_cycle_atom;
      }
    }
    if (this_atom_id.find("O") != std::string::npos)
    {
      O5prime = this_cycle_atom;
    }
  }
  if ((C5prime == NULL)||(O5prime == NULL))
  {
    gmml::log(__LINE__, __FILE__, gmml::ERR, "O5prime or C5prime is null in calculate Omega Angle");
    return -9999;
  }

  //Get C6'
  if (local_debug > 0)
  {
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Get C6'");
  }
  std::vector<MolecularModeling::Atom*> C5prime_neighbors = C5prime->GetNode()->GetNodeNeighbors();
  for (std::vector<MolecularModeling::Atom*>::iterator it = C5prime_neighbors.begin(); it != C5prime_neighbors.end(); ++it)
  {
    MolecularModeling::Atom* this_neighbor = (*it);
    if ( this_neighbor == reducing_mono_carbon_)
    {
      C6prime = this_neighbor;
    }
  }
  if(C6prime == NULL)
  {
    gmml::log(__LINE__, __FILE__, gmml::ERR, "C6prime is null in calculate Omega Angle");
    return -9999;
  }
  //Get glycosidic oxygen
  if (local_debug > 0)
  {
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Get glycosidic oxygen");
  }
  std::vector<MolecularModeling::Atom*> C6prime_neighbors = C6prime->GetNode()->GetNodeNeighbors();
  for (std::vector<MolecularModeling::Atom*>::iterator it = C6prime_neighbors.begin(); it != C6prime_neighbors.end(); ++it)
  {
    MolecularModeling::Atom* this_neighbor = (*it);
    std::string this_atom_id = this_neighbor->GetId();
    if (this_atom_id.find("O") != std::string::npos && !this_neighbor->GetIsCycle())
    {
      glycosidicO = this_neighbor;
    }
  }
  if(glycosidicO == NULL)
  {
    gmml::log(__LINE__, __FILE__, gmml::ERR, "glycosidicO is null in calculate Omega Angle");
    return -9999;
  }
  if (local_debug > 0)
  {
    gmml::log(__LINE__, __FILE__,  gmml::INF, "About to calcuate torsion");
  }
  double Omega = non_reducing_mono_->assembly_->CalculateTorsionAngleByAtoms(glycosidicO, C6prime, C5prime, O5prime);
  Omega = gmml::ConvertRadian2Degree(Omega);
  return Omega;
}

