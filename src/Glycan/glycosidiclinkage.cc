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
      std::vector<MolecularModeling::Atom*> targetMonoAtoms = targetMono->cycle_atoms_[0]->GetResidue()->GetAtoms();
      for(std::vector<MolecularModeling::Atom*>::iterator it = targetMonoAtoms.begin(); it != targetMonoAtoms.end(); it++)
      {
        MolecularModeling::Atom* thisMonoAtom = (*it);
        if(thisMonoAtom->GetId() == target_carbon_ID)
        {
          reducing_mono_carbon_ = thisMonoAtom;
        }
      }
    }
    else if(non_reducing_mono_ == targetMono)
    {
      reducing_mono_ = sourceMono;
      std::vector<MolecularModeling::Atom*> sourceMonoAtoms = sourceMono->cycle_atoms_[0]->GetResidue()->GetAtoms();
      for(std::vector<MolecularModeling::Atom*>::iterator it = sourceMonoAtoms.begin(); it != sourceMonoAtoms.end(); it++)
      {
        MolecularModeling::Atom* thisMonoAtom = (*it);
        if(thisMonoAtom->GetId() == source_carbon_ID)
        {
          reducing_mono_carbon_ = thisMonoAtom;
        }
      }
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

  //GET THE ATOMS NEEDED FOR TORSION ANGLES AND THEN PASS THEM TO THE BELOW FUNCTIONS.
  //This vector will have the non_reducing sugar's ring oxygen, then anomeric carbon,
  //then the glycosidic O (Ox of the reducing sugar), Cx of the reducing sugar, Cx-1, and for omega (1 or 2 - 6 linked) Cx-2
  std::vector<MolecularModeling::Atom*> linkageAtoms;
  for(std::vector<MolecularModeling::Atom*>::iterator it = non_reducing_mono_->cycle_atoms_.begin(); it != non_reducing_mono_->cycle_atoms_.end(); it++)
  {
    MolecularModeling::Atom* thisCycleAtom = (*it);
    if(thisCycleAtom->GetElementSymbol() == "O")
    {//Ring Oxygen
      linkageAtoms.push_back(thisCycleAtom);
    }
  }
  linkageAtoms.push_back(non_reducing_mono_carbon_);
  std::vector<MolecularModeling::Atom*> anomericNeighbors = non_reducing_mono_carbon_->GetNode()->GetNodeNeighbors();
  for(std::vector<MolecularModeling::Atom*>::iterator it = anomericNeighbors.begin(); it != anomericNeighbors.end(); it++)
  {
    MolecularModeling::Atom* thisAnomericNeighbor = (*it);
    if(reducing_mono_!= NULL)
    {
      if(reducing_mono_carbon_->GetResidue() == thisAnomericNeighbor->GetResidue())
      {
        if(thisAnomericNeighbor->GetElementSymbol() == "O")
        {//Then this is the glycosidic oxygen
          linkageAtoms.push_back(thisAnomericNeighbor);
          glycosidic_oxygen_ = thisAnomericNeighbor;
        }
        else if (local_debug > 0)
        {//TODO make this handle non oxygen glycosidic atoms
          ss.str("");
          ss << "Glycosidic Oxygen not present. The atom present is " << thisAnomericNeighbor->GetElementSymbol() << ".";
          gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
          ss.str("");
        }
      }
    }
    else if (non_reducing_mono_2_ != NULL)
    { //The oxygen can belong to either sugar
      // if(non_reducing_mono_2_carbon_->GetResidue() == thisAnomericNeighbor->GetResidue())
      // {
      if((thisAnomericNeighbor->GetElementSymbol() == "O") && thisAnomericNeighbor->GetIsRing() == false)
      {//Then this is the glycosidic oxygen
        linkageAtoms.push_back(thisAnomericNeighbor);
        glycosidic_oxygen_ = thisAnomericNeighbor;
      }
      // else if (local_debug > 0)
      // {//TODO make this handle non oxygen glycosidic atoms
      //   ss.str("");
      //   ss << "Glycosidic Oxygen not present. The atom present is " << thisAnomericNeighbor->GetElementSymbol() << ".";
      //   gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
      //   ss.str("");
      // }
      // }
      // else if (local_debug > 0)
      // {
      //   ss.str("");
      //   ss << "Unable to find Glycosidic O.";
      //   gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
      //   ss.str("");
      // }
    }
    else
    {
      gmml::log(__LINE__, __FILE__, gmml::ERR, "Cannot identify both sugars for the linkage, something went horribly wrong!");
    }
  }
  if(!anomeric_anomeric_linkage_)
  {
    linkageAtoms.push_back(reducing_mono_carbon_);
    int x = std::atoi(&reducing_mono_carbon_->GetId()[1]);
    std::vector<MolecularModeling::Atom*> CxNeighbors = reducing_mono_carbon_->GetNode()->GetNodeNeighbors();
    for(std::vector<MolecularModeling::Atom*>::iterator it = CxNeighbors.begin(); it != CxNeighbors.end(); it++)
    {
      MolecularModeling::Atom* thisCxNeighbor = (*it);
      int thisNeighborNum = std::atoi(&thisCxNeighbor->GetId()[1]);
      if((thisNeighborNum = (x - 1)) && (thisCxNeighbor->GetElementSymbol() == "C"))
      {//This is Cx-1
        linkageAtoms.push_back(thisCxNeighbor);
        if((linkage_type_ == "1-6") || (inverse_linkage_type_ == "1-6")||(linkage_type_ == "2-6") || (inverse_linkage_type_ == "2-6"))
        {//Get Cx-2
          std::vector<MolecularModeling::Atom*> CxMinusOneNeighbors = thisCxNeighbor->GetNode()->GetNodeNeighbors();
          for(std::vector<MolecularModeling::Atom*>::iterator it2 = CxMinusOneNeighbors.begin(); it2 != CxMinusOneNeighbors.end(); it2++)
          {
            MolecularModeling::Atom* thisCxMinusOneNeighbor = (*it2);
            int thisCxMinusOneNeighborNum = std::atoi(&thisCxMinusOneNeighbor->GetId()[1]);
            if((thisCxMinusOneNeighborNum = (x - 2)) && (thisCxMinusOneNeighbor->GetElementSymbol() == "C"))
            {//This is Cx-2
              linkageAtoms.push_back(thisCxMinusOneNeighbor);
            }
          }
        }
      }
    }
    phi_angle_ = CalculatePhiAngle(linkageAtoms);
    psi_angle_ = CalculatePsiAngle(linkageAtoms);
    if((linkage_type_ == "1-6") || (inverse_linkage_type_ == "1-6")||(linkage_type_ == "2-6") || (inverse_linkage_type_ == "2-6"))
    {
      omega_angle_ = CalculateOmegaAngle(linkageAtoms);
    }
  }
  else
  {//Anomeric-anomeric linkage uses different atoms
    //Get anomeric C and ring O of other sugar
    linkageAtoms.push_back(non_reducing_mono_2_carbon_);
    for(std::vector<MolecularModeling::Atom*>::iterator it = non_reducing_mono_2_->cycle_atoms_.begin(); it != non_reducing_mono_2_->cycle_atoms_.end(); it++)
    {
      MolecularModeling::Atom* thisAtom = (*it);
      if(thisAtom->GetElementSymbol() == "O")
      {//Then this is the ring oxygen
        linkageAtoms.push_back(thisAtom);
      }
    }
    phi_angle_ = CalculatePhiAngle(linkageAtoms);
    //This is going to calculate phi' even though its using the Psi function
    //because we changed what atoms are being fed into it (ring O instead of Cx-1)
    phi_prime_angle_ = CalculatePsiAngle(linkageAtoms);
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

double Glycan::GlycosidicLinkage::CalculatePhiAngle(std::vector<MolecularModeling::Atom*> linkageAtoms)
{
  double Phi;
  int local_debug = -1;
  std::stringstream ss;
  if(local_debug > 0)
  {
    gmml::log(__LINE__, __FILE__, gmml::INF, "CalculatingPhiAngle");
    if(linkageAtoms.size() > 3)
    {
      ss << "Atoms being used are: ";
      ss << linkageAtoms[0]->GetId() << ", " << linkageAtoms[1]->GetId() << ", ";
      ss << linkageAtoms[2]->GetId() << ", and " << linkageAtoms[3]->GetId();
      gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
      ss.str("");
    }
  }
  if(linkageAtoms.size() < 4)
  {
    gmml::log(__LINE__, __FILE__, gmml::ERR, "Not enough atoms in vector to calculate Phi");
    Phi = -9999;
    return Phi;
  }
  // MolecularModeling::Atom* O5 = NULL;
  // MolecularModeling::Atom* C1 = non_reducing_mono_carbon_;
  // MolecularModeling::Atom* glycosidicO = NULL;
  // MolecularModeling::Atom* Cx = NULL;
  //
  // //Get Anomeric Carbon and ring oxygen
  // if(non_reducing_mono_ != NULL)
  // {
  //   for(std::vector<MolecularModeling::Atom*>::iterator it = non_reducing_mono_->cycle_atoms_.begin();
  //       it != non_reducing_mono_->cycle_atoms_.end(); it++)
  //   {
  //     MolecularModeling::Atom* this_cycle_atom = (*it);
  //     std::string this_element = this_cycle_atom->GetElementSymbol();
  //     std::string this_atom_id = this_cycle_atom->GetId();
  //
  //     if (this_atom_id.find("O") != std::string::npos)
  //     {
  //       O5 = this_cycle_atom;
  //     }
  //   }
  // }
  // else
  // {
  //   gmml::log(__LINE__, __FILE__, gmml::ERR, "C1 is null in calculate Phi Angle");
  //   return -9999;
  // }
  // if((O5 == NULL))
  // {
  //   gmml::log(__LINE__, __FILE__, gmml::ERR, "O5 is null in calculate Phi Angle");
  //   return -9999;
  // }
  // if(local_debug > 0)
  // {
  //   ss << "O5: " << O5->GetId();
  //   gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
  //   ss.str("");
  //   ss << "C1: " << C1->GetId();
  //   gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
  //   ss.str("");
  // }
  // //Find glycosidicO
  // if(local_debug > 0)
  // {
  //   gmml::log(__LINE__, __FILE__,  gmml::INF, "Get glycosidic O";
  // }
  // std::vector<MolecularModeling::Atom*> C1_neighbors = C1->GetNode()->GetNodeNeighbors();
  // for (std::vector<MolecularModeling::Atom*>::iterator it = C1_neighbors.begin(); it != C1_neighbors.end(); it++)
  // {
  //   MolecularModeling::Atom* this_neighbor = (*it);
  //   std::string this_element = this_neighbor->GetElementSymbol();
  //   if(local_debug > 0)
  //   {
  //     ss << "Is Ring: " << this_neighbor->GetIsRing();
  //     gmml::log(__LINE__, __FILE__,  gmml::INF, ss.str());
  //     ss.str("");
  //     ss << this_element << " Neighbor ID: " << this_neighbor->GetId();
  //     gmml::log(__LINE__, __FILE__,  gmml::INF, ss.str());
  //     ss.str("");
  //   }
  //   if (this_element =="O" && !this_neighbor->GetIsRing())
  //   {
  //     glycosidicO = this_neighbor;
  //   }
  // }
  // if(glycosidicO == NULL)
  // {
  //   // gmml::log(__LINE__, __FILE__, gmml::INF, glycosidicO->GetId());
  //   gmml::log(__LINE__, __FILE__, gmml::ERR, "glycosidicO is null in calculate Phi Angle");
  //   return -9999;
  // }
  // else
  // {
  //   glycosidic_oxygen_ = glycosidicO;
  // }
  //
  // //Get Cx
  // std::vector<MolecularModeling::Atom*> glycosidicO_neighbors = glycosidicO->GetNode()->GetNodeNeighbors();
  // for (std::vector<MolecularModeling::Atom*>::iterator it = glycosidicO_neighbors.begin(); it != glycosidicO_neighbors.end(); it++)
  // {
  //   MolecularModeling::Atom* this_neighbor = (*it);
  //   if (local_debug > 0)
  //   {
  //     ss << "Cx: " << this_neighbor->GetId();
  //     gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
  //     ss.str("");
  //   }
  //   if (this_neighbor != C1)
  //   {
  //     Cx = this_neighbor;
  //   }
  // }
  // if(Cx == NULL)
  // {
  //   // gmml::log(__LINE__, __FILE__, gmml::INF, Cx->GetId());
  //   gmml::log(__LINE__, __FILE__, gmml::ERR, "Cx is null in calculate Phi Angle");
  //   return -9999;
  // }
  //
  // if(local_debug > 0)
  // {
  //   gmml::log(__LINE__, __FILE__, gmml::INF, "glycosidicO");
  //   gmml::log(__LINE__, __FILE__, gmml::INF, glycosidicO->GetId());
  //   gmml::log(__LINE__, __FILE__, gmml::INF, "Cx");
  //   gmml::log(__LINE__, __FILE__, gmml::INF, Cx->GetId());
  // }

  //(O5-C1-O-Cx') {Ring oxygen of child_oligo}-{child_atom_id}-{glycosidic_atom_id}-{parent_atom_id}
  Phi = non_reducing_mono_->assembly_->CalculateTorsionAngleByAtoms(linkageAtoms[0],linkageAtoms[1],linkageAtoms[2],linkageAtoms[3]);
  Phi = gmml::ConvertRadian2Degree(Phi);
  if (local_debug > 0)
  {
    std::stringstream ss;
    ss << "Phi: " << Phi;
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
  }
  return Phi;
}

double Glycan::GlycosidicLinkage::CalculatePsiAngle(std::vector<MolecularModeling::Atom*> linkageAtoms)
{
  //(C1-O-Cx'-C[x-1]') {child_atom_id}-{glycosidic_atom_id}-{parent_atom_id}-{parent_atom_id - 1}
  double Psi;
  int local_debug = -1;
  std::stringstream ss;
  if(local_debug > 0)
  {
    gmml::log(__LINE__, __FILE__, gmml::INF, "CalculatingPsiAngle");
    if(linkageAtoms.size() > 4)
    {
      ss << "Atoms being used are: ";
      ss << linkageAtoms[1]->GetId() << ", " << linkageAtoms[2]->GetId() << ", ";
      ss << linkageAtoms[3]->GetId() << ", and " << linkageAtoms[4]->GetId();
      gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
      ss.str("");
    }
  }
  if(linkageAtoms.size() < 5)
  {
    gmml::log(__LINE__, __FILE__, gmml::ERR, "Not enough atoms in vector to calculate Psi");
    Psi = -9999;
    return Psi;
  }
  // MolecularModeling::Atom* C1 = non_reducing_mono_carbon_;
  // MolecularModeling::Atom* glycosidicO = NULL;
  // MolecularModeling::Atom* Cx = NULL;
  // MolecularModeling::Atom* Cx_1 = NULL;
  //
  // //Get Anomeric Carbon
  // if(C1 == NULL)
  // {
  //   gmml::log(__LINE__, __FILE__, gmml::ERR, "C1 is null in calculate Psi Angle");
  //   return -9999;
  // }
  //
  // //Find glycosidicO
  // std::vector<MolecularModeling::Atom*> C1_neighbors = C1->GetNode()->GetNodeNeighbors();
  // for (std::vector<MolecularModeling::Atom*>::iterator it = C1_neighbors.begin(); it != C1_neighbors.end(); it++)
  // {
  //   MolecularModeling::Atom* this_neighbor = (*it);
  //   std::string this_atom_id = this_neighbor->GetId();
  //   if (this_atom_id.find("O") != std::string::npos && !this_neighbor->GetIsCycle())
  //   {
  //     glycosidicO = this_neighbor;
  //   }
  // }
  // if(glycosidicO == NULL)
  // {
  //   gmml::log(__LINE__, __FILE__, gmml::ERR, "glycosidicO is null in calculate Psi Angle");
  //   return -9999;
  // }
  //
  // //Get Cx
  // std::vector<MolecularModeling::Atom*> glycosidicO_neighbors = glycosidicO->GetNode()->GetNodeNeighbors();
  // for (std::vector<MolecularModeling::Atom*>::iterator it = glycosidicO_neighbors.begin(); it != glycosidicO_neighbors.end(); it++)
  // {
  //   MolecularModeling::Atom* this_neighbor = (*it);
  //   if (this_neighbor != C1)
  //   {
  //     Cx = this_neighbor;
  //   }
  // }
  // if(Cx == NULL)
  // {
  //   gmml::log(__LINE__, __FILE__, gmml::ERR, "Cx is null in calculate Psi Angle");
  //   return -9999;
  // }
  //
  // //Get Cx-1
  // std::vector<MolecularModeling::Atom*> Cx_neighbors = Cx->GetNode()->GetNodeNeighbors();
  // for (std::vector<MolecularModeling::Atom*>::iterator it = Cx_neighbors.begin(); it != Cx_neighbors.end(); it++)
  // {
  //   MolecularModeling::Atom* this_neighbor = (*it);
  //   std::string this_atom_id = this_neighbor->GetId();
  //   int Cx_atom_number = Cx->GetId().at(1) - '0';
  //   int this_atom_number = this_atom_id.at(1) - '0';
  //   if ((Cx_atom_number != 1) && ((Cx_atom_number - 1) == this_atom_number))
  //   {
  //     Cx_1 = this_neighbor;
  //   }
  //   else if((Cx_atom_number == 1) && ((Cx_atom_number + 1) == this_atom_number))
  //   {
  //     Cx_1 = this_neighbor;
  //   }
  // }
  // if(Cx_1 == NULL)
  // {
  //   // gmml::log(__LINE__, __FILE__, gmml::INF, Cx->GetId());
  //   gmml::log(__LINE__, __FILE__, gmml::ERR, "Cx_1 is null in calculate Psi Angle");
  //   return -9999;
  // }

  Psi = non_reducing_mono_->assembly_->CalculateTorsionAngleByAtoms(linkageAtoms[1],linkageAtoms[2],linkageAtoms[3],linkageAtoms[4]);
  Psi = gmml::ConvertRadian2Degree(Psi);
  if (local_debug > 0)
  {
    std::stringstream ss;
    ss << "Psi: " << Psi;
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
  }
  return Psi;
}

double Glycan::GlycosidicLinkage::CalculateOmegaAngle(std::vector<MolecularModeling::Atom*> linkageAtoms)
{
  double Omega;
  int local_debug = -1;
  std::stringstream ss;
  //(O-C6'-C5'-O5') {glycosidic_atom_id}-{parent_atom_id}-{Carbon 5 in parent oligo}-{Ring oxygen in parent_oligo}
  if (local_debug > 0)
  {
    gmml::log(__LINE__, __FILE__,  gmml::INF, "CalculateOmegaAngle()");
    if(linkageAtoms.size() > 5)
    {
      ss << "Atoms being used are: ";
      ss << linkageAtoms[2]->GetId() << ", " << linkageAtoms[3]->GetId() << ", ";
      ss << linkageAtoms[4]->GetId() << ", and " << linkageAtoms[5]->GetId();
      gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
      ss.str("");
    }
  }
  if(linkageAtoms.size() < 6)
  {
    gmml::log(__LINE__, __FILE__, gmml::ERR, "Not enough atoms in vector to calculate Omega");
    Omega = -9999;
    return Omega;
  }
  // MolecularModeling::Atom* glycosidicO = NULL;
  // MolecularModeling::Atom* C6prime = NULL;
  // MolecularModeling::Atom* C5prime = NULL;
  // MolecularModeling::Atom* O5prime = NULL;
  //
  // //Get C5' and O5' first
  // if (local_debug > 0)
  // {
  //   gmml::log(__LINE__, __FILE__,  gmml::INF, "Get C5' and O5' first");
  // }
  // for(std::vector<MolecularModeling::Atom*>::iterator it = reducing_mono_->cycle_atoms_.begin();
  //     it != reducing_mono_->cycle_atoms_.end(); it++)
  // {
  //   MolecularModeling::Atom* this_cycle_atom = (*it);
  //   std::string this_atom_id = this_cycle_atom->GetId();
  //   int this_atom_number = this_atom_id.at(1) - '0';
  //   std::string this_element = this_cycle_atom->GetElementSymbol();
  //   if (local_debug > 0)
  //   {
  //     std::stringstream ss;
  //     ss << "Is Ring Atom: " << this_cycle_atom->GetIsRing();
  //     gmml::log(__LINE__, __FILE__,  gmml::INF, ss.str());
  //     ss.clear();
  //     ss << " | " << this_cycle_atom->GetElementSymbol() << "_" << this_atom_number << " Neighbor ID: " << this_cycle_atom->GetId();
  //     gmml::log(__LINE__, __FILE__,  gmml::INF, ss.str());
  //     ss.clear();
  //   }
  //   if (this_atom_number == 5)
  //   {
  //     if (this_element == "C")
  //     {
  //       C5prime = this_cycle_atom;
  //     }
  //     else if (this_element == "O")
  //     {
  //       O5prime = this_cycle_atom;
  //     }
  //   }
  // }
  // if ((C5prime == NULL)||(O5prime == NULL))
  // {
  //   gmml::log(__LINE__, __FILE__, gmml::ERR, "O5prime or C5prime is null in calculate Omega Angle");
  //   return -9999;
  // }
  //
  // //Get C6'
  // if (local_debug > 0)
  // {
  //   gmml::log(__LINE__, __FILE__,  gmml::INF, "Get C6'");
  // }
  // std::vector<MolecularModeling::Atom*> C5prime_neighbors = C5prime->GetNode()->GetNodeNeighbors();
  // for (std::vector<MolecularModeling::Atom*>::iterator it = C5prime_neighbors.begin(); it != C5prime_neighbors.end(); it++)
  // {
  //   MolecularModeling::Atom* this_neighbor = (*it);
  //   if (local_debug > 0)
  //   {
  //     std::stringstream ss;
  //     ss << "Is Ring Atom: " << this_neighbor->GetIsRing();
  //     gmml::log(__LINE__, __FILE__,  gmml::INF, ss.str());
  //     ss.clear();
  //     ss << " | " << this_neighbor->GetElementSymbol() << " Neighbor ID: " << this_neighbor->GetId();
  //     gmml::log(__LINE__, __FILE__,  gmml::INF, ss.str());
  //     ss.clear();
  //   }
  //   if ( !this_neighbor->GetIsRing() &&  this_neighbor->GetElementSymbol() == "C")
  //   {
  //     C6prime = this_neighbor;
  //   }
  // }
  // if(C6prime == NULL)
  // {
  //   gmml::log(__LINE__, __FILE__, gmml::ERR, "C6prime is null in calculate Omega Angle");
  //   return -9999;
  // }
  // //Get glycosidic oxygen
  // if (local_debug > 0)
  // {
  //   gmml::log(__LINE__, __FILE__,  gmml::INF, "Get glycosidic oxygen");
  // }
  // std::vector<MolecularModeling::Atom*> C6prime_neighbors = C6prime->GetNode()->GetNodeNeighbors();
  // for (std::vector<MolecularModeling::Atom*>::iterator it = C6prime_neighbors.begin(); it != C6prime_neighbors.end(); it++)
  // {
  //   MolecularModeling::Atom* this_neighbor = (*it);
  //   if (local_debug > 0)
  //   {
  //     std::stringstream ss;
  //     ss << "Is Cycle: " << this_neighbor->GetIsCycle();
  //     gmml::log(__LINE__, __FILE__,  gmml::INF, ss.str());
  //     ss.str("");
  //     ss << " | " << this_neighbor->GetElementSymbol() << " Neighbor ID: " << this_neighbor->GetId();
  //     gmml::log(__LINE__, __FILE__,  gmml::INF, ss.str());
  //     ss.str("");
  //   }
  //   if (this_neighbor->GetElementSymbol() == "O" && !this_neighbor->GetIsCycle())
  //   {
  //     glycosidicO = this_neighbor;
  //   }
  // }
  // if(glycosidicO == NULL)
  // {
  //   gmml::log(__LINE__, __FILE__, gmml::ERR, "glycosidicO is null in calculate Omega Angle");
  //   return -9999;
  // }
  // if (local_debug > 0)
  // {
  //   gmml::log(__LINE__, __FILE__,  gmml::INF, "About to calcuate torsion");
  // }
  Omega = non_reducing_mono_->assembly_->CalculateTorsionAngleByAtoms(linkageAtoms[2],linkageAtoms[3],linkageAtoms[4],linkageAtoms[5]);
  Omega = gmml::ConvertRadian2Degree(Omega);
  return Omega;
}
