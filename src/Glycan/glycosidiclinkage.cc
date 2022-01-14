//Created 12/11/18
//Dave Montgomery

#include <vector>
#include <string>
#include <sstream>

#include "../../includes/MolecularModeling/assembly.hpp"
//#include "../../includes/utils.hpp"
#include "../../includes/MolecularModeling/residue.hpp"
#include "../../includes/Glycan/oligosaccharide.hpp"
#include "../../includes/Glycan/monosaccharide.hpp"
#include "../../includes/Glycan/glycosidiclinkage.hpp"
#include "../../includes/GeometryTopology/geometrytopology.hpp"
#include "includes/CodeUtils/logging.hpp"

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
    {
      gmml::log(__LINE__, __FILE__, gmml::INF, "Source Mono Anomeric Carbon: " + sourceMono->anomeric_carbon_pointer_->GetId());
    }

    gmml::log(__LINE__, __FILE__, gmml::INF, "Source Carbon ID " + source_carbon_ID);

    if(targetMono->anomeric_carbon_pointer_ != NULL)
    {
      gmml::log(__LINE__, __FILE__, gmml::INF, "Target Mono Anomeric Carbon: " + targetMono->anomeric_carbon_pointer_->GetId());
    }
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
          //TODO: This should be able to identify nonstandard glycosidic atoms, like C and N
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
        //TODO: This should be able to identify nonstandard glycosidic atoms, like C and N
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
  Phi = GeometryTopology::CalculateDihedralAngle(linkageAtoms[0]->GetCoordinate(),linkageAtoms[1]->GetCoordinate(),
                                                 linkageAtoms[2]->GetCoordinate(),linkageAtoms[3]->GetCoordinate());

  //Below being replaced as we are trying to use one function for torsions
  //non_reducing_mono_->assembly_->CalculateTorsionAngleByAtoms(linkageAtoms[0],linkageAtoms[1],linkageAtoms[2],linkageAtoms[3]);
  //Phi = gmml::ConvertRadian2Degree(Phi);
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

  Psi = GeometryTopology::CalculateDihedralAngle(linkageAtoms[1]->GetCoordinate(),linkageAtoms[2]->GetCoordinate(),
                                                 linkageAtoms[3]->GetCoordinate(),linkageAtoms[4]->GetCoordinate());
  //Psi = gmml::ConvertRadian2Degree(Psi);
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
  Omega = GeometryTopology::CalculateDihedralAngle(linkageAtoms[2]->GetCoordinate(),linkageAtoms[3]->GetCoordinate(),
                                                   linkageAtoms[4]->GetCoordinate(),linkageAtoms[5]->GetCoordinate());
  //Omega = gmml::ConvertRadian2Degree(Omega);
  return Omega;
}

double Glycan::GlycosidicLinkage::CalculatePhiChiEnergy()
{
  /////////////////////////////////////////////////////////////////////
  //CHI_Energy function taken from https://doi.org/10.1002/jcc.23517 //
  /////////////////////////////////////////////////////////////////////

  //       N               (x - b_i)^2
  //f(x)=  Σ  a_i * e * -(――――――――――――――) + Offset
  //      i=1                  c_i

  //There are 4 different energy functions for CHI Phi Energy
  bool phi_alpha = false;
  bool phi_beta = false;
  bool phi_2_3 = false;
  bool phi_2_6 = false;

  //All depend on the phi_angle,
  //-9999 means something went wrong or it wasn't set
  if(phi_angle_ != -9999)
  {
    ////////////////////////////////////////////
    //Logic to determine which function to use//
    ////////////////////////////////////////////

    //CHI Energy functions have only been created for certain linkages
    //The following list applies for alpha and beta linkages
    std::list<std::string> ChiLinkTypes = { "1-1", "1-2", "1-3",
                                                  "1-4", "1-6" };
    if(find(ChiLinkTypes.begin(), ChiLinkTypes.end(), linkage_type_) != ChiLinkTypes.end())
    {
      if(anomeric_configuration_ == "alpha")
      {
        bool phi_alpha = true;
      }
      else if(anomeric_configuration_ == "beta")
      {
        bool phi_beta = true;
      }
      else
      {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Unable to determine the configuration (α/β) of the linkage");
      }
    }
    else if(find(betaVector.begin(), betaVector.end(), linkage_type_) != betaVector.end())
    {

    }

    ////////////////////////
    //CHI Energy Functions//
    ////////////////////////
    if (phi_alpha)
    {
      double LH = 2.97696467271672,
             Lc = -199.494365163839,
             LW = 677.808323900125;
      double RH = 102.253303636096,
             Rc = 170.599580473404,
             RW = 1696.7844369942;
      double aH = 10.7448005875571,
             ac = -105.313553566706,
             aW = 4724.58364072706;
      double bH = 3.67344580413578,
             bc = 6.20116671874232,
             bW = 1347.7205625156;
      double cH = 2.06094652655659,
             cc = 91.6553021324274,
             cW = 1500.02002601097;
      double dH = 6.19388683252667,
             dc = -22.9786969888816,
             dW = 2122.27783139301;
      double eH = -2.11153017593601,
             ec = 83.6019123356148,
             eW = 1254.13371108961;
      double fH = -98.0013005657107,
             fc = 170.012289132741,
             fW = 1598.73272567307;

      double Off = 1.00501e-30;

      double Leftx, Rightx, ax, bx, cx, dx, ex, fx, Totx;

      Leftx = LH * exp(-pow((phi_angle_-(Lc)),2.0)/LW);
      Rightx = RH * exp(-pow((phi_angle_-(Rc)),2.0)/RW);
      ax = aH * exp(-pow((phi_angle_-(ac)),2.0)/aW);
      bx = bH * exp(-pow((phi_angle_-(bc)),2.0)/bW);
      cx = cH * exp(-pow((phi_angle_-(cc)),2.0)/cW);
      dx = dH * exp(-pow((phi_angle_-(dc)),2.0)/dW);
      ex = eH * exp(-pow((phi_angle_-(ec)),2.0)/eW);
      fx = fH * exp(-pow((phi_angle_-(fc)),2.0)/fW);
      Totx = Rightx + Leftx + ax + bx + cx + dx + ex + fx;
      return Totx;
    }
    else if (phi_beta)
    {
      double Lc = -330.769995527134, aH = 5.93533323829663, ac = -152.080139620062, aW = 6049.77220005964, bH = 22.467372096061, bc = -23.5159916173247, bW = 606.89715970453, cH = 10.0360057033439, cc = 120.962836525241, cW = 4037.89330459304, LH = 450.540038600828, LW = 4449.7622241787, RH = 23.7118506901333, Rc = 304.624980492529, RW = 8375.1929028027, /*Off = -2.27829251796721,*/ Off = -2.1283, dH = -18.1406478565247, dc = -24.2677756921736, dW = 543.050986049266, eH = 5.88226333077368, ec = 19.6321032903376, eW = 897.92664572344, Leftx, Rightx, ax, bx, cx, dx, ex, x, Totx;
      x=phi_angle;
      Leftx = LH * exp(-pow((x-(Lc)),2.0)/LW);
      Rightx = RH * exp(-pow((x-(Rc)),2.0)/RW);
      ax = aH * exp(-pow((x-(ac)),2.0)/aW);
      bx = bH * exp(-pow((x-(bc)),2.0)/bW);
      cx = cH * exp(-pow((x-(cc)),2.0)/cW);
      dx = dH * exp(-pow((x-(dc)),2.0)/dW);
      ex = eH * exp(-pow((x-(ec)),2.0)/eW);
      Totx = Rightx + Leftx + ax + bx + cx + dx + ex + Off;
      return Totx;
    }

    //In the SI for Vina Carb
    //(https://pubs.acs.org/doi/suppl/10.1021/acs.jctc.5b00834/suppl_file/ct5b00834_si_002.pdf)
    //The ranges for Phi of 2-3 and 2-6 linkages are defined as:
    //(The range values had typos but I checked the figures)
    //
    //               0.0018(x-60)^2,          if x ϵ [0,120]
    //ΔE_Sia23 (x) = 0.0018(x-180)^2 + 1.95,  if x ϵ [120,240]
    //               0.0018(x-300)^2 + 0.385, if x ϵ [240,360]
    //
    //               0.0018(x-60)^2,          if x ϵ [0,120]
    //ΔE_Sia26 (x) = 0.0018(x-180)^2 + 1.19,  if x ϵ [120,240]
    //               0.0018(x-300)^2 + 1.41,  if x ϵ [240,360]
    //
    //However, the ranges overlap at 120, 240, and 360/0.
    //The values are pretty close, so I'm just going to use
    //if x ϵ [0,120]
    //if x ϵ (120,240]
    //if x ϵ (240,360]

    else if (phi_2_3)
    {
      //if x ϵ [0,120]
      if((0 <= phi_angle_ ) && (phi_angle_ <= 120))
      {
        phi_CHI_Energy_ = 0.0018 * pow( (phi_angle - 60), 2);
      }
      //if x ϵ (120,240]
      else if((120 < phi_angle_) && (phi_angle_ <= 240))
      {
        phi_CHI_Energy_ = 0.0018 * pow( (phi_angle - 180), 2) + 1.95;
      }
      //if x ϵ (240,360]
      else if((240 < phi_angle_) && (phi_angle_ <= 360))
      {
        phi_CHI_Energy_ = 0.0018 * pow( (phi_angle - 300), 2) + 0.385;
      }
      else
      {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Phi was not calculated properly. Unable to run CHI Energy analysis");
      }
    }
    else if (phi_2_6)
    {
      //if x ϵ [0,120]
      if((0 <= phi_angle_ ) && (phi_angle_ <= 120))
      {
        phi_CHI_Energy_ = 0.0018 * pow( (phi_angle - 60), 2);
      }
      //if x ϵ (120,240]
      else if((120 < phi_angle_) && (phi_angle_ <= 240))
      {
        phi_CHI_Energy_ = 0.0018 * pow( (phi_angle - 180), 2) + 1.19;
      }
      //if x ϵ (240,360]
      else if((240 < phi_angle_) && (phi_angle_ <= 360))
      {
        phi_CHI_Energy_ = 0.0018 * pow( (phi_angle - 300), 2) + 1.41;
      }
      else
      {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Phi was not calculated properly. Unable to run CHI Energy analysis");
      }
    }
  }
  else
  {
    gmml::log(__LINE__, __FILE__, gmml::ERR, "Phi was not calculated properly.  Unable to run CHI Energy analysis");
  }
}

double Glycan::GlycosidicLinkage::CalculatePsiChiEnergy()
{
  //////////////////////////////////////////////////////////////////////
  // CHI_Energy function taken from https://doi.org/10.1002/jcc.23517 //
  //////////////////////////////////////////////////////////////////////

  //       N               (x - b_i)^2
  //f(x)=  Σ  a_i * e * -(――――――――――――――) + Offset
  //      i=1                  c_i

  //////////////////////////////////////////////////////////////////////
  //    For 1-6, 2-3, & 2-6 energy function details see the SI in     //
  //          https://doi.org/10.1021/acs.jctc.5b00834                //
  //////////////////////////////////////////////////////////////////////



  //There are 4 different energy functions for CHI Psi Energy
  //each with their own criteria
  bool psi_CHI_1 = false; //1-2 & 1-4 axial, 1-3, 2-3, and 2-6 equatorial
  bool psi_CHI_2 = false; //
  bool psi_CHI_3 = false; //
  bool psi_CHI_4 = false; //

  double model::psi_2A3E_energy(double psi_angle)
  {
  double LH = 4.62366277694845, Lc = 5.045583934308, LW = 5005.75236060956, RH = 4.61387450239844, Rc = 362.487847702007, RW = 2090.63190217702, aH = 4.94191727813274, ac = 121.202321824468, aW = 2093.75214491931, bH = 0.402901504292045, bc = 241.428583877882, bW = 456.828754790442, cH = 0.798876573705798, cc = 68.425080241155, cW = 678.807178379645, Off = -0.125645118474882, dc = 192.925748017071, dW = 347.244734136509, dH = 0.222992242737354, Leftx, Rightx, ax, bx, cx, dx, x, Totx;
  if(psi_angle<0)
  {
  psi_angle=360+psi_angle;}
  x=psi_angle;
  Leftx = LH * exp(-pow((x-(Lc)),2.0)/LW);
  Rightx = RH * exp(-pow((x-(Rc)),2.0)/RW);
  ax = aH * exp(-pow((x-(ac)),2.0)/aW);
  bx = bH * exp(-pow((x-(bc)),2.0)/bW);
  cx = cH * exp(-pow((x-(cc)),2.0)/cW);
  dx = dH * exp(-pow((x-(dc)),2.0)/dW);
  Totx = Rightx + Leftx + ax + bx + cx + dx + Off;
  return Totx;
  }

  double model::psi_2E3A_energy(double psi_angle)
  {
  double LH = 4.46811874171788, Lc = 1e-30, LW = 1279.58772056199, RH = 4.38204018882823, Rc = 357.770654336205, RW = 6050.14162479438, aH = 284.944948778136, ac = 146.644068129462, aW = 1551.75673776163, bH = 4.76134025362478, bc = 220.683118921686, bW = 5892.94143218231, cH = -169.197666368856, cc = 147.370828680332, cW = 1742.47541063603, Off = 1.0219924486158, dc = 146.05660843428, dW = 1359.82873591396, dH = -118.440552792375, Leftx, Rightx, ax, bx, cx, dx, x, Totx;
  if(psi_angle<0){
  psi_angle=360+psi_angle;}
  x=psi_angle;
  Leftx = LH * exp(-pow((x-(Lc)),2.0)/LW);
  Rightx = RH * exp(-pow((x-(Rc)),2.0)/RW);
  ax = aH * exp(-pow((x-(ac)),2.0)/aW);
  bx = bH * exp(-pow((x-(bc)),2.0)/bW);
  cx = cH * exp(-pow((x-(cc)),2.0)/cW);
  dx = dH * exp(-pow((x-(dc)),2.0)/dW);
  Totx = Rightx + Leftx + ax + bx + cx + dx + Off;
  return Totx;
  }

  double model::psi_6A_energy(double psi_angle)
  {
  double aH = 67.9431348410598, ac = -59.5393395706705, aW = 993.323581145538, bH = 6.13421142432396, bc = 10.4786088782815, bW = 945.770771330812, cH = 3.27628727235978, cc = 54.2960678151208, cW = 851.528141801851, dH = 0.727486729062442, dc = 131.067737803489, dW = 1037.41211378392, eH = 2.57362265878937, ec = 245.102891425541, eW = 2012.99451568206, fH = 5.75995973448166, fc = 359.999988549478, fW = 1153.3974275673, gH = 3.47492643928157, gc = 321.677942414686, gW = 2080.97053159226, hH = -0.741000462200939, hc = 199.106903524814, hW = 522.180434119001, ax, bx, cx, dx, ex, fx, gx, hx, x, Totx;

  if(psi_angle<0)
  {
  psi_angle=360+psi_angle;
  }
  x=psi_angle;
  ax = aH * exp(-pow((x-(ac)),2.0)/aW);
  bx = bH * exp(-pow((x-(bc)),2.0)/bW);
  cx = cH * exp(-pow((x-(cc)),2.0)/cW);
  dx = dH * exp(-pow((x-(dc)),2.0)/dW);
  ex = eH * exp(-pow((x-(ec)),2.0)/eW);
  fx = fH * exp(-pow((x-(fc)),2.0)/fW);
  gx = gH * exp(-pow((x-(gc)),2.0)/gW);
  hx = hH * exp(-pow((x-(hc)),2.0)/hW);
  Totx = ax + bx + cx + dx + ex + fx + gx + hx;
  return Totx;
  }

  double model::psi_6E_energy(double psi_angle)
  {
  double aH = 7.24858655753829, ac = 3.60600554520403, aW = 2459.23916629141, bH = 1.9, bc = 96.5930821702371, bW = 2683.88656516991, cH = 0.741022592342903, cc = 141.663521919709, cW = 1150.04756181103, dH = 0.2, dc = 162, dW = 400, eH = 0.287090039932611, ec = 228.171314273305, eW = 272.201363844744, fH = 1.22591971967808, fc = 292.206221787048, fW = 1134.52455512381, gH = 7.41063235334191, gc = 369.867701147817, gW = 3499.15994772992, hH = -0.61489499584011, hc = 271.319024293053, hW = 532.437194483944, iH = -0.35, ic = 183, iW = 100, ax, bx, cx, dx, ex, fx, gx, hx, ix, x, Totx;

  if(psi_angle<0)
  {
  psi_angle=360+psi_angle;
  }
  x=psi_angle;
  ax = aH * exp(-pow((x-(ac)),2.0)/aW);
  bx = bH * exp(-pow((x-(bc)),2.0)/bW);
  cx = cH * exp(-pow((x-(cc)),2.0)/cW);
  dx = dH * exp(-pow((x-(dc)),2.0)/dW);
  ex = eH * exp(-pow((x-(ec)),2.0)/eW);
  fx = fH * exp(-pow((x-(fc)),2.0)/fW);
  gx = gH * exp(-pow((x-(gc)),2.0)/gW);
  hx = hH * exp(-pow((x-(hc)),2.0)/hW);
  ix = iH * exp(-pow((x-(ic)),2.0)/iW);
  Totx = ax + bx + cx + dx + ex + fx + gx + hx + ix;
  return Totx;
  }
}

double Glycan::GlycosidicLinkage::CalculateOmegaChiEnergy()
{
  double model::omega_6A_energy(double omega_angle)
  {
  double x, energy, b, k=0.0025;
          if(omega_angle<0)
          {
          x=360+omega_angle;
          }
          else
          {
          x=omega_angle;
          }
          if((x>=0.0 && x<120.0)||(x>=360.0 && x<120.0))
          {
          b=0.0;
          energy=k*pow((x-60),2)+b;
          }
          else if(x>=120.0 && x<240.0)
          {
          b=0.3;
          energy=k*pow((x-180),2)+b;
          }
          else if(x>=240.0 && x<360.0)
          {
          b=1.0;
          energy=k*pow((x-300),2)+b;
          }
  return energy;
  }


double model::omega_6E_energy(double omega_angle)
  {
  double x, energy, b, k=0.0025;
          if(omega_angle<0)
          {
          x=360+omega_angle;
          }
          else
          {
          x=omega_angle;
          }
          if((x>=0.0 && x<120.0)||(x>=360.0 && x<120.0))
          {
          b=0.21;
          energy=k*pow((x-60),2)+b;
          }
          else if(x>=120.0 && x<240.0)
          {
          b=1.39;
          energy=k*pow((x-180),2)+b;
          }
          else if(x>=240.0 && x<360.0)
          {
          b=0.0;
          energy=k*pow((x-300),2)+b;
          }
  return energy;
  }
}
