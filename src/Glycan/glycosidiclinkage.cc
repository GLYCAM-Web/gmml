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
  std::stringstream ss;//for gmml logs
  if(local_debug > 0)
  {
    if(sourceMono->anomeric_carbon_pointer_ != NULL)
    {//Print sourceMono anomeric C info
      ss.str("");
      ss << "Source Mono Anomeric Carbon: ";
      ss << sourceMono->anomeric_carbon_pointer_->GetId();
      gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
      ss.str("");
    }
    else
    {//Print error
      ss.str("");
      ss << "Source Mono Anomeric Carbon is NULL";
      gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
      ss.str("");
    }
    //Print source_carbon_ID info
    ss.str("");
    ss << "Source Carbon ID: " << source_carbon_ID;
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
    ss.str("");

    if(targetMono->anomeric_carbon_pointer_ != NULL)
    {//Print targetMono anomeric C info
      ss.str("");
      ss << "Target Mono Anomeric Carbon: ";
      ss << targetMono->anomeric_carbon_pointer_->GetId();
      gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
      ss.str("");
    }
    //Print target_carbon_ID info
    ss.str("");
    ss << "Target Carbon ID: " << target_carbon_ID;
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
    ss.str("");
  }

  /*Determine Anomeric Carbon(s) in Linkage
      Knowing which monosaccharide has the anomeric carbon lets us know
      which direction the linkage is going, unless it is anomeric-anomeric.
      This lets us pick the correct atoms for calculating torsion angles.

      Just looking at mono->anomeric_carbon_pointer_ didn't work for all
      sugars, so I added a check for carbon number: C1 is anomeric unless
      there is no C1, then C2 is.  This checking should really happen at
      the monosaccharide level, and doesn't work if the ring atoms are
      numbered incorrectly
  */
  for(std::vector<MolecularModeling::Atom*>::iterator it = sourceMono->cycle_atoms_.begin(); it != sourceMono->cycle_atoms_.end(); it++)
  {
    MolecularModeling::Atom* thisAtom = *it;
    if((thisAtom->GetId() == source_carbon_ID) &&
      (sourceMono->anomeric_carbon_pointer_ != NULL))
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
    if((thisAtom->GetId() == target_carbon_ID) &&
      (targetMono->anomeric_carbon_pointer_ != NULL))
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

  if(local_debug > 0)
  {
    if(non_reducing_mono_ == NULL)
    {//Print error
      gmml::log(__LINE__, __FILE__, gmml::ERR, "Unable to identify non_reducing_mono_");
    }
    if((reducing_mono_ == NULL) && (non_reducing_mono_2_ == NULL))
    {//Print error
      gmml::log(__LINE__, __FILE__, gmml::ERR, "Unable to identify reducing_mono_ or 2nd non_reducing_mono_");
    }
  }

  /*Determine the anomeric_configuration_
  */
  if(non_reducing_mono_ != NULL)
  {
    anomeric_configuration_ = non_reducing_mono_->sugar_name_.configuration_;
    if(local_debug > 0)
    {
      ss.str("");
      ss << "The anomeric configuration is: ";
      ss << anomeric_configuration_;
      gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
      ss.str("");
    }
  }

  /*Determine the Linkage Type (linkage_type_)
      This also determines inverse_linkage_type_
      "1-4", "2-6", etc.
  */
  /*TODO Fix how linkage type is determined
      Currently this returns weird things like "-", "A-B", etc.
      We should have/use a function for the ring and exocyclic
      atom numbering for monosaccharides

      We should also split the linkage type determination into another function
  */
  if(non_reducing_mono_carbon_ != NULL)
  {
    std::stringstream LinkSs, linkNameSs, linkResNameSs;
    if(non_reducing_mono_->sugar_name_.monosaccharide_short_name_ != "")
    {
      linkNameSs << non_reducing_mono_->sugar_name_.monosaccharide_short_name_;
    }
    else if(non_reducing_mono_->sugar_name_.monosaccharide_stereochemistry_short_name_ != "")
    {
      linkNameSs << non_reducing_mono_->sugar_name_.monosaccharide_stereochemistry_short_name_;
    }
    else
    {//PRINT ERROR
      if(local_debug > 0)
        {
          ss.str("");
          ss << "Unable to determine Linkage name";
          gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
          ss.str("");
        }
    }
    linkResNameSs << non_reducing_mono_carbon_->GetResidue()->GetId();
    linkResNameSs << " (";
    linkResNameSs << non_reducing_mono_carbon_->GetId()[1] << "-";
    linkNameSs << non_reducing_mono_carbon_->GetId()[1];
    linkNameSs << "-";
    LinkSs << non_reducing_mono_carbon_->GetId()[1];
    LinkSs  << "-";
    if(anomeric_anomeric_linkage_)
    {
      if(non_reducing_mono_2_carbon_ != NULL)
      {
        LinkSs << non_reducing_mono_2_carbon_->GetId()[1];
        linkNameSs << non_reducing_mono_2_carbon_->GetId()[1];
        linkResNameSs << non_reducing_mono_2_carbon_->GetId()[1];
        linkResNameSs << ") ";
        linkResNameSs << non_reducing_mono_2_carbon_->GetResidue()->GetId();
        if(non_reducing_mono_2_->sugar_name_.monosaccharide_short_name_ != "")
        {
          linkNameSs << non_reducing_mono_2_->sugar_name_.monosaccharide_short_name_;
        }
        else if(non_reducing_mono_2_->sugar_name_.monosaccharide_stereochemistry_short_name_ != "")
        {
          linkNameSs << non_reducing_mono_2_->sugar_name_.monosaccharide_stereochemistry_short_name_;
        }
        else
        {//PRINT ERROR
          if(local_debug > 0)
            {
              ss.str("");
              ss << "Unable to determine Linkage name";
              gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
              ss.str("");
            }
        }
      }
      else
      {
        if(local_debug > 0)
        {
          ss.str("");
          ss << "Unable to determine Linkage, non-reducing ";
          ss << "monosaccharide 2 carbon not set.";
          gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
          ss.str("");
        }
      }
    }
    else
    {
      if(reducing_mono_carbon_ != NULL)
      {
        LinkSs << reducing_mono_carbon_->GetId()[1];
        linkNameSs << reducing_mono_carbon_->GetId()[1];
        linkResNameSs << reducing_mono_carbon_->GetId()[1];
        linkResNameSs << ") ";
        linkResNameSs << reducing_mono_carbon_->GetResidue()->GetId();
        if(reducing_mono_->sugar_name_.monosaccharide_short_name_ != "")
        {
          linkNameSs << reducing_mono_->sugar_name_.monosaccharide_short_name_;
        }
        else if(reducing_mono_->sugar_name_.monosaccharide_stereochemistry_short_name_ != "")
        {
          linkNameSs << reducing_mono_->sugar_name_.monosaccharide_stereochemistry_short_name_;
        }
        else
        {//PRINT ERROR
          if(local_debug > 0)
            {
              ss.str("");
              ss << "Unable to determine Linkage name";
              gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
              ss.str("");
            }
        }
      }
      else
      {
        if(local_debug > 0)
        {
          ss.str("");
          ss << "Unable to determine Linkage, reducing ";
          ss << "monosaccharide carbon not set.";
          gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
          ss.str("");
        }
      }
    }
    if(local_debug > 0)
    {//Print info on LinkSs
      ss.str("");
      ss << "Linkage stringstream is: " << LinkSs.str();
      gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
      ss.str("");
    }
    residue_linkage_name_ = linkResNameSs.str();
    linkage_name_ = linkNameSs.str();
    linkage_type_ = LinkSs.str();
    inverse_linkage_type_ = linkage_type_;
    std::reverse(inverse_linkage_type_.begin(), inverse_linkage_type_.end());
    if(local_debug > 0)
    {//Print info on Linkage
      ss.str("");
      ss << "Linkage is " << linkage_type_;
      ss << ", and inverse linkage is " << inverse_linkage_type_;
      gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
      ss.str("");
      ss << "Linkage name is " << linkage_name_;
      gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
      ss.str("");
      ss << "Linkage residue name is " << residue_linkage_name_;
      gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
      ss.str("");
    }
  }
  else
  {//something went wrong
    if(local_debug > 0)
    {
      ss.str("");
      ss << "Unable to determine Linkage, non-reducing ";
      ss << "monosaccharide carbon not set.";
      gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
      ss.str("");
    }
  }

  /*@TODO BugFix in Atom Selection
  Fix the following code to determine what types of sugars are involved
  ie Aldopyranose or Ulofuranose; it will change which atoms are needed
  for torsion angle calculations, and should reduce the number of bugs
  This will need to change the logic for calculating omega too
  as a 2-8 linkage would have an omega as well
  This is an example for Phi taken from https://iupac.qmul.ac.uk/misc/psac.html
  Where i is the residue number

  Aldopyranose: O5(i )-C(i )-OX(i-1)-CX(i-1)
  Aldofuranose: O4(i )-C1(i )-OX(i-1)-CX(i-1)
  Ulopyranose: O6(i )-C2(i )-OX(i-1)-CX(i-1)
  Ulofuranose: O5(i )-C2(i )-OX(i-1)-CX(i-1)

  Also, a lot of angles are missing when the archive is run
  */
  /*TODO Add other torsion angle definitions
      IE NMR uses different atoms, including Hydrogens
      which are often left out of crystallographic models

      This is not essential and is mostly a wishlist for
      completeness item
  */
  /*GET THE ATOMS NEEDED FOR TORSION ANGLE FUNCTIONS.
      This vector will have:
        The non_reducing sugar's ring oxygen,
        The non_reducing sugar's anomeric carbon,
        The reducing sugar's glycosidic O (Ox),
        The reducing sugar's Cx,
        The reducing sugar's Cx-1,
      and for omega (1/2-6 linked for now),
        The reducing sugar's Cx-2
  */

  std::vector<MolecularModeling::Atom*> linkageAtoms;

  bool setRingO = false;
  bool setAnomericC = false;
  bool setGlycosidicO = false;
  bool setCx = false;
  bool setCx1 = false;//Cx-1
  bool setCx2 = false;//Cx-2

  //Add Ring Oxygen
  /*TODO Add Functionality for ring Nitrogen
      This is not a priority but it's a wishlist item from Rob
  */
  int testingInt = 0;
  if(non_reducing_mono_ != NULL)
  {
    for(std::vector<MolecularModeling::Atom*>::iterator it = non_reducing_mono_->cycle_atoms_.begin(); it != non_reducing_mono_->cycle_atoms_.end(); it++)
    {
      testingInt++;
      // if(local_debug > 0)
      // {
      //   ss.str("");
      //   ss << "Cycle Atom: " << testingInt;
      //   gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
      //   ss.str("");
      //   ss << "Atom: " << (*it)->GetElementSymbol();
      // }
      MolecularModeling::Atom* thisCycleAtom = (*it);
      if(thisCycleAtom->GetElementSymbol() == "O")
      {//Then it's the ring Oxygen
        linkageAtoms.push_back(thisCycleAtom);
        setRingO = true;
      }
    }
    if(setRingO != true)//should this break/exit?
    {//Could not ID ring O
      //TODO Update once GLOBAL_DEBUG is in place
      if(local_debug > 0)
      {//Print Error
        ss.str("");
        ss << "Unable to find the ring O of the non-reducing monosaccharide: ";
        ss << non_reducing_mono_->residue_name_ << ".";
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        ss.str("");
      }
    }

    //Add Anomeric Carbon
    if(non_reducing_mono_carbon_ != NULL)
    {
      linkageAtoms.push_back(non_reducing_mono_carbon_);
      setAnomericC = true;
    }
    else//should this break/exit?
    {//could not ID anomeric C
      if(local_debug > 0)
      {//Print Error
        ss.str("");
        ss << "Unable to find the anomeric C of the non-reducing monosaccharide: ";
        ss << non_reducing_mono_->residue_name_ << ".";
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        ss.str("");
      }
    }

    //Add Glycosidic Oxygen
    /*TODO Add functionality for any element.
        Other than Oxygen, this would typically
        be C, N, or S, but can also be Se
        (See 2-Carb-33.4. Selenoglycosides at
        https://iupac.qmul.ac.uk/2carb/33.html#334),
        and I'm sure others.
        Halides can exist here, but as they can only form 1 bond,
        they will be aglycones and not glycosidic linkage atoms
    */
    std::vector<MolecularModeling::Atom*> anomericNeighbors = non_reducing_mono_carbon_->GetNode()->GetNodeNeighbors();
    for(std::vector<MolecularModeling::Atom*>::iterator it = anomericNeighbors.begin(); it != anomericNeighbors.end(); it++)
    {
      MolecularModeling::Atom* thisAnomericNeighbor = (*it);
      if(reducing_mono_!= NULL)
      {
        //Below should work; SetIsCycle is called on all atoms in the monosaccharide's ring
        if((reducing_mono_carbon_->GetResidue() ==
            thisAnomericNeighbor->GetResidue()) &&
          (thisAnomericNeighbor->GetIsCycle() == false))
        {
          if(thisAnomericNeighbor->GetElementSymbol() == "O")
          {//Then this is the glycosidic oxygen
            linkageAtoms.push_back(thisAnomericNeighbor);
            glycosidic_oxygen_ = thisAnomericNeighbor;
            setGlycosidicO = true;
          }
          else if (local_debug > 0)
          {//Print Warning for now
            ss.str("");
            ss << "Glycosidic Oxygen not present. The atom present is " << thisAnomericNeighbor->GetElementSymbol() << ".";
            gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
            ss.str("");
          }
        }
      }
      else if (non_reducing_mono_2_ != NULL)
      { //The oxygen can belong to either sugar
        if((thisAnomericNeighbor->GetElementSymbol() == "O") &&
          (thisAnomericNeighbor->GetIsCycle() == false))
        {//Then this is the glycosidic oxygen
          linkageAtoms.push_back(thisAnomericNeighbor);
          glycosidic_oxygen_ = thisAnomericNeighbor;
          setGlycosidicO = true;
        }
        else if (local_debug > 0)
        {//Print Warning for now
          ss.str("");
          ss << "Glycosidic Oxygen not present. The atom present is " << thisAnomericNeighbor->GetElementSymbol() << ".";
          gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
          ss.str("");
        }
      }
      else
      {
        gmml::log(__LINE__, __FILE__, gmml::ERR, "Cannot identify both sugars for the linkage, something went horribly wrong!");
      }
    }
    if(setGlycosidicO != true)//should this break/exit?
    {//Could not ID Glycosidic O
      //TODO Update once GLOBAL_DEBUG is in place
      if(local_debug > 0)
      {//Print Error
        ss.str("");
        ss << "Unable to find the glycosidic O of the reducing monosaccharide: ";
        ss << reducing_mono_->residue_name_ << ".";
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        ss.str("");
      }
    }

    //Add Cx
    if(!anomeric_anomeric_linkage_ && reducing_mono_carbon_ != NULL)
    {
      linkageAtoms.push_back(reducing_mono_carbon_);
      setCx = true;
    }
    else if(anomeric_anomeric_linkage_ && non_reducing_mono_2_carbon_ != NULL)
    {
      linkageAtoms.push_back(non_reducing_mono_2_carbon_);
      setCx = true;
    }
    else//should this break/exit?
    {//Could not ID Cx
      //TODO Update once GLOBAL_DEBUG is in place
      if(local_debug > 0)
      {//Print Error
        ss.str("");
        ss << "Unable to find Cx of ";
        if(!anomeric_anomeric_linkage_)
        {
          ss << "the reducing monosaccharide: ";
          ss << reducing_mono_->residue_name_ << ".";
        }
        else
        {//anomeric-anomeric
          ss << "non-reducing monosaccharide 2: ";
          ss << non_reducing_mono_2_->residue_name_ << ".";
        }
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        ss.str("");
      }
    }

    //Get Cx-1 and Cx-2 if applicable
    if(!anomeric_anomeric_linkage_)
    {//Get atoms for calculations and calculate

      /*TODO Write smarter logic
          For determining Cx-1, not all files or molecules have
          ring atoms that are numbered correctly, so this would
          ideally use ring atom position information that needs
          to be added to the monosacharide class
      */
      int x = std::atoi(&reducing_mono_carbon_->GetId()[1]);
      std::vector<MolecularModeling::Atom*> CxNeighbors = reducing_mono_carbon_->GetNode()->GetNodeNeighbors();
      for(std::vector<MolecularModeling::Atom*>::iterator it = CxNeighbors.begin(); it != CxNeighbors.end(); it++)
      {
        MolecularModeling::Atom* thisCxNeighbor = (*it);
        int thisNeighborNum = std::atoi(&thisCxNeighbor->GetId()[1]);
        if((thisNeighborNum = (x - 1)) && (thisCxNeighbor->GetElementSymbol() == "C"))
        {//This is Cx-1
          //Add Cx-1
          linkageAtoms.push_back(thisCxNeighbor);
          setCx1 = true;

          /*TODO Update logic for omega
              All glycosidic linkages with more than 1 bridge node
                IE 1-6, where O6 and C6 are bridge nodes connecting
                  the anomeric carbon with C5
              should have omega angles
          */
          if((linkage_type_ == "1-6") ||
            (inverse_linkage_type_ == "1-6")||
            (linkage_type_ == "2-6") ||
            (inverse_linkage_type_ == "2-6"))
          {//Get Cx-2
            std::vector<MolecularModeling::Atom*> CxMinusOneNeighbors = thisCxNeighbor->GetNode()->GetNodeNeighbors();
            for(std::vector<MolecularModeling::Atom*>::iterator it2 = CxMinusOneNeighbors.begin(); it2 != CxMinusOneNeighbors.end(); it2++)
            {
              MolecularModeling::Atom* thisCxMinusOneNeighbor = (*it2);
              int thisCxMinusOneNeighborNum = std::atoi(&thisCxMinusOneNeighbor->GetId()[1]);
              if((thisCxMinusOneNeighborNum = (x - 2)) && (thisCxMinusOneNeighbor->GetElementSymbol() == "C"))
              {//This is Cx-2
                //Add Cx-2
                linkageAtoms.push_back(thisCxMinusOneNeighbor);
                setCx2 = true;
              }
            }
          }
        }
      }
      if(setCx1 != true)//should this break/exit?
      {//Could not ID Cx-1
        //TODO Update once GLOBAL_DEBUG is in place
        if(local_debug > 0)
        {//Print Error
          ss.str("");
          ss << "Unable to find Cx-1 of the reducing monosaccharide: ";
          ss << reducing_mono_->residue_name_ << ".";
          gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
          ss.str("");
        }
      }
      if(setCx2 !=true)//should this break/exit?
      {//Could not ID Cx-2
        //TODO Update logic for omega
        if((linkage_type_ == "1-6") ||
          (inverse_linkage_type_ == "1-6")||
          (linkage_type_ == "2-6") ||
          (inverse_linkage_type_ == "2-6"))
        {
          //TODO Update once GLOBAL_DEBUG is in place
          if(local_debug > 0)
          {//Print Error
            ss.str("");
            ss << "Unable to find Cx-2 of the reducing monosaccharide: ";
            ss << reducing_mono_->residue_name_ << ".";
            gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
            ss.str("");
          }
        }
      }
    }
    else
    {//Anomeric-anomeric linkage uses different atoms
      /*Using angle definitions provided by Dr. Al French.
        The Ring O of non_reducing_mono_2_ is used instead of Cx-1
      */
      for(std::vector<MolecularModeling::Atom*>::iterator it = non_reducing_mono_2_->cycle_atoms_.begin(); it != non_reducing_mono_2_->cycle_atoms_.end(); it++)
      {
        MolecularModeling::Atom* thisAtom = (*it);
        if(thisAtom->GetElementSymbol() == "O")
        {//Then this is the ring oxygen
          //Add Cx-1 (I know this is poor naming)
          linkageAtoms.push_back(thisAtom);
          setCx1 = true;
        }
      }
      phi_angle_ = CalculatePhiAngle(linkageAtoms);
      //This is going to calculate phi' even though its using the Psi function
      //because we changed what atoms are being fed into it (ring O instead of Cx-1)
      phi_prime_angle_ = CalculatePsiAngle(linkageAtoms);
    }

    //Calculate Torsion Angles
    if(setRingO && setAnomericC && setGlycosidicO && setCx && setCx1)
    {
      phi_angle_ = CalculatePhiAngle(linkageAtoms);

      if(!anomeric_anomeric_linkage_)
      {
        psi_angle_ = CalculatePsiAngle(linkageAtoms);
      }
      else
      {//anomeric-anomeric
        phi_prime_angle_ = CalculatePsiAngle(linkageAtoms);
      }

      //TODO Update logic for omega
      if(((linkage_type_ == "1-6") ||
        (inverse_linkage_type_ == "1-6")||
        (linkage_type_ == "2-6") ||
        (inverse_linkage_type_ == "2-6")))
      {
        if(setCx2)
        {
          omega_angle_ = CalculateOmegaAngle(linkageAtoms);
        }
        else
        {
          //TODO Update once GLOBAL_DEBUG is in place
          if(local_debug > 0)
          {//Print Error
            ss.str("");
            ss << "Unable to calculate the omega angle. Missing Cx-2 of: ";
            ss << reducing_mono_->residue_name_ << ".";
            gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
            ss.str("");
          }
        }
      }
    }
    else
    {//Problem setting linkage atoms
      //TODO Update once GLOBAL_DEBUG is in place
      if(local_debug > 0)
      {//Print Error
        ss.str("");
        ss << "Unable to calculate the glycosidic torsion angles. ";
        ss << "Atoms needed for calculations could not be identified.";
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        ss.str("");
      }
    }
  

    /*Checking to see if angles are negative/weird
    Needed for CHI Energy Functions
    For now they are expected not to be negative
    (unless they are the default of -9999, which
    means something went wrong),
    but if that changes this logic needs to be updated.
    This will also add 360 to any value where
    -360 <= angle < 0,
    and output an error if angle < -360 or angle > 360
    */
    if(phi_angle_ != -9999)
    {
      if((0 <= phi_angle_) && (phi_angle_ < 360))
      {//normal
        if(local_debug > 0)
        {
          ss.str("");//clear stringstream just in case
          ss << "Phi angle of " << phi_angle_ << " falls into expected range (0-360 degrees)";
          gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
          ss.str("");//clean up stringstream
        }
      }
      else if((-360 <= phi_angle_) && (phi_angle_ < 0))
      {//add 360 degrees
        phi_angle_ = phi_angle_ + 360;
      }

      else if(phi_angle_ < -360)
      {//Warning
        if(local_debug > 0)
        {
          ss.str("");//clear stringstream just in case
          ss << "Phi angle of " << phi_angle_ << " is less than -360 degrees";
          gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
          ss.str("");//clean up stringstream
        }
      }
      else if(phi_angle_ > 360)
      {//Warning
        if(local_debug > 0)
        {
          ss.str("");//clear stringstream just in case
          ss << "Phi angle of " << phi_angle_ << " is greater than 360 degrees";
          gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
          ss.str("");//clean up stringstream
        }
      }
    }
    else if(local_debug > 0)
    {//Print error
      ss.str("");//clear stringstream just in case
      ss << "Phi angle was not calculated properly";
      gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
      ss.str("");//clean up stringstream
    }

    if(psi_angle_ != -9999)
    {
      if((0 <= psi_angle_) && (psi_angle_ < 360))
      {//normal
        if(local_debug > 0)
        {
          ss.str("");//clear stringstream just in case
          ss << "Psi angle of " << psi_angle_ << " falls into expected range (0-360 degrees)";
          gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
          ss.str("");//clean up stringstream
        }
      }
      else if((-360 <= psi_angle_) && (psi_angle_ < 0))
      {//add 360 degrees
        psi_angle_ = psi_angle_ + 360;
      }

      else if(psi_angle_ < -360)
      {//Warning
        if(local_debug > 0)
        {
          ss.str("");//clear stringstream just in case
          ss << "Psi angle of " << psi_angle_ << " is less than -360 degrees";
          gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
          ss.str("");//clean up stringstream
        }
      }
      else if(psi_angle_ > 360)
      {//Warning
        if(local_debug > 0)
        {
          ss.str("");//clear stringstream just in case
          ss << "Psi angle of " << psi_angle_ << " is greater than 360 degrees";
          gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
          ss.str("");//clean up stringstream
        }
      }
    }
    else if(local_debug > 0)
    {//Print error
      ss.str("");//clear stringstream just in case
      ss << "Psi angle was not calculated properly";
      gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
      ss.str("");//clean up stringstream
    }

    if(omega_angle_ != -9999)
    {
      if((0 <= omega_angle_) && (omega_angle_ < 360))
      {//normal
        if(local_debug > 0)
        {
          ss.str("");//clear stringstream just in case
          ss << "Omega angle of " << omega_angle_ << " falls into expected range (0-360 degrees)";
          gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
          ss.str("");//clean up stringstream
        }
      }
      else if((-360 <= omega_angle_) && (omega_angle_ < 0))
      {//add 360 degrees
        omega_angle_ = omega_angle_ + 360;
      }

      else if(omega_angle_ < -360)
      {//Warning
        if(local_debug > 0)
        {
          ss.str("");//clear stringstream just in case
          ss << "Omega angle of " << omega_angle_ << " is less than -360 degrees";
          gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
          ss.str("");//clean up stringstream
        }
      }
      else if(omega_angle_ > 360)
      {//Warning
        if(local_debug > 0)
        {
          ss.str("");//clear stringstream just in case
          ss << "Omega angle of " << omega_angle_ << " is greater than 360 degrees";
          gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
          ss.str("");//clean up stringstream
        }
      }
    }
    else if((local_debug > 0))
    {//Print error
      ss.str("");//clear stringstream just in case
      ss << "Omega angle was not calculated properly";
      gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
      ss.str("");//clean up stringstream
    }


    /* TODO Add Aglycone Torsions
    Add code to determine the glycosidic linkage and torsion angles between
    the sugar and aglycone, especially if it is a protein sidechain.

    For now, we don't consider this a glycosidic linkage at all, as the
    constructor takes monosaccharide objects

    Also, determine the correct definitions (atom selection) of these angles for
    N-linked, O-linked, etc.  It will likely change by amino acid.
    */

    //Get the orientation of the hydroxy involved in the linkage
    hydroxyl_configuration_ = determineLinkageConfiguration();
    if(hydroxyl_configuration_ == "Unknown")
    {
      if(local_debug > 0)
      {//Print error
        ss.str("");//clear stringstream just in case
        ss << "Linkage hydroxyl configuration unknown.";
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        ss.str("");//clean up stringstream
      }
    }
    else if((hydroxyl_configuration_ != "axial") && (hydroxyl_configuration_ != "equatorial"))
    {
      if(local_debug > 0)
      {//Print error
        ss.str("");//clear stringstream just in case
        ss << "Something went horribly wrong determining the linkage hydroxyl configuration.";
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        ss.str("");//clean up stringstream
      }
    }

    /////////////////////////////
    // CHI Energy Calculations //
    /////////////////////////////
    /*CHI Energy function caveats
      •Developed for rings in 1C4 or 4C1 chair conformations
          As I don't fully understand/trust the BFMP code, and I'm
        curious about how ring shape affects CHI Energy, I will be
        runnning CHI Energy functions regardless of ring shape.
          When Cremer-Pople ring puckering is added to GMML, this may
        change to only calculate on chair forms (defined by either BFMP
        or Cremer-Pople) depending on the results I get for non-ring
        CHI Energies and what Rob thinks.

      •Developed for a subset of possible linkage types
        These are checked in the CHI Energy functions, so this is more
      of an FYI.
    */

    if(local_debug > 0)
    {//check ring conformation and write logs
      if(((non_reducing_mono_->bfmp_ring_conformation_.find("1C4")
            != std::string::npos) ||
          (non_reducing_mono_->bfmp_ring_conformation_.find("4C1")
            != std::string::npos)) &&
        ((reducing_mono_->bfmp_ring_conformation_.find("1C4")
            != std::string::npos) ||
          (reducing_mono_->bfmp_ring_conformation_.find("4C1")
            != std::string::npos)))
      {//Both are in chair form
        ss.str("");//clear stringstream just in case
        ss << "Both monosaccharide rings are in chair form.";
        gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
        ss.str("");//clean up stringstream
      }
      else
      {//Print warning
        ss.str("");//clear stringstream just in case
        ss << "CHI Energy function not developed for non chair rings.";
        gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
        ss.str("");//clean up stringstream
      }
    }

    if((!anomeric_anomeric_linkage_) &&
      (non_reducing_mono_->sugar_name_.ring_type_ == "P") &&
      (reducing_mono_->sugar_name_.ring_type_ == "P"))
    {
      if((phi_angle_ >= 0) && (phi_angle_ < 360))
      {
        phi_CHI_Energy_ = CalculatePhiChiEnergy();
      }
      else if(local_debug > 0)
      {//Warning
        ss.str("");
        ss << "phi_angle_ value of " << phi_angle_ << ".  ";
        ss << "Something went wrong";
        gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
        ss.str("");
      }

      if((psi_angle_ >= 0) && (psi_angle_ < 360))
      {
        psi_CHI_Energy_   = CalculatePsiChiEnergy();
      }
      else if(local_debug > 0)
      {//Warning
        ss.str("");
        ss << "psi_angle_ value of " << psi_angle_ << ".  ";
        ss << "Something went wrong";
        gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
        ss.str("");
      }

      //TODO Update logic for omega
      if(((linkage_type_ == "1-6") ||
        (inverse_linkage_type_ == "1-6")||
        (linkage_type_ == "2-6") ||
        (inverse_linkage_type_ == "2-6")))
      {
        if((omega_angle_ >= 0) && (omega_angle_ < 360))
        {
          omega_CHI_Energy_ = CalculateOmegaChiEnergy();
        }
        else if(local_debug > 0)
        {//Warning
          ss.str("");
          ss << "omega_angle_ value of " << omega_angle_ << ".  ";
          ss << "Something went wrong";
          gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
          ss.str("");
        }
      }
    }
    else if(local_debug > 0)
    {//Print warning
      ss.str("");//clear stringstream just in case
      ss << "CHI Energy function not developed for ";
      ss << "this linkage.";
      gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
      ss.str("");//clean up stringstream
    }

    if(local_debug > 0)
    {//Print info
      ss.str("");//clear stringstream just in case
      ss << "CHI Energies: Phi = " << phi_CHI_Energy_ << ", ";
      ss << "Psi = " << psi_CHI_Energy_ << ", ";
      ss << "Omega = " << omega_CHI_Energy_ << ". ";
      gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
      ss.str("");//clean up stringstream
      ss << "CHI Energy functions: Phi = " << phi_CHI_function_ << ", ";
      ss << "Psi = " << psi_CHI_function_ << ", ";
      ss << "Omega = " << omega_CHI_function_ << ". ";
      gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
      ss.str("");//clean up stringstream
    }
    if((phi_CHI_Energy_ == -1)||(psi_CHI_Energy_ == -1)||
      ((omega_CHI_Energy_ == -1) && (omega_angle_ != -9999)))
    {
      if(local_debug > 0)
      {//Print Warnings
        if(phi_CHI_Energy_ == -1)
        {//Print Warning
          ss << "Something went wrong calculating Phi CHI Energy.";
          gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
          ss.str("");//clean up stringstream
        }
        if(psi_CHI_Energy_ == -1)
        {//Print Warning
          ss << "Something went wrong calculating Psi CHI Energy.";
          gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
          ss.str("");//clean up stringstream
        }
        if((omega_CHI_Energy_ == -1) && (omega_angle_ != -9999))
        {//Print Warning
          ss << "Something went wrong calculating Omega CHI Energy.";
          gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
          ss.str("");//clean up stringstream
        }
        //one of them is wrong so total is invalid
        //Print error
        ss << "Something went wrong calculating the total CHI Energy.";
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        ss.str("");//clean up stringstream
      }
    }
    else
    {
      total_CHI_Energy_ = phi_CHI_Energy_ + psi_CHI_Energy_;
      if((omega_angle_ != -9999))
      {
        total_CHI_Energy_ = total_CHI_Energy_ + omega_CHI_Energy_;
      }
      if(local_debug > 0)
      {//Print info
      ss.str("");//clear stringstream just in case
      ss << "Total CHI Energy for " << linkage_name_;
      ss << " is: " <<total_CHI_Energy_;
      gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
      ss.str("");//clean up stringstream
      }
    }
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
      ss.str("");
      ss << "Atoms being used are: ";
      ss << linkageAtoms[0]->GetId() << ", " << linkageAtoms[1]->GetId() << ", ";
      ss << linkageAtoms[2]->GetId() << ", and " << linkageAtoms[3]->GetId();
      gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
      ss.str("");
    }
  }
  if(linkageAtoms.size() < 4)
  {//Print error
    gmml::log(__LINE__, __FILE__, gmml::ERR, "Not enough atoms in vector to calculate Phi");
    Phi = -9999;
    return Phi;
  }

  /* Replaced to determine atoms seperately from calculating Phi
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
  */
  /*Replaced as we are trying to use one function for torsions
  //non_reducing_mono_->assembly_->CalculateTorsionAngleByAtoms(linkageAtoms[0],linkageAtoms[1],linkageAtoms[2],linkageAtoms[3]);
  //Phi = gmml::ConvertRadian2Degree(Phi);
  */

  Phi = GeometryTopology::CalculateDihedralAngle(
    linkageAtoms[0]->GetCoordinate(),
    linkageAtoms[1]->GetCoordinate(),
    linkageAtoms[2]->GetCoordinate(),
    linkageAtoms[3]->GetCoordinate());


  if (local_debug > 0)
  {//Output Phi value
    ss.str("");
    ss << "Phi: " << Phi;
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
    ss.str("");
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
      ss.str("");
      ss << "Atoms being used are: ";
      ss << linkageAtoms[1]->GetId() << ", " << linkageAtoms[2]->GetId() << ", ";
      ss << linkageAtoms[3]->GetId() << ", and " << linkageAtoms[4]->GetId();
      gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
      ss.str("");
    }
  }
  if(linkageAtoms.size() < 5)
  {//Print error
    gmml::log(__LINE__, __FILE__, gmml::ERR, "Not enough atoms in vector to calculate Psi");
    Psi = -9999;
    return Psi;
  }

  /* Replaced to determine atoms seperately from calculating Psi
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
  */

  Psi = GeometryTopology::CalculateDihedralAngle(
    linkageAtoms[1]->GetCoordinate(),
    linkageAtoms[2]->GetCoordinate(),
    linkageAtoms[3]->GetCoordinate(),
    linkageAtoms[4]->GetCoordinate());

  if (local_debug > 0)
  {//Output Psi value
    ss.str("");
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
      ss.str("");
      ss << "Atoms being used are: ";
      ss << linkageAtoms[2]->GetId() << ", " << linkageAtoms[3]->GetId() << ", ";
      ss << linkageAtoms[4]->GetId() << ", and " << linkageAtoms[5]->GetId();
      gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
      ss.str("");
    }
  }
  if(linkageAtoms.size() < 6)
  {//Print error
    gmml::log(__LINE__, __FILE__, gmml::ERR, "Not enough atoms in vector to calculate Omega");
    Omega = -9999;
    return Omega;
  }

  /* Replaced to determine atoms seperately from calculating Omega
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
  */

  Omega = GeometryTopology::CalculateDihedralAngle(
    linkageAtoms[2]->GetCoordinate(),
    linkageAtoms[3]->GetCoordinate(),
    linkageAtoms[4]->GetCoordinate(),
    linkageAtoms[5]->GetCoordinate());

    if (local_debug > 0)
    {//Print Omega value
      ss.str("");
      ss << "Omega: " << Omega;
      gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
      ss.str("");
    }
  return Omega;
}

/*TODO documentation & warning/error reporting
    Instead of (or in addition to) writing to the gmml::log,
    warnings and errors should be generated that have a sensible
    way to identify the linkage to be able to show which linkages
    in an oligosaccharide have issues on GlyFinder's Web interface

    The code would be something like
    if(parent_oligo_ != NULL)
    {
      Glycan::Note* linkage_note = new Glycan::Note();
      linkage_note->type_ = Glycan::ERROR;
      linkage_note->category_ = Glycan::GLYCOSIDIC;
      ss.str("");
      ss << linkage_name_ << " (";//linkage_name_ not yet determined
      ss << detailed_linkage_name_ << ") ";//detailed_linkage_name_ not yet determined
      ss << "has some horrific error here";
      linkage_note->description_ = ss.str();
      parent_oligo_->AddNote(thisNote);
      ss.str("");
    }
    Or really this class should have std::vector<Note*> glycosidic_linkage_notes_;
*/

double Glycan::GlycosidicLinkage::CalculatePhiChiEnergy()
{
  /////////////////////////////////////////////////////////////////////
  //CHI_Energy function taken from https://doi.org/10.1002/jcc.23517 //
  /////////////////////////////////////////////////////////////////////
  /*CHI_Energy function equation
        N               (x - b_i)^2
  f(x)=  Σ  a_i * e * -(――――――――――――――) + Offset
       i=1                  c_i
  */

  int local_debug = -1;
  std::stringstream ss;//for debugging & logs
  double phiCHI;

  //There are 4 different energy functions for CHI Phi Energy
  bool phi_alpha = false;
  bool phi_beta = false;
  bool phi_2_3 = false;
  bool phi_2_6 = false;

  /*L-sugars note from paper
      "To apply the energy curves shown in Figure 5 to linkages containing
    L-sugars, it is simply necessary to use the mirror images of the
    relevant energy curve."

      According to Rob: "It only applies to the Phi angle, and so only applies to the non-reducing sugar in the linkage."
  */

  /*Phi CHI Energy uses -180 to 180
    For whatever reason, instead of 0 to 360 degrees
    like the other energy funtions, the CHI Energy
    function for Phi uses -180 to 180.

    Since there are checks for Phi angles before this
    function is called, it should already be 0 to 360 degrees
  */
  double adjustedPhi;
  if(phi_angle_ > 180.0)
  {//Make it negative
    adjustedPhi = phi_angle_ - 360;
    if(local_debug > 0)
    {//Print info
      ss.str("");
      ss << "Phi angle is: ";
      ss << phi_angle_;
      ss << ". As CHI Energy from Phi goes from -180 to 180 degrees,";
      ss << "phi was adjusted to: ";
      ss << adjustedPhi << " for this function only.";
      gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
      ss.str("");//clear stringstream for other messages
    }
  }
  else
  {//don't adjust phi
    adjustedPhi = phi_angle_;
  }

  double phi;
  if(non_reducing_mono_->sugar_name_.isomer_ == "D")
  {//D non-reducing Sugar
    phi = adjustedPhi;
  }
  else if(non_reducing_mono_->sugar_name_.isomer_ == "L")
  {//L non-reducing Sugar, use inverse function
    phi = adjustedPhi * -1;
  }
  else
  {//Isomer unknown
    gmml::log(__LINE__, __FILE__, gmml::ERR, "Unable to determine the isomer (D/L) of the non-reducing monosaccharide involved in the linkage.  Exiting CalculatePhiChiEnergy()");
    return 0;
  }

  ////////////////////////////////////////////
  //Logic to determine which function to use//
  ////////////////////////////////////////////

  /*TODO Energy Function Weighting
      Rob wants the energy function to be weighted based on ring shape
    and isomer.  IE 4C1 alpha thats 75% 4C1 and 25% 1C4 (method to
    determine this TBA), then CHI = 0.75 alpha and 0.25 inverse beta
    because 1C4 D alpha is like 4C1 L beta

    This means the logic for determining which function to use will get
    complicated.

    For now, I will only worry about 4C1 and 1C4
    DM 1-31-22
  */

  //CHI Energy functions have only been created for certain linkages
  //The following list applies for alpha and beta linkages
  std::list<std::string> ChiLinkTypes = { "1-1", "1-2", "1-3",
                                                "1-4", "1-6" };
  if((find(ChiLinkTypes.begin(), ChiLinkTypes.end(),
          linkage_type_) != ChiLinkTypes.end()) ||
     (find(ChiLinkTypes.begin(), ChiLinkTypes.end(),
          inverse_linkage_type_) != ChiLinkTypes.end()))
  {
    if(anomeric_configuration_ == "a")
    {
      phi_alpha = true;
      phi_CHI_function_ = 1;
      if(local_debug > 0)
      {//Print Info
        ss.str("");
        ss << "The Phi CHI Energy function for alpha linkages (#1) was used.";
        gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
        ss.str("");
      }
    }
    else if(anomeric_configuration_ == "b")
    {
      phi_beta = true;
      phi_CHI_function_ = 2;
      if(local_debug > 0)
      {//Print Info
        ss.str("");
        ss << "The Phi CHI Energy function for beta linkages (#2) was used.";
        gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
        ss.str("");
      }
    }
    else
    {
      gmml::log(__LINE__, __FILE__, gmml::ERR, "Unable to determine the configuration (α/β) of the linkage. Exiting CalculatePhiChiEnergy()");
      return 0;
    }
  }
  //TODO: Figure out if inverse_linkage_type_ lookup is needed.
  else if((linkage_type_ == "2-3") || (inverse_linkage_type_ == "2-3"))
  {
    phi_2_3 = true;
    phi_CHI_function_ = 3;
    if(local_debug > 0)
    {//Print Info
      ss.str("");
      ss << "The Phi CHI Energy function for 2-3 linkages (#3) was used.";
      gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
      ss.str("");
    }
  }
  else if((linkage_type_ == "2-6") || (inverse_linkage_type_ == "2-6"))
  {
    phi_2_6 = true;
    phi_CHI_function_ = 4;
    if(local_debug > 0)
    {//Print Info
      ss.str("");
      ss << "The Phi CHI Energy function for 2-6 linkages (#4) was used.";
      gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
      ss.str("");
    }
  }
  //TODO: Add other Phi CHI Energy functions if/when they are developed
  else
  {//Warning
    if(local_debug > 0)
    {//Print warning
      ss.str("");
      ss << "No CHI Energy Function exists for the Phi angle of ";
      ss << linkage_type_;
      ss << " linkages. Exiting CalculatePhiChiEnergy()";
      gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
      ss.str("");//clear stringstream for other messages
    }
    return 0;
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
    /*Offset so small it's not even in Totx formula
    double Off = 1.00501e-30;
    */
    double Leftx, Rightx, ax, bx, cx, dx, ex, fx;

    Leftx  = LH * exp(-pow((phi-(Lc)),2.0)/LW);
    Rightx = RH * exp(-pow((phi-(Rc)),2.0)/RW);
    ax     = aH * exp(-pow((phi-(ac)),2.0)/aW);
    bx     = bH * exp(-pow((phi-(bc)),2.0)/bW);
    cx     = cH * exp(-pow((phi-(cc)),2.0)/cW);
    dx     = dH * exp(-pow((phi-(dc)),2.0)/dW);
    ex     = eH * exp(-pow((phi-(ec)),2.0)/eW);
    fx     = fH * exp(-pow((phi-(fc)),2.0)/fW);
    phiCHI = Rightx + Leftx + ax + bx + cx + dx + ex + fx;
    return phiCHI;
  }
  else if (phi_beta)
  {
    double LH = 450.540038600828,
           Lc = -330.769995527134,
           LW = 4449.7622241787;
    double RH = 23.7118506901333,
           Rc = 304.624980492529,
           RW = 8375.1929028027;
    double aH = 5.93533323829663,
           ac = -152.080139620062,
           aW = 6049.77220005964;
    double bH = 22.467372096061,
           bc = -23.5159916173247,
           bW = 606.89715970453;
    double cH = 10.0360057033439,
           cc = 120.962836525241,
           cW = 4037.89330459304;
    double dH = -18.1406478565247,
           dc = -24.2677756921736,
           dW = 543.050986049266;
    double eH = 5.88226333077368,
           ec = 19.6321032903376,
           eW = 897.92664572344;
    /*Off = -2.27829251796721; Leaving this here just in case.
    Not sure why it was in the original CHI_Energy code from VinaCarb
    */
    double Off = -2.1283;
    double Leftx, Rightx, ax, bx, cx, dx, ex;

    Leftx  = LH * exp(-pow((phi-(Lc)),2.0)/LW);
    Rightx = RH * exp(-pow((phi-(Rc)),2.0)/RW);
    ax     = aH * exp(-pow((phi-(ac)),2.0)/aW);
    bx     = bH * exp(-pow((phi-(bc)),2.0)/bW);
    cx     = cH * exp(-pow((phi-(cc)),2.0)/cW);
    dx     = dH * exp(-pow((phi-(dc)),2.0)/dW);
    ex     = eH * exp(-pow((phi-(ec)),2.0)/eW);
    phiCHI = Rightx + Leftx + ax + bx + cx + dx + ex + Off;
    return phiCHI;
  }

  /*NOTE for 2-3/6 linkages
  In the SI for Vina Carb
  (https://pubs.acs.org/doi/suppl/10.1021/acs.jctc.5b00834/suppl_file/ct5b00834_si_002.pdf)
  The ranges for Phi of 2-3 and 2-6 linkages are defined as:
  (The range values had typos but I checked the figures)

                0.0018(x-60)^2,          if x ϵ [0,120]
  ΔE_Sia23 (x) = 0.0018(x-180)^2 + 1.95,  if x ϵ [120,240]
                0.0018(x-300)^2 + 0.385, if x ϵ [240,360]

                0.0018(x-60)^2,          if x ϵ [0,120]
  ΔE_Sia26 (x) = 0.0018(x-180)^2 + 1.19,  if x ϵ [120,240]
                0.0018(x-300)^2 + 1.41,  if x ϵ [240,360]

  However, the ranges overlap at 120, 240, and 360/0.
  The values are pretty close, so I'm just going to use
  if x ϵ [0,120]
  if x ϵ (120,240]
  if x ϵ (240,360]
  */

  else if (phi_2_3)
  {
    //2-3 & 2-6 functions use 0-360, not -180 - 180;
    if(non_reducing_mono_->sugar_name_.isomer_ == "L")
    {//L non-reducing Sugar, use inverse function
      phi = phi_angle_ * -1;
    }
    else
    {
      phi = phi_angle_;
    }
    if((0 <= phi) && (phi <= 120))
    {//if x ϵ [0,120]
      phiCHI = 0.0018 * pow( (phi - 60), 2);
      return phiCHI;
    }
    else if((120 < phi) && (phi <= 240))
    {//if x ϵ (120,240]
      phiCHI = 0.0018 * pow( (phi - 180), 2) + 1.95;
      return phiCHI;
    }
    else if((240 < phi) && (phi <= 360))
    {//if x ϵ (240,360]
      phiCHI = 0.0018 * pow( (phi - 300), 2) + 0.385;
      return phiCHI;
    }
    else
    {
      if(local_debug > 0)
      {//Print Error
        ss.str("");
        ss << "Phi angle of " << phi;
        ss << " outside of normal range (0-360).";
        ss << "Unable to run CHI Energy analysis";
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        ss.str("");
      }
      return 0;
    }
  }
  else if (phi_2_6)
  {
    //2-3 & 2-6 functions use 0-360, not -180 - 180;
    if(non_reducing_mono_->sugar_name_.isomer_ == "L")
    {//L non-reducing Sugar, use inverse function
      phi = phi_angle_ * -1;
    }
    else
    {
      phi = phi_angle_;
    }

    if((0 <= phi ) && (phi <= 120))
    {//if x ϵ [0,120]
      phiCHI = 0.0018 * pow( (phi - 60), 2);
      return phiCHI;
    }
    else if((120 < phi) && (phi <= 240))
    {//if x ϵ (120,240]
      phiCHI = 0.0018 * pow( (phi - 180), 2) + 1.19;
      return phiCHI;
    }
    else if((240 < phi) && (phi <= 360))
    {//if x ϵ (240,360]
      phiCHI = 0.0018 * pow( (phi - 300), 2) + 1.41;
      return phiCHI;
    }
    else
    {
      if(local_debug > 0)
      {//Print Error
        ss.str("");
        ss << "Phi angle of " << phi;
        ss << " outside of normal range (0-360).";
        ss << "Unable to run CHI Energy analysis";
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        ss.str("");
      }
      return 0;
    }
  }
  else
  {//something went horribly wrong
    //there is a check for this after setting booleans for
    //which function to use, so it should never get here.
    if(local_debug > 0)
    {
      ss.str("");
      ss << "No CHI Energy Function exists for the Phi angle of ";
      ss << linkage_type_;
      ss << " linkages.";
      gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
      ss.str("");//clear stringstream for other messages
    }
    return 0;
  }
  return 0;
}

double Glycan::GlycosidicLinkage::CalculatePsiChiEnergy()
{
  //////////////////////////////////////////////////////////////////////
  // CHI_Energy function taken from https://doi.org/10.1002/jcc.23517 //
  //////////////////////////////////////////////////////////////////////
  //           N               (x - b_i)^2                            //
  //    f(x)=  Σ  a_i * e * -(――――――――――――――) + Offset                //
  //          i=1                  c_i                                //
  //////////////////////////////////////////////////////////////////////
  //    For 1-6, 2-3, & 2-6 energy function details see the SI in     //
  //          https://doi.org/10.1021/acs.jctc.5b00834                //
  //////////////////////////////////////////////////////////////////////

  int local_debug = -1;
  std::stringstream ss;//for debugging & logs
  double psiCHI;
  /* There are 4 different energy functions for CHI Psi Energy
  //each with their own criteria
  bool psi_CHI_1 = false; //1-2 & 1-4 axial, 1-3 & 2-3 equatorial
  bool psi_CHI_2 = false; //1-2 & 1-4 equatorial, 1-3 & 2-3 axial
  bool psi_CHI_3 = false; //1-6 & 2-6 alpha
  bool psi_CHI_4 = false; //1-6 & 2-6 beta
  */

  psi_CHI_function_ = pickPsiChiEnergyFunction();

  ////////////////////////
  //CHI Energy Functions//
  ////////////////////////
  if(psi_CHI_function_ == -1)
  {//Warning & Exit
    if(local_debug > 0)
    {
      ss.str("");
      ss << "Unable to identify the correct Psi CHI function.";
      gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
    }
    return -1.0;
  }
  else if(psi_CHI_function_ == 1)//2A3E in VinaCarb
  {
    double  LH = 4.62366277694845,
            Lc = 5.045583934308,
            LW = 5005.75236060956;
    double  RH = 4.61387450239844,
            Rc = 362.487847702007,
            RW = 2090.63190217702;
    double  aH = 4.94191727813274,
            ac = 121.202321824468,
            aW = 2093.75214491931;
    double  bH = 0.402901504292045,
            bc = 241.428583877882,
            bW = 456.828754790442;
    double  cH = 0.798876573705798,
            cc = 68.425080241155,
            cW = 678.807178379645;
    double  dH = 0.222992242737354,
            dc = 192.925748017071,
            dW = 347.244734136509;
    double  Off = -0.125645118474882;
    double Leftx, Rightx, ax, bx, cx, dx;

    Leftx  = LH * exp(-pow((psi_angle_-(Lc)),2.0)/LW);
    Rightx = RH * exp(-pow((psi_angle_-(Rc)),2.0)/RW);
    ax     = aH * exp(-pow((psi_angle_-(ac)),2.0)/aW);
    bx     = bH * exp(-pow((psi_angle_-(bc)),2.0)/bW);
    cx     = cH * exp(-pow((psi_angle_-(cc)),2.0)/cW);
    dx     = dH * exp(-pow((psi_angle_-(dc)),2.0)/dW);
    psiCHI = Rightx + Leftx + ax + bx + cx + dx + Off;
    return psiCHI;
  }
  else if(psi_CHI_function_ == 2)//2E3A in VinaCarb
  {
    double LH = 4.46811874171788,
           Lc = 1e-30,
           LW = 1279.58772056199;
    double RH = 4.38204018882823,
           Rc = 357.770654336205,
           RW = 6050.14162479438;
    double aH = 284.944948778136,
           ac = 146.644068129462,
           aW = 1551.75673776163;
    double bH = 4.76134025362478,
           bc = 220.683118921686,
           bW = 5892.94143218231;
    double cH = -169.197666368856,
           cc = 147.370828680332,
           cW = 1742.47541063603;
    double dH = -118.440552792375,
           dc = 146.05660843428,
           dW = 1359.82873591396;
    double Off = 1.0219924486158;
    double Leftx, Rightx, ax, bx, cx, dx;

    Leftx  = LH * exp(-pow((psi_angle_-(Lc)),2.0)/LW);
    Rightx = RH * exp(-pow((psi_angle_-(Rc)),2.0)/RW);
    ax     = aH * exp(-pow((psi_angle_-(ac)),2.0)/aW);
    bx     = bH * exp(-pow((psi_angle_-(bc)),2.0)/bW);
    cx     = cH * exp(-pow((psi_angle_-(cc)),2.0)/cW);
    dx     = dH * exp(-pow((psi_angle_-(dc)),2.0)/dW);
    psiCHI = Rightx + Leftx + ax + bx + cx + dx + Off;
    return psiCHI;
  }
  else if(psi_CHI_function_ == 3)//6A in VinaCarb
  {
    double aH = 67.9431348410598,
           ac = -59.5393395706705,
           aW = 993.323581145538;
    double bH = 6.13421142432396,
           bc = 10.4786088782815,
           bW = 945.770771330812;
    double cH = 3.27628727235978,
           cc = 54.2960678151208,
           cW = 851.528141801851;
    double dH = 0.727486729062442,
           dc = 131.067737803489,
           dW = 1037.41211378392;
    double eH = 2.57362265878937,
           ec = 245.102891425541,
           eW = 2012.99451568206;
    double fH = 5.75995973448166,
           fc = 359.999988549478,
           fW = 1153.3974275673;
    double gH = 3.47492643928157,
           gc = 321.677942414686,
           gW = 2080.97053159226;
    double hH = -0.741000462200939,
           hc = 199.106903524814,
           hW = 522.180434119001;
    double ax, bx, cx, dx, ex, fx, gx, hx;

    ax = aH * exp(-pow((psi_angle_-(ac)),2.0)/aW);
    bx = bH * exp(-pow((psi_angle_-(bc)),2.0)/bW);
    cx = cH * exp(-pow((psi_angle_-(cc)),2.0)/cW);
    dx = dH * exp(-pow((psi_angle_-(dc)),2.0)/dW);
    ex = eH * exp(-pow((psi_angle_-(ec)),2.0)/eW);
    fx = fH * exp(-pow((psi_angle_-(fc)),2.0)/fW);
    gx = gH * exp(-pow((psi_angle_-(gc)),2.0)/gW);
    hx = hH * exp(-pow((psi_angle_-(hc)),2.0)/hW);
    psiCHI = ax + bx + cx + dx + ex + fx + gx + hx;
    return psiCHI;
  }
  else if(psi_CHI_function_ == 4)//6E in VinaCarb
  {
    double aH = 7.24858655753829,
           ac = 3.60600554520403,
           aW = 2459.23916629141;
    double bH = 1.9,
           bc = 96.5930821702371,
           bW = 2683.88656516991;
    double cH = 0.741022592342903,
           cc = 141.663521919709,
           cW = 1150.04756181103;
    double dH = 0.2,
           dc = 162,
           dW = 400;
    double eH = 0.287090039932611,
           ec = 228.171314273305,
           eW = 272.201363844744;
    double fH = 1.22591971967808,
           fc = 292.206221787048,
           fW = 1134.52455512381;
    double gH = 7.41063235334191,
           gc = 369.867701147817,
           gW = 3499.15994772992;
    double hH = -0.61489499584011,
           hc = 271.319024293053,
           hW = 532.437194483944;
    double iH = -0.35,
           ic = 183,
           iW = 100;
    double ax, bx, cx, dx, ex, fx, gx, hx, ix;

    ax = aH * exp(-pow((psi_angle_-(ac)),2.0)/aW);
    bx = bH * exp(-pow((psi_angle_-(bc)),2.0)/bW);
    cx = cH * exp(-pow((psi_angle_-(cc)),2.0)/cW);
    dx = dH * exp(-pow((psi_angle_-(dc)),2.0)/dW);
    ex = eH * exp(-pow((psi_angle_-(ec)),2.0)/eW);
    fx = fH * exp(-pow((psi_angle_-(fc)),2.0)/fW);
    gx = gH * exp(-pow((psi_angle_-(gc)),2.0)/gW);
    hx = hH * exp(-pow((psi_angle_-(hc)),2.0)/hW);
    ix = iH * exp(-pow((psi_angle_-(ic)),2.0)/iW);
    psiCHI = ax + bx + cx + dx + ex + fx + gx + hx + ix;
    return psiCHI;
  }
  else
  {//ERROR
    ss.str("");
    ss << "Something went horribly wrong and the code should ";
    ss << "not have reached this point.";
    gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
    return 0;
  }
  return 0;
}

int Glycan::GlycosidicLinkage::pickPsiChiEnergyFunction()
{
  /* Return value for CHI Energy Functions
   -1:    Something went horribly wrong or this linkage
            does not have an energy function.
    1:    1-2 & 1-4 axial, 1-3 & 2-3 equatorial
    2:    1-2 & 1-4 equatorial, 1-3 & 2-3 axial
    3:    1-6 & 2-6 alpha
    4:    1-6 & 2-6 beta
  */

  int local_debug = -1;
  std::stringstream ss;//for debugging & logs
  std::string thisOrientation;

  thisOrientation = determineLinkageConfiguration();

  if(local_debug > 0)
  {//Log Information
    ss.str("");
    ss << "The reducing monosachharide's hydroxyl oxygen (O";
    if(linkage_type_[2] != '6')
    {
      ss << linkage_type_[2];//Number for Reducing Mono's O
      ss << ") involved in the linkage has an ";
    }
    else
    {
      ss << "4) [whose orientation is used for CHI energy of omega angles]";
    }
    ss << thisOrientation << " orientation.";
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());

    if(glycosidic_oxygen_ != NULL)
    {//Print Info
      ss.str("");//clear stringstream
      ss << "The glycosidic_oxygen_oxygen atom is ";
      ss << glycosidic_oxygen_->GetId();
      gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
      ss.str("");//clear stringstream
    }
    else
    {//Print Error
      ss.str("");//clear stringstream
      ss << "The glycosidic oxygen atom is NULL";
      gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
      ss.str("");//clear stringstream
    }
  }

  if(thisOrientation == "Unknown")
  {//Error
    if(local_debug > 0)
    {//Print Error
      ss.str("");//clear stringstream
      ss << "Something went wrong determining the orientation (";
      ss << thisOrientation;
      ss << ") of the reducing monosaccharide's (glycosidic) oxygen.";
      gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
      ss.str("");//clear stringstream
    }
    return -1;
  }
  else
  {
    std::list<std::string> ChiLinkTypes = {"1-2","1-3","1-4",
                                           "1-6","2-3","2-6"};
    if(find(ChiLinkTypes.begin(), ChiLinkTypes.end(),
            linkage_type_) != ChiLinkTypes.end())
    {
      if((linkage_type_ == "1-2") || (linkage_type_ == "1-4"))
      {
        if(thisOrientation == "axial")
        {//1-2 or 1-4 axial
          return 1;
        }
        else if(thisOrientation == "equatorial")
        {//1-2 or 1-4 equatorial
          return 2;
        }
        else
        {//ERROR
          ss.str("");
          ss << "Something went horribly wrong and the code should ";
          ss << "not have reached this point.";
          gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
          ss.str("");//clear stringstream for other messages
          return -1;
        }
      }
      else if((linkage_type_ == "1-3") || (linkage_type_ == "2-3"))
      {
        if(thisOrientation == "axial")
        {//1-3 or 2-3 axial
          return 2;
        }
        else if(thisOrientation == "equatorial")
        {//1-3 or 2-3 equatorial
          return 1;
        }
        else
        {//ERROR
          ss.str("");
          ss << "Something went horribly wrong and the code should ";
          ss << "not have reached this point.";
          gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
          ss.str("");//clear stringstream for other messages
          return -1;
        }
      }
      else if((linkage_type_ == "1-6") || (linkage_type_ == "2-6"))
      {
        if(anomeric_configuration_ == "a")
        {
          return 3;
        }
        else if(anomeric_configuration_ == "b")
        {
          return 4;
        }
        else
        {//ERROR
          ss.str("");
          ss << "Something went horribly wrong and the code should ";
          ss << "not have reached this point.";
          gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
          ss.str("");//clear stringstream for other messages
          return -1;
        }
      }
    }
    else
    {//No CHI Energy function developed
      if(local_debug > 0)
      {//Print Warning
        ss.str("");
        ss << linkage_type_;
        ss << " linkages do not have a Psi CHI energy function.";
        gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
        ss.str("");//clear stringstream for other messages
      }
      return -1;
    }
  }

  //code shouldn't reach here
  return -1;
}

double Glycan::GlycosidicLinkage::CalculateOmegaChiEnergy()
{
  //////////////////////////////////////////////////////////////////////
  //    For 1-6 & 2-6 CHI energy function details see the SI in       //
  //          https://doi.org/10.1021/acs.jctc.5b00834                //
  //////////////////////////////////////////////////////////////////////

  int local_debug = -1;
  std::stringstream ss;//for debugging & logs
  double omegaCHI;

  //There are 2 different energy functions for Omega CHI Energy
  bool omega_axial = false;
  bool omega_equatorial = false;

  ////////////////////////////////////////////
  //Logic to determine which function to use//
  ////////////////////////////////////////////
  std::string thisOrientation = determineLinkageConfiguration();//will get orientation of O4
  if(thisOrientation == "Unknown")
  {
    if(local_debug > 0)
    {//Print Error
      ss.str("");
      ss << "Could not get orientation of O4; Exiting CalculateOmegaChiEnergy().";
      gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
      ss.str("");//clear stringstream for other messages
    }
    return 0;
  }
  else if(thisOrientation == "axial")
  {
    omega_axial = true;
    omega_CHI_function_ = 1;
  }
  else if(thisOrientation == "equatorial")
  {
    omega_equatorial =  true;
    omega_CHI_function_ = 2;
  }
  else
  {
    if(local_debug > 0)
    {//Print Error
      ss.str("");
      ss << "Something went horribly wrong";
      gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
      ss.str("");//clear stringstream for other messages
    }
    return 0;
  }

  ////////////////////////
  //CHI Energy Functions//
  ////////////////////////

  if(omega_axial)//omega_6A in VinaCarb
  {
    double k=0.0025;
    if((0.0 <= omega_angle_) && (omega_angle_ < 120.0))
    {
      omegaCHI = k * pow((omega_angle_-60),2);
      return omegaCHI;
    }
    else if(omega_angle_>=120.0 && omega_angle_<240.0)
    {
      omegaCHI = k * pow((omega_angle_-180),2) + 0.3;
      return omegaCHI;
    }
    else if(omega_angle_>=240.0 && omega_angle_<360.0)
    {
      omegaCHI = k * pow((omega_angle_-300),2) + 1.0;
      return omegaCHI;
    }
    else
    {//Something went wrong
      if(local_debug > 0)
      {//Print Error
        ss.str("");
        ss << "Omega angle of " << omega_angle_;
        ss << " outside of normal range (0-360).";
        ss << "Unable to run CHI Energy analysis";
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        ss.str("");
      }
      return 0;
    }
  }
  else if(omega_equatorial)//omega_6E in VinaCarb
  {
    double k=0.0025;
    if(omega_angle_>=0.0 && omega_angle_<120.0)
    {
      omegaCHI = k * pow((omega_angle_-60),2) + 0.21;
      return omegaCHI;
    }
    else if(omega_angle_>=120.0 && omega_angle_<240.0)
    {
      omegaCHI = k * pow((omega_angle_-180),2) + 1.39;
      return omegaCHI;
    }
    else if(omega_angle_>=240.0 && omega_angle_<360.0)
    {
      omegaCHI = k * pow((omega_angle_-300),2);
      return omegaCHI;
    }
    else
    {//Something went wrong
      if(local_debug > 0)
      {//Print Error
        ss.str("");
        ss << "Omega angle of " << omega_angle_;
        ss << " outside of normal range (0-360).";
        ss << "Unable to run CHI Energy analysis";
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        ss.str("");
      }
      return 0;
    }
  }
  else
  {//something went horribly wrong
    //there is a check for this after setting booleans for
    //which function to use, so the code should never get here.
    if(local_debug > 0)
    {
      ss.str("");
      ss << "Unable to calculate the Omega angle of ";
      ss << non_reducing_mono_->cycle_atoms_[0]->GetResidue()->GetName();
      ss << " ";
      ss << linkage_type_;
      ss << " ";
      ss << reducing_mono_->cycle_atoms_[0]->GetResidue()->GetName();
      gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
      ss.str("");//clear stringstream for other messages
    }
    return 0;
  }
  return 0;
}

std::string Glycan::GlycosidicLinkage::determineLinkageConfiguration()
{
  int local_debug = -1;
  std::stringstream ss;//for debugging & logs
  // std::string orientation = "Unknown";/**<- Do I even need this? */
  int anomericCarbonNum;
  int otherCarbonNum;/**<- the carbon of the reducing (or other non reducing) monosaccharide*/

  Monosaccharide* otherMono;
  if(anomeric_anomeric_linkage_)
  {
    otherMono = non_reducing_mono_2_;
  }
  else
  {
    otherMono = reducing_mono_;
  }
  /*TODO Update all code to use this instead of repeated anomeric linkage logic
  */

  //Check for NULL atoms, return Unknown if found
  if(non_reducing_mono_carbon_ == NULL)
  {
    if(local_debug > 0)
    {//Print Error
      ss.str("");
      ss << "non_reducing_mono_carbon_ is NULL";
      gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
    }
    return "Unknown";
  }
  if((reducing_mono_carbon_ == NULL) && !anomeric_anomeric_linkage_)
  {
    if(local_debug > 0)
    {//Print Error
      ss.str("");
      ss << "reducing_mono_carbon_ is NULL";
      gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
    }
    return "Unknown";
  }
  else if((non_reducing_mono_2_carbon_ == NULL) && anomeric_anomeric_linkage_)
  {
    if(local_debug > 0)
    {//Print Error
      ss.str("");
      ss << "non_reducing_mono_2_carbon_ is NULL";
      gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
    }
    return "Unknown";
  }


  //Get non_reducing_mono_carbon_ number
  if(non_reducing_mono_->sugar_name_.chemical_code_string_.find("-1") != std::string::npos)
  {//Then this monosaccharide's anomeric Carbon should be C2
    anomericCarbonNum = 2;
    if(non_reducing_mono_carbon_->GetId()[1] != '2')
    {//Something went wrong
      if(local_debug > 0)
      {//Print Error
        ss.str("");
        ss << "The carbon of " << non_reducing_mono_->residue_name_;
        ss<< "is expected to be C2 but is C" << non_reducing_mono_carbon_->GetId()[1];
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        ss.str("");
      }
    }
  }
  else
  {//Anomeric Carbon should be C1
    anomericCarbonNum = 1;
    if(non_reducing_mono_carbon_->GetId()[1] != '1')
    {//Something went wrong
      if(local_debug > 0)
      {//Print Error
        ss.str("");
        ss << "The carbon of " << non_reducing_mono_->residue_name_;
        ss<< "is expected to be C1 but is C" << non_reducing_mono_carbon_->GetId()[1];
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
        ss.str("");
      }
    }
  }

  //Get other carbon number
  if(anomeric_anomeric_linkage_)
  {
    if(non_reducing_mono_2_->sugar_name_.chemical_code_string_.find("-1") != std::string::npos)
    {//Then this monosaccharide's anomeric Carbon should be C2
      otherCarbonNum = 2;
      if(non_reducing_mono_2_carbon_->GetId()[1] != '2')
      {//Something went wrong
        if(local_debug > 0)
        {//Print Error
          ss.str("");
          ss << "The carbon of " << non_reducing_mono_2_->residue_name_;
          ss<< "is expected to be C2 but is C" << non_reducing_mono_2_carbon_->GetId()[1];
          gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
          ss.str("");
        }
      }
    }
    else
    {//Anomeric Carbon should be C1
      otherCarbonNum = 1;
      if(non_reducing_mono_2_carbon_->GetId()[1] != '1')
      {//Something went wrong
        if(local_debug > 0)
        {//Print Error
          ss.str("");
          ss << "The carbon of " << non_reducing_mono_2_->residue_name_;
          ss<< "is expected to be C1 but is C" << non_reducing_mono_2_carbon_->GetId()[1];
          gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
          ss.str("");
        }
      }
    }
  }
  else
  {//not anomeric_anomeric_linkage_
    if(local_debug > 0)
    {//Print Debugging
      ss.str("");
      ss << "reducing_mono_->cycle_atoms_str_:" << reducing_mono_->cycle_atoms_str_;
      gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
      ss.str("");
      ss << "reducing_mono_carbon_->GetId()[1]: "<< reducing_mono_carbon_->GetId()[1];
      gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
      ss.str("");
    }

    if(isdigit(reducing_mono_carbon_->GetId()[1]))
    {
      otherCarbonNum = std::stoi(reducing_mono_carbon_->GetId().substr(1,1));
    }
    else
    {
      otherCarbonNum = 0;
      if(local_debug > 0)
      {//Print warning
        ss.str("");
        ss << "The carbon of " << reducing_mono_->residue_name_;
        ss<< " is labeled C" << reducing_mono_carbon_->GetId()[1];
        gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
        ss.str("");
      }
    }

    if(1 < otherCarbonNum && otherCarbonNum < 10)
    {//Other carbon's number is between 2 and 9
      //Everything is fine?
      if(local_debug > 0)
      {
        ss.str("");
        ss << "otherCarbonNum: " << otherCarbonNum;
        gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
        ss.str("");
      }
    }
    else
    {//something is wrong with the numbering or weird with the ring
      if(local_debug > 0)
      {
        ss.str("");
        ss << "The reducing monosaccharide's carbon numbering is incorrect.";
        ss << " The connecting atom is labeled: ";
        ss << reducing_mono_carbon_->GetId();
        gmml::log(__LINE__,__FILE__,gmml::WAR, ss.str());
      }

      if(reducing_mono_carbon_->GetIsExocyclicCarbon())
      {//atom not in ring and numbered wrong

        /*TODO Refactor
            This is a really ugly way to do it. It should be a BFS algorithm.
            It should really be done elsewhere, so all of the monosaccharide's
            carbons and side groups are numbered correctly and easy to identify
        */

        //If neighbor is in ring, then this is C6 (or 7 if -1 Carbon exists)
        std::vector<MolecularModeling::Atom*> neighborAtoms = reducing_mono_carbon_->GetNode()->GetNodeNeighbors();
        for(std::vector<MolecularModeling::Atom*>::iterator it = neighborAtoms.begin(); it != neighborAtoms.end(); it++)
        {
          MolecularModeling::Atom* thisNeighbor = (*it);
          if(thisNeighbor->GetResidue()->GetId() != reducing_mono_carbon_->GetResidue()->GetId())
          {//skip
            continue;
          }
          if(reducing_mono_->cycle_atoms_str_.find(thisNeighbor->GetId()) != std::string::npos)
          {//This neighbor is in the ring
            otherCarbonNum = 6;
            if(reducing_mono_->sugar_name_.chemical_code_string_.find("-1") != std::string::npos)
            {
              otherCarbonNum++;
            }
            break;
          }
          else
          {//go to neighbor's neighbors
            std::vector<MolecularModeling::Atom*> neighbor2Atoms = thisNeighbor->GetNode()->GetNodeNeighbors();
            for(std::vector<MolecularModeling::Atom*>::iterator it = neighbor2Atoms.begin(); it != neighbor2Atoms.end(); it++)
            {
              MolecularModeling::Atom* thisN2neighbor = (*it);
              if(thisN2neighbor->GetResidue()->GetId() != reducing_mono_carbon_->GetResidue()->GetId())
              {//skip
                continue;
              }
              if(reducing_mono_->cycle_atoms_str_.find(thisN2neighbor->GetId()) != std::string::npos)
              {//This neighbor is in the ring
                otherCarbonNum = 7;
                if(reducing_mono_->sugar_name_.chemical_code_string_.find("-1") != std::string::npos)
                {
                  otherCarbonNum++;
                }
                break;
              }
              else
              {
                std::vector<MolecularModeling::Atom*> neighbor3Atoms = thisN2neighbor->GetNode()->GetNodeNeighbors();
                for(std::vector<MolecularModeling::Atom*>::iterator it = neighbor3Atoms.begin(); it != neighbor3Atoms.end(); it++)
                {
                  MolecularModeling::Atom* thisN3neighbor = (*it);
                  if(thisN3neighbor->GetResidue()->GetId() != reducing_mono_carbon_->GetResidue()->GetId())
                  {//skip
                    continue;
                  }
                  if(reducing_mono_->cycle_atoms_str_.find(thisN3neighbor->GetId()) != std::string::npos)
                  {//This neighbor's neighbor's neighbor (Ugh I know) is in the ring
                    otherCarbonNum = 8;
                    if(reducing_mono_->sugar_name_.chemical_code_string_.find("-1") != std::string::npos)
                    {
                      otherCarbonNum++;
                    }
                    break;
                  }
                  else
                  {//In case of C9
                    std::vector<MolecularModeling::Atom*> neighbor3Atoms = thisN2neighbor->GetNode()->GetNodeNeighbors();
                    for(std::vector<MolecularModeling::Atom*>::iterator it = neighbor3Atoms.begin(); it != neighbor3Atoms.end(); it++)
                    {
                      MolecularModeling::Atom* thisN3neighbor = (*it);
                      if(thisN3neighbor->GetResidue()->GetId() != reducing_mono_carbon_->GetResidue()->GetId())
                      {//skip
                        continue;
                      }
                      if(reducing_mono_->cycle_atoms_str_.find(thisN3neighbor->GetId()) != std::string::npos)
                      {//This is in the ring
                        otherCarbonNum = 9;
                        if(reducing_mono_->sugar_name_.chemical_code_string_.find("-1") != std::string::npos)
                        {//Warning as none of the code handles c10 in sugars
                          ss.str("");
                          ss << "The code thinks there is a C10 carbon, but no ";
                          ss << "code exists to handle them, including the lookup";
                          ss << " table.";
                          gmml::log(__LINE__,__FILE__, gmml::WAR, ss.str());
                          ss.str("");
                        }
                        break;
                      }
                      else
                      {
                        if(local_debug > 0)
                        {
                          ss.str("");
                          ss << "Could not determine exocyclic Carbon number";
                          gmml::log(__LINE__,__FILE__, gmml::ERR, ss.str());
                          ss.str("");
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      else
      {
        int position = 1; //start at anomeric carbon
        if(reducing_mono_->sugar_name_.chemical_code_string_.find("-1") != std::string::npos)
        {//Anomeric Carbon is C2
          position = 2;
        }
        for(std::vector<MolecularModeling::Atom*>::iterator it = reducing_mono_->cycle_atoms_.begin(); it != reducing_mono_->cycle_atoms_.end() - 1; it++)
        {//Find out the corerct numbering
          MolecularModeling::Atom* thisAtom = (*it);
          if(reducing_mono_carbon_->GetId() == thisAtom->GetId())
          {
            otherCarbonNum = position;
            break; //exit loop
          }
          position++;
        }
      }
    }
  }

  if(anomericCarbonNum != 1 && anomericCarbonNum != 2)
  {//Print error, return Unknown
    ss.str("");
    ss << "Unable to determine the non_reducing_mono_carbon_ number";
    gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
    return "Unknown";
  }
  if(otherCarbonNum < 0 || otherCarbonNum > 9)
  {//Print error, return Unknown
    ss.str("");
    ss << "Unable to determine the other Carbon's number";
    gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
    return "Unknown";
  }

  //Get the position of the linkage in the chemical_code_string_
  if(local_debug > 0)
  {//Print linkage info
    ss.str("");
    ss << "The linkage is " << anomericCarbonNum << "-" << otherCarbonNum;
    gmml::log(__LINE__,__FILE__,gmml::INF, ss.str());
    ss.str("");
  }
  std::string upOrDown;
  std::size_t position;
  std::string* chemicalCode = &otherMono->sugar_name_.chemical_code_string_;
  /*x-6 linkages are affected by O4's orientation.
  //Also, as these carbons are not in the ring, they can't really have eq or
  //ax hydroxyls.
  if((chemicalCode.find("+") != std::string::npos) && otherCarbonNum > 5)
  {//only + atoms (exocyclic carbons) get renumbered in chemical_code_
    int plusNum = otherCarbonNum - 5;
    if(chemicalCode.find("-1") != std::string::npos)
    {//Need to subtract 1
      plusNum--;
    }
    if(otherMono->sugar_name_.ring_type_ == "F")
    {//Furanose has 1 less carbon in the ring, so need to add 1
      plusNum++;
    }
    std::stringstream posSs;//to find the position of the orientation
    posSs << "+" << plusNum;
    position = chemicalCode.find(posSs.str());
  }*/
  if(otherCarbonNum == 6)
  {
    position = chemicalCode->find("4");
  }
  else
  {//numbering is the same
    position = chemicalCode->find(std::to_string(otherCarbonNum));
  }
  if(position == std::string::npos)
  {//Error
    return "Unknown";
  }
  if(position == 0)
  {//Error
    return "Unknown";
  }
  //Figure out orientation of linkage
  //The orientation is the character before the carbon number
  if(local_debug > 0)
  {//Printing info to debug
    ss.str("");
    ss << "otherMono's chemicalCode: " << *chemicalCode;
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());

    ss.str("");
    ss << "position [" << position << "] has value ";
    ss << chemicalCode->at(position);
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());

    ss.str("");
    ss << "position [" << position << "-1] has value ";
    ss << chemicalCode->at(position-1);
    gmml::log(__LINE__, __FILE__, gmml::INF, ss.str());
  }
  upOrDown = chemicalCode->at(position - 1);
  if((upOrDown != "_") && (upOrDown != "^"))
  {//Error
    if(local_debug > 0)
    {
      ss.str("");
      ss << "Could not determine the orientation of the hydroxyl. ";
      ss << "The chemical code orientation for " << position;
      ss << " is: " << upOrDown;
      ss << " Exiting determineLinkageConfiguration()";
      gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
    }
    return "Unknown";
  }

  /*convert upOrDown to axial or equatorial
      Assuming 4C1 chair for D and 1C4 for L isomers
        This will all change when weighting function is added

      Not considering furanose rings for now as they do not have CHI
      Energy functions.  Also they're kind of a pain.

      This should really be looking at if it is perpendicular to the
      plane of the ring which would be axial.  For now this will have to do
  */
  if(otherMono->sugar_name_.ring_type_ == "F")
  {//Unknown for now
    return "Unknown";
  }

  if(otherMono->sugar_name_.isomer_ == "D")
  {//and 4C1
    if(otherCarbonNum % 2 == 0)
    {//carbon is even
      if(chemicalCode->find("-1") != std::string::npos)
      {//reverses even/odd rules
        //up is equatorial, down is axial
        if(upOrDown == "^")
        {
          return "equatorial";
        }
        else//"-"
        {//already checked to make sure it's only ^ or _
          return "axial";
        }
      }
      else
      {//up is axial, down is equatorial
        if(upOrDown == "^")
        {
          return "axial";
        }
        else//"-"
        {//already checked to make sure it's only ^ or _
          return "equatorial";
        }
      }
    }
    else
    {//carbon is odd
      if(chemicalCode->find("-1") != std::string::npos)
      {//reverses even/odd rules
        //up is axial, down is equatorial
        if(upOrDown == "^")
        {
          return "axial";
        }
        else//"-"
        {//already checked to make sure it's only ^ or _
          return "equatorial";
        }
      }
      else
      {//up is equatorial, down is axial
        if(upOrDown == "^")
        {
          return "equatorial";
        }
        else//"-"
        {//already checked to make sure it's only ^ or _
          return "axial";
        }
      }
    }
  }
  else if(otherMono->sugar_name_.isomer_ == "L")
  {//again, assuming 1C4
    if(otherCarbonNum % 2 == 0)
    {//carbon is even
      if(chemicalCode->find("-1") != std::string::npos)
      {//reverses even/odd rules
        //up is axial, down is equatorial
        if(upOrDown == "^")
        {
          return "axial";
        }
        else//"-"
        {//already checked to make sure it's only ^ or _
          return "equatorial";
        }
      }
      else
      {//up is equatorial, down is axial
        if(upOrDown == "^")
        {
          return "equatorial";
        }
        else//"-"
        {//already checked to make sure it's only ^ or _
          return "axial";
        }
      }
    }
    else
    {//carbon is odd
      if(chemicalCode->find("-1") != std::string::npos)
      {//reverses even/odd rules
        //up is equatorial, down is axial
        if(upOrDown == "^")
        {
          return "equatorial";
        }
        else//"-"
        {//already checked to make sure it's only ^ or _
          return "axial";
        }
      }
      else
      {//up is axial, down is equatorial
        if(upOrDown == "^")
        {
          return "axial";
        }
        else//"-"
        {//already checked to make sure it's only ^ or _
          return "equatorial";
        }
      }
    }
  }
  else
  {//Error
    if(local_debug > 0)
    {
      ss.str("");
      ss << "Could not determine which isomer the sugar is.";
      ss << "Exiting determineLinkageConfiguration()";
      gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
    }
    return "Unknown";
  }

  //If it reaches here without returning, something went wrong
  if(local_debug > 0)
  {//Print Error
    ss.str("");
    ss << "Unable to determine the orientation of the linkage at the";
    ss << "reducing monosaccharide";
    gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
  }
  return "Unknown";
}
