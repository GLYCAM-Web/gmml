
#include "../../includes/MolecularModeling/assembly.hpp"
#include "../../includes/MolecularModeling/atom.hpp"
#include "../../includes/MolecularModeling/atomnode.hpp"

bool MolecularModeling::Assembly::guessIfC_CDoubleBond(MolecularModeling::Atom* carbon1, MolecularModeling::Atom* carbon2)
{
  bool bothCarbon = false;
  bool areNeighbors = false;
  bool areCloseEnough = false;
  //bool haveAll120angles = false;  // commented out because set but not used
  bool areTrigonalPlanar = false;
  int local_debug = -1;
  std::stringstream debugStr;
  if(local_debug > 0)
  {
    gmml::log(__LINE__, __FILE__,  gmml::INF, "In Guess Function");
    gmml::log(__LINE__, __FILE__,  gmml::INF, carbon1->GetName());
    gmml::log(__LINE__, __FILE__,  gmml::INF, carbon2->GetName());
  }
  std::string carbon1name = carbon1->GetName();
  std::string carbon2name = carbon2->GetName();
  //Check if both Carbon
  if((carbon1name.find("C") != std::string::npos)&&(carbon2name.find("C") != std::string::npos))
  {
    bothCarbon = true;
    if(local_debug > 0)
    {
      debugStr << carbon1->GetId() << " and " << carbon2->GetId() << " are carbon\n";
      gmml::log(__LINE__, __FILE__, gmml::INF, debugStr.str());
      debugStr.str("");
    }
  }
  else
  {
    return false;
  }

  //See if they are neighbors
  MolecularModeling::Assembly::AtomVector neighbors = carbon1->GetNode()->GetNodeNeighbors();
  for(MolecularModeling::Assembly::AtomVector::iterator it = neighbors.begin(); it != neighbors.end(); it++)
  {
    MolecularModeling::Atom* thisNeighbor = *it;
    if(thisNeighbor->GetId() == carbon2->GetId())
    {
      areNeighbors = true;
      if(local_debug > 0)
      {
        debugStr << carbon1->GetId() << " and " << carbon2->GetId() << " are neigbors\n";
        gmml::log(__LINE__, __FILE__, gmml::INF, debugStr.str());
        debugStr.str("");
      }
    }
  }
  if (areNeighbors == false)
  {
    return false;
  }

  //Get Distance; if < 1.48 (my guess from double bond meta data) it may be a double bond.
  if ((carbon1->GetDistanceToAtom(carbon2) < 1.33))
  {
    areCloseEnough = true;
    if(local_debug > 0)
    {
      debugStr << carbon1->GetId() << " and " << carbon2->GetId() << " are less than 1.48 angstroms apart\n";
      gmml::log(__LINE__, __FILE__, gmml::INF, debugStr.str());
      debugStr.str("");
    }
  }
  else
  {
    return false;
  }

  //Check geometry;
  MolecularModeling::Assembly::AtomVector carbon1neighbors = carbon1->GetNode()->GetNodeNeighbors();
  MolecularModeling::Assembly::AtomVector carbon2neighbors = carbon2->GetNode()->GetNodeNeighbors();

  if((carbon1neighbors.size() < 4) && (carbon2neighbors.size() < 4))//If both have three or less neighbors
  {
    if(local_debug > 0)
    {
      debugStr << carbon1->GetId() << " and " << carbon2->GetId() << " both have three or fewer neighbors\n";
      gmml::log(__LINE__, __FILE__, gmml::INF, debugStr.str());
      debugStr.str("");
    }
    for(MolecularModeling::Assembly::AtomVector::iterator it = carbon1neighbors.begin(); it != carbon1neighbors.end(); it++)
    {
      MolecularModeling::Atom* thisC1Neighbor = *it;
      if(thisC1Neighbor->GetId() != carbon2->GetId())
      {
        for(MolecularModeling::Assembly::AtomVector::iterator it2 = carbon2neighbors.begin(); it2 != carbon2neighbors.end(); it2++)
        {
          MolecularModeling::Atom* thisC2Neighbor = *it2;
          if(thisC2Neighbor->GetId() != carbon1->GetId())
          {
            if(local_debug > 0)
            {
              debugStr << thisC2Neighbor->GetId() << ", " << carbon1->GetId() << ", " << carbon2->GetId() << " have an angle of " << MolecularModeling::Assembly::CalculateBondAngleByAtoms(thisC2Neighbor, carbon2, carbon1);
              gmml::log(__LINE__, __FILE__, gmml::INF, debugStr.str());
              debugStr.str("");
              debugStr << thisC1Neighbor->GetId() << ", " << carbon1->GetId() << ", " << carbon2->GetId() << " have an angle of " << MolecularModeling::Assembly::CalculateBondAngleByAtoms(thisC1Neighbor, carbon2, carbon1);
              gmml::log(__LINE__, __FILE__, gmml::INF, debugStr.str());
              debugStr.str("");
            }
            if((MolecularModeling::Assembly::CalculateBondAngleByAtoms(thisC2Neighbor, carbon2, carbon1) < 130) && (MolecularModeling::Assembly::CalculateBondAngleByAtoms(thisC2Neighbor, carbon2, carbon1) > 110) && (MolecularModeling::Assembly::CalculateBondAngleByAtoms(thisC1Neighbor, carbon2, carbon1) < 130) && (MolecularModeling::Assembly::CalculateBondAngleByAtoms(thisC1Neighbor, carbon2, carbon1) > 110))
            {
            //  haveAll120angles = true;  // commenting out because not used
              if(local_debug > 0)
              {
                debugStr << carbon1->GetId() << " and " << carbon2->GetId() << " and " << thisC2Neighbor->GetId() << " form an angle between 115 and 125 degrees\n";
                gmml::log(__LINE__, __FILE__, gmml::INF, debugStr.str());
                debugStr.str("");
              }
            }
            else
            {
          //    haveAll120angles = false;  // commenting out because not used
            }
          
            //If all 4 torsion angles are about 0 (+-5?)
            if(local_debug > 0)
            {
              debugStr << carbon1->GetId() << ", " << carbon2->GetId() << ", " << thisC2Neighbor->GetId() << " and " << thisC1Neighbor->GetId() << " have a torsion angle of " << MolecularModeling::Assembly::CalculateTorsionAngleByAtoms(thisC1Neighbor, carbon1, carbon2, thisC2Neighbor);
              gmml::log(__LINE__, __FILE__, gmml::INF, debugStr.str());
              debugStr.str("");
            }
            if((MolecularModeling::Assembly::CalculateTorsionAngleByAtoms(thisC1Neighbor, carbon1, carbon2, thisC2Neighbor) > 6.2 /*greater than about 355 degrees as the function returns radians*/) && (MolecularModeling::Assembly::CalculateTorsionAngleByAtoms(thisC1Neighbor, carbon1, carbon2, thisC2Neighbor) < 0.09 /*less than about 5 degrees*/))
            {
              areTrigonalPlanar = true;
              if(local_debug > 0)
              {
                debugStr << thisC1Neighbor->GetId() << ", " << carbon1->GetId() << ", " << carbon2->GetId() << ", and " << thisC2Neighbor->GetId() << " are planar\n";
                gmml::log(__LINE__, __FILE__, gmml::INF, debugStr.str());
                debugStr.str("");
              }
            }
            else if((MolecularModeling::Assembly::CalculateTorsionAngleByAtoms(thisC1Neighbor, carbon1, carbon2, thisC2Neighbor) < -6.2 /*less than about -355 degrees as the function returns radians*/) && (MolecularModeling::Assembly::CalculateTorsionAngleByAtoms(thisC1Neighbor, carbon1, carbon2, thisC2Neighbor) > -0.09 /*greater than about -5 degrees*/))
            {
              areTrigonalPlanar = true;
              if(local_debug > 0)
              {
                debugStr << thisC1Neighbor->GetId() << ", " << carbon1->GetId() << ", " << carbon2->GetId() << ", and " << thisC2Neighbor->GetId() << " are planar\n";
                gmml::log(__LINE__, __FILE__, gmml::INF, debugStr.str());
                debugStr.str("");
              }
            }
            else if ((MolecularModeling::Assembly::CalculateTorsionAngleByAtoms(thisC1Neighbor, carbon1, carbon2, thisC2Neighbor) > 3.05 /*greater than about 175 degrees as the function returns radians*/) && (MolecularModeling::Assembly::CalculateTorsionAngleByAtoms(thisC1Neighbor, carbon1, carbon2, thisC2Neighbor) < 3.22 /*less than about 185 degrees*/))
            {
              areTrigonalPlanar = true;
              if(local_debug > 0)
              {
                debugStr << thisC1Neighbor->GetId() << ", " << carbon1->GetId() << ", " << carbon2->GetId() << ", and " << thisC2Neighbor->GetId() << " are planar\n";
                gmml::log(__LINE__, __FILE__, gmml::INF, debugStr.str());
                debugStr.str("");
              }
            }
            else if ((MolecularModeling::Assembly::CalculateTorsionAngleByAtoms(thisC1Neighbor, carbon1, carbon2, thisC2Neighbor) < -3.05 /*less than about -175 degrees as the function returns radians*/) && (MolecularModeling::Assembly::CalculateTorsionAngleByAtoms(thisC1Neighbor, carbon1, carbon2, thisC2Neighbor) > -3.22 /*greater than about -185 degrees*/))
            {
              areTrigonalPlanar = true;
              if(local_debug > 0)
              {
                debugStr << thisC1Neighbor->GetId() << ", " << carbon1->GetId() << ", " << carbon2->GetId() << ", and " << thisC2Neighbor->GetId() << " are planar\n";
                gmml::log(__LINE__, __FILE__, gmml::INF, debugStr.str());
                debugStr.str("");
              }
            }
            else
            {
              areTrigonalPlanar = false;
            }
          }
        }
      }
    }


  }

  if((bothCarbon == true) && (areNeighbors ==true) && (areCloseEnough == true) && (areTrigonalPlanar == true))
  {
    return true;
  }

  return false;
}
