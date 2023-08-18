#include "includes/CentralDataStructure/Selections/shaperSelections.hpp"

/* Below has following flaws:
  1) fail to find inner rotatable bonds here:
   __    __
__/  \__/  \__
  \__/  \__/
  If the above was within one residue.

  2) fail to find both connections of Res4 here:
 Res1 __
        \__ Res4
 Res2 __/

  3) It's ok with fused rings!

*/
std::vector<cds::Atom*> cdsSelections::FindCyclePoints(cds::Atom* atom, cds::Residue* residue)
{
    //    std::cout << "Entered FindCyclePoints with " << atom->getName() << std::endl;
    std::vector<cds::Atom*> rotation_points;
    std::vector<cds::Atom*> atom_path;
    bool found = false;
    if (residue->GetType() == cds::ResidueType::Protein)
    {
        cds::Atom* caAtom = residue->FindAtom("CA");
        // Find any cycle points. Note: starting from connection point to other residue
        cds::Atom* cycle_point;
        found = false;
        atom_path.clear();
        // This should only find cycles in Tyr, Trp etc. Not Asn, Ser, Thr as there aren't any unless bonding is messed
        // up.
        if (cdsSelections::FindCyclePoint(atom, residue, atom, &atom_path, &found, cycle_point))
        {
            rotation_points.push_back(cycle_point);
            found = false;
            atom_path.clear();
            cdsSelections::ClearAtomLabels(residue);
            // std::cout << "       >............." << std::endl;
            cdsSelections::FindCyclePoint(caAtom, residue, caAtom, &atom_path, &found, cycle_point);
            rotation_points.push_back(cycle_point);
        }
        // Always want this at the end of the vector
        rotation_points.push_back(caAtom);
    }
    else
    {
        cds::Atom* rotation_point;
        atom_path.clear();
        found = false;
        // Find path to first cycle atom, i.e. anomeric carbon
        //       std::cout << "Non-protein, checking for cycles..." << std::endl;
        if (cdsSelections::FindCyclePoint(atom, residue, atom, &atom_path, &found, rotation_point))
        {
            rotation_points.push_back(rotation_point);
        }
        // Ok, deal with non-protein non-cycles
        else
        { // Look for atom(s) with neighbors with no other neighbors within residue
            //           std::cout << "Dealing with non protein non-cycle" << std::endl;
            // cdsSelections::FindRotationPointsForNonCycles(atom, atom, &rotation_points);
            // I can't deal with non protein non cycles yet. How would I incode all the generic metadata?
            // Just set it all as rigid, and use the connecting atom as the "cycle point".
            rotation_points.push_back(atom);
        }
    }
    return rotation_points;
}

// Will not ignore fused rings. Explores everything to find all cycle points. Looks for cycle point closest to start
// atom.
bool cdsSelections::FindCyclePoint(cds::Atom* previous_atom, cds::Residue* residue, cds::Atom* current_atom,
                                   std::vector<cds::Atom*>* atom_path, bool* found_cycle_point, cds::Atom*& cycle_point)
{ // I definitely don't want cycles involving hydrogens (they only every form one bond, unless
    // bond by distance has bonded them).
    if (current_atom->getElement() != "H")
    {
        // Need this to explore everything. It will find same cycle point more than once, but that doesn't matter.
        current_atom->setLabels("VisitedByFindCyclePoint");
        //   std::cout << "Checking neighbors of " << current_atom->getName() << "\n";
        //   std::cout << "Found cycle points is currently: " << std::boolalpha << *found_cycle_point << std::endl;
        atom_path->push_back(current_atom);
        std::vector<cds::Atom*> neighbors = current_atom->getNeighbors();
        for (std::vector<cds::Atom*>::iterator it1 = neighbors.begin(); it1 != neighbors.end(); ++it1)
        {
            cds::Atom* neighbor = *it1;
            // If not previous atom and not from a different residue
            if ((neighbor->getIndex() != previous_atom->getIndex()) && (residue->contains(neighbor)))
            // if ( neighbor->getIndex()() != previous_atom->getIndex()() ) // Good for testing multiple cycles
            {
                //   std::cout << "Coming from previous atom " << previous_atom->getName() << " we see bonding: " <<
                //   current_atom->getName() << "->" << neighbor->getName() << "\n";
                if (std::find(atom_path->begin(), atom_path->end(), neighbor) !=
                    atom_path->end()) // If we've been at this atom before
                {
                    //                    std::cout << "Found a potential cycle point! Found cycle point already is: "
                    //                    << std::boolalpha << *found_cycle_point << std::endl;
                    if (*found_cycle_point) // If there are more than one cycle points
                    {
                        // Finds position of atoms in atom_path. Want earliest possible cycle point i.e. closest to
                        // start atom
                        std::ptrdiff_t new_cycle_position = std::distance(
                            atom_path->begin(), std::find(atom_path->begin(), atom_path->end(), neighbor));
                        std::ptrdiff_t current_cycle_position = std::distance(
                            atom_path->begin(), std::find(atom_path->begin(), atom_path->end(), cycle_point));
                        if (new_cycle_position < current_cycle_position)
                        {
                            //                            std::cout << "Updating cycle point to be: " <<
                            //                            neighbor->GetId() << std::endl;
                            cycle_point = neighbor;
                        }
                    }
                    else
                    {
                        *found_cycle_point = true;
                        cycle_point        = neighbor;
                        //                        std::cout << "Found the cycle point to be: " << cycle_point->GetId()
                        //                        << "\n";
                    }
                }
                if (neighbor->getLabel().compare("VisitedByFindCyclePoint") != 0) // Don't look back!
                {
                    // std::cout << "DEEPER STILL\n";
                    //                    std::cout << "Going one deeper with found cycle set as " << std::boolalpha <<
                    //                    *found_cycle_point << std::endl;
                    cdsSelections::FindCyclePoint(current_atom, residue, neighbor, atom_path, found_cycle_point,
                                                  cycle_point);
                }
            }
        }
    }
    return *found_cycle_point;
}

bool cdsSelections::FindPathBetweenTwoAtoms(cds::Atom* current_atom, cds::Residue* currentResidue,
                                            cds::Atom* target_atom, cds::Residue* targetResidue,
                                            std::vector<cds::Atom*>* atom_path, bool* found)
{
    // atom_path->push_back(current_atom);
    current_atom->setLabels("VistedByFindPathBetweenTwoAtoms");
    std::vector<cds::Atom*> neighbors = current_atom->getNeighbors();
    for (std::vector<cds::Atom*>::iterator it1 = neighbors.begin(); it1 != neighbors.end(); ++it1)
    {
        cds::Atom* neighbor = *it1;
        // std::cout << "neighbor: " << neighbor->getId() << " with " << neighbor->getNeighbors().size() <<
        // "connections.\n";
        if (neighbor->getIndex() == target_atom->getIndex())
        {
            *found = true;
            atom_path->push_back(neighbor);
        }
        // If not found && not previously visited atom && ( if neighbor residue is current residue || target_atom
        // residue)
        if ((*found == false) && (neighbor->getLabel() != "VistedByFindPathBetweenTwoAtoms") &&
            ((currentResidue->contains(neighbor)) || (targetResidue->contains(neighbor))))
        {
            // std::cout << "STEEEPER\n";
            cdsSelections::FindPathBetweenTwoAtoms(neighbor, currentResidue, target_atom, targetResidue, atom_path,
                                                   found);
        }
    }
    if (*found) // As you fall back out from the found target, create a list of the atoms.
    {
        // std::cout << "path atom: " << current_atom->getId() << std::endl;
        atom_path->push_back(current_atom);
    }
    return *found;
}

// I want a generic recursive function, where I can pass in the condition(s). Lots of Repeating code here.
// This one was written before the others. Could update with previous atom being passed in, though that makes the
// initial call confusing...
void cdsSelections::FindAtomsConnectingResidues(cds::Atom* current_atom, const cds::Residue* currentResidue,
                                                const cds::Residue* otherResidue,
                                                std::vector<cds::Atom*>* connecting_atoms, bool* found_neighbor)
{
    current_atom->setLabels("VisitedByFindAtomsConnectingResidues");
    // std::cout << "Checking neighbors of " << current_atom->getName() << " in " << currentResidue->getId() << "\n";
    for (auto& neighbor : current_atom->getNeighbors())
    {
        //        std::cout << "    " << neighbor->getName() << "_" << neighbor->getIndex() << ", is it in " <<
        //        otherResidue->getId() << ":\n"; for (auto & otherResidueAtom : otherResidue->getAtoms())
        //        {
        //            std::cout << ", " << otherResidueAtom->getName() << "_" << otherResidueAtom->getIndex();
        //        }
        //        std::cout << "\n";
        if (otherResidue->contains(neighbor))
        {
            *found_neighbor = true;
            connecting_atoms->push_back(current_atom);
            connecting_atoms->push_back(neighbor);
            //    std::cout << "Found the connection point: " << current_atom->getId() << " - " << neighbor->getId() <<
            //    "\n";
        }
        // If haven't visited this atom already AND don't move onto other residues
        else if ((neighbor->getLabel() != "VisitedByFindAtomsConnectingResidues") &&
                 (currentResidue->contains(neighbor)))
        {
            cdsSelections::FindAtomsConnectingResidues(neighbor, currentResidue, otherResidue, connecting_atoms,
                                                       found_neighbor);
        }
    }
    return;
}

// Find a neighbor of a cycle point to define the dihedral. Must not be in path used to come to this cycle_point.
// Oh gawd this code is horrible.
// The logic for selecting the atom to define a dihedral is messy.
// Comparing strings is messy when I really care about the number of e.g. C2 Vs O5, but need to factor in a C2 vs O
// comparision
cds::Atom* cdsSelections::FindCyclePointNeighbor(const std::vector<cds::Atom*> atom_path, cds::Atom* cycle_point,
                                                 cds::Residue* cyclePointResidue)
{
    cds::Atom* selected_neighbor;
    if (cycle_point->getName().compare("CA") == 0) // This is a protein, and we always want the N atom.
    {
        selected_neighbor = cyclePointResidue->FindAtom("N");
    } // If this is a C2 like in Sia, then we always want the C1 atom unless that atom is in the linkage path (like
      // fructose 1-1)
    else if ((cycle_point->getName().compare("C2") == 0) &&
             (std::find(atom_path.begin(), atom_path.end(), cyclePointResidue->FindAtom("C1")) == atom_path.end()))
    {
        selected_neighbor = cyclePointResidue->FindAtom("C1");
    }
    else if (cycle_point->getName().compare("C1") == 0)
    {
        selected_neighbor =
            cyclePointResidue->FindAtom("C2"); // We can always set phi to 180 this way, regardless of alpha/beta
    }
    else // This bit is overdone now, as I was looking for higher numbered atoms of C1, but now I know I always want C2,
         // so put that above.
    {
        std::vector<cds::Atom*> neighbors = cycle_point->getNeighbors();
        // Ok must first get a list of neighbors that weren't in the connection path
        std::vector<cds::Atom*> good_neighbors; // Couldn't think of a better name. Everybody needs these.
        for (auto& neighbor : neighbors)
        {
            if (!(std::find(atom_path.begin(), atom_path.end(), neighbor) !=
                  atom_path.end())) // If we've NOT been at this atom on way to cycle point
            {
                if (neighbor->getName().at(0) != 'H') // Don't find hydrogens. Later we swap out to use a hydrogen to
                                                      // define a dihedral, but that's a very specific one.
                {
                    good_neighbors.push_back(neighbor);
                }
            }
        }
        if (good_neighbors.size() == 0) // Ok take hydrogens then.
        {
            for (auto& neighbor : neighbors)
            {
                if (!(std::find(atom_path.begin(), atom_path.end(), neighbor) !=
                      atom_path.end())) // If we've NOT been at this atom on way to cycle point
                {
                    good_neighbors.push_back(neighbor);
                }
            }
        }
        if (good_neighbors.size() == 0)
        {
            //            std::cout << "About to segfault in MolecularModeling/Selections/selections.cpp
            //            cdsSelections::FindCyclePointNeighbor" << std::endl;
        }
        selected_neighbor = good_neighbors.at(0); // Set to any to start. If there are not good_neighbors then you
                                                  // deserve to crash and burn std::cout << "Good neighbors are: ";
        for (std::vector<cds::Atom*>::iterator it1 = good_neighbors.begin(); it1 != good_neighbors.end(); ++it1)
        {
            cds::Atom* neighbor = *it1;
            //  std::cout << neighbor->getName() << " ,";
            if (selected_neighbor->getName().size() >= 2)
            {
                if (neighbor->getName().size() >=
                    2) // This is the only time I want to compare and select the larger number
                {
                    if (neighbor->getName().at(1) > selected_neighbor->getName().at(1))
                    {
                        selected_neighbor = neighbor;
                    }
                } // Otherwise any neighbor is ok. Yes I'm comparing char's, but that is fine unless C9 Vs C10, but in
                  // that case I don't care again.
            }
        }
        //  std::cout << "\n";
    }
    //    std::cout << "Returning with neighbor: " << selected_neighbor->getName() << "\n";
    return selected_neighbor;
}

// After getting connection atoms between two residues, we have the linear path between them e.g.:
// (Residue1 C2 - O7 - C7 - O6 - C6 Residue2), where C2 and C6 are the cycle points. In cases like 2-7 or 2-8 linkages,
// we have significant branches from this linear path
//                    /
//                   C8-O8
//                  /
//                 C9-O9
// And we need to set reasonable values for the C9-C8, and C8-C7 dihedral angles.
// This is further complicated by the possibility of "DNeup5Aca2-7[DNeup5Aca2-8]DNeup5Aca2-OH", where one of the
// branches is part of another linkage. This will work for 7/8/9 linked sialic acids, but if we get more heavily
// branched linkages this will need to change to iteratively find branches from branches
//

void cdsSelections::FindEndsOfBranchesFromLinkageAtom(cds::Atom* currentAtom, cds::Atom* previousAtom,
                                                      cds::Residue* residue, Branch* branch)
{
    branch->ChangeDepth(1);
    currentAtom->setLabels("VistedByFindEndsOfBranchesFromLinkageAtom");
    bool deadEndAtom              = true;
    bool connectsToAnotherResidue = false;
    for (auto& neighbor : currentAtom->getNeighbors())
    {
        if (neighbor->getLabel() != "VistedByFindEndsOfBranchesFromLinkageAtom" && *neighbor != *previousAtom &&
            residue->contains(neighbor)) // Don't explore across residues.
        {
            if (neighbor->getNeighbors().size() > 1)
            {
                deadEndAtom = false;
                // std::cout << "At depth " << branch->GetDepth() << " going deeper from " << currentAtom->GetId() << "
                // to " << neighbor->GetId() << "\n";
                FindEndsOfBranchesFromLinkageAtom(neighbor, currentAtom, residue, branch);
                branch->ChangeDepth(-1);
            }
        }
        if (!residue->contains(neighbor))
        {
            connectsToAnotherResidue = true;
        }
    }
    // std::cout << "  Status at " << currentAtom->GetId() << " is deadEndAtom:" << std::boolalpha << deadEndAtom << ",
    // depth: " << branch->GetDepth() << ", connectsToOther: " << connectsToAnotherResidue << std::endl;
    if (deadEndAtom && !connectsToAnotherResidue && branch->GetDepth() > 1 && branch->AtMaxDepth())
    {
        // std::cout << "      Found a dead end: " << currentAtom->GetId() << " at depth " << branch->GetDepth() <<
        // std::endl;
        branch->SetEnd(currentAtom);
    }
    return;
}

void cdsSelections::ClearAtomLabels(cds::Residue* residue)
{
    for (auto& atom : residue->getAtoms())
    {
        atom->clearLabels();
    }
    return;
}

cds::ResidueLinkage* cdsSelections::selectLinkageWithIndex(std::vector<cds::ResidueLinkage>& inputLinkages,
                                                           const long long unsigned int indexQuery)
{
    for (auto& linkage : inputLinkages)
    {
        if (linkage.GetIndex() == indexQuery)
        {
            return &linkage;
        }
    }
    // Error
    std::stringstream ss;
    ss << "Linkage numbered " << indexQuery << " not found in linkages for this carbohydrate\n";
    gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
    throw std::runtime_error(ss.str());
}

// Gonna choke on cyclic glycans. Add a check for IsVisited when that is required.
/* This is a straight copy from glycoprotein_builder. I need a high level class that deals with both
 * cds::ResidueLinkages, ring shapes etc. That way I can create X shapes of a molecule. For now this will do to figure
 * out some implementation details like file naming.
 */
std::vector<cds::ResidueLinkage>
cdsSelections::SplitLinkagesIntoPermutants(std::vector<cds::ResidueLinkage>& inputLinkages)
{
    std::vector<cds::ResidueLinkage> sortedLinkages;
    for (auto& linkage : inputLinkages)
    {
        if (linkage.CheckIfConformer())
        {
            sortedLinkages.push_back(linkage);
        }
        else // if not a conformer
        {
            std::vector<cds::RotatableDihedral> rotatableDihedrals =
                linkage
                    .GetRotatableDihedralsWithMultipleRotamers(); // only want the rotatabe dihedrals within a linkage
                                                                  // that have multiple rotamers. Some bonds won't.
            for (auto& rotatableDihedral : rotatableDihedrals)
            {
                cds::ResidueLinkage splitLinkage         = linkage; // Copy it to get correct info into class
                std::vector<cds::RotatableDihedral> temp = {rotatableDihedral};
                splitLinkage.SetRotatableDihedrals(temp);
                sortedLinkages.push_back(splitLinkage);
                // std::cout << "Split out " << splitLinkage.GetFromThisResidue1()->GetId() << "-" <<
                // splitLinkage.GetToThisResidue2()->GetId() << " rotamer with number of shapes: " <<
                // rotatableDihedral.GetNumberOfRotamers() << "\n";
            }
        }
    }
    return sortedLinkages;
}
