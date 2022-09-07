#include "../../../includes/MolecularModeling/Selections/selections.hpp"
#include "../../../includes/MolecularModeling/assembly.hpp"
#include <regex>
#include <iostream>
MolecularModeling::AtomVector selection::AtomsWithinDistanceOf(MolecularModeling::Atom *query_atom, double distance, MolecularModeling::AtomVector atoms)
{
    MolecularModeling::AtomVector atoms_within_distance;
    for(MolecularModeling::AtomVector::iterator it1 = atoms.begin(); it1 != atoms.end(); ++it1)
    {
        MolecularModeling::Atom *atom1 = (*it1);
        if (atom1->GetDistanceToAtom(query_atom) < distance )
        {
            atoms_within_distance.push_back(atom1);
        }
    }
    return atoms_within_distance;
}

// I want a generic recursive function, where I can pass in the condition(s). Lots of Repeating code here.
// This one was written before the others. Could update with previous atom being passed in, though that makes the initial call confusing...
void selection::FindAtomsConnectingResidues(MolecularModeling::Atom *current_atom, MolecularModeling::Residue *second_residue, MolecularModeling::AtomVector *connecting_atoms, bool *found_neighbor)
{
    current_atom->SetDescription("VisitedByFindAtomsConnectingResidues");
    MolecularModeling::AtomVector neighbors = current_atom->GetNode()->GetNodeNeighbors();
    for(MolecularModeling::AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); ++it1)
    {
        MolecularModeling::Atom *neighbor = *it1;
        if(neighbor->GetResidue()->GetId().compare(second_residue->GetId())==0) // Should be using indices for Residues too.
        {
            *found_neighbor = true;
            connecting_atoms->push_back(current_atom);
            connecting_atoms->push_back(neighbor);
          //  std::cout << "Found the connection point: " << current_atom->GetId() << " - " << neighbor->GetId() << "\n";
        }
        // If haven't visited this atom already AND don't move onto other residues
        else if ((neighbor->GetDescription().compare("VisitedByFindAtomsConnectingResidues")!=0) && (current_atom->GetResidue()->GetId().compare(neighbor->GetResidue()->GetId())==0))
        {
            selection::FindAtomsConnectingResidues(neighbor, second_residue, connecting_atoms, found_neighbor);
        }
    }
    return;
}

// Will not ignore fused rings. Explores everything to find all cycle points. Looks for cycle point closest to start atom.
bool selection::FindCyclePoint(MolecularModeling::Atom *previous_atom, MolecularModeling::Atom *current_atom, MolecularModeling::AtomVector *atom_path, bool *found_cycle_point, MolecularModeling::Atom *&cycle_point)
{
    // I wish there was a more solid way to check element type. But here we go, I definitely don't want cycles involving hydrogens (they only every form one bond, unless
    // bond by distance has bonded them).
    if ( current_atom->GetName().at(0) != 'H')
    {
        // Need this to explore everything. It will find same cycle point more than once, but that doesn't matter.
        current_atom->SetDescription("VisitedByFindCyclePoint");
     //   std::cout << "Checking neighbors of " << current_atom->GetName() << "\n";
     //   std::cout << "Found cycle points is currently: " << std::boolalpha << *found_cycle_point << std::endl;

        atom_path->push_back(current_atom);
        MolecularModeling::AtomVector neighbors = current_atom->GetNode()->GetNodeNeighbors();
        for(MolecularModeling::AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); ++it1)
        {
            MolecularModeling::Atom *neighbor = *it1;
            // If not previous atom and not from a different residue
            if ( (neighbor->GetIndex() != previous_atom->GetIndex()) && (current_atom->GetResidue()->GetId().compare(neighbor->GetResidue()->GetId())==0))
                //if ( neighbor->GetIndex() != previous_atom->GetIndex() ) // Good for testing multiple cycles
            {
            //   std::cout << "Coming from previous atom " << previous_atom->GetName() << " we see bonding: " << current_atom->GetName() << "->" << neighbor->GetName() << "\n";
                if ( std::find(atom_path->begin(), atom_path->end(), neighbor) != atom_path->end() ) // If we've been at this atom before
                {
//                    std::cout << "Found a potential cycle point! Found cycle point already is: " << std::boolalpha << *found_cycle_point << std::endl;
                    if(*found_cycle_point) // If there are more than one cycle points
                    {
                        // Finds position of atoms in atom_path. Want earliest possible cycle point i.e. closest to start atom
                        std::ptrdiff_t new_cycle_position = std::distance(atom_path->begin(), std::find(atom_path->begin(), atom_path->end(), neighbor));
                        std::ptrdiff_t current_cycle_position = std::distance(atom_path->begin(), std::find(atom_path->begin(), atom_path->end(), cycle_point));
                        if (new_cycle_position < current_cycle_position)
                        {
//                            std::cout << "Updating cycle point to be: " << neighbor->GetId() << std::endl;
                            cycle_point = neighbor;
                        }
                    }
                    else
                    {
                        *found_cycle_point = true;
                        cycle_point = neighbor;
//                        std::cout << "Found the cycle point to be: " << cycle_point->GetId() << "\n";
                    }
                }
                if(neighbor->GetDescription().compare("VisitedByFindCyclePoint")!=0) // Don't look back!
                {
                    //std::cout << "DEEPER STILL\n";
//                    std::cout << "Going one deeper with found cycle set as " << std::boolalpha << *found_cycle_point << std::endl;
                    selection::FindCyclePoint(current_atom, neighbor, atom_path, found_cycle_point, cycle_point);
                }
            }
        }
    }
    return *found_cycle_point;
}

// If current_atom is connected only to atoms that have zero other intra residue connections, it is a rotation point
bool selection::FindRotationPointsForNonCycles(MolecularModeling::Atom *previous_atom, MolecularModeling::Atom *current_atom, MolecularModeling::AtomVector *rotation_points)
{
    // My problem with the below is we are now dealing with non-protein, non-cycles that we don't have any metadata for.
    // Rather than figure out how each bond would rotate, deal with branching etc, I'll just leave it all rigid for now.

//    MolecularModeling::AtomVector neighbors = current_atom->GetNode()->GetNodeNeighbors();
//    MolecularModeling::AtomVector intra_node_neighbors;
//    for (auto &neighbor : neighbors)
//    { // If not previous atom and not from a different residue
//        if ( (neighbor->GetIndex() != previous_atom->GetIndex()) && (current_atom->GetResidue()->GetId().compare(neighbor->GetResidue()->GetId())==0))
//        {
//            intra_node_neighbors.push_back(neighbor);
//        }
//    }
//    bool found_rotation_point = true;
//    for (auto &neighbor : intra_node_neighbors)
//    {
//        MolecularModeling::AtomVector neighbor_neighbors = neighbor->GetNode()->GetNodeNeighbors();
//        MolecularModeling::AtomVector intra_node_neighbor_neighbors;
//        for (auto &neighbor_neighbor : neighbor_neighbors)
//        {
//            if ( (neighbor_neighbor->GetIndex() != current_atom->GetIndex()) && (current_atom->GetResidue()->GetId().compare(neighbor_neighbor->GetResidue()->GetId())==0))
//            {
//                intra_node_neighbor_neighbors.push_back(neighbor);
//            }
//        }
//        if (intra_node_neighbor_neighbors.size() > 0 )
//        {
//            found_rotation_point = false;
//        }
//    }
//    if ( found_rotation_point ) // Deadend Chitchatchoes. And we aren't right at the beginning looking down a deadend branch.
//    {
//        rotation_points->push_back(current_atom);
//        //std::cout << "Pushed back: " << current_atom->GetId();
//    }
//    else // Keep looking for the ends
//    {
//        for (auto &neighbor : intra_node_neighbors)
//        {
//            // Only want one dead end for now, can't handle mulitple ones yet. Also, don't land on a dead end. Should not happen except at very start of recursion e.g(otherRes-C(=O)-CH3)
//            //std::cout << neighbor->GetId() << "neighbor of\n" << current_atom->GetId() << "has " << neighbor->GetNode()->GetIntraNodeNeighbors().size()  << " intranode neih\n";
//            if ( (rotation_points->empty()) && (neighbor->GetNode()->GetNodeNeighbors().size() > 1) )
//                selection::FindRotationPointsForNonCycles(current_atom, neighbor, rotation_points);
//        }
//    }
//    return (!rotation_points->empty()); // return false if rotation_points is empty
}


bool selection::FindPathBetweenTwoAtoms(MolecularModeling::Atom *current_atom, MolecularModeling::Atom *target_atom, MolecularModeling::AtomVector *atom_path, bool *found)
{
    //atom_path->push_back(current_atom);
    current_atom->SetDescription("VistedByFindPathBetweenTwoAtoms");
    MolecularModeling::AtomVector neighbors = current_atom->GetNode()->GetNodeNeighbors();
    for(MolecularModeling::AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); ++it1)
    {
        MolecularModeling::Atom *neighbor = *it1;
        if(neighbor->GetIndex() == target_atom->GetIndex())
        {
            *found = true;
            atom_path->push_back(neighbor);
        }
        // If not found && not previously visited atom && ( if neighbor residue is current residue || target_atom residue)
        if ( (*found == false) 
            && (neighbor->GetDescription() != "VistedByFindPathBetweenTwoAtoms" ) 
            //(neighbor->GetDescription() != "VistedByFindPathBetweenTwoAtoms" ) 
            && ( (current_atom->GetResidue() == neighbor->GetResidue()) || (target_atom->GetResidue() == neighbor->GetResidue()) )
            //&&  ( (current_atom->GetResidue()->GetId() == neighbor->GetResidue()->GetId()) || (target_atom->GetResidue()->GetId() == neighbor->GetResidue()->GetId()) )
            )
        {
            //std::cout << "STEEEPER\n";
            selection::FindPathBetweenTwoAtoms(neighbor, target_atom, atom_path, found);
        }
    }
    if(*found) // As you fall back out from the found target, create a list of the atoms.
    {
        //std::cout << "path atom: " << current_atom->GetId() << std::endl;
        atom_path->push_back(current_atom);
    }
    return *found;
}

void selection::ClearAtomDescriptions(MolecularModeling::Residue *residue)
{
    MolecularModeling::AtomVector atoms = residue->GetAtoms();
    for(MolecularModeling::AtomVector::iterator it = atoms.begin(); it != atoms.end(); ++it)
    {
        MolecularModeling::Atom *atom = *it;
        atom->SetDescription("");
    }
    return;
}

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
MolecularModeling::AtomVector selection::FindCyclePoints(MolecularModeling::Atom *atom)
{
//    std::cout << "Entered FindCyclePoints with " << atom->GetId() << std::endl;
    MolecularModeling::AtomVector rotation_points;
    MolecularModeling::AtomVector atom_path;
    bool found = false;
    if(atom->GetResidue()->CheckIfProtein())
    {
        MolecularModeling::Atom *caAtom = atom->GetResidue()->GetAtom("CA");
        // Find any cycle points. Note: starting from connection point to other residue
        MolecularModeling::Atom *cycle_point;
        found = false;
        atom_path.clear();
        // This should only find cycles in Tyr, Trp etc. Not Asn, Ser, Thr as there aren't any unless bonding is messed up.
        if(selection::FindCyclePoint(atom, atom, &atom_path, &found, cycle_point))
        {
            rotation_points.push_back(cycle_point);
            found = false;
            atom_path.clear();
            //Want residue.clearAtomDescriptions
            selection::ClearAtomDescriptions(atom->GetResidue());
           // std::cout << "       >............." << std::endl;
            selection::FindCyclePoint(caAtom, caAtom, &atom_path, &found, cycle_point);
            rotation_points.push_back(cycle_point);
        }
        // Always want this at the end of the vector
        rotation_points.push_back(caAtom);
    }
    else 
    {
        MolecularModeling::Atom *rotation_point;
        atom_path.clear();
        found = false;
        //Find path to first cycle atom, i.e. anomeric carbon
 //       std::cout << "Non-protein, checking for cycles..." << std::endl;
        if(selection::FindCyclePoint(atom, atom, &atom_path, &found, rotation_point))
        {
            rotation_points.push_back(rotation_point);
        }
        // Ok, deal with non-protein non-cycles
        else
        { // Look for atom(s) with neighbors with no other neighbors within residue
 //           std::cout << "Dealing with non protein non-cycle" << std::endl;
            //selection::FindRotationPointsForNonCycles(atom, atom, &rotation_points);
            // I can't deal with non protein non cycles yet. How would I incode all the generic metadata?
            // Just set it all as rigid, and use the connecting atom as the "cycle point".
            rotation_points.push_back(atom);
        }
    }
    return rotation_points;
}

// Find a neighbor of a cycle point to define the dihedral. Must not be in path used to come to this cycle_point.
// Oh gawd this code is horrible.
// The logic for selecting the atom to define a dihedral is messy.
// Comparing strings is messy when I really care about the number of e.g. C2 Vs O5, but need to factor in a C2 vs O comparision
MolecularModeling::Atom* selection::FindCyclePointNeighbor(const MolecularModeling::AtomVector atom_path, MolecularModeling::Atom *cycle_point)
{
    MolecularModeling::Atom *selected_neighbor;
    if (cycle_point->GetName().compare("CA")==0) // This is a protein, and we always want the N atom.
    {
        selected_neighbor = cycle_point->GetResidue()->GetAtom("N");
    } // If this is a C2 like in Sia, then we always want the C1 atom unless that atom is in the linkage path (like fructose 1-1) 
    else if ( (cycle_point->GetName().compare("C2")==0) && (std::find(atom_path.begin(), atom_path.end(), cycle_point->GetResidue()->GetAtom("C1")) == atom_path.end()) )
    {
        selected_neighbor = cycle_point->GetResidue()->GetAtom("C1");
    }
    else if (cycle_point->GetName().compare("C1")==0)
    {
        selected_neighbor = cycle_point->GetResidue()->GetAtom("C2"); // We can always set phi to 180 this way, regardless of alpha/beta
    }
    else // This bit is overdone now, as I was looking for higher numbered atoms of C1, but now I know I always want C2, so put that above.
    {
        MolecularModeling::AtomVector neighbors = cycle_point->GetNode()->GetNodeNeighbors();
        // Ok must first get a list of neighbors that weren't in the connection path
        MolecularModeling::AtomVector good_neighbors; // Couldn't think of a better name. Everybody needs these.
        for(auto &neighbor : neighbors)
        {
            if ( ! (std::find(atom_path.begin(), atom_path.end(), neighbor) != atom_path.end()) ) // If we've NOT been at this atom on way to cycle point
            {
                if ( neighbor->GetName().at(0) != 'H' ) // Don't find hydrogens. Later we swap out to use a hydrogen to define a dihedral, but that's a very specific one.
                {
                    good_neighbors.push_back(neighbor);
                }
            }
        }
        if(good_neighbors.size()==0) // Ok take hydrogens then.
        {
            for(auto &neighbor : neighbors)
            {
                if ( ! (std::find(atom_path.begin(), atom_path.end(), neighbor) != atom_path.end()) ) // If we've NOT been at this atom on way to cycle point
                {
                    good_neighbors.push_back(neighbor);
                }

            }
        }
        if(good_neighbors.size()==0)
        {
//            std::cout << "About to segfault in MolecularModeling/Selections/selections.cpp selection::FindCyclePointNeighbor" << std::endl;
        }
        selected_neighbor = good_neighbors.at(0); // Set to any to start. If there are not good_neighbors then you deserve to crash and burn
       // std::cout << "Good neighbors are: ";
        for(MolecularModeling::AtomVector::iterator it1 = good_neighbors.begin(); it1 != good_neighbors.end(); ++it1)
        {
            MolecularModeling::Atom *neighbor = *it1;
          //  std::cout << neighbor->GetName() << " ,";
            if(selected_neighbor->GetName().size() >= 2)
            {
                if(neighbor->GetName().size() >= 2) // This is the only time I want to compare and select the larger number
                {
                    if(neighbor->GetName().at(1) > selected_neighbor->GetName().at(1))
                    {
                        selected_neighbor = neighbor;
                    }
                } // Otherwise any neighbor is ok. Yes I'm comparing char's, but that is fine unless C9 Vs C10, but in that case I don't care again.
            }
        }
      //  std::cout << "\n";
    }
//    std::cout << "Returning with neighbor: " << selected_neighbor->GetName() << "\n";
    return selected_neighbor;
}

// Ok not really a selection, but like, chill for a minute ok?
double selection::GetMaxDistanceBetweenAtoms(MolecularModeling::AtomVector atoms)
{
    double max_distance = 0.0;
    for(MolecularModeling::AtomVector::iterator it1 = atoms.begin(); it1 != atoms.end(); ++it1)
    {
        MolecularModeling::Atom *atom1 = (*it1);
        for(MolecularModeling::AtomVector::iterator it2 = it1; it2 != atoms.end(); ++it2)
        {
            MolecularModeling::Atom *atom2 = (*it2);
            if (atom1->GetDistanceToAtom(atom2) > max_distance)
            {
                max_distance = atom1->GetDistanceToAtom(atom2);
            }
        }
    }
    return max_distance;
}

MolecularModeling::AtomVector selection::GetAtomsCommonToBothAtomVectors(MolecularModeling::AtomVector a, MolecularModeling::AtomVector b)
{
    MolecularModeling::AtomVector returned_selection;
    for(MolecularModeling::AtomVector::iterator it = a.begin(); it != a.end(); ++it)
    {
        MolecularModeling::Atom *atom_a = *it;
        for(MolecularModeling::AtomVector::iterator it1 = b.begin(); it1 != b.end(); ++it1)
        {
            MolecularModeling::Atom *atom_b = *it1;
            if (atom_a->GetIndex() == atom_b->GetIndex())
            {
                returned_selection.push_back(atom_a);
            }
        }
    }
    return returned_selection;
}

MolecularModeling::AtomVector selection::GetAtomsin_a_Notin_b_AtomVectors(MolecularModeling::AtomVector a, MolecularModeling::AtomVector b)
{
    MolecularModeling::AtomVector returned_selection;
    bool isCommon = false;
    for(MolecularModeling::AtomVector::iterator it = a.begin(); it != a.end(); ++it)
    {
        isCommon = false;
        MolecularModeling::Atom *atom_a = *it;
        for(MolecularModeling::AtomVector::iterator it1 = b.begin(); it1 != b.end(); ++it1)
        {
            MolecularModeling::Atom *atom_b = *it1;
            if (atom_a->GetIndex() == atom_b->GetIndex())
            {
                isCommon = true;
            }
        }
        if (isCommon==false)
        {
            returned_selection.push_back(atom_a);
        }
    }
    return returned_selection;
}

//MolecularModeling::Atom* selection::FindAtomNeighborThatMatchesQuery(MolecularModeling::Atom *atom, std::string query)
//{
//    MolecularModeling::Atom *selected_atom;
//    std::regex regex(query, std::regex_constants::ECMAScript);
//    MolecularModeling::AtomVector neighbors = atom->GetNode()->GetNodeNeighbors();
//    for (const auto& neighbor : neighbors)
//    {
//        if (std::regex_search(neighbor->GetName(), regex))
//        {
//            std::cout << "Matching neighbor for " << atom->GetName() << " is " << neighbor->GetName() << "\n";
//            selected_atom = &(*neighbor);
//        }
//    }
//    return selected_atom;
//}

//MolecularModeling::AtomVector selection::FindBondedAtomsThatMatchQuery(MolecularModeling::Atom *anomeric_carbon, std::vector<std::regex> atom_name_queries)
//{
//    // Jesus.
//    // Get every atom in Residue1 and Residue2.
//    // Get every atom that matches first atom in query

//    MolecularModeling::AtomVector found_atom_path;
//    bool success = false;
//    int depth = 0;
//     selection::seek_neighbors_with_regex_query(anomeric_carbon, &found_atom_path, depth, atom_name_queries, &success);
////    std::regex first_atom = atom_name_queries.at(0);
////    if (std::regex_search(anomeric_carbon->GetName(), first_atom))
////    {
////        found_atom_path.push_back(anomeric_carbon);
////        ++depth;
////        selection::seek_neighbors_with_regex_query(anomeric_carbon, &found_atom_path, depth, atom_name_queries, &success);
////    }
////    else
////    {
////        depth = 0;
// //       selection::seek_neighbors_with_regex_query(anomeric_carbon, &found_atom_path, depth, atom_name_queries, &success);
////    }
////    std::ptrdiff_t pos = std::distance(atom_name_queries.begin(), std::find(atom_name_queries.begin(), atom_name_queries.end(), old_name_));
////    if(pos >= Names.size()) {
////        //old_name_ not found
////    }
//}

//MolecularModeling::Atom* selection::FindClosestNeighbor

//void selection::seek_neighbors_with_regex_query(MolecularModeling::Atom *current_atom, MolecularModeling::AtomVector *found_atom_path, int depth, std::vector<std::regex> atom_name_queries, bool &success)
//{
//    std::regex current_regex = atom_name_queries.at(depth);
//    MolecularModeling::AtomVector neighbors = current_atom->GetNode()->GetNodeNeighbors();
//    for (const auto& neighbor : neighbors)
//    {
//        if (std::regex_search(neighbor->GetName(), current_regex))
//        {
//            if (depth == 3) // we have made it
//            {
//                *success = true;
//            }
//            else
//            {
//                ++depth;
//                selection::seek_neighbors_with_regex_query(current_atom, depth, atom_name_queries, &success);
//            }
//        }
//    }
//}

// Assumes query is in the format ?_222 or A_222, where 222 is the residue number and ?/A is the chain ID. ? if not specified.
MolecularModeling::Residue* selection::FindResidue(MolecularModeling::Assembly &assembly, const std::string query)
{
    for(auto &residue : assembly.GetAllResiduesOfAssembly())
    {
        std::string id = residue->GetId();
        std::string formatted_query = "_" + query + "_";
        if( id.compare(3, formatted_query.size(), formatted_query) == 0)
        {
            return residue;
        }
    }
    std::cerr << "Residue " << query << " not found in assembly!" << std::endl;
    return nullptr;
}

MolecularModeling::AtomVector selection::FindOtherAtomsWithinMolecule(MolecularModeling::Atom *queryAtom)
{
    MolecularModeling::AtomVector foundAtoms;
    queryAtom->FindConnectedAtoms(foundAtoms);
    return foundAtoms;
}

 // A function that compares atom numbers to see which is higher. Used for sorting residue neighbors.
bool selection::compareAtomNumbers(MolecularModeling::Atom *a1, MolecularModeling::Atom *a2)
{
    int a1Number = (a1->GetName().at(1) - '0'); // The atom "number" is a char. e.g. the 2 in C2.
    int a2Number = (a2->GetName().at(1) - '0');
    if (a1Number < a2Number)
    {
        return true;
    }
    return false;
}

MolecularModeling::ResidueVector selection::SortResidueNeighborsByAcendingConnectionAtomNumber(MolecularModeling::AtomVector atomsConnectedToOtherResidues)
{   // Sort neighbors by lowest index so Residue_linkage labeling matches front end. 
    // Index refers to a 1-6 vs a 1-3 branch, the 1-3 should be labeled with a lower number.

    // First we sort the atoms in the current residue by their atom number. I.e. C2, C6, C4, becomes C2, C4, C6
    // These should be atoms that are connected to other residues

    std::sort(atomsConnectedToOtherResidues.begin(), atomsConnectedToOtherResidues.end(), selection::compareAtomNumbers);

    // Now we want to get the residue of each atom connected to these atoms, but some of these will be in the same residue and should be ignored.
    MolecularModeling::ResidueVector sortedNeighbors; 
    for(auto & sortedAtom : atomsConnectedToOtherResidues)
    {
       // std::cout << "Atoms connected to other residues are " << sortedAtom->GetId() << "\n";
        for(auto & neighboringAtom : sortedAtom->GetNode()->GetNodeNeighbors())
        {   // want to find residue of atoms that are connected to sortedAtom
            if (neighboringAtom->GetResidue()->GetIndex() != sortedAtom->GetResidue()->GetIndex())
            {
            //    std::cout << "Neighbor of that atom which is in another resiude is " << neighboringAtom->GetResidue()->GetId() << "@" << neighboringAtom->GetId() << "\n";
                sortedNeighbors.push_back(neighboringAtom->GetResidue());
            }
        }
    }
    return sortedNeighbors;
}


// After getting connection atoms between two residues, we have the linear path between them e.g.:
// (Residue1 C2 - O7 - C7 - O6 - C6 Residue2), where C2 and C6 are the cycle points. In cases like 2-7 or 2-8 linkages, we have significant branches from this linear path
//                    /
//                   C8-O8
//                  /
//                 C9-O9 
// And we need to set reasonable values for the C9-C8, and C8-C7 dihedral angles.
// This is further complicated by the possibility of "DNeup5Aca2-7[DNeup5Aca2-8]DNeup5Aca2-OH", where one of the branches is part of another linkage.
// This will work for 7/8/9 linked sialic acids, but if we get more heavily branched linkages this will need to change to iteratively find branches from branches
//

void selection::FindEndsOfBranchesFromLinkageAtom(MolecularModeling::Atom* currentAtom, MolecularModeling::Atom* previousAtom, Branch *branch)
{
    branch->ChangeDepth(1);
    currentAtom->SetDescription("VistedByFindEndsOfBranchesFromLinkageAtom");
    bool deadEndAtom = true;
    bool connectsToAnotherResidue = false;
    for (auto &neighbor : currentAtom->GetNode()->GetNodeNeighbors())
    {
        if (neighbor->GetDescription() != "VistedByFindEndsOfBranchesFromLinkageAtom" 
            && *neighbor != *previousAtom 
            && neighbor->GetResidue() == previousAtom->GetResidue()) // Don't explore across residues.
        { 
            if (neighbor->GetNode()->GetNodeNeighbors().size() > 1)
            {
                deadEndAtom = false;
                //std::cout << "At depth " << branch->GetDepth() << " going deeper from " << currentAtom->GetId() << " to " << neighbor->GetId() << "\n";
                FindEndsOfBranchesFromLinkageAtom(neighbor, currentAtom, branch);
                branch->ChangeDepth(-1);
            }
        }
        if (neighbor->GetResidue() != previousAtom->GetResidue())
        {
            connectsToAnotherResidue = true;
        }
    }
    //std::cout << "  Status at " << currentAtom->GetId() << " is deadEndAtom:" << std::boolalpha << deadEndAtom << ", depth: " << branch->GetDepth() << ", connectsToOther: " << connectsToAnotherResidue << std::endl; 
    if (deadEndAtom && !connectsToAnotherResidue && branch->GetDepth() > 1 && branch->AtMaxDepth())
    {
        //std::cout << "      Found a dead end: " << currentAtom->GetId() << " at depth " << branch->GetDepth() << std::endl;
        branch->SetEnd(currentAtom);
    }
    return;
}

std::string selection::GetNonCarbonHeavyAtomNumbered(AtomVector atoms, std::string queryNumber)
{
    for (auto &atom : atoms)
    {    // Assuming atoms like C2, O2 or N2. Nothing else should match.
        std::string number = atom->GetName().substr(1);
        if ( (number == queryNumber) && (atom->GetName().at(0) != 'C') && (atom->GetName().at(0) != 'H'))
        {
            return atom->GetName();
        }
    }
    return "";    
}

AtomVector selection::FindHeavyAtoms(AtomVector query)
{
	AtomVector foundAtoms;
	for (auto &atom : query)
	{
		char element = atom->GetName().at(0);
	    if ( ( element == 'S') || ( element == 'C') || ( element == 'O') || ( element == 'N') || ( element == 'P') )
	    {
	    	foundAtoms.push_back(atom);
	    }
	}
	return foundAtoms;
}

int selection::FindHighestResidueNumber(MolecularModeling::Assembly &ass, std::string inputChainID)
{
    int highest = 0;
    for (auto &residue : FindResiduesWithChainID(ass, inputChainID))
    {
        int resNumInt = std::stoi(residue->GetNumber());
        if (resNumInt > highest)
        {
            highest = resNumInt;
        }
    }
    return highest;
}

std::vector<MolecularModeling::Residue*> selection::FindResiduesWithChainID(MolecularModeling::Assembly &ass, std::string inputChainID)
{
    std::vector<MolecularModeling::Residue*> foundResidues;
    for (auto &residue : ass.GetResidues())
    {
        if (inputChainID == residue->GetChainID())
        {
            foundResidues.push_back(residue);
        }
    }
    return foundResidues;
}
