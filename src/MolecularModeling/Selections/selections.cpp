#include "../../../includes/MolecularModeling/Selections/selections.h"

AtomVector selection::AtomsWithinDistanceOf(MolecularModeling::Atom *query_atom, double distance, AtomVector atoms)
{
    AtomVector atoms_within_distance;
    for(AtomVector::iterator it1 = atoms.begin(); it1 != atoms.end(); ++it1)
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
void selection::FindAtomsConnectingResidues(Atom *current_atom, Residue *second_residue, AtomVector *connecting_atoms, bool *found_neighbor)
{
    current_atom->SetDescription("VisitedByFindAtomsConnectingResidues");
    AtomVector neighbors = current_atom->GetNode()->GetNodeNeighbors();
    for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); ++it1)
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
bool selection::FindCyclePoint(Atom *previous_atom, Atom *current_atom, AtomVector *atom_path, bool *found_cycle_point, Atom *&cycle_point)
{
    // I wish there was a more solid way to check element type. But here we go, I definitely don't want cycles involving hydrogens (they only every form one bond, unless
    // bond by distance has bonded them).
    if ( current_atom->GetName().at(0) != 'H')
    {
        // Need this to explore everything. It will find same cycle point more than once, but that doesn't matter.
        current_atom->SetDescription("VisitedByFindCyclePoint");
//        std::cout << "Checking neighbors of " << current_atom->GetName() << "\n";
//        std::cout << "Found cycle points is currently: " << std::boolalpha << *found_cycle_point << std::endl;

        atom_path->push_back(current_atom);
        AtomVector neighbors = current_atom->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); ++it1)
        {
            MolecularModeling::Atom *neighbor = *it1;
            // If not previous atom and not from a different residue
            if ( (neighbor->GetIndex() != previous_atom->GetIndex()) && (current_atom->GetResidue()->GetId().compare(neighbor->GetResidue()->GetId())==0))
                //if ( neighbor->GetIndex() != previous_atom->GetIndex() ) // Good for testing multiple cycles
            {
//                std::cout << "Coming from previous atom " << previous_atom->GetName() << " we see bonding: " << current_atom->GetName() << "->" << neighbor->GetName() << "\n";
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
bool selection::FindRotationPointsForNonCycles(Atom *previous_atom, Atom *current_atom, AtomVector *rotation_points)
{
    AtomVector neighbors = current_atom->GetNode()->GetNodeNeighbors();
    AtomVector intra_node_neighbors;
    for (auto &neighbor : neighbors)
    { // If not previous atom and not from a different residue
        if ( (neighbor->GetIndex() != previous_atom->GetIndex()) && (current_atom->GetResidue()->GetId().compare(neighbor->GetResidue()->GetId())==0))
        {
            intra_node_neighbors.push_back(neighbor);
        }
    }
    bool found_rotation_point = true;
    for (auto &neighbor : intra_node_neighbors)
    {
        AtomVector neighbor_neighbors = neighbor->GetNode()->GetNodeNeighbors();
        AtomVector intra_node_neighbor_neighbors;
        for (auto &neighbor_neighbor : neighbor_neighbors)
        {
            if ( (neighbor_neighbor->GetIndex() != current_atom->GetIndex()) && (current_atom->GetResidue()->GetId().compare(neighbor_neighbor->GetResidue()->GetId())==0))
            {
                intra_node_neighbor_neighbors.push_back(neighbor);
            }
        }
        if (intra_node_neighbor_neighbors.size() > 0 )
        {
            found_rotation_point = false;
        }
    }
    if ( found_rotation_point ) // Deadend Chitchatchoes
    {
        rotation_points->push_back(current_atom);
    }
    else // Keep looking for the ends
    {
        for (auto &neighbor : intra_node_neighbors)
        {
            selection::FindRotationPointsForNonCycles(current_atom, neighbor, rotation_points);
        }
    }
    return (!rotation_points->empty()); // return false if rotation_points is empty
}


bool selection::FindPathBetweenTwoAtoms(Atom *current_atom, Atom *target_atom, AtomVector *atom_path, bool *found)
{
    //atom_path->push_back(current_atom);
    current_atom->SetDescription("VistedByFindPathBetweenTwoAtoms");
    AtomVector neighbors = current_atom->GetNode()->GetNodeNeighbors();
    for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); ++it1)
    {
        MolecularModeling::Atom *neighbor = *it1;
        if(neighbor->GetIndex() == target_atom->GetIndex())
        {
            *found = true;
            atom_path->push_back(neighbor);
        }
        // If not found && not previously visited atom && ( if neighbor residue is current residue || target_atom residue)
        if ( (*found == false) && (neighbor->GetDescription().compare("VistedByFindPathBetweenTwoAtoms")!=0) && ((current_atom->GetResidue()->GetId().compare(neighbor->GetResidue()->GetId())==0) || (target_atom->GetResidue()->GetId().compare(neighbor->GetResidue()->GetId())==0)) )
        {
          //  std::cout << "STEEEPER\n";
            selection::FindPathBetweenTwoAtoms(neighbor, target_atom, atom_path, found);
        }
    }
    if(*found)
    {
        atom_path->push_back(current_atom);
    }
    return *found;
}

void selection::ClearAtomDescriptions(Residue *residue)
{
    AtomVector atoms = residue->GetAtoms();
    for(AtomVector::iterator it = atoms.begin(); it != atoms.end(); ++it)
    {
        Atom *atom = *it;
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

  2) fail to find both connections of the middle Res here:
 Res1 __
        \__ Res4
 Res2 __/
*/
AtomVector selection::FindCyclePoints(Atom *atom)
{
//    std::cout << "Entered FindCyclePoints with " << atom->GetId() << std::endl;
    AtomVector rotation_points;
    AtomVector atom_path;
    bool found = false;
    if(atom->GetResidue()->CheckIfProtein())
    {
        Atom *caAtom = atom->GetResidue()->GetAtom("CA");
        // Find any cycle points. Note: starting from connection point to other residue
        Atom *cycle_point;
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
    else // Non protein.
    {
        Atom *rotation_point;
        atom_path.clear();
        found = false;
        //Find path to first cycle atom, i.e. anomeric carbon
        if(selection::FindCyclePoint(atom, atom, &atom_path, &found, rotation_point))
        {
            rotation_points.push_back(rotation_point);
        }
        // Ok, deal with non-protein non-cycles
        else
        { // Look for atom(s) with neighbors with no other neighbors within residue
            selection::FindRotationPointsForNonCycles(atom, atom, &rotation_points);
        }
    }
    return rotation_points;
}

// Find a neighbor of a cycle point to define the dihedral. Must not be in path used to come to this cycle_point.
// Oh gawd this code is horrible.
// The logic for selecting the atom to define a dihedral is messy.
// Comparing strings is messy when I really care about the number of e.g. C2 Vs O5, but need to factor in a C2 vs O comparision
Atom* selection::FindCyclePointNeighbor(const AtomVector atom_path, Atom *cycle_point)
{
    Atom *selected_neighbor;
    if (cycle_point->GetName().compare("CA")==0) // This is a protein, and we always want the N atom.
    {
        selected_neighbor = cycle_point->GetResidue()->GetAtom("N");
    }
    else if (cycle_point->GetName().compare("C2")==0) // This is an ulose (e.g. Sia), and we always want the C1 atom.
    {
        selected_neighbor = cycle_point->GetResidue()->GetAtom("C1");
    }
    else if (cycle_point->GetName().compare("C1")==0)
    {
        selected_neighbor = cycle_point->GetResidue()->GetAtom("C2"); // We can always set phi to 180 this way, regardless of alpha/beta
    }
    else // This bit is overdone now, as I was looking for higher numbered atoms of C1, but now I know I always want C2, so put that above.
    {
        AtomVector neighbors = cycle_point->GetNode()->GetNodeNeighbors();
        // Ok must first get a list of neighbors that weren't in the connection path
        AtomVector good_neighbors; // Couldn't think of a better name. Everybody needs these.
        for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); ++it1)
        {
            Atom *neighbor = *it1;
            if ( ! (std::find(atom_path.begin(), atom_path.end(), neighbor) != atom_path.end()) ) // If we've NOT been at this atom on way to cycle point
            {
                if ( neighbor->GetName().at(0) != 'H' ) // Don't find hydrogens. Later we swap out to use a hydrogen to define a dihedral, but that's a very specific one.
                {
                    good_neighbors.push_back(neighbor);
                }
            }
        }
        selected_neighbor = good_neighbors.at(0); // Set to any to start. If there are not good_neighbors then you deserve to crash and burn
       // std::cout << "Good neighbors are: ";
        for(AtomVector::iterator it1 = good_neighbors.begin(); it1 != good_neighbors.end(); ++it1)
        {
            Atom *neighbor = *it1;
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
   // std::cout << "Returning with neighbor: " << selected_neighbor->GetName() << "\n";
    return selected_neighbor;
}

//Atom* selection::FindAtomNeighborThatMatchesQuery(Atom *atom, std::string query)
//{
//    Atom *selected_atom;
//    std::regex regex(query, std::regex_constants::ECMAScript);
//    AtomVector neighbors = atom->GetNode()->GetNodeNeighbors();
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

//AtomVector selection::FindBondedAtomsThatMatchQuery(Atom *anomeric_carbon, std::vector<std::regex> atom_name_queries)
//{
//    // Jesus.
//    // Get every atom in Residue1 and Residue2.
//    // Get every atom that matches first atom in query

//    AtomVector found_atom_path;
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

//Atom* selection::FindClosestNeighbor

//void selection::seek_neighbors_with_regex_query(Atom *current_atom, AtomVector *found_atom_path, int depth, std::vector<std::regex> atom_name_queries, bool &success)
//{
//    std::regex current_regex = atom_name_queries.at(depth);
//    AtomVector neighbors = current_atom->GetNode()->GetNodeNeighbors();
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


