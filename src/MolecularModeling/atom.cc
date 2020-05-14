#include "../../includes/MolecularModeling/atom.hpp"
#include "../../includes/MolecularModeling/quantommechanicatom.hpp"
#include "../../includes/MolecularModeling/moleculardynamicatom.hpp"
#include "../../includes/MolecularModeling/dockingatom.hpp"
#include "../../includes/MolecularModeling/oligosaccharidedetectionatom.hpp"
#include "../../includes/MolecularModeling/atomnode.hpp"
#include "../../includes/MolecularModeling/residue.hpp"
#include "../../includes/GeometryTopology/coordinate.hpp"
#include "cmath"
#include "algorithm"

using MolecularModeling::Atom;
//////////////////////////////////////////////////////////
//                       CONSTRUCTORS                   //
//////////////////////////////////////////////////////////
Atom::Atom()
{
	this->index_ = this->generateAtomIndex();
	// Call to private helper function.
    this->SetAttributes(NULL, "", GeometryTopology::CoordinateVector(), "", "", "", std::vector<AtomNode*>(), "", false, "");
	this->SetBFactor(0);
} // end Default Constructor

Atom::Atom(MolecularModeling::Residue* residue, std::string name, GeometryTopology::CoordinateVector coordinates)
{
	this->index_ = this->generateAtomIndex();
	std::stringstream ss;
	ss << name << "_" << this->GetIndex() << "_" << residue->GetName() << "_?_1_?_?_1";
	// Call to private helper function.
	this->SetAttributes(residue, name, coordinates, "", "", "", std::vector<AtomNode*>(), ss.str(), false, "");
	//this->SetBFactor(residue[atom]->GetBFactor())
} // end Constructor

Atom::Atom(MolecularModeling::Residue* residue, std::string name, GeometryTopology::Coordinate coordinate)
{
	this->index_ = this->generateAtomIndex();
	std::stringstream ss;
	ss << name << "_" << this->GetIndex() << "_" << residue->GetName() << "_?_1_?_?_1";
    GeometryTopology::CoordinateVector coordinates;
	coordinates.push_back(new GeometryTopology::Coordinate(coordinate.GetX(), coordinate.GetY(), coordinate.GetZ()));
	this->SetAttributes(residue, name, coordinates, "", "", "", std::vector<AtomNode*>(), ss.str(), false, "");
} // end Constructor

Atom::Atom(const Atom* atom)
{
    this->index_ = this->generateAtomIndex();
    this->Copy(atom);
} // end Copy Constructor(*)


/*! \todo  Figure out why the constructor below gets the error below and fix it.

src/MolecularModeling/atom.cc: In copy constructor 'MolecularModeling::Atom::Atom(const MolecularModeling::Atom&)':
src/MolecularModeling/atom.cc:51:1: warning: base class 'class MolecularModeling::MolecularDynamicAtom' should be explicitly initialized in the copy constructor [-Wextra]
 Atom::Atom(const Atom& atom)
 ^
*/
// Added by Ayush
Atom::Atom(const Atom& atom)
{
	this->index_ = this->generateAtomIndex();
	this->Copy(&atom);
} // end Copy Constructor(&)

//////////////////////////////////////////////////////////
//                       DESTRUCTOR                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                         ACCESSORS                    //
//////////////////////////////////////////////////////////
MolecularModeling::Residue* Atom::GetResidue() const
{
	return this->residue_;
} // end GetResidue

std::string Atom::GetName() const
{
    return this->name_;
} // end GetName

GeometryTopology::CoordinateVector Atom::GetCoordinates() const
{
	return this->coordinates_;
} // end GetCoordinates

// Most of the time, you just want the first coordinate
GeometryTopology::Coordinate* Atom::GetCoordinate()
{
    return this->coordinates_.at(0);
} // end GetCoordinate


std::string Atom::GetChemicalType() const
{
 	return this->chemical_type_;
} // end GetChemicalType

std::string Atom::GetDescription() const
{
	return this->description_;
} // end GetDescription

std::string Atom::GetElementSymbol() const
{
	return this->element_symbol_;
} // end GetElementSymbol

MolecularModeling::AtomNode* Atom::GetNode(int index) const
{
	int max_index = this->nodes_.size() -1;
	if (max_index >= index){
	    return this->nodes_.at(index);
	}
        return NULL;
} // end GetNode

std::vector<MolecularModeling::AtomNode*> Atom::GetNodes() const
{
        return this->nodes_;
}


std::string Atom::GetId() const
{
	return this->id_;
} // end GetId

bool Atom::GetIsRing() const
{
	return this->is_ring_;
} // end GetIsRing

bool Atom::GetIsExocyclicCarbon() const
{
	return this->is_exocyclic_C_;
} // end GetIsExocyclicCarbon

unsigned long long Atom::GetIndex() const
{
	return this->index_;
} // end GetIndex

//Added by ayush on 13/11/17 for molecules in assembly to set the atom type as an attribute like O,H, etc.

std::string Atom::GetAtomType() const
{
	return this->atom_type_;
}

//Added by Dave on 03/23/18 for adding B Factor to ontology

float Atom::GetBFactor() const
{
	return b_factor_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void Atom::SetResidue(MolecularModeling::Residue* residue)
{
	this->residue_ = residue;
} // end SetResidue

void Atom::SetName(std::string name)
{
	this->name_ = name;
} // end SetName

void Atom::SetCoordinates(GeometryTopology::CoordinateVector coordinates)
{
	// First need to delete any previous Coordinates, so we don't have any memory leaks.
    for(GeometryTopology::CoordinateVector::iterator it = this->coordinates_.begin(); it != this->coordinates_.end(); it++ )
	{
		GeometryTopology::Coordinate* coordinate = (*it);
		if(coordinate != NULL)
		{
			delete coordinate;
			coordinate = NULL;
		}
	}
	this->coordinates_.clear();
    for(GeometryTopology::CoordinateVector::iterator it = coordinates.begin(); it != coordinates.end(); it++) {
		this->coordinates_.push_back(*it);
	}
}// end SetCoordinates

void Atom::AddCoordinate(GeometryTopology::Coordinate* coordinate)
{
	this->coordinates_.push_back(coordinate);
} // end AddCoordinate

void Atom::AddNode(MolecularModeling::AtomNode* node)
{
        this->nodes_.push_back(node);
}


void Atom::SetChemicalType(std::string chemical_type)
{
	this->chemical_type_ = chemical_type;
} // end SetChemicalType

void Atom::SetDescription(std::string description)
{
	this->description_ = description;
} // end SetDescription

void Atom::SetElementSymbol(std::string element_symbol)
{
	this->element_symbol_ = element_symbol;
} // end SetElementSymbol

void Atom::SetNodes(std::vector<MolecularModeling::AtomNode*> nodes)
{
	this->nodes_ = nodes;
} // end SetNode

void Atom::SetNode(MolecularModeling::AtomNode* node)
{
        if (this->nodes_.empty()){
	    this->AddNode(node);
	}
}

void Atom::SetId(std::string id)
{
	this->id_ = id;
} // end SetId

void Atom::SetIsRing(bool is_ring)
{
	this->is_ring_ = is_ring;
}
//Added by ayush on 13/11/17 for molecules in assembly

void Atom::SetIsExocyclicCarbon(bool is_exocyclic_C)
{
	this->is_exocyclic_C_ = is_exocyclic_C;
}

void Atom::SetAtomType(std::string atom_type)
{
	this->atom_type_ = atom_type;
}

//Added by Dave on 03/23/18 for adding B Factor to ontology

void Atom::SetBFactor(float b_factor)
{
	this->b_factor_ = b_factor;
}


unsigned long long Atom::generateAtomIndex() {
	static unsigned long long s_AtomIndex = 0; // static keyword means it is created only once and persists beyond scope of code block.
	return s_AtomIndex++; // makes copy of s_AtomIndex, increments the real s_AtomIndex, then returns the value in the copy
} // end generateAtomIndex

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
void Atom::FindConnectedAtoms(AtomVector &visitedAtoms, int coord_index)
{
    visitedAtoms.push_back(this);
	AtomVector neighbors = this->GetNode(coord_index)->GetNodeNeighbors();
	bool alreadyVisited = false;
	for(AtomVector::iterator neighbor = neighbors.begin(); neighbor != neighbors.end(); neighbor++)
	{
		alreadyVisited = false; // reset for each neighbor
        for(AtomVector::iterator visitedAtom = visitedAtoms.begin(); visitedAtom != visitedAtoms.end(); visitedAtom++)
		{
			if((*neighbor)->GetIndex() == (*visitedAtom)->GetIndex())
			alreadyVisited = true;
		}
		if(!alreadyVisited)
		{
            (*neighbor)->FindConnectedAtoms(visitedAtoms, coord_index); // recursive function call
		}
	}
} // end FindConnectedAtoms

double Atom::GetDistanceToAtom(Atom* otherAtom)
{
    return GetDistanceToCoordinate(otherAtom->GetCoordinate());
} // end GetDistanceToAtom

double Atom::GetDistanceToCoordinate(GeometryTopology::Coordinate* coordinate)
{
    double x = (this->GetCoordinate()->GetX() - coordinate->GetX());
    double y = (this->GetCoordinate()->GetY() - coordinate->GetY());
    double z = (this->GetCoordinate()->GetZ() - coordinate->GetZ());
    return std::abs(sqrt((x * x) + (y * y) + (z * z)));
} // end GetDistanceToCoordinate

bool Atom::CheckIfOtherAtomIsWithinBondingDistance(Atom* otherAtom)
{
    if (this->GetIndex() == otherAtom->GetIndex())
    {
        std::cerr << "Warning have just checked distance between an atom and itself!" << std::endl;
        return true;
    }
    bool withinDistance = false;
    if (std::abs(this->GetCoordinate()->GetX() - otherAtom->GetCoordinate()->GetX()) < gmml::maxCutOff)
    {
        if (std::abs(this->GetCoordinate()->GetY() - otherAtom->GetCoordinate()->GetY()) < gmml::maxCutOff)
        {
            if (std::abs(this->GetCoordinate()->GetZ() - otherAtom->GetCoordinate()->GetZ()) < gmml::maxCutOff)
            {
                //If each dimension is within cutoff, then calculate 3D distance
                if (this->GetDistanceToAtom(otherAtom) < gmml::maxCutOff)
                {
                    withinDistance = true;
                }
            }
        }
    }
    return withinDistance;
}

std::string Atom::DetermineChirality() //Added by Yao 08/26/3019 Return values are R,S,A. A = achiral
{
    AtomVector primary_neighbors = this->GetNode()->GetNodeNeighbors();
    for (unsigned int i = 0; i < primary_neighbors.size(); i++){
	//std::cout << "Primary neighbor: " << primary_neighbors[i]->GetName() << std::endl;
    }
    if (primary_neighbors.size() >= 5){
	//std::cout << "Cannot determine chirality for atoms with 5 or more bonds." << std::endl;
	return "A";
    }

    else if (primary_neighbors.size() == 4){
        std::vector<int> ranks(primary_neighbors.size(), 1);
	AtomVector ordered_primary_neighbors = this->GetRankedPrimaryNeighbors(ranks);

	std::set<int> unique_ranks;
        for (std::vector<int>::iterator int_it = ranks.begin(); int_it != ranks.end(); int_it++){
	    //std::cout << "Atom " << primary_neighbors[std::distance(ranks.begin(), int_it)]->GetName() << " has rank " << *int_it << std::endl;
            unique_ranks.insert(*int_it);
        }

	if (unique_ranks.size() == ordered_primary_neighbors.size()){
	    return this->DetermineRSAssignment(ordered_primary_neighbors, NULL);
	}
	else {
	    //std::cout << "Primary neighbors not unique" << std::endl;
	    return "A";
	}
    }

    else if (primary_neighbors.size() == 3){
	bool sp2 = false;
	//Allow up to 11.47 deg from parallel. This cutoff is just slightly above side group N2 in prep file, which is 0.98491, or 9.97 deg from ideal trigonal planar.
	if (std::abs(std::cos(GetDihedral(primary_neighbors[0], this, primary_neighbors[1], primary_neighbors[2]))) > 0.98){
            sp2 = true; //if improper torsion is close to 0 or 180 deg, then carbon neighbor is planar, hence sp2 carbon.
        }
	//std::cout << "Sp2 judgment torsion is: " << primary_neighbors[0]->GetName() << "-" <<  this->GetName() << "-" <<  primary_neighbors[1]->GetName() << "-" <<  primary_neighbors[2]->GetName()
		//<< ": " << std::abs(std::cos(GetDihedral(primary_neighbors[0], this, primary_neighbors[1], primary_neighbors[2]))) << std::endl;
	if (sp2){
	    //std::cout << "This atom is trigonal planar/sp2" << std::endl;
	    return "A";
	}
	else{

            std::vector<int> ranks(primary_neighbors.size(), 1);
	    AtomVector ordered_primary_neighbors = this-> GetRankedPrimaryNeighbors(ranks);

            std::set<int> unique_ranks;
            for (std::vector<int>::iterator int_it = ranks.begin(); int_it != ranks.end(); int_it++){
                unique_ranks.insert(*int_it);
            }

            if (unique_ranks.size() == ordered_primary_neighbors.size()){
	        Atom fake_hydrogen = this->PlaceFakeHydrogen();
		fake_hydrogen.SetName("FakeH");
		//std::cout << "Fake hydrogen placed." << std::endl;
                return this->DetermineRSAssignment(ordered_primary_neighbors, &fake_hydrogen);
            }
            else {
                //std::cout << "Primary neighbors not unique" << std::endl;
                return "A";
            }

	}
    }

    else {
	return "A";
    }
}

Atom Atom::PlaceFakeHydrogen()
{
    Atom fake_hydrogen = Atom();
    GeometryTopology::CoordinateVector coordinates;
    GeometryTopology::Coordinate* coord = new GeometryTopology::Coordinate();

    double total_x=0.0;
    double total_y=0.0;
    double total_z=0.0;

    MolecularModeling::AtomVector primary_neighbors = this->GetNode()->GetNodeNeighbors();
    for (AtomVector::iterator atom_it = primary_neighbors.begin(); atom_it != primary_neighbors.end(); atom_it++){
	//std::cout << "Bonded neighbor: " << (*atom_it)->GetName() << std::endl;
	//std::cout << "Coordinate" << (*atom_it)->GetCoordinates().at(0)->GetX() << "," << (*atom_it)->GetCoordinates().at(0)->GetY() << "," << (*atom_it)->GetCoordinates().at(0)->GetZ() << std::endl;
	GeometryTopology::Coordinate bond;
	bond.SetX((*atom_it)->GetCoordinates().at(0)->GetX() - this->GetCoordinates().at(0)->GetX());
	bond.SetY((*atom_it)->GetCoordinates().at(0)->GetY() - this->GetCoordinates().at(0)->GetY());
	bond.SetZ((*atom_it)->GetCoordinates().at(0)->GetZ() - this->GetCoordinates().at(0)->GetZ());
	bond.Normalize();

        total_x += (bond.GetX());
        total_y += (bond.GetY());
        total_z += (bond.GetZ());
    }

    coord->SetX(this->GetCoordinates().at(0)->GetX() + (-1 * total_x));
    coord->SetY(this->GetCoordinates().at(0)->GetY() + (-1 * total_y));
    coord->SetZ(this->GetCoordinates().at(0)->GetZ() + (-1 * total_z));

    //std::cout << "FakeH X: " << coord->GetX() << std::endl;
    //std::cout << "FakeH Y: " << coord->GetY() << std::endl;
    //std::cout << "FakeH Z: " << coord->GetZ() << std::endl;

    coordinates.push_back(coord);
    fake_hydrogen.SetCoordinates(coordinates);

    return fake_hydrogen;
}

MolecularModeling::AtomVector Atom::GetRankedPrimaryNeighbors(std::vector<int>& ranks)
{
    std::map<std::vector<int>, std::vector<int> > duplicate_value_indices_versus_higher_rank_indices = this-> ComparePrimaryNeighbors(ranks);
    AtomVector primary_neighbors = this->GetNode()->GetNodeNeighbors();
    AtomVector ordered_primary_neighbors;

    if (duplicate_value_indices_versus_higher_rank_indices.empty()){
        std::set<int> unique_ranks;
        for (std::vector<int>::iterator int_it = ranks.begin(); int_it != ranks.end(); int_it++){
            unique_ranks.insert(*int_it);
        }
        for (std::set<int>::iterator set_it = unique_ranks.begin(); set_it != unique_ranks.end(); set_it++){
            int rank = *set_it;
            for (std::vector<int>::iterator int_it = ranks.begin(); int_it != ranks.end(); int_it++){
                if (rank == *int_it){
                    int index = std::distance(ranks.begin(), int_it);
                    AtomVector::iterator atom_it = primary_neighbors.begin();
                    std::advance(atom_it, index);
                    Atom* primary_neighbor = *atom_it;
                    ordered_primary_neighbors.push_back(primary_neighbor);
                }
            }
        }

	return ordered_primary_neighbors;
    }

    else {
        bool dead_end = true;
	std::vector<unsigned int> empty_count (duplicate_value_indices_versus_higher_rank_indices.size(), 0);
        for (std::map<std::vector<int>, std::vector<int> >::iterator mapit = duplicate_value_indices_versus_higher_rank_indices.begin(); mapit !=
            duplicate_value_indices_versus_higher_rank_indices.end(); mapit++){

	    int map_index = std::distance(duplicate_value_indices_versus_higher_rank_indices.begin(), mapit);
	    //std::cout << "Map index : " << map_index << std::endl;
	    std::vector<int> duplicate_indices = mapit->first;
            std::vector<int> higher_rank_indices = mapit->second;
            for (std::vector<int>::iterator int_it = duplicate_indices.begin(); int_it != duplicate_indices.end(); int_it++){
	        int index = *int_it;
                AtomVector::iterator atomit = primary_neighbors.begin();
                std::advance(atomit, index);
                Atom* atom = *(atomit);
                AtomVector neighbors = atom->GetNode()->GetNodeNeighbors();
                neighbors.erase(std::find(neighbors.begin(), neighbors.end(),this));

	        if (neighbors.empty()){
		    empty_count[map_index]++;
	        }
	    }
	    //std::cout << "Empty count is: " << empty_count[map_index] << std::endl;

	}

        for (std::map<std::vector<int>, std::vector<int> >::iterator mapit = duplicate_value_indices_versus_higher_rank_indices.begin(); mapit !=
          duplicate_value_indices_versus_higher_rank_indices.end(); mapit++){

	    int map_index = std::distance(duplicate_value_indices_versus_higher_rank_indices.begin(), mapit);
	    std::vector<int> duplicate_indices = mapit->first;
	    std::vector<int> higher_rank_indices = mapit->second;
            for (std::vector<int>::iterator int_it = duplicate_indices.begin(); int_it != duplicate_indices.end(); int_it++){
                int index = *int_it;
                AtomVector::iterator atomit = primary_neighbors.begin();
                std::advance(atomit, index);
                Atom* atom = *(atomit);
                AtomVector neighbors = atom->GetNode()->GetNodeNeighbors();
                neighbors.erase(std::find(neighbors.begin(), neighbors.end(),this));
                if (neighbors.empty() && empty_count[map_index] < duplicate_indices.size()){
		    //Although neighbors have identical ranking at the moment, if a neighbor has no downstream neighbor, it is smaller than other neighbors.
		    //std::cout << "Adjusting rank, which shouldn't happen." << std::endl;
		    ranks[index]++;
		    for (std::vector<int>::iterator int_it2 = higher_rank_indices.begin(); int_it2 != higher_rank_indices.end(); int_it2++){
			ranks[*int_it2]++;
		    }
                    //dead_end = true;
                }
		else {
		    dead_end = false;
		}

            }
        }

	for (unsigned int s = 0; s < ranks.size(); s++){
	    //std::cout << "Test before recursion, corresponding atom is: " << primary_neighbors[s]->GetName() << std::endl;
	    //std::cout << "Test before recusion rank: " << ranks[s] << std::endl;
	}
	std::set<int> unique_ranks2;
	for (std::vector<int>::iterator int_it = ranks.begin(); int_it != ranks.end(); int_it++){
	    unique_ranks2.insert(*int_it);
	}
        if (!dead_end && unique_ranks2.size() != ranks.size()){ //If not at a dead end yet and ranking is not unique, need to start recursive comparison.
            AtomVector visited_atoms;
            visited_atoms.insert(visited_atoms.end(), primary_neighbors.begin(), primary_neighbors.end());
            visited_atoms.push_back(this);

            std::map<Atom*, std::vector<AtomVector> > comparison_progress_tracker;
            this->InitializeComparisonTracker(duplicate_value_indices_versus_higher_rank_indices, visited_atoms, comparison_progress_tracker);
            this->RecursivelyCompareBranches(duplicate_value_indices_versus_higher_rank_indices, ranks, visited_atoms, comparison_progress_tracker);
        }

        for (unsigned int a = 0; a < primary_neighbors.size(); a++){
            //std::cout << "Name: " << primary_neighbors[a]->GetName() << std::endl;
            //std::cout << "Rank: " << ranks[a] << std::endl;
        }
        //std::cout << "Done chirality" << std::endl;

        std::set<int> unique_ranks;
        for (std::vector<int>::iterator int_it = ranks.begin(); int_it != ranks.end(); int_it++){
            unique_ranks.insert(*int_it);
        }

        AtomVector ordered_primary_neighbors;
        for (std::set<int>::iterator set_it = unique_ranks.begin(); set_it != unique_ranks.end(); set_it++){
            int rank = *set_it;
            for (std::vector<int>::iterator int_it = ranks.begin(); int_it != ranks.end(); int_it++){
                if (rank == *int_it){
                    int index = std::distance(ranks.begin(), int_it);
                    AtomVector::iterator atom_it = primary_neighbors.begin();
                    std::advance(atom_it, index);
                    Atom* primary_neighbor = *atom_it;
                    ordered_primary_neighbors.push_back(primary_neighbor);
                }
            }
        }

        return ordered_primary_neighbors;
    }
}

std::string Atom::DetermineRSAssignment(MolecularModeling::AtomVector& ordered_primary_neighbors, Atom* fake_hydrogen)
{
    Atom* smallest_atom = NULL;
    Atom* second_smallest = NULL;
    Atom* third_smallest = NULL;
    Atom* largest = NULL;

    if (fake_hydrogen != NULL){
	smallest_atom = fake_hydrogen;
	second_smallest = ordered_primary_neighbors[2];
	third_smallest = ordered_primary_neighbors[1];
	largest = ordered_primary_neighbors[0];
    }
    else{
	smallest_atom = ordered_primary_neighbors[3];
	second_smallest = ordered_primary_neighbors[2];
	third_smallest = ordered_primary_neighbors[1];
	largest = ordered_primary_neighbors[0];
    }

    //std::cout << "Largest: " << largest->GetName() << std::endl;
    //std::cout << "Third smallest: " << third_smallest->GetName() << std::endl;
    //std::cout << "Second smallest: " << second_smallest->GetName() << std::endl;
    //std::cout << "smalleest: " << smallest_atom->GetName() << std::endl;
    double torsion1 = this->GetDihedral(largest, this, smallest_atom, third_smallest);
    double torsion2 = this->GetDihedral(third_smallest, this, smallest_atom, second_smallest);

    if (torsion1 > 0 && torsion2 > 0){
	//std::cout << "Clockwise" << std::endl;
	return "R";
    }
    else if (torsion1 < 0 && torsion2 < 0){
	//std::cout << "Counterclockwise" << std::endl;
	return "S";
    }
    else {
	//std::cout << "Don't know what happened." << std::endl;
	//std::cout << "Torsion 1: " << torsion1 << " Torsion 2: " << torsion2 << std::endl;
	return "A";
    }

}

double Atom::GetDihedral(MolecularModeling::Atom *atom1, MolecularModeling::Atom *atom2, MolecularModeling::Atom *atom3, MolecularModeling::Atom *atom4)
{
    double current_dihedral = 0.0;
    GeometryTopology::Coordinate* a1 = atom1->GetCoordinates().at(0);
    GeometryTopology::Coordinate* a2 = atom2->GetCoordinates().at(0);
    GeometryTopology::Coordinate* a3 = atom3->GetCoordinates().at(0);
    GeometryTopology::Coordinate* a4 = atom4->GetCoordinates().at(0);

    GeometryTopology::Coordinate b1 = a2;
    b1.operator -(*a1);
    GeometryTopology::Coordinate b2 = a3;
    b2.operator -(*a2);
    GeometryTopology::Coordinate b3 = a4;
    b3.operator -(*a3);
    GeometryTopology::Coordinate b4 = b2;
    b4.operator *(-1);

    GeometryTopology::Coordinate b2xb3 = b2;
    b2xb3.CrossProduct(b3);

    GeometryTopology::Coordinate b1_m_b2n = b1;
    b1_m_b2n.operator *(b2.length());

    GeometryTopology::Coordinate b1xb2 = b1;
    b1xb2.CrossProduct(b2);

    current_dihedral = atan2(b1_m_b2n.DotProduct(b2xb3), b1xb2.DotProduct(b2xb3));
    //return current_dihedral /3.1415 * 180;
    return current_dihedral;
}

void Atom::InitializeComparisonTracker(std::map<std::vector<int>, std::vector<int> >& duplicate_value_indices_versus_higher_rank_indices, MolecularModeling::AtomVector& visited_atoms,
  std::map<Atom*, std::vector<MolecularModeling::AtomVector> >& comparison_progress_tracker)
{
    AtomVector primary_neighbors = this->GetNode()->GetNodeNeighbors();
    for (std::map<std::vector<int>, std::vector<int> >::iterator mapit = duplicate_value_indices_versus_higher_rank_indices.begin(); mapit !=
        duplicate_value_indices_versus_higher_rank_indices.end(); mapit++){
        std::vector<int> duplicate_indices = mapit->first;
        for (std::vector<int>::iterator int_it = duplicate_indices.begin(); int_it != duplicate_indices.end(); int_it++){
            int index = *int_it;
            AtomVector::iterator primary_neighbor_it = primary_neighbors.begin();
            std::advance(primary_neighbor_it, index);
            Atom* primary_neighbor = *primary_neighbor_it;

            std::vector<AtomVector> initial_branches;
            AtomVector downstream_neighbors;
            AtomVector secondary_neighbors = primary_neighbor->GetNode()->GetNodeNeighbors();

            for (AtomVector::iterator atomit = secondary_neighbors.begin(); atomit != secondary_neighbors.end(); atomit++){
                Atom* secondary_neighbor = *atomit;
                if (std::find(visited_atoms.begin(), visited_atoms.end(), secondary_neighbor) == visited_atoms.end()){
                    downstream_neighbors.push_back(secondary_neighbor);
                }
            }

            initial_branches.push_back(downstream_neighbors);
            comparison_progress_tracker[primary_neighbor] = initial_branches;
        }
    }
}

std::map<std::vector<int>, std::vector<int> > Atom::ComparePrimaryNeighbors(std::vector<int>& ranks)
{
    AtomVector primary_neighbors = this->GetNode()->GetNodeNeighbors();

    AtomVector one_time_use_vector = AtomVector(); //No visited atoms to block. We're at the root.
    std::multimap<int, Atom*> primary_neighbor_atomic_number_map = GetAtomicNumbersOfNeighbors(primary_neighbors, one_time_use_vector);

    std::set<int> unique_element_numbers; //Utilize the unique property of a set to sort unique atomic numbers(key) in multimap.
    for (std::multimap<int, Atom*>::iterator mapit = primary_neighbor_atomic_number_map.begin(); mapit != primary_neighbor_atomic_number_map.end(); mapit++){
        int element_number = mapit->first;
        unique_element_numbers.insert(element_number);
    }

    for (std::multimap<int, Atom*>::iterator mapit = primary_neighbor_atomic_number_map.begin(); mapit != primary_neighbor_atomic_number_map.end(); mapit++){
        int atomic_number = mapit->first;
        Atom* primary_neighbor = mapit->second;
	//std::cout << "Test comparing primary neighbor: " << primary_neighbor->GetName() << std::endl;
        int rank = std::distance(unique_element_numbers.rbegin(), std::find(unique_element_numbers.rbegin(), unique_element_numbers.rend(), atomic_number)) + 1;
	//std::cout << "Rank is: " << rank << std::endl;
        int index = std::distance(primary_neighbors.begin(), std::find(primary_neighbors.begin(), primary_neighbors.end(), primary_neighbor));
        ranks[index] = rank;
    }

    std::map<std::vector<int>, std::vector<int> > duplicate_value_indices_versus_higher_rank_indices = this-> MakeDuplicateRanksHigherRankIndicesMap(ranks);
    return duplicate_value_indices_versus_higher_rank_indices;
}

void Atom::RecursivelyCompareBranches(std::map<std::vector<int>, std::vector<int> >& duplicate_higher_indices_map, std::vector<int>& ranks,
  MolecularModeling::AtomVector visited_atoms, std::map<Atom*, std::vector<AtomVector> >& comparison_progress_tracker)
{
    AtomVector primary_neighbors = this->GetNode()->GetNodeNeighbors();
    for (std::map<std::vector<int>, std::vector<int> >::iterator mapit = duplicate_higher_indices_map.begin(); mapit !=
      duplicate_higher_indices_map.end(); mapit++){

	//std::cout << "segfault test 1" << std::endl;
        std::vector<int> indices_with_identical_ranking = mapit->first;
        std::vector<int> indices_with_higher_ranking = mapit->second;
        AtomVector primary_neighbors_with_identical_ranking;

	bool branch_empty = false;
        for (AtomVector::iterator atomit = primary_neighbors.begin(); atomit != primary_neighbors.end(); atomit++){
            int index = std::distance(primary_neighbors.begin(), atomit);
            if (std::find(indices_with_identical_ranking.begin(), indices_with_identical_ranking.end(), index) != indices_with_identical_ranking.end()){
                Atom* neighbor = *atomit;
		//std::cout << "Primary neighbors with identical ranking: " << neighbor->GetName() << std::endl;
                primary_neighbors_with_identical_ranking.push_back(neighbor);
		std::vector<AtomVector> branches = comparison_progress_tracker[neighbor];
		for (unsigned int i = 0; i < branches.size(); i++){
		    for (unsigned int j = 0; j < branches[i].size(); j++){
			//std::cout << "Comparison tracker atoms: " << branches[i][j]->GetName() << std::endl;
		    }
		}
		if (branches.empty()){
		    branch_empty = true;
		}
            }
        }
	if (branch_empty){
	    continue;
	}
	//std::cout << "segfault test 2" << std::endl;

	std::map<Atom*, std::vector<std::vector<int> > > comparison_result_tracker;

        //Compare downstream neighbors and determine priority.
        for (AtomVector::iterator atom_it = primary_neighbors_with_identical_ranking.begin(); atom_it != primary_neighbors_with_identical_ranking.end(); atom_it++){
            Atom* primary_neighbor = *atom_it;
            std::vector<AtomVector> branches = comparison_progress_tracker[primary_neighbor];
	    if (branches.empty()){
		std::cout << primary_neighbor->GetResidue()->GetName() << "-" << primary_neighbor->GetName() << "is empty. Prevent previous recursion from visiting it." << std::endl;
		std::cout << "Primary neighbor neighbors size: "<< primary_neighbor->GetNode()->GetNodeNeighbors().size() << std::endl;
	    }

            std::vector<std::vector<int> > branch_info = this-> ObtainBranchInfo (branches, visited_atoms);
            //Sort branches.
            std::vector<std::vector<int> > descent_order_sorted_branch_info = this->SortBranchesInDescendingOrder(branch_info);
            //Add descent sorted branch info to comparison tracker. Descent sorted so that later the bigger brances are compared first.
            comparison_result_tracker[primary_neighbor] = descent_order_sorted_branch_info;



        }//Done comparing one set of neighbors
	//std::cout << "segfault test 3" << std::endl;

	//Append to visisted atoms.
	for (AtomVector::iterator atom_it = primary_neighbors_with_identical_ranking.begin(); atom_it != primary_neighbors_with_identical_ranking.end(); atom_it++){
	    Atom* primary_neighbor = *atom_it;
	    std::vector<AtomVector> branches = comparison_progress_tracker[primary_neighbor];
            for (std::vector<AtomVector>::iterator branch_it = branches.begin(); branch_it != branches.end(); branch_it++){
                AtomVector individual_branch = *branch_it;
		if (!individual_branch.empty()){
                    visited_atoms.insert(visited_atoms.end(), individual_branch.begin(), individual_branch.end());
		}
            }

            //Make new branches for next round of comparison. Branch expansion rule: each atom in the current branch is replaced by its neighbors, which forms a new branch.
            std::vector<AtomVector> new_branches = this-> MakeNextLevelOfBranches(branches,visited_atoms);
	    //if (new_branches.empty()) std::cout << "New branches is empty" << std::endl;
	    for (unsigned int i = 0; i < new_branches.size(); i++){
		for (unsigned int j = 0; j < new_branches[i].size(); j++){
		    //std::cout << "Next level of branch atom: " << new_branches[i][j]->GetName() << std::endl;
		}
	    }
            comparison_progress_tracker[primary_neighbor] = new_branches;
	}

        //Access the comparison result tracker to determine the relative ranking of this set of primary neighbors;
        std::vector<std::pair<Atom*, int> > descent_sorted_primary_neighbors = this->SortPrimaryNeighborBranchesInDescendingOrder(comparison_result_tracker);
        for (std::vector<std::pair<Atom*, int> >::iterator pair_it = descent_sorted_primary_neighbors.begin(); pair_it != descent_sorted_primary_neighbors.end(); pair_it++){
	    //std::cout << "Testing relative rank after sorting: "  << pair_it->first->GetName() << "-" << pair_it->second << std::endl;
	}

        //Adjust overall ranking based on relative ranking
	std::set<int> unique_relative_ranking;
        for (std::vector<std::pair<Atom*, int> >::iterator pair_it = descent_sorted_primary_neighbors.begin(); pair_it != descent_sorted_primary_neighbors.end(); pair_it++){
	    unique_relative_ranking.insert(pair_it->second);
	}

	//Check if certain identically ranked primary neighbors and empty. If some, but not all, are empty, increase their relative rank by 1.
	std::vector<std::pair<AtomVector, AtomVector> > duplicate_neighbors_higher_ranked_neighbors;
	for (std::set<int>::iterator setit = unique_relative_ranking.begin(); setit != unique_relative_ranking.end(); setit++){
	    int unique_rel_rank = *setit;
	    AtomVector duplicate_neighbors = AtomVector();
	    AtomVector higher_ranked_neighbors = AtomVector();

	    for (std::vector<std::pair<Atom*, int> >::iterator pair_it = descent_sorted_primary_neighbors.begin(); pair_it != descent_sorted_primary_neighbors.end(); pair_it++){
		Atom* primary_neighbor = pair_it->first;
		int rel_rank = pair_it->second;
		if (rel_rank == unique_rel_rank){
		    duplicate_neighbors.push_back(primary_neighbor);
		}
		else if (rel_rank > unique_rel_rank){
		    higher_ranked_neighbors.push_back(primary_neighbor);
		}
	    }

	    if (duplicate_neighbors.size() > 1){
		duplicate_neighbors_higher_ranked_neighbors.push_back(std::make_pair(duplicate_neighbors, higher_ranked_neighbors));
	    }
	}

	for (std::vector<std::pair<AtomVector, AtomVector> >::iterator vec_it = duplicate_neighbors_higher_ranked_neighbors.begin(); vec_it !=
	  duplicate_neighbors_higher_ranked_neighbors.end(); vec_it++){
	    AtomVector duplicate_neighbors = vec_it->first;
	    AtomVector higher_ranked_neighbors = vec_it->second;
	    bool all_empty = true;
	    for (AtomVector::iterator atomit = duplicate_neighbors.begin(); atomit != duplicate_neighbors.end(); atomit++){
		Atom* pri_neighbor = *atomit;
		std::vector<AtomVector> new_branches = comparison_progress_tracker[pri_neighbor];
                for (std::vector<AtomVector>::iterator branchit = new_branches.begin(); branchit != new_branches.end(); branchit++){
                    if (!branchit->empty()){
                        all_empty = false;
                    }
                }
	    }

	    for (AtomVector::iterator atomit = duplicate_neighbors.begin(); atomit != duplicate_neighbors.end(); atomit++){
		Atom* pri_neighbor = *atomit;
		bool empty = true;
		std::vector<AtomVector> new_branches = comparison_progress_tracker[pri_neighbor];
		for (std::vector<AtomVector>::iterator branchit = new_branches.begin(); branchit != new_branches.end(); branchit++){
		    if (!branchit->empty()){
			empty = false;
			all_empty = false;
		    }
		}
		if (empty && !all_empty){
		    //std::cout << pri_neighbor->GetName() << "is a dead end. " << std::endl;
		    for (std::vector<std::pair<Atom*, int> >::iterator pair_it = descent_sorted_primary_neighbors.begin(); pair_it != descent_sorted_primary_neighbors.end(); pair_it++){
			if (pair_it->first == pri_neighbor){
			    pair_it->second++; //Adjust relative ranking if empty.
			}
		    }
		}
	    }

	    if (!all_empty){
	        for (AtomVector::iterator atomit = higher_ranked_neighbors.begin(); atomit != higher_ranked_neighbors.end(); atomit++){
		    Atom* pri_neighbor = *atomit;
		    for (std::vector<std::pair<Atom*, int> >::iterator pair_it = descent_sorted_primary_neighbors.begin(); pair_it != descent_sorted_primary_neighbors.end(); pair_it++){
                        if (pair_it->first == pri_neighbor){
                            pair_it->second++; //Adjust higher ranked primary neighbors.
                        }
                    }
	        }
	    }
	}

	//std::cout << "segfault test 4" << std::endl;

	int min_relative_rank = 99999;
        for (std::vector<std::pair<Atom*, int> >::iterator pair_it = descent_sorted_primary_neighbors.begin(); pair_it != descent_sorted_primary_neighbors.end(); pair_it++){
            int relative_rank = pair_it->second;
	    if (relative_rank < min_relative_rank){
		min_relative_rank = relative_rank;
	    }
	}


        for (std::vector<std::pair<Atom*, int> >::iterator pair_it = descent_sorted_primary_neighbors.begin(); pair_it != descent_sorted_primary_neighbors.end(); pair_it++){
            Atom* primary_neighbor = pair_it->first;
	    //std::cout << "Primary neighbor is: " << primary_neighbor->GetName() << std::endl;
            int primary_neighbor_index = std::distance(primary_neighbors.begin(), std::find(primary_neighbors.begin(), primary_neighbors.end(), primary_neighbor));
	    //std::cout << "But neighbor according to index is: " << primary_neighbors[primary_neighbor_index]->GetName() << std::endl;
            int relative_rank = pair_it->second;
	    //std::cout << "relative rank is: " << relative_rank << " ," << "min rank is: " << min_relative_rank << std::endl;
	    //std::cout << "Rank before adjustment is: " << ranks[primary_neighbor_index] << std::endl;
            ranks[primary_neighbor_index] += (relative_rank - min_relative_rank);
	    //std::cout << "Rank after adjustment is: " << ranks[primary_neighbor_index] << std::endl;
	}

        int maximum_relative_rank = descent_sorted_primary_neighbors.back().second;
        for (std::vector<int>::iterator int_it = indices_with_higher_ranking.begin(); int_it != indices_with_higher_ranking.end(); int_it++){
            int higher_rank_index = *int_it;
            ranks[higher_rank_index] += (maximum_relative_rank -1);
        }
	//std::cout << "segfault test 5" << std::endl;

	for (unsigned int k = 0; k < ranks.size() ; k ++){
	    //std::cout << "Corresponding primary neighbor: " << primary_neighbors[k]->GetName() << std::endl;
	    //std::cout << "Checking rank: " << ranks[k] << std::endl;
	}
    }

    //Check if primary_rankings are unique.If not, start new recursive call
    std::map<std::vector<int>, std::vector<int> > child_duplicate_higher_indices_map = this-> MakeDuplicateRanksHigherRankIndicesMap(ranks);
    int non_empty_primary_neighbor_count = 0;
    if (!child_duplicate_higher_indices_map.empty()){
        for (std::map<std::vector<int>, std::vector<int> >::iterator map_it = child_duplicate_higher_indices_map.begin(); map_it != child_duplicate_higher_indices_map.end();
            map_it++){
            //std::cout << "Still non-unique" << map_it->first[0] << "," << map_it->first[1] << std::endl;
	    std::vector<int> indices_with_identical_ranking = map_it->first;
            std::vector<int> indices_with_higher_ranking = map_it->second;
            for (AtomVector::iterator atomit = primary_neighbors.begin(); atomit != primary_neighbors.end(); atomit++){
                int index = std::distance(primary_neighbors.begin(), atomit);
                if (std::find(indices_with_identical_ranking.begin(), indices_with_identical_ranking.end(), index) != indices_with_identical_ranking.end()){
                    Atom* neighbor = *atomit;
                    if (!comparison_progress_tracker[neighbor].empty()){
                        non_empty_primary_neighbor_count++;
                    }    
		}
            }
        }

    }
    if (non_empty_primary_neighbor_count > 1){
        this->RecursivelyCompareBranches(child_duplicate_higher_indices_map, ranks, visited_atoms, comparison_progress_tracker);
    }

}

std::vector<MolecularModeling::AtomVector> Atom::MakeNextLevelOfBranches(std::vector<MolecularModeling::AtomVector>& current_branches, MolecularModeling::AtomVector& visited_atoms)
{
    std::vector<AtomVector> new_branches = std::vector<AtomVector>();
    for (std::vector<AtomVector>::iterator branch_it = current_branches.begin(); branch_it != current_branches.end(); branch_it++){
        AtomVector individual_branch = *branch_it;
        for (AtomVector::iterator atom_it = individual_branch.begin(); atom_it != individual_branch.end(); atom_it++){
            Atom* individual_atom = *atom_it;
            AtomVector individual_new_branch;
	    AtomVector neighbors = individual_atom->GetNode()->GetNodeNeighbors();
	    for (AtomVector::iterator atom_it = neighbors.begin(); atom_it != neighbors.end(); atom_it++){
		Atom* neighbor = *atom_it;
		//std::cout << "Neighbor: " << neighbor->GetName() << std::endl;
		if (std::find(visited_atoms.begin(), visited_atoms.end(), neighbor) == visited_atoms.end()){
		    //std::cout << "This neighbor is unvisited." << std::endl;
		    //std::cout << "Next level neighbor " << neighbor->GetName() << " is the offspring of " << individual_atom->GetName() << std::endl;
		    individual_new_branch.push_back(neighbor);
		}
	    }

	    if (!individual_new_branch.empty()){
                new_branches.push_back(individual_new_branch);
	    }
        }
    }
    //std::cout << "Next lvl of branches size is: " << new_branches.size() << std::endl;
    return new_branches;
}

std::vector<std::vector<int> > Atom::ObtainBranchInfo (std::vector<MolecularModeling::AtomVector>& branches, MolecularModeling::AtomVector& visited_atoms)
{
    std::vector<std::vector<int> > branch_info;
    for (std::vector<AtomVector>::iterator branch_it = branches.begin(); branch_it != branches.end(); branch_it++){
        AtomVector individual_branch = *branch_it;
	//TODO:Right now each neighbor is counted only once, regardless of the bond type. Need to find double/triple bonded neighbors and count
	//that neighbor multiple times.
        std::multimap<int, Atom*> neighbor_atomic_number_map = this->GetAtomicNumbersOfNeighbors(individual_branch, visited_atoms);

        for (std::multimap<int, Atom*>::iterator nanm_it = neighbor_atomic_number_map.begin(); nanm_it != neighbor_atomic_number_map.end(); nanm_it ++){
	    //std::cout << "atomic number map: " << nanm_it->first << "-" << nanm_it->second->GetName() << std::endl;
	}

        std::set<int> unique_atomic_numbers;
        for (std::multimap<int, Atom*>::iterator nanm_it = neighbor_atomic_number_map.begin(); nanm_it != neighbor_atomic_number_map.end(); nanm_it ++){
            unique_atomic_numbers.insert(nanm_it->first);
        }

        std::vector<int> individual_branch_description;
        for (std::set<int>::reverse_iterator uan_rit = unique_atomic_numbers.rbegin(); uan_rit != unique_atomic_numbers.rend(); uan_rit++){
            int atomic_number = *uan_rit;
	    //std::cout << "Atomic number: " << atomic_number << std::endl;
            std::pair <std::multimap<int, Atom*>::iterator, std::multimap<int, Atom*>::iterator > range = neighbor_atomic_number_map.equal_range(atomic_number);
            int count = std::distance(range.first, range.second);
	    //std::cout << "Count: " << count << std::endl;
            individual_branch_description.push_back(atomic_number);
            individual_branch_description.push_back(count);
        }

        branch_info.push_back(individual_branch_description);
    }////////
    return branch_info;
}

std::map<std::vector<int>, std::vector<int> > Atom::MakeDuplicateRanksHigherRankIndicesMap(std::vector<int> &ranks)
{
    int min_rank = *std::min_element(ranks.begin(), ranks.end());
    int max_rank = *std::max_element(ranks.begin(), ranks.end());
    std::map<std::vector<int>, std::vector<int> > duplicate_value_indices_versus_higher_rank_indices;

    for (int rank_value = min_rank; rank_value <= max_rank; rank_value++){
        std::vector<int> rank_value_indices;
        std::vector<int> indices_with_higher_rank;
        for (std::vector<int>::iterator it = ranks.begin(); it != ranks.end(); it++){
            int index = std::distance(ranks.begin(), it);
            if (*it == rank_value){
                rank_value_indices.push_back(index);
            }
            else if (*it > rank_value){
                indices_with_higher_rank.push_back(index);
            }
        }
        //Filter out unique values
        if (rank_value_indices.size() > 1){
            duplicate_value_indices_versus_higher_rank_indices[rank_value_indices] = indices_with_higher_rank;
        }
    }/////
    return duplicate_value_indices_versus_higher_rank_indices;
}


std::vector<std::pair<Atom*, int> > Atom::SortPrimaryNeighborBranchesInDescendingOrder(std::map<Atom*, std::vector<std::vector<int> > >& comparison_result_tracker)
{
    //std::cout << "Compare primary neighbors: " << std::endl;
    /*for (std::map<Atom*, std::vector<std::vector<int> > >::iterator m = comparison_result_tracker.begin(); m != comparison_result_tracker.end(); m++ ){
	std::cout << m->first->GetName() << ",";
    }
    std::cout << std::endl;*/

    AtomVector descent_sorted_primary_neighbor;
    std::map<Atom*, std::vector<std::vector<int> > >::iterator begin_it = comparison_result_tracker.begin();
    Atom* first_atom = begin_it->first;
    descent_sorted_primary_neighbor.push_back(first_atom);

    std::map<Atom*, std::vector<std::vector<int> > >::iterator begin_pos = comparison_result_tracker.begin();
    std::advance(begin_pos ,1);
    for (std::map<Atom*, std::vector<std::vector<int> > >::iterator mapit = begin_pos; mapit != comparison_result_tracker.end(); mapit++){
	Atom* primary_neighbor = mapit->first;
	//std::cout << "Primary neighbor: " << primary_neighbor->GetName() << std::endl << std::flush;
	std::vector<std::vector<int> > secondary_branches = mapit->second;
	if (this->CompareTwoSetsOfBranches(secondary_branches, comparison_result_tracker[descent_sorted_primary_neighbor.front()]) != 1){ //incoming set >= existing set
	    //std::cout << "S1 >= S2" << std::endl;
	    //std::cout << "Enter if" << std::endl << std::flush;
	    descent_sorted_primary_neighbor.insert(descent_sorted_primary_neighbor.begin(), primary_neighbor);
	}
	else if (this->CompareTwoSetsOfBranches(secondary_branches, comparison_result_tracker[descent_sorted_primary_neighbor.back()]) != 0){ //incoming set <= existing set
	    //std::cout << "S1 <= S2" << std::endl;
	    //std::cout << "Enter elif" << std::endl << std::flush;
	    descent_sorted_primary_neighbor.push_back(primary_neighbor);
	}
	else {
	    //std::cout << "Enter else" << std::endl << std::flush;
	    for (AtomVector::iterator atom_it = descent_sorted_primary_neighbor.begin(); atom_it != descent_sorted_primary_neighbor.end() -1; atom_it++){
		Atom* existing_atom = *atom_it;
		AtomVector::iterator next_atom_it = atom_it +1;
		Atom* next_existing_atom = *next_atom_it;
		if (this->CompareTwoSetsOfBranches(secondary_branches, comparison_result_tracker[existing_atom]) != 0 &&
		  this->CompareTwoSetsOfBranches(secondary_branches,comparison_result_tracker[next_existing_atom]) != 1){
		    descent_sorted_primary_neighbor.insert(next_atom_it, primary_neighbor);
		}
	    }
	}
    }
    std::vector<std::pair<Atom*, int> > descent_sorted_descent_sorted_primary_neighbors_and_ranking;
    int relative_rank = 1;
    descent_sorted_descent_sorted_primary_neighbors_and_ranking.push_back (std::make_pair(descent_sorted_primary_neighbor.front(),relative_rank));

    for (AtomVector::iterator atom_it = descent_sorted_primary_neighbor.begin(); atom_it != descent_sorted_primary_neighbor.end() -1; atom_it++){
	Atom* existing_atom = *atom_it;
	AtomVector::iterator next_atom_it = atom_it +1;
	Atom* next_existing_atom = *next_atom_it;
	if (this->CompareTwoSetsOfBranches(comparison_result_tracker[next_existing_atom], comparison_result_tracker[existing_atom]) == 1){
	    relative_rank++;
	}
	descent_sorted_descent_sorted_primary_neighbors_and_ranking.push_back(std::make_pair(next_existing_atom, relative_rank));
    }
    return descent_sorted_descent_sorted_primary_neighbors_and_ranking;
}

int Atom::CompareTwoSetsOfBranches(std::vector<std::vector<int> > set1, std::vector<std::vector<int> > set2)
{
    //std::cout << "Compare set1: " << std::endl;
    for (unsigned int a = 0; a < set1.size(); a++){
	for (unsigned int a1 = 0 ; a1 < set1[a].size(); a1++){
	    //std::cout << set1[a][a1] << ",";
	}
	//std::cout << std::endl;
    }

    //std::cout << "Compare set2: " << std::endl;
    for (unsigned int a = 0; a < set2.size(); a++){
	for (unsigned int a1 = 0 ; a1 < set2[a].size(); a1++){
	    //std::cout << set2[a][a1] << ",";
	}
	//std::cout << std::endl;
    }

    unsigned int smaller_size = std::min(set1.size(), set2.size());
    for (unsigned int i = 0; i < smaller_size; i++){
	if (this->CompareBranches(set1[i], set2[i]) == 0){
            return 0; //se1 greater than set2
        }
        else if (this->CompareBranches(set1[i], set2[i]) == 1){
            return 1; //v1 smaller than v2
        }
        else{
            if (i == smaller_size - 1){
                if (set1.size() > set2.size()){
                    return 0;
                }
                else if (set1.size() < set2.size()){
                    return 1;
                }
                else {
		    break; //Could've returned 2 here.But this will cause compilation warning that control reaches non-void function.So break out this for loop and return at the end.
                }
            }
        }
    }
    return 2; //Set and set2 are equal
}

std::vector<std::vector<int> > Atom::SortBranchesInDescendingOrder(std::vector<std::vector<int> >& originally_ordered_branches)
{
    std::vector<std::vector<int> > descent_order_sorted_branch_info;
    std::vector<int> first_branch = originally_ordered_branches[0];
    descent_order_sorted_branch_info.push_back(first_branch);

    for (std::vector<std::vector<int> >::iterator bi_it = originally_ordered_branches.begin() +1; bi_it != originally_ordered_branches.end(); bi_it++){
        std::vector<int> branch_info = *bi_it;
        if (this->CompareBranches(branch_info, descent_order_sorted_branch_info.front()) != 1){
            descent_order_sorted_branch_info.insert(descent_order_sorted_branch_info.begin(), branch_info);
        }
        else if (this->CompareBranches(branch_info, descent_order_sorted_branch_info.back()) != 0){
            descent_order_sorted_branch_info.push_back(branch_info);
        }
        else{
            for (std::vector<std::vector<int> >::iterator dosbi_it = descent_order_sorted_branch_info.begin(); dosbi_it != descent_order_sorted_branch_info.end() -1;
              dosbi_it++){
                std::vector<std::vector<int> >::iterator next_it = dosbi_it +1;
                if (this->CompareBranches(branch_info, *dosbi_it) != 0 && this->CompareBranches(branch_info, *next_it) != 1){
                    descent_order_sorted_branch_info.insert(next_it, branch_info);
                }
            }
        }
    }
    return descent_order_sorted_branch_info;
}

int Atom::CompareBranches(std::vector<int>& vec1, std::vector<int>& vec2){
    unsigned int smaller_size = std::min(vec1.size(), vec2.size());
    for (unsigned int i = 0; i < smaller_size; i++){
        if (vec1[i] > vec2[i]){
            return 0; //v1 greater than v2
        }
        else if (vec1[i] < vec2[i]){
            return 1; //v1 smaller than v2
        }
        else{
            if (i == smaller_size - 1){
                if (vec1.size() > vec2.size()){
                    return 0;
                }
                else if (vec1.size() < vec2.size()){
                    return 1;
                }
                else {
                    break; //Could've returned 2 here,which cause compilation warning "control reaches end of non-void function". So break out of for loop and return at the end.
                }
            }
        }
    }
    return 2; //v1 and v2 are equal;
}

std::multimap<int, Atom*> Atom::GetAtomicNumbersOfNeighbors(AtomVector& neighbors, AtomVector& visited_atoms)
{
    std::map<std::string, int> atomic_number_lookup;
    atomic_number_lookup["H"] = 1;
    atomic_number_lookup["He"] = 2;
    atomic_number_lookup["Li"] = 3;
    atomic_number_lookup["Be"] = 4;
    atomic_number_lookup["B"] = 5;
    atomic_number_lookup["C"] = 6;
    atomic_number_lookup["N"] = 7;
    atomic_number_lookup["O"] = 8;
    atomic_number_lookup["F"] = 9;
    atomic_number_lookup["Ne"] = 10;
    atomic_number_lookup["Na"] = 11;
    atomic_number_lookup["Mg"] = 12;
    atomic_number_lookup["Al"] = 13;
    atomic_number_lookup["Si"] = 14;
    atomic_number_lookup["P"] = 15;
    atomic_number_lookup["S"] = 16;
    atomic_number_lookup["Cl"] = 17;
    atomic_number_lookup["Ar"] = 18;
    std::multimap<int, Atom*> neighbor_atomic_number_map;
    //AtomVector center_atom_neighbors = center_atom->GetNode()->GetNodeNeighbors();
    for (AtomVector::iterator atom_it = neighbors.begin(); atom_it != neighbors.end(); atom_it++){
        Atom* atom = *atom_it;
	if (std::find(visited_atoms.begin(), visited_atoms.end(), atom) == visited_atoms.end()){
            std::string element_symbol = atom->GetElementSymbol();
            if (element_symbol.empty()){
                element_symbol = atom->GetName().at(0);
            }

            int atomic_number = atomic_number_lookup[element_symbol];
            neighbor_atomic_number_map.insert(std::make_pair(atomic_number, atom));
	}
    }
    return neighbor_atomic_number_map;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void Atom::Print(std::ostream& out)
{
	out << "Atom name: "    << this->GetName() << std::endl;
	out << "Element: "      << this->GetElementSymbol() << std::endl;
	out << "Atom Type: "    << this->GetAtomType() << std::endl;
	out << "Coordinates: "  << std::endl;
    GeometryTopology::CoordinateVector coordinates = this->GetCoordinates();
    for(GeometryTopology::CoordinateVector::iterator it = coordinates.begin(); it != coordinates.end(); it++)
	{
		GeometryTopology::Coordinate* coordinate = (*it);
		if(coordinate != NULL)
		{
			out << "\t";
			coordinate->Print(out);
			out << std::endl;
		}
	}
	out << "**************** Structure *****************" << std::endl;
	if(this->GetNode() != NULL)
	{
		this->GetNode()->Print(out);
	}
	out << std::endl;
} // end Print

//////////////////////////////////////////////////////////
//                   OVERLOADED OPERATORS               //
//////////////////////////////////////////////////////////
bool Atom::operator== (const Atom &otherAtom)
{
    return (this->GetIndex() == otherAtom.GetIndex());
}

bool Atom::operator!= (const Atom &otherAtom)
{
    return (this->GetIndex() != otherAtom.GetIndex());
}
//////////////////////////////////////////////////////////
//                   HELPER FUNCTIONS                   //
//////////////////////////////////////////////////////////
void Atom::Copy(const Atom* atom)
{
    // Copy the easy stuff.
    this->SetName(atom->GetName());
    this->SetChemicalType(atom->GetChemicalType());
    this->SetDescription(atom->GetDescription());
    this->SetElementSymbol(atom->GetElementSymbol());
    this->SetId(atom->GetId());
    this->SetIsRing(atom->GetIsRing());
    this->SetAtomType(atom->GetAtomType());
    // Deep Copy objects
    // this->residue_ = new MolecularModeling::Residue(atom->GetResidue());
    this->SetResidue(atom->GetResidue());
	if(this->nodes_.at(0) != NULL)
	{
	    delete this->nodes_.at(0);
	}
    this->nodes_.emplace_back (new MolecularModeling::AtomNode(atom->GetNode()));
    GeometryTopology::CoordinateVector atomCoordinates = atom->GetCoordinates();
    for(GeometryTopology::CoordinateVector::iterator it = atomCoordinates.begin(); it != atomCoordinates.end(); it++ )
    {
        GeometryTopology::Coordinate* tempCoordinate = (*it);
        this->coordinates_.push_back(new GeometryTopology::Coordinate(tempCoordinate));
    }
} // end Copy

void Atom::SetAttributes(	MolecularModeling::Residue* residue, std::string name, GeometryTopology::CoordinateVector coordinates,
							std::string chemical_type, std::string description, std::string element_symbol,
							std::vector<MolecularModeling::AtomNode*> nodes, std::string id, bool is_ring, std::string atom_type)
{
	// Having this function call the Setter functions for everything allows for
	//	simple error handling and debugging because the setting of the variables
	//	only happen in one place. DT
	this->SetResidue(residue);
	this->SetName(name);
	this->SetCoordinates(coordinates);
	this->SetChemicalType(chemical_type);
	this->SetDescription(description);
	this->SetElementSymbol(element_symbol);
	this->SetNodes(nodes);
	this->SetId(id);
	this->SetIsRing(is_ring);
	// This function doesn't set index because of the attributes uniqueness, it should
	//	only be called from the Constructors. DT
	this->SetAtomType(atom_type);
} // end SetAttributes
// @TODO DT - The below operators are left for future development. See atom.hpp for
//			details on why they could be useful.
// void Atom::operator=(const Atom&)
// {
//
// } // end Atom::operator=(&)

// void Atom::operator=(const Atom*)
// {
//
// } // end operator=(*)

// std::ostream& operator<<(std::ostream& out, const Atom& atom)
// {
// 	out << "Atom name: "    << atom.GetName() << std::endl;
// 	out << "Element: "      << atom.GetElementSymbol() << std::endl;
// 	out << "Atom Type: "    << atom.MolecularDynamicAtom::GetAtomType()	<< std::endl;
// 	out << "Coordinates: "  << std::endl;
// 	for( CoordinateVector::iterator it = atom.GetCoordinates().begin(); it != atom.GetCoordinates().end(); it++ )
// 	{
// 		GeometryTopology::Coordinate * coordinate = (*it);
// 		out << "\t";
// 		coordinate->Print(out);
// 		out << std::endl;
// 	}
// 	if(atom.GetNode() != NULL)
// 	{
// 		out << "**************** Structure *****************" << std::endl;
// 		atom.GetNode()->Print(out);
// 	}
// 	return out;
// } // end opeartor<<
