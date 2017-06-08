#include <iostream>
#include "../includes/gmml.hpp"
#include "cmath"

using namespace std;
using namespace MolecularModeling;
typedef std::vector<Atom*> AtomVector;
typedef std::vector<GeometryTopology::Coordinate*> CoordinateVector;
typedef std::vector<Residue*> ResidueVector;
typedef std::vector<Assembly*> AssemblyVector;
typedef GeometryTopology::Coordinate Vector;
typedef std::vector<GeometryTopology::Coordinate> CoordinateSet;




void GetResidueRingCenter (Residue* residue, GeometryTopology::Coordinate* center);
double GetDistanceToAromaticRingCenter(GeometryTopology::Coordinate* centroid);
//GeometryTopology::Coordinate FindVector(GeometryTopology::Coordinate coord_1, GeometryTopology::Coordinate coord_2);
void FindVector(Atom* atom_1, Atom* atom_2, GeometryTopology::Coordinate*v1);
//void FindConnectedAliphaticC(Atom* h_ptr, AtomVector &list, Atom* c_ptr);
Atom* FindConnectedAliphaticC(Atom* h_ptr);
double FindAngleToNormalVector(GeometryTopology::Coordinate* normal_vector, GeometryTopology::Coordinate v1);
void CreateCHpiEnergyList(double r, double theta, std::vector<double> *cHpiEnergy);
double TotalCHpiEnergy(std::vector<double> energyList);
GeometryTopology::Coordinate FindNormalVector(MolecularModeling::Atom *centroid_ptr, GeometryTopology::Coordinate* coord_1, GeometryTopology::Coordinate* coord_2);

//measure distance of H coord to ring. If < 6 A, note that it is a plausible CH_pi pair so I can do math on it later.

//pseudocode
//store center of rings in a vector
//iterate through center of ring coordinates, and find distances to aliphatic H atoms. H atom coordinates by themselves do not need to be stored, and can use pointers to access their coordinates.
//Once distance is found, find the angle to the normal vector.
//write a function that will input distance and angle to get the CH-pi energies. You do not need to save distance and angles. 
//write a function that will save the energies.

//pseudocode 2
//create yohAssembly, and from there extract all residues of assembly.
//create CHAssembly, which will contain values from pointers to important CH's. This is how I will access information about atoms. Get Center of Ring (done).
//within CHAssembly, iterate through the C atoms and get distances to each H atom. Lachele says there is a function already that can get distances to each atom.
//calculate the angle from a normal vector to the plane of the ring.
//input distance and angle into energy function to get CH-pi energy. Store it in a vector.

int main()
{
    std::vector<std::string> amino_libs, glycam_libs, other_libs, prep;
    amino_libs.push_back("/home/yohanna/gems/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/amino12.lib");
    amino_libs.push_back("/home/yohanna/gems/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/aminoct12.lib");
    amino_libs.push_back("/home/yohanna/gems/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/aminont12.lib");

    glycam_libs.push_back("/home/yohanna/gems/gmml/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_amino_06j_12SB.lib");
    glycam_libs.push_back("/home/yohanna/gems/gmml/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_aminoct_06j_12SB.lib");
    glycam_libs.push_back("/home/yohanna/gems/gmml/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_aminont_06j_12SB.lib");

    other_libs.push_back("/home/yohanna/gems/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/nucleic12.lib");
    other_libs.push_back("/home/yohanna/gems/gmml/dat/CurrentParams/leaprc.ff12SB_2014-04-24/nucleic12.lib");
    other_libs.push_back("/home/yohanna/gems/gmml/dat/CurrentParams/other/solvents.lib");

    prep.push_back("/home/yohanna/gems/gmml/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_06j-1.prep");

    string parameter_file_path = "/home/yohanna/gems/gmml/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_06j.dat";
    string ion_parameter_file_path = "/home/yohanna/gems/gmml/dat/CurrentParams/other/atomic_ions.lib";

    //my example
    Assembly yohAssembly;

    string pdb_file_path = "/home/yohanna/work/3.CH-pi_interactions/test_cases_for_chpi_code/TYR-0SA-test2-aligned_to_3.8.pdb";
    yohAssembly.BuildAssemblyFromPdbFile(pdb_file_path, amino_libs, glycam_libs, other_libs, prep, parameter_file_path);
    yohAssembly.BuildStructureByDistance();
    //puts all residues in yohResidues. Creates vector of pointers to Residue
    ResidueVector yohResidues = yohAssembly.GetAllResiduesOfAssembly();
    ResidueVector aromaticResidues;
    ResidueVector sugarResidues;
    GeometryTopology::Coordinate* centerOfRing = new GeometryTopology::Coordinate();
    Atom* aliphaticC;
//    GeometryTopology::Coordinate *AliphaticH = new GeometryTopology::Coordinate();
    double distanceToCentroid;
//    GeometryTopology::Coordinate* cHVector = new GeometryTopology::Coordinate();
    GeometryTopology::Coordinate cHVector;
    double angle;
    std::vector<double> cHpiEnergies;
    GeometryTopology::Coordinate* coordinateCG_Ptr = NULL;
//    GeometryTopology::Coordinate coordinateCD1;
    GeometryTopology::Coordinate* coordinateCD1_Ptr = NULL;

//    std::vector<pair<double,double>> cHPair_r_theta;

//    Residue aromatic;

    //create aromatic and sugar residue vectors
//    yohAssembly.PrintHetResidues();
    for(ResidueVector::iterator it = yohResidues.begin(); it != yohResidues.end(); ++it)
    {
        Residue* current_residue = *it;
        if(current_residue->GetName().compare("TYR")==0 || current_residue->GetName().compare("PHE")==0 ||
           current_residue->GetName().compare("TRP")==0 || current_residue->GetName().compare("HIS")==0)
        {
          // Want it to look like this after moding Assembly in GMML:
         //centerOfRing = current_residue->GetResidueRingCenter;

//          std::cout << "Residue is " << current_residue->GetId() << std::endl;
//          centerOfRing->Print();
//          std::cout << std::endl;
          //add this Residue to CHAssembly
            aromaticResidues.push_back(current_residue);            
        }
    }
    for(ResidueVector::iterator it = yohResidues.begin(); it != yohResidues.end(); ++it)
    {
        Residue* current_residue = *it;
        if(current_residue->GetName().compare("0SA")==0)
        {
            sugarResidues.push_back(current_residue);
//            std::cout << "Residue is: " << current_residue->GetId() << std::endl;
        }
    }
//    for (ResidueVector::iterator it = aromaticResidues.begin(); it != aromaticResidues.end(); ++it)
//    {
//       (*it)->Print();
//    }
//    for (ResidueVector::iterator it = sugarResidues.begin(); it != sugarResidues.end(); ++it)
//    {
//        (*it)->Print();
//    }


    //get ring center and calculate distance to H
    for (ResidueVector::iterator it = aromaticResidues.begin(); it !=aromaticResidues.end(); ++it)
    {
        Residue* current_residue = *it;
        if (current_residue->GetName().compare("TYR")==0 || current_residue->GetName().compare("PHE")==0 ||
            current_residue->GetName().compare("TRP")==0 || current_residue->GetName().compare("HIS")==0)
        {
            GetResidueRingCenter(current_residue, centerOfRing);
            CoordinateVector centroid_ptr;
            centroid_ptr.push_back(centerOfRing);
            Atom centroidPseudoAtom(current_residue, "Cp", centroid_ptr);
            Atom* centroidPseudoAtom_ptr = &centroidPseudoAtom;
            std::cout << "centroid pseudo atom: " << std::endl;
            centroidPseudoAtom.Print();

            AtomVector aromatic_atoms = current_residue->GetAtoms();
            for(Assembly::AtomVector::iterator it = aromatic_atoms.begin(); it != aromatic_atoms.end(); ++it)
            {
                Atom* current_aromatic_atom = *it;
                //assume that v1 and v2 are with respect to origin
                if (current_aromatic_atom->GetName().compare("CG")==0 || current_aromatic_atom->GetName().compare("CD1")==0)
                {
                    //save coordinates of CG and CD1 for getting plane vectors
                    if (current_aromatic_atom->GetName().compare("CG")==0)
                    {
                        std::cout << "inside CG" << std::endl;
                        GeometryTopology::Coordinate coordinateCG(current_aromatic_atom->GetCoordinates().at(0));
                        coordinateCG.Print();
                        std::cout << "" << std::endl;
                        coordinateCG_Ptr = &coordinateCG;

                    }
                    if (current_aromatic_atom->GetName().compare("CD1")==0)
                    {
//                        std::cout << "inside CD1" << std::endl;
                        GeometryTopology::Coordinate coordinateCD1(current_aromatic_atom->GetCoordinates().at(0));
                        coordinateCD1_Ptr = &coordinateCD1;
                    }
                    if (coordinateCG_Ptr != NULL && coordinateCD1_Ptr != NULL)
                    {
                        std::cout << "CG...not a null pointer: " << std::endl;
                        coordinateCG_Ptr->Print();
                        std::cout << " " << std::endl;
                         //calculate distance from ring center to H and CH Vector
                        for (ResidueVector::iterator it = sugarResidues.begin(); it !=sugarResidues.end(); ++it)
                        {
                            Residue* sugar_residue = *it;
                            if (sugar_residue->GetName().compare("0SA")==0)
                            {
                                AtomVector glycan_atoms = sugar_residue->GetAtoms();
                                for(Assembly::AtomVector::iterator it = glycan_atoms.begin(); it != glycan_atoms.end(); ++it)
                                {
                                    Atom* current_atom = *it;
                                    if (current_atom->GetName().compare("H3A")==0)
                                    {
                                        Atom* aliphaticH = *it;
    //                                    std::cout <<"atom is: " << std::endl;
    //                                    aliphaticH->Print();
    //                                    std::cout << " " << std::endl;
                                        aliphaticC = FindConnectedAliphaticC(aliphaticH);
                                        FindVector(aliphaticH, aliphaticC, &cHVector);
                                        std::cout << "CH Vector" << std::endl;
                                        cHVector.Print();
                                        std::cout << "coordinates of CG pointer are: " << std::endl;
                                        coordinateCD1_Ptr->Print();
                                        GeometryTopology::Coordinate normal_v = FindNormalVector(centroidPseudoAtom_ptr, coordinateCG_Ptr, coordinateCD1_Ptr);
                                        GeometryTopology::Coordinate* normal_v_Ptr = &normal_v;
                                        std::cout << "Normal Vector" << std::endl;
                                        normal_v_Ptr->Print();
                                        distanceToCentroid = aliphaticC->GetDistanceToAtom(centroidPseudoAtom_ptr);
    //                                    std::cout << "distance to Centroid: " << distanceToCentroid  << std::endl;
    //                                    std::cout << "outside function: " <<std::endl;
    //                                    cHVector.Print();
    //                                    std::cout << " " << std::endl;
                                        angle = FindAngleToNormalVector(normal_v_Ptr, cHVector);
    //                                    angle = 1.0;
                                        CreateCHpiEnergyList(distanceToCentroid, angle, &cHpiEnergies);
                                        std::cout << "angle outside function: " << std::endl;
                                        std::cout << angle << std::endl;
                                        std::cout << " " << std::endl;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    double finalEnergy = TotalCHpiEnergy(cHpiEnergies);
    std::cout << "Final Energy is: " << finalEnergy << std::endl;
//    PdbFileSpace::PdbFile* pdbOutput = yohAssembly.BuildPdbFileStructureFromAssembly(-1,0);
//    pdbOutput->Write("/home/yohanna/work/3.CH-pi_interactions/test_cases_for_chpi_code/output_TYR-0SA-test2-aligned_to_3.8.pdb");
    return 0;
}
//            centerOfRing->Print();
        //calculate distance from ring center to H
//        for (ResidueVector::iterator it = sugarResidues.begin(); it !=sugarResidues.end(); ++it)
//        {
//            Residue* sugar_residue = *it;
//            if (sugar_residue->GetName().compare("0SA")==0)
//            {
//                AtomVector glycan_atoms = sugar_residue->GetAtoms();
//                for(Assembly::AtomVector::iterator it = glycan_atoms.begin(); it != glycan_atoms.end(); ++it)
//                {
//                    if (((*it)->GetName()).compare("H3E")==0 || ((*it)->GetName()).compare("H3A")==0)
//                    {
//                    //get coordinate of aliphatic H
//                    // write function similar to GetDistanceToAtom(Atom *otherAtom) but change to GetDistanceToAromaticRingCenter. I can either write my own function, or make center of ring as an atom type.
//                        Atom* aliphaticH = *it;
//                        distanceToCentroid = aliphaticH->GetDistanceToAtom(centroidPseudoAtom_ptr);
//                        std::cout << "distance to Centroid: " << distanceToCentroid  << std::endl;
//                        AtomVector ch_list;
//                        aliphaticH->FindConnectedAtoms(h_list);//I think this is finding the whole sugar connected to it, not just the first? What I think based on Print().

//                        aliphaticC = FindConnectedAliphaticC(aliphaticH);
//                        aliphaticC->Print();
//                        GeometryTopology::Coordinate cHVector = FindVector(aliphaticH, aliphaticC);
//                        std::cout << "CH Vector: " << std::endl;
//                        cHVector.Print();
//                        std::cout << std::endl;
//                        cHVectorSet.push_back(cHVector);
//                        //check to see if adding to cHVectorSet
//                        for(CoordinateSet::iterator it = cHVectorSet.begin(); it !=cHVectorSet.end(); ++it)
//                        {
//                            GeometryTopology::Coordinate current_vector = *it;
//                            current_vector.Print();
//                            double cHAngle = FindAngleToNormalVector(current_vector, )
//                        }

//                        double cHAngle =FindAngleToNormalVector(normal_v,cHVector); //create a getcHAngle() later

//                    }
//                }
//            }
//        }

//        }
//   }
//  return 0;
//}

//FUNCTIONS///////////////////////////////////////////////////////////////////////////////

void GetResidueRingCenter (Residue *residue, GeometryTopology::Coordinate *center)
{
    double sumX = 0.0, sumY = 0.0, sumZ = 0.0;
    int numberOfRingAtoms = 0;
    AtomVector atoms = residue->GetAtoms();

    for(Assembly::AtomVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        Atom* atom = *it;
        if (atom->GetName().compare("CG")==0 || atom->GetName().compare("CD2")==0 ||
              atom->GetName().compare("CE2")==0 || atom->GetName().compare("CZ")==0 ||
              atom->GetName().compare("CE1")==0 || atom->GetName().compare("CD1")==0 ||
              atom->GetName().compare("NE1")==0 || atom->GetName().compare("CZ2")==0 ||
              atom->GetName().compare("CH2")==0 || atom->GetName().compare("CZ3")==0 ||
              atom->GetName().compare("CE3")==0 || atom->GetName().compare("ND1")==0 ||
              atom->GetName().compare("NE2")==0)
            {
            numberOfRingAtoms++;
//            std::cout << "Atom is ring: " << atom->GetName() << endl;
            sumX += atom->GetCoordinates().at(0)->GetX();
            sumY += atom->GetCoordinates().at(0)->GetY();
            sumZ += atom->GetCoordinates().at(0)->GetZ();
            }
    }
    center->SetX( sumX / numberOfRingAtoms  );
    center->SetY( sumY / numberOfRingAtoms  );
    center->SetZ( sumZ / numberOfRingAtoms  );
    return;
}

Atom* FindConnectedAliphaticC(Atom* h_ptr)
{
//    visitedAtoms.push_back(h_ptr);
    AtomVector neighbors = h_ptr->GetNode()->GetNodeNeighbors();
//    bool alreadyVisited = false;

    for(AtomVector::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
    {
        Atom* neighbor = *it;
//        alreadyVisited = false;
        Atom* c_ptr = neighbor;
        std::cout << "c_ptr in FindConnectedAliphaticC" << std::endl;
        c_ptr->Print();
        return c_ptr;
//        for(AtomVector::iterator visitedAtom = visitedAtoms.begin(); visitedAtom != visitedAtoms.end(); visitedAtom++)
//        {
//            if ( (*neighbor)->GetIndex() == (*visitedAtom)->GetIndex() );
//            alreadyVisited = true;
//            unsigned long long test_index = (*visitedAtom)->Get

//        }

    }
}


//    for(AtomVector::iterator neighbor = neighbors.begin(); neighbor != neighbors.end(); neighbor++){
//        alreadyVisited = false; // reset for each neighbor
//        for(AtomVector::iterator visitedAtom = visitedAtoms.begin(); visitedAtom != visitedAtoms.end(); visitedAtom++){
//            if ( (*neighbor)->GetIndex() == (*visitedAtom)->GetIndex() )
//                alreadyVisited = true;


//void PlaneOfAromatic(Atom* center, Atom* CG)

 void FindVector(Atom* atom_1, Atom* atom_2, GeometryTopology::Coordinate* v1)
{
    GeometryTopology::Coordinate coord_1(atom_1->GetCoordinates().at(0));
//    std::cout << "coord_1 is: " << std::endl;
//    coord_1.Print();
    GeometryTopology::Coordinate coord_2(atom_2->GetCoordinates().at(0));
    coord_1.operator -(coord_2);
    GeometryTopology::Coordinate coord_out(coord_1);
//    std::cout << "inside function: " << std::endl;
//    coord_out.Print();
    std::cout << std::endl;
//    v1 = &coord_out;
//    double x = v1->GetX();
//    double y = v1->GetY();
//    double z = v1->GetZ();
    double x = coord_out.GetX();
    double y = coord_out.GetY();
    double z = coord_out.GetZ();
    v1->SetX(x);
    v1->SetY(y);
    v1->SetZ(z);
    return;
}

 double FindAngleToNormalVector(GeometryTopology::Coordinate* normal_vector, GeometryTopology::Coordinate v1)
 {
     GeometryTopology::Coordinate norm_coord(normal_vector);
     double dot = norm_coord.DotProduct(v1);
     double lenSq1 = v1.length();
     double lenSq2 = norm_coord.length();
     double angle = acos(dot/(lenSq1 * lenSq2));
     std::cout << "angle inside function: " << std::endl;
     std::cout << angle << std::endl;
     return angle;
 }

void CreateCHpiEnergyList(double r, double theta, std::vector<double>* cHpiEnergy)
{
    double answer = (-14064.0/pow(r,8))*(1 + cos(2*theta));
    cHpiEnergy->push_back(answer);
    return;
}

double TotalCHpiEnergy(std::vector<double> energyList)
{
    double sum_of_elems = 0.0;
    for(std::vector<double>:: iterator it = energyList.begin(); it != energyList.end(); ++it)
    {
        sum_of_elems += *it;
    }
    return sum_of_elems;
}

GeometryTopology::Coordinate FindNormalVector(MolecularModeling::Atom* centroid_ptr, GeometryTopology::Coordinate* coord_1, GeometryTopology::Coordinate* coord_2)
{
    //something wrong with getting coordinates from pointers. Fix later
    GeometryTopology::Plane plane;
    GeometryTopology::Coordinate centroid_coord(centroid_ptr->GetCoordinates().at(0));
    double x_coord_1 = coord_1->GetX();
    double y_coord_1 = coord_1->GetY();
    double z_coord_1 = coord_1->GetZ();
    std::cout << "x_coord is " << x_coord_1 << std::endl;
    GeometryTopology::Coordinate coord1(x_coord_1, y_coord_1, z_coord_1);
    double x_coord_2 = coord_2->GetX();
    double y_coord_2 = coord_2->GetY();
    double z_coord_2 = coord_2->GetZ();
    GeometryTopology::Coordinate coord2(x_coord_2, y_coord_2, z_coord_2);
//    GeometryTopology::Coordinate reference(reference_coord->GetCoordinates().at(0));
    std::cout << "Print coord_1 before" << std::endl;
    coord1.Print();
    std::cout << " " << std::endl;
    coord1.operator-(centroid_coord);
    coord2.operator-(centroid_coord);
    std::cout << "Print coord_1 after" << std::endl;
    coord1.Print();
    plane.SetV1(coord1);
    plane.SetV2(coord2);
    GeometryTopology::Coordinate normal_v = plane.GetUnitNormalVector();
    std::cout << "normal v. inside function" << std::endl;
    std::cout << " " << std::endl;
    return normal_v;
}


// double Coordinate::DotProduct(Coordinate coordinate)
// {
//     return ((x_ * coordinate.x_) + (y_ * coordinate.y_) + (z_ * coordinate.z_));

//Coordinate prev_atom_coord = Coordinate(*prev_atom->GetCoordinates().at(model_index_));
//Coordinate current_atom_coord = Coordinate(*current_atom->GetCoordinates().at(model_index_));
//Coordinate next_atom_coord = Coordinate(*next_atom->GetCoordinates().at(model_index_));
//prev_atom_coord.operator -(current_atom_coord) ;
//next_atom_coord.operator -(current_atom_coord) ;
//Plane plane = Plane();
//plane.SetV1(prev_atom_coord);
//plane.SetV2(next_atom_coord);
//Coordinate normal_v = plane.GetUnitNormalVector();
