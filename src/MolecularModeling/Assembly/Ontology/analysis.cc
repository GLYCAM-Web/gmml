#include <math.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <queue>
#include <stack>

#include "../../../../includes/MolecularModeling/assembly.hpp"
#include "../../../../includes/MolecularModeling/residue.hpp"
#include "../../../../includes/MolecularModeling/atom.hpp"
#include "../../../../includes/MolecularModeling/atomnode.hpp"
#include "../../../../includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp"
#include "../../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp"
#include "../../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceamberprepresidue.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyfile.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyassembly.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyresidue.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyatom.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologybond.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologybondtype.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyangle.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyangletype.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologydihedral.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyatompair.hpp"
#include "../../../../includes/InputSet/CoordinateFileSpace/coordinatefile.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbfile.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbtitlecard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbmodelcard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbmodel.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbheterogenatomcard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbatom.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbconnectcard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdblinkcard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdblink.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdblinkresidue.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp"
#include "../../../../includes/ParameterSet/LibraryFileSpace/libraryfile.hpp"
#include "../../../../includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp"
#include "../../../../includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfile.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../../../includes/utils.hpp"
#include "../../../../includes/common.hpp"
#include "../../../../includes/GeometryTopology/grid.hpp"
#include "../../../../includes/GeometryTopology/cell.hpp"

#include <unistd.h>
#include <errno.h>
#include <string.h>

using namespace std;
using namespace MolecularModeling;
using namespace TopologyFileSpace;
using namespace CoordinateFileSpace;
using namespace PrepFileSpace;
using namespace PdbFileSpace;
using namespace PdbqtFileSpace;
using namespace ParameterFileSpace;
using namespace GeometryTopology;
using namespace LibraryFileSpace;
using namespace gmml;
using namespace Glycan;
using namespace CondensedSequenceSpace;

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
vector<double> Assembly::CalculateBondlengthsStatisticsBasedOnOntologyInfo(string atom_name1, string atom_name2, bool is_atom2_ring, string mono_name)
{
    stringstream query;
    query << "sparql PREFIX : <http://gmmo.uga.edu/#> " <<
             "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> " <<
             "PREFIX owl: <http://www.w3.org/2002/07/owl#> " <<
             "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> " <<
             "PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>" <<
             "SELECT ?atom1_crd ?atom2_crd WHERE {" ;

    query <<  "?mono          :hasSugarName   ?sn.";
    query <<  "?sn            :monosaccharideShortName   \"" << mono_name << "\".\n";
    query <<  "?mono          :hasRingAtom   ?atom1.";
    if(is_atom2_ring)
    {
        query <<  "?mono          :hasRingAtom   ?atom2.";
        query <<  "?atom1         :hasNeighbor    ?atom2.";
    }
    else
        query <<  "?atom1          :hasSideAtom   ?atom2.";
    query <<  "?atom1         :identifier    ?atom1_id.";
    query <<  "?atom2         :identifier    ?atom2_id.";
    query << "FILTER regex(?atom1_id, \"" << atom_name1 << "_" << "\", \"i\")";
    query << "FILTER regex(?atom2_id, \"" << atom_name2 << "_" << "\", \"i\")";

    query <<  "?atom1         :coordinate    ?atom1_crd.";
    query <<  "?atom2         :coordinate    ?atom2_crd.";
    query << "};";

    std::ofstream sparql;
    sparql.open("bonds.sparql", fstream::app);
    sparql << query.str();
    sparql.close();

    system("/home/delaram/virtuoso-7.2.4/bin/isql 1111 dba dba \< bonds.sparql \>  bond_results.txt");
    remove("bonds.sparql");

    ///Read query result file
    string line;
    ifstream in("bond_results.txt");
    while (getline (in, line))///skip the first lines until the coordinates
    {
        if(line.find("____") != string::npos)
        {
            getline (in, line);
            break;
        }
    }
    Coordinate* atom1_crd = new Coordinate();
    Coordinate* atom2_crd = new Coordinate();
    double distance = 0.0;
    double sum_of_bond_lengths = 0.0;
    double sum_of_bond_lengths_squared = 0.0;
    int number_of_bond_lengths = 0;

    ///Reading the coordinates from the result file, calculating the distance
    while (getline (in, line) && !line.empty())
    {
        vector<string> line_tokens = Split(line," ");
        atom1_crd->SetX(ConvertString<double>(Split(line_tokens.at(0), ",").at(0)));
        atom1_crd->SetY(ConvertString<double>(Split(line_tokens.at(1), ",").at(0)));
        atom1_crd->SetZ(ConvertString<double>(line_tokens.at(2)));

        atom2_crd->SetX(ConvertString<double>(Split(line_tokens.at(3), ",").at(0)));
        atom2_crd->SetY(ConvertString<double>(Split(line_tokens.at(4), ",").at(0)));
        atom2_crd->SetZ(ConvertString<double>(line_tokens.at(5)));

        distance = atom1_crd->Distance(*(atom2_crd));
        sum_of_bond_lengths += distance;
        sum_of_bond_lengths_squared += (distance*distance);
        number_of_bond_lengths++;
    }
    in.close();
    remove("bond_results.txt");

    double mean  = sum_of_bond_lengths/number_of_bond_lengths;
    double standard_deviation = sqrt(((sum_of_bond_lengths_squared - ((sum_of_bond_lengths*sum_of_bond_lengths)/number_of_bond_lengths))/number_of_bond_lengths));

    vector<double> statistics = vector<double>();
    statistics.push_back(mean);
    statistics.push_back(standard_deviation);

    return statistics;
}

vector<double> Assembly::CalculateBondAnglesStatisticsBasedOnOntologyInfo(string atom_name1, string atom_name2, string atom_name3, bool is_atom3_ring, string mono_name)
{
    stringstream query;
    query << "sparql PREFIX : <http://gmmo.uga.edu/#> " <<
             "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> " <<
             "PREFIX owl: <http://www.w3.org/2002/07/owl#> " <<
             "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> " <<
             "PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>" <<
             "SELECT ?atom1_crd ?atom2_crd ?atom3_crd WHERE {" ;

    query <<  "?mono          :hasSugarName   ?sn.";
    query <<  "?sn            :monosaccharideShortName   \"" << mono_name << "\".\n";
    query <<  "?mono          :hasRingAtom   ?atom1.";
    query <<  "?mono          :hasRingAtom   ?atom2.";
    query <<  "?atom1         :hasNeighbor    ?atom2.";
    if(is_atom3_ring)
    {
        query <<  "?mono          :hasRingAtom   ?atom3.";
        query <<  "?atom2         :hasNeighbor    ?atom3.";
        query <<  "FILTER (?atom1 != ?atom3)";
    }
    else
        query <<  "?atom2          :hasSideAtom   ?atom3.";
    query <<  "?atom1         :identifier    ?atom1_id.";
    query <<  "?atom2         :identifier    ?atom2_id.";
    query <<  "?atom3         :identifier    ?atom3_id.";
    query << "FILTER regex(?atom1_id, \"" << atom_name1 << "_" << "\", \"i\")";
    query << "FILTER regex(?atom2_id, \"" << atom_name2 << "_" << "\", \"i\")";
    query << "FILTER regex(?atom3_id, \"" << atom_name3 << "_" << "\", \"i\")";

    query <<  "?atom1         :coordinate    ?atom1_crd.";
    query <<  "?atom2         :coordinate    ?atom2_crd.";
    query <<  "?atom3         :coordinate    ?atom3_crd.";
    query << "};";

    std::ofstream sparql;
    sparql.open("bond_angles.sparql", fstream::app);
    sparql << query.str();
    sparql.close();

    system("/home/delaram/virtuoso-7.2.4/bin/isql 1111 dba dba \< bond_angles.sparql \>  bond_angle_results.txt");
    remove("bond_angles.sparql");

    ///Read query result file
    string line;
    ifstream in("bond_angle_results.txt");
    while (getline (in, line))///skip the first lines until the coordinates
    {
        if(line.find("____") != string::npos)
        {
            getline (in, line);
            break;
        }
    }
    Coordinate* atom1_crd = new Coordinate();
    Coordinate* atom2_crd = new Coordinate();
    Coordinate* atom3_crd = new Coordinate();
    double bond_angle = 0.0;
    double sum_of_bond_angles = 0.0;
    double sum_of_bond_angles_squared = 0.0;
    int number_of_bond_angles = 0;

    ///Reading the coordinates from the result file, calculating the distance
    while (getline (in, line) && !line.empty())
    {
        vector<string> line_tokens = Split(line," ");
        atom1_crd->SetX(ConvertString<double>(Split(line_tokens.at(0), ",").at(0)));
        atom1_crd->SetY(ConvertString<double>(Split(line_tokens.at(1), ",").at(0)));
        atom1_crd->SetZ(ConvertString<double>(line_tokens.at(2)));

        atom2_crd->SetX(ConvertString<double>(Split(line_tokens.at(3), ",").at(0)));
        atom2_crd->SetY(ConvertString<double>(Split(line_tokens.at(4), ",").at(0)));
        atom2_crd->SetZ(ConvertString<double>(line_tokens.at(5)));

        atom3_crd->SetX(ConvertString<double>(Split(line_tokens.at(6), ",").at(0)));
        atom3_crd->SetY(ConvertString<double>(Split(line_tokens.at(7), ",").at(0)));
        atom3_crd->SetZ(ConvertString<double>(line_tokens.at(8)));

        bond_angle = CalculateBondAngleByCoordinates(atom1_crd, atom2_crd, atom3_crd);
        sum_of_bond_angles += bond_angle;
        sum_of_bond_angles_squared += (bond_angle*bond_angle);
        number_of_bond_angles++;
    }

    in.close();
    remove("bond_angle_results.txt");

    double mean  = sum_of_bond_angles/number_of_bond_angles;
    double standard_deviation = sqrt(((sum_of_bond_angles_squared - ((sum_of_bond_angles*sum_of_bond_angles)/number_of_bond_angles))/number_of_bond_angles));

    vector<double> statistics = vector<double>();
    statistics.push_back(mean);
    statistics.push_back(standard_deviation);

    return statistics;
}

void Assembly::ExtractTorsionAnglesFromPDB(vector<string> amino_lib_files, string disaccharide)
{
    string pdb_file_path = this->GetSourceFile();
    string pdb = pdb_file_path.substr(pdb_file_path.find_last_of("/") + 1, pdb_file_path.size());
    OligosaccharideVector oligos = this->ExtractSugars(amino_lib_files);
    double phi_angle = 0.0;
    double psi_angle = 0.0;
    int link_index = disaccharide.find_first_of("-");
    char mono1_carbon_index = disaccharide.at(link_index - 1);
    char mono2_carbon_index = disaccharide.at(link_index + 1);
    string first_mono = disaccharide.substr(0, link_index - 1);
    string second_mono = disaccharide.substr(link_index + 2, disaccharide.size());

    std::ofstream out_file;
    out_file.open("torsions.txt", fstream::app);

    ifstream in("torsions.txt");///Checking if the file is empty
    size_t out_file_size = 0;
    in.seekg(0,ios_base::end);
    out_file_size = in.tellg();
    in.close();

    if(out_file_size == 0)
        out_file << left << setw(15) << "PDB" << setw(15) << "Phi Angle" << setw(15) << "Psi Angle" << setw(25) << "Pattern" << "Oligosaccharide" << endl;

    for(OligosaccharideVector::iterator it = oligos.begin(); it != oligos.end(); it++)
    {
        Oligosaccharide* oligo = (*it);
        string oligo_name = oligo->oligosaccharide_name_;
        ///The root of the sequence, last monosaccharide in the oligo name. sequence created backward from root to first mono
        if(oligo_name.compare("") != 0 && oligo_name.find(disaccharide) != string::npos)///found root and the oligo contains the disaccharide.
        {
            queue<Oligosaccharide*> oligo_queue;
            oligo_queue.push(oligo);
            if(MatchDisaccharide(oligo_queue, phi_angle, psi_angle, first_mono, mono1_carbon_index, second_mono, mono2_carbon_index))
                out_file << left << setw(15) << pdb << setw(15) << ConvertRadian2Degree(phi_angle) << setw(15) << ConvertRadian2Degree(psi_angle) << setw(25) << disaccharide << oligo_name << endl;
        }
    }
    out_file.close();
}

void Assembly::CalculateTorsionStatistics(string torsion_file, int low_range, int high_range)
{
    /* Sample output of the function
   |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  1 |  0 |  0 |  0 |  0 |  0 |  0 |  0 | 180
   |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |
   |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  1 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |
   |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  3 |  0 |  1 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |
   |  0 |  0 |  0 |  0 |  0 |  0 |  1 |  0 |  0 |  1 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |
   |  0 |  0 |  0 |  0 |  0 |  1 |  1 |  0 |  1 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |
   |  0 |  0 |  0 |  0 |  0 |  1 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |
   |  0 |  0 |  0 |  0 |  1 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |
Psi|  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |
   |  0 |  0 |  0 |  1 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |
   |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |
   |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |
   |  2 |  0 |  0 |  0 |  0 |  1 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |
   |  0 |  0 |  0 |  0 |  0 |  1 |  1 |  0 |  0 |  5 |  1 |  0 |  0 |  1 |  0 |  1 |  0 |  0 |
   |  0 |  0 |  1 |  0 | 15 | 24 | 15 |  0 |  0 |  0 |  0 |  2 |  3 |  5 |  0 |  0 |  1 |  0 |
   |  0 |  0 |  0 |  1 | 18 | 41 | 43 |  0 |  0 |  0 |  1 | 14 |  6 |  0 |  0 |  0 |  1 |  0 |
   |  0 |  0 |  1 |  0 |  1 | 23 |  9 |  0 |  1 |  2 |  4 | 31 |  1 |  1 |  0 |  0 |  0 |  0 |
   |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  0 |  2 |  0 |  5 |  1 |  0 |  0 |  0 |  0 |  0 |  0 | -180
-180                                         Phi                                          180
     */
    vector<vector<int> > psi_phi_matrix(18, vector<int> (18, 0)); ///matrix 18x18
    vector<vector<int> > psi_omega_matrix(18, vector<int> (18, 0)); ///matrix 18x18
    vector<vector<int> > phi_omega_matrix(18, vector<int> (18, 0)); ///matrix 18x18
    int scale = (high_range - low_range) / 18; ///range of each cell in the matrix
    int x = 0;
    int y = 0;
    bool has_omega = false;

    string line;
    string input_file = "";
    if(torsion_file.compare("") == 0)
        input_file = "torsions.txt";
    else
        input_file = torsion_file;

    ifstream in(input_file.c_str());
    if (!in.is_open())
        cout << "Error in reading the torsion results file" << endl;

    ///Read the head of the file untill the data about torsions
    while (getline (in, line) )
    {
        if(line.find("PDB") != string::npos)
            break;
    }

    while (getline (in, line) && !line.empty())
    {
        x = 0;
        y = 0;
        vector<string> line_tokens = Split(line, " ");
        x = (abs(low_range) + ConvertString<double>(line_tokens.at(1))) / scale; ///x axis in psi_phi_matrix
        y = (abs(low_range) + ConvertString<double>(line_tokens.at(2))) / scale; ///y axis in psi_phi_matrix
        psi_phi_matrix.at(y).at(x) += 1;

        if(line_tokens.size() > 3 && line_tokens.at(3).compare("") != 0) ///x axis in psi_omega_matrix
        {
            has_omega = true;

            x = (abs(low_range) + ConvertString<double>(line_tokens.at(3))) / scale; ///x axis in psi_omega_matrix and phi_omega_matrix
            psi_omega_matrix.at(y).at(x) += 1;

            y = (abs(low_range) + ConvertString<double>(line_tokens.at(1))) / scale; ///y axis in phi_omega_matrix
            phi_omega_matrix.at(y).at(x) += 1;
        }
    }
    in.close();

    for(int j = 17; j >= 0; j-- )
    {
        if(j == 9)
            cout << "Psi| ";
        else
            cout << setw(5) << right << "| ";
        for(int k = 0; k <= 17; k++ )
        {
            cout << setw(2) << psi_phi_matrix.at(j).at(k) << " | ";
        }
        if(j == 17)
            cout << high_range;
        if(j == 0)
            cout << low_range;
        cout << endl;
    }
    cout << setw(5) << left << low_range << setw(43) << right << "Phi" << setw(45) << right << high_range << endl << endl;

    if(has_omega)
    {
        for(int j = 17; j >= 0; j-- )
        {
            if(j == 9)
                cout << "Psi| ";
            else
                cout << setw(5) << right << "| ";
            for(int k = 0; k <= 17; k++ )
            {
                cout << setw(2) << psi_omega_matrix.at(j).at(k) << " | ";
            }
            if(j == 17)
                cout << high_range;
            if(j == 0)
                cout << low_range;
            cout << endl;
        }
        cout << setw(5) << left << low_range << setw(43) << right << "Omega" << setw(45) << right << high_range << endl << endl;

        for(int j = 17; j >= 0; j-- )
        {
            if(j == 9)
                cout << "Phi| ";
            else
                cout << setw(5) << right << "| ";
            for(int k = 0; k <= 17; k++ )
            {
                cout << setw(2) << phi_omega_matrix.at(j).at(k) << " | ";
            }
            if(j == 17)
                cout << high_range;
            if(j == 0)
                cout << low_range;
            cout << endl;
        }
        cout << setw(5) << left << low_range << setw(43) << right << "Omega" << setw(45) << right << high_range << endl;
    }
}

void Assembly::CalculateGlyprobityGeometryOutliers(Monosaccharide* mono)
{
    ///EXTRACTING BOND LENGTHS OF MONOSACCHARIDE ATOM PAIRS AND BOND ANGLES
    vector<string> visited_bonds = vector<string>();
    vector<string> visited_angles = vector<string>();
    vector<double> bond_statistics = vector<double>();
    vector<double> bond_angle_statistics = vector<double>();
    stringstream bond_lengths_stream;
    stringstream bond_angles_stream;
    double bond_angle = 0.0;
    bool is_atom2_ring = false;
    bool is_atom3_ring = false;

    bond_lengths_stream << "GLYPROBITY REPORT" << endl <<
                           "<--Geometry Outliers-->" << endl <<
                           "Bond lengths (current bond length, ontology mean, ontology standard deviation)" << endl;

    bond_angles_stream << "Bond angles (current bond angle, ontology mean, ontology standard deviation)" << endl;

    for(AtomVector::iterator ring_atom_it = mono->cycle_atoms_.begin(); ring_atom_it != mono->cycle_atoms_.end(); ring_atom_it++)
    {
        Atom* atom1 = (*ring_atom_it);
        AtomVector atom1_neighbors = atom1->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator atom1_neighbors_it = atom1_neighbors.begin(); atom1_neighbors_it != atom1_neighbors.end(); atom1_neighbors_it++)
        {
            Atom* atom2 = (*atom1_neighbors_it);
            stringstream check_bond;
            stringstream check_bond_reverse;

            check_bond << atom1->GetId() << "-" << atom2->GetId();
            check_bond_reverse << atom2->GetId() << "-" << atom1->GetId();
            if(find(visited_bonds.begin(), visited_bonds.end(), check_bond.str()) == visited_bonds.end() &&
                    find(visited_bonds.begin(), visited_bonds.end(), check_bond_reverse.str()) == visited_bonds.end())///if the bond has not been visited before
            {
                visited_bonds.push_back(check_bond.str());
                visited_bonds.push_back(check_bond_reverse.str());
                bond_lengths_stream << atom1->GetName() << "-" << atom2->GetName() << ": " <<
                                       atom1->GetCoordinates().at(model_index_)->Distance(*(atom2->GetCoordinates().at(model_index_)));

                ///if atom 2 is not a ring atom set the flag to true
                if(mono->cycle_atoms_str_.find(atom2->GetId()) == string::npos)
                    is_atom2_ring = false;
                else
                    is_atom2_ring = true;

                ///Find same bonds in the same monosaccharide in the ontology, calculate mean and standard deviation from ontology
                bond_statistics = CalculateBondlengthsStatisticsBasedOnOntologyInfo(atom1->GetName(), atom2->GetName(), is_atom2_ring, mono->sugar_name_.monosaccharide_short_name_);
                bond_lengths_stream << ", " << bond_statistics.at(0) << ", " << bond_statistics.at(1) << endl;
            }

            if(is_atom2_ring)///Calculating bond angles only for ring atoms and combination of ring atoms anf their immediate neighbors
            {
                ///EXTRACTING BOND ANGLES
                AtomVector atom2_neighbors = atom2->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator atom2_neighbors_it = atom2_neighbors.begin(); atom2_neighbors_it != atom2_neighbors.end(); atom2_neighbors_it++)
                {
                    if((*ring_atom_it) != (*atom2_neighbors_it))///if the neighbor of the second atom is not the first atom we chose
                    {
                        Atom* atom3 = (*atom2_neighbors_it);
                        stringstream check_angle;
                        stringstream check_angle_reverse;
                        check_angle << atom1->GetId() << "-" << atom2->GetId() << "-" << atom3->GetId();
                        check_angle_reverse  << atom3->GetId() << "-" << atom2->GetId() << "-" << atom1->GetId();
                        if(find(visited_angles.begin(), visited_angles.end(), check_angle.str()) == visited_angles.end() &&
                                find(visited_angles.begin(), visited_angles.end(), check_angle_reverse.str()) == visited_angles.end())///if the bond angle has not been visited before
                        {
                            visited_angles.push_back(check_angle.str());
                            visited_angles.push_back(check_angle_reverse.str());

                            bond_angle = CalculateBondAngleByAtoms(atom1, atom2, atom3);
                            bond_angles_stream << atom1->GetName() << "-" << atom2->GetName() << "-" << atom3->GetName() << ": " << bond_angle;

                            ///if atom 3 is not a ring atom set the flag to true
                            if(mono->cycle_atoms_str_.find(atom3->GetId()) == string::npos)
                                is_atom3_ring = false;
                            else
                                is_atom3_ring = true;

                            ///Find same bond angles in the same monosaccharide in the ontology, calculate mean and standard deviation from ontology
                            bond_angle_statistics = CalculateBondAnglesStatisticsBasedOnOntologyInfo(atom1->GetName(), atom2->GetName(), atom3->GetName(), is_atom3_ring,
                                                                                                     mono->sugar_name_.monosaccharide_short_name_);
                            bond_angles_stream << ", " << bond_angle_statistics.at(0) << ", " << bond_angle_statistics.at(1) << endl;
                        }
                    }
                }
            }
        }
    }
    cout << endl << bond_lengths_stream.str() << endl << bond_angles_stream.str() << "<--------------------->" << endl << endl;
}
