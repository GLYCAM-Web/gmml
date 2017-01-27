#include <math.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <queue>
#include <stack>

#include "../../../includes/MolecularModeling/assembly.hpp"
#include "../../../includes/MolecularModeling/residue.hpp"
#include "../../../includes/MolecularModeling/atom.hpp"
#include "../../../includes/MolecularModeling/atomnode.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceamberprepresidue.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyfile.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyassembly.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyresidue.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyatom.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologybond.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologybondtype.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyangle.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyangletype.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologydihedral.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp"
#include "../../../includes/InputSet/TopologyFileSpace/topologyatompair.hpp"
#include "../../../includes/InputSet/CoordinateFileSpace/coordinatefile.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbfile.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbtitlecard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodelcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodel.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogenatomcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatom.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbconnectcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblinkcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblink.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdblinkresidue.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfile.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfile.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/GeometryTopology/grid.hpp"
#include "../../../includes/GeometryTopology/cell.hpp"

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
Assembly::AtomVector Assembly::Select(string pattern)
{
    AtomVector selection = AtomVector();

    HierarchicalContainmentMap hierarchical_map;
    stringstream index;
    index << this->GetSequenceNumber();
    this->GetHierarchicalMapOfAssembly(hierarchical_map, index);

    for(HierarchicalContainmentMap::iterator it = hierarchical_map.begin(); it != hierarchical_map.end(); it++)
        cout << (*it).first << " " << (*it).second.size() << endl;

    SelectPatternMap select_pattern_map = ParsePatternString(pattern);

    //    for(SelectPatternMap::iterator it = select_pattern_map.begin(); it != select_pattern_map.end(); it++)
    //    {
    //        cout << "================================= " << (*it).first << " ================================" << endl;
    //        map<string, vector<string> > value = (*it).second;
    //        for(map<string, vector<string> >::iterator it1 = value.begin(); it1 != value.end(); it1++)
    //        {
    //            cout << "-------------------------------- " << (*it1).first << " --------------------------" << endl;
    //            vector<string> atom_names = (*it1).second;
    //            for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
    //            {
    //                cout << (*it3) << endl;
    //            }
    //        }
    //    }

    for(SelectPatternMap::iterator it = select_pattern_map.begin(); it != select_pattern_map.end(); it++)
    {
        string key = (*it).first;
        if(key.find("*") == string::npos)
        {
            ResidueVector residues_of_assembly = hierarchical_map.at(key);
            map<string, vector<string> > value = (*it).second;
            for(map<string, vector<string> >::iterator it1 = value.begin(); it1 != value.end(); it1++)
            {
                string residue_name = (*it1).first;
                int residue_name_search_type = 0;
                if(residue_name.at(0) == '^')
                {
                    residue_name = residue_name.substr(1);
                    residue_name_search_type = -1;
                }
                else if(residue_name.at(residue_name.size() - 1) == '$')
                {
                    residue_name = residue_name.substr(0, residue_name.size() - 2);
                    residue_name_search_type = 1;
                }
                else if(residue_name.at(0) == '#')
                {
                    residue_name = residue_name.substr(1);
                    residue_name_search_type = -2;
                }
                else
                    residue_name_search_type = 0;
                if(residue_name.find("*") == string::npos)
                {
                    for(ResidueVector::iterator it2 = residues_of_assembly.begin(); it2 != residues_of_assembly.end(); it2++)
                    {
                        Residue* residue = *it2;
                        switch(residue_name_search_type)
                        {
                            case 0:  /// Search in residue set by matching the whole name of the residue
                            {
                                if(residue->GetName().compare(residue_name) == 0)
                                {
                                    vector<string> atom_names = (*it1).second;
                                    if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                    {
                                        AtomVector atoms = residue->GetAtoms();
                                        for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                        {
                                            string atom_name = *it3;
                                            int atom_name_search_type = 0;
                                            if(atom_name.at(0) == '^')
                                            {
                                                atom_name = atom_name.substr(1);
                                                atom_name_search_type = -1;
                                            }
                                            else if(atom_name.at(atom_name.size() - 1) == '$')
                                            {
                                                atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                atom_name_search_type = 1;
                                            }
                                            else if(atom_name.at(0) == '#')
                                            {
                                                atom_name = atom_name.substr(1);
                                                atom_name_search_type = -2;
                                            }
                                            else
                                                atom_name_search_type = 0;
                                            for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                            {
                                                Atom* atom = *it4;
                                                switch(atom_name_search_type)
                                                {
                                                    case 0:
                                                    {
                                                        if(atom->GetName().compare(atom_name) == 0)
                                                            selection.push_back(atom);
                                                        break;
                                                    }
                                                    case 1:
                                                    {
                                                        if(atom->GetName().find(atom_name) != string::npos &&
                                                                atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                            selection.push_back(atom);
                                                        break;
                                                    }
                                                    case -1:
                                                    {
                                                        if(atom->GetName().find(atom_name) != string::npos &&
                                                                atom->GetName().find(atom_name) == 0)
                                                            selection.push_back(atom);
                                                        break;
                                                    }
                                                    case -2:
                                                    {
                                                        string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                        int range_selection = 0;
                                                        if(atom_name.find("-") != string::npos)
                                                        {
                                                            range_selection = 1;

                                                        }
                                                        else
                                                            range_selection = 0;
                                                        switch(range_selection)
                                                        {
                                                            case 0:
                                                            {
                                                                if(atom_serial_number.find(atom_name) != string::npos &&
                                                                        atom_serial_number.find(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case 1:
                                                            {
                                                                string start_index = Split(atom_name, "-").at(0);
                                                                string end_index = Split(atom_name, "-").at(1);
                                                                int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                int start_index_int = ConvertString<int>(start_index);
                                                                int end_index_int = ConvertString<int>(end_index);
                                                                if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                        }
                                                        break;
                                                    }
                                                }

                                            }
                                        }
                                    }
                                    else
                                    {
                                        AtomVector atoms = residue->GetAtoms();
                                        for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                        {
                                            Atom* atom = *it4;
                                            selection.push_back(atom);
                                        }
                                    }
                                }
                                break;
                            }
                            case 1:   /// Searching the residue set by maching the end of the residue names
                            {
                                if(residue->GetName().find(residue_name) != string::npos &&
                                        residue->GetName().find(residue_name) == residue->GetName().size() - residue_name.size())
                                {
                                    vector<string> atom_names = (*it1).second;
                                    if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                    {
                                        AtomVector atoms = residue->GetAtoms();
                                        for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                        {
                                            string atom_name = *it3;
                                            int atom_name_search_type = 0;
                                            if(atom_name.at(0) == '^')
                                            {
                                                atom_name = atom_name.substr(1);
                                                atom_name_search_type = -1;
                                            }
                                            else if(atom_name.at(atom_name.size() - 1) == '$')
                                            {
                                                atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                atom_name_search_type = 1;
                                            }
                                            else if(atom_name.at(0) == '#')
                                            {
                                                atom_name = atom_name.substr(1);
                                                atom_name_search_type = -2;
                                            }
                                            else
                                                atom_name_search_type = 0;
                                            for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                            {
                                                Atom* atom = *it4;
                                                switch(atom_name_search_type)
                                                {
                                                    case 0:
                                                    {
                                                        if(atom->GetName().compare(atom_name) == 0)
                                                            selection.push_back(atom);
                                                        break;
                                                    }
                                                    case 1:
                                                    {
                                                        if(atom->GetName().find(atom_name) != string::npos &&
                                                                atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                            selection.push_back(atom);
                                                        break;
                                                    }
                                                    case -1:
                                                    {
                                                        if(atom->GetName().find(atom_name) != string::npos &&
                                                                atom->GetName().find(atom_name) == 0)
                                                            selection.push_back(atom);
                                                        break;
                                                    }
                                                    case -2:
                                                    {
                                                        string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                        int range_selection = 0;
                                                        if(atom_name.find("-") != string::npos)
                                                        {
                                                            range_selection = 1;

                                                        }
                                                        else
                                                            range_selection = 0;
                                                        switch(range_selection)
                                                        {
                                                            case 0:
                                                            {
                                                                if(atom_serial_number.find(atom_name) != string::npos &&
                                                                        atom_serial_number.find(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case 1:
                                                            {
                                                                string start_index = Split(atom_name, "-").at(0);
                                                                string end_index = Split(atom_name, "-").at(1);
                                                                int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                int start_index_int = ConvertString<int>(start_index);
                                                                int end_index_int = ConvertString<int>(end_index);
                                                                if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                        }
                                                        break;
                                                    }
                                                }

                                            }
                                        }
                                    }
                                    else
                                    {
                                        AtomVector atoms = residue->GetAtoms();
                                        for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                        {
                                            Atom* atom = *it4;
                                            selection.push_back(atom);
                                        }
                                    }
                                }
                                break;
                            }
                            case -1:  /// Searching the residue set by matching the begining of the residue names
                            {
                                if(residue->GetName().find(residue_name) != string::npos &&
                                        residue->GetName().find(residue_name) == 0)
                                {
                                    vector<string> atom_names = (*it1).second;
                                    if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                    {
                                        AtomVector atoms = residue->GetAtoms();
                                        for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                        {
                                            string atom_name = *it3;
                                            int atom_name_search_type = 0;
                                            if(atom_name.at(0) == '^')
                                            {
                                                atom_name = atom_name.substr(1);
                                                atom_name_search_type = -1;
                                            }
                                            else if(atom_name.at(atom_name.size() - 1) == '$')
                                            {
                                                atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                atom_name_search_type = 1;
                                            }
                                            else if(atom_name.at(0) == '#')
                                            {
                                                atom_name = atom_name.substr(1);
                                                atom_name_search_type = -2;
                                            }
                                            else
                                                atom_name_search_type = 0;
                                            for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                            {
                                                Atom* atom = *it4;
                                                switch(atom_name_search_type)
                                                {
                                                    case 0:
                                                    {
                                                        if(atom->GetName().compare(atom_name) == 0)
                                                            selection.push_back(atom);
                                                        break;
                                                    }
                                                    case 1:
                                                    {
                                                        if(atom->GetName().find(atom_name) != string::npos &&
                                                                atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                            selection.push_back(atom);
                                                        break;
                                                    }
                                                    case -1:
                                                    {
                                                        if(atom->GetName().find(atom_name) != string::npos &&
                                                                atom->GetName().find(atom_name) == 0)
                                                            selection.push_back(atom);
                                                        break;
                                                    }
                                                    case -2:
                                                    {
                                                        string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                        int range_selection = 0;
                                                        if(atom_name.find("-") != string::npos)
                                                        {
                                                            range_selection = 1;

                                                        }
                                                        else
                                                            range_selection = 0;
                                                        switch(range_selection)
                                                        {
                                                            case 0:
                                                            {
                                                                if(atom_serial_number.find(atom_name) != string::npos &&
                                                                        atom_serial_number.find(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case 1:
                                                            {
                                                                string start_index = Split(atom_name, "-").at(0);
                                                                string end_index = Split(atom_name, "-").at(1);
                                                                int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                int start_index_int = ConvertString<int>(start_index);
                                                                int end_index_int = ConvertString<int>(end_index);
                                                                if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                        }
                                                        break;
                                                    }
                                                }

                                            }
                                        }
                                    }
                                    else
                                    {
                                        AtomVector atoms = residue->GetAtoms();
                                        for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                        {
                                            Atom* atom = *it4;
                                            selection.push_back(atom);
                                        }
                                    }
                                }
                                break;
                            }
                            case -2:  /// Searching the residue set by matching the residue sequence number
                            {
                                string residue_sequence_number = Split(residue->GetId(), "_").at(2);
                                int range_residue_selection = 0;
                                if(residue_sequence_number.find("-") != string::npos)
                                {
                                    range_residue_selection = 1;
                                }
                                else
                                    range_residue_selection = 0;
                                switch(range_residue_selection)
                                {
                                    case 0:
                                    {
                                        if(residue_sequence_number.find(residue_name) != string::npos &&
                                                residue_sequence_number.find(residue_name) == 0)
                                        {
                                            vector<string> atom_names = (*it1).second;
                                            if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                            {
                                                AtomVector atoms = residue->GetAtoms();
                                                for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                                {
                                                    string atom_name = *it3;
                                                    int atom_name_search_type = 0;
                                                    if(atom_name.at(0) == '^')
                                                    {
                                                        atom_name = atom_name.substr(1);
                                                        atom_name_search_type = -1;
                                                    }
                                                    else if(atom_name.at(atom_name.size() - 1) == '$')
                                                    {
                                                        atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                        atom_name_search_type = 1;
                                                    }
                                                    else if(atom_name.at(0) == '#')
                                                    {
                                                        atom_name = atom_name.substr(1);
                                                        atom_name_search_type = -2;
                                                    }
                                                    else
                                                        atom_name_search_type = 0;
                                                    for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                    {
                                                        Atom* atom = *it4;
                                                        switch(atom_name_search_type)
                                                        {
                                                            case 0:
                                                            {
                                                                if(atom->GetName().compare(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case 1:
                                                            {
                                                                if(atom->GetName().find(atom_name) != string::npos &&
                                                                        atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case -1:
                                                            {
                                                                if(atom->GetName().find(atom_name) != string::npos &&
                                                                        atom->GetName().find(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case -2:
                                                            {
                                                                string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                                int range_selection = 0;
                                                                if(atom_name.find("-") != string::npos)
                                                                {
                                                                    range_selection = 1;

                                                                }
                                                                else
                                                                    range_selection = 0;
                                                                switch(range_selection)
                                                                {
                                                                    case 0:
                                                                    {
                                                                        if(atom_serial_number.find(atom_name) != string::npos &&
                                                                                atom_serial_number.find(atom_name) == 0)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                    case 1:
                                                                    {
                                                                        string start_index = Split(atom_name, "-").at(0);
                                                                        string end_index = Split(atom_name, "-").at(1);
                                                                        int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                        int start_index_int = ConvertString<int>(start_index);
                                                                        int end_index_int = ConvertString<int>(end_index);
                                                                        if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                }
                                                                break;
                                                            }
                                                        }

                                                    }
                                                }
                                            }
                                            else
                                            {
                                                AtomVector atoms = residue->GetAtoms();
                                                for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                {
                                                    Atom* atom = *it4;
                                                    selection.push_back(atom);
                                                }
                                            }
                                        }
                                        break;
                                    }
                                    case 1:
                                    {
                                        string residue_start_index = Split(residue_name, "-").at(0);
                                        string residue_end_index = Split(residue_name, "-").at(1);
                                        int residue_sequence_number_int = ConvertString<int>(residue_sequence_number);
                                        int residue_start_index_int = ConvertString<int>(residue_start_index);
                                        int residue_end_index_int = ConvertString<int>(residue_end_index);
                                        if(residue_sequence_number_int >= residue_start_index_int && residue_sequence_number_int <= residue_end_index_int)
                                        {
                                            vector<string> atom_names = (*it1).second;
                                            if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                            {
                                                AtomVector atoms = residue->GetAtoms();
                                                for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                                {
                                                    string atom_name = *it3;
                                                    int atom_name_search_type = 0;
                                                    if(atom_name.at(0) == '^')
                                                    {
                                                        atom_name = atom_name.substr(1);
                                                        atom_name_search_type = -1;
                                                    }
                                                    else if(atom_name.at(atom_name.size() - 1) == '$')
                                                    {
                                                        atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                        atom_name_search_type = 1;
                                                    }
                                                    else if(atom_name.at(0) == '#')
                                                    {
                                                        atom_name = atom_name.substr(1);
                                                        atom_name_search_type = -2;
                                                    }
                                                    else
                                                        atom_name_search_type = 0;
                                                    for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                    {
                                                        Atom* atom = *it4;
                                                        switch(atom_name_search_type)
                                                        {
                                                            case 0:
                                                            {
                                                                if(atom->GetName().compare(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case 1:
                                                            {
                                                                if(atom->GetName().find(atom_name) != string::npos &&
                                                                        atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case -1:
                                                            {
                                                                if(atom->GetName().find(atom_name) != string::npos &&
                                                                        atom->GetName().find(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case -2:
                                                            {
                                                                string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                                int range_selection = 0;
                                                                if(atom_name.find("-") != string::npos)
                                                                {
                                                                    range_selection = 1;

                                                                }
                                                                else
                                                                    range_selection = 0;
                                                                switch(range_selection)
                                                                {
                                                                    case 0:
                                                                    {
                                                                        if(atom_serial_number.find(atom_name) != string::npos &&
                                                                                atom_serial_number.find(atom_name) == 0)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                    case 1:
                                                                    {
                                                                        string start_index = Split(atom_name, "-").at(0);
                                                                        string end_index = Split(atom_name, "-").at(1);
                                                                        int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                        int start_index_int = ConvertString<int>(start_index);
                                                                        int end_index_int = ConvertString<int>(end_index);
                                                                        if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                }
                                                                break;
                                                            }
                                                        }

                                                    }
                                                }
                                            }
                                            else
                                            {
                                                AtomVector atoms = residue->GetAtoms();
                                                for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                {
                                                    Atom* atom = *it4;
                                                    selection.push_back(atom);
                                                }
                                            }
                                            break;
                                        }
                                    }
                                }
                                break;
                            }
                        }
                    }
                }
                else
                {
                    for(ResidueVector::iterator it2 = residues_of_assembly.begin(); it2 != residues_of_assembly.end(); it2++)
                    {
                        Residue* residue = *it2;
                        vector<string> atom_names = (*it1).second;
                        if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                        {
                            AtomVector atoms = residue->GetAtoms();
                            for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                            {
                                string atom_name = *it3;
                                int atom_name_search_type = 0;
                                if(atom_name.at(0) == '^')
                                {
                                    atom_name = atom_name.substr(1);
                                    atom_name_search_type = -1;
                                }
                                else if(atom_name.at(atom_name.size() - 1) == '$')
                                {
                                    atom_name = atom_name.substr(0, atom_name.size() - 2);
                                    atom_name_search_type = 1;
                                }
                                else if(atom_name.at(0) == '#')
                                {
                                    atom_name = atom_name.substr(1);
                                    atom_name_search_type = -2;
                                }
                                else
                                    atom_name_search_type = 0;
                                for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                {
                                    Atom* atom = *it4;
                                    switch(atom_name_search_type)
                                    {
                                        case 0:
                                        {
                                            if(atom->GetName().compare(atom_name) == 0)
                                                selection.push_back(atom);
                                            break;
                                        }
                                        case 1:
                                        {
                                            if(atom->GetName().find(atom_name) != string::npos &&
                                                    atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                selection.push_back(atom);
                                            break;
                                        }
                                        case -1:
                                        {
                                            if(atom->GetName().find(atom_name) != string::npos &&
                                                    atom->GetName().find(atom_name) == 0)
                                                selection.push_back(atom);
                                            break;
                                        }
                                        case -2:
                                        {
                                            string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                            int range_selection = 0;
                                            if(atom_name.find("-") != string::npos)
                                            {
                                                range_selection = 1;

                                            }
                                            else
                                                range_selection = 0;
                                            switch(range_selection)
                                            {
                                                case 0:
                                                {
                                                    if(atom_serial_number.find(atom_name) != string::npos &&
                                                            atom_serial_number.find(atom_name) == 0)
                                                        selection.push_back(atom);
                                                    break;
                                                }
                                                case 1:
                                                {
                                                    string start_index = Split(atom_name, "-").at(0);
                                                    string end_index = Split(atom_name, "-").at(1);
                                                    int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                    int start_index_int = ConvertString<int>(start_index);
                                                    int end_index_int = ConvertString<int>(end_index);
                                                    if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                        selection.push_back(atom);
                                                    break;
                                                }
                                            }
                                            break;
                                        }
                                    }

                                }
                            }
                        }
                        else
                        {
                            AtomVector atoms = residue->GetAtoms();
                            for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                            {
                                Atom* atom = *it4;
                                selection.push_back(atom);
                            }
                        }
                    }
                }

            }

        }
        else
        {
            int start_index = key.find("*");
            if(start_index == 0)
            {
                for(HierarchicalContainmentMap::iterator it5 = hierarchical_map.begin(); it5 != hierarchical_map.end(); it5++)
                {
                    ResidueVector residues_of_assembly = (*it5).second;
                    map<string, vector<string> > value = (*it).second;
                    for(map<string, vector<string> >::iterator it1 = value.begin(); it1 != value.end(); it1++)
                    {
                        string residue_name = (*it1).first;
                        int residue_name_search_type = 0;
                        if(residue_name.at(0) == '^')
                        {
                            residue_name = residue_name.substr(1);
                            residue_name_search_type = -1;
                        }
                        else if(residue_name.at(residue_name.size() - 1) == '$')
                        {
                            residue_name = residue_name.substr(0, residue_name.size() - 2);
                            residue_name_search_type = 1;
                        }
                        else if(residue_name.at(0) == '#')
                        {
                            residue_name = residue_name.substr(1);
                            residue_name_search_type = -2;
                        }
                        else
                            residue_name_search_type = 0;
                        if(residue_name.find("*") == string::npos)
                        {
                            for(ResidueVector::iterator it2 = residues_of_assembly.begin(); it2 != residues_of_assembly.end(); it2++)
                            {
                                Residue* residue = *it2;
                                switch(residue_name_search_type)
                                {
                                    case 0:  /// Search in residue set by matching the whole name of the residue
                                    {
                                        if(residue->GetName().compare(residue_name) == 0)
                                        {
                                            vector<string> atom_names = (*it1).second;
                                            if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                            {
                                                AtomVector atoms = residue->GetAtoms();
                                                for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                                {
                                                    string atom_name = *it3;
                                                    int atom_name_search_type = 0;
                                                    if(atom_name.at(0) == '^')
                                                    {
                                                        atom_name = atom_name.substr(1);
                                                        atom_name_search_type = -1;
                                                    }
                                                    else if(atom_name.at(atom_name.size() - 1) == '$')
                                                    {
                                                        atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                        atom_name_search_type = 1;
                                                    }
                                                    else if(atom_name.at(0) == '#')
                                                    {
                                                        atom_name = atom_name.substr(1);
                                                        atom_name_search_type = -2;
                                                    }
                                                    else
                                                        atom_name_search_type = 0;
                                                    for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                    {
                                                        Atom* atom = *it4;
                                                        switch(atom_name_search_type)
                                                        {
                                                            case 0:
                                                            {
                                                                if(atom->GetName().compare(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case 1:
                                                            {
                                                                if(atom->GetName().find(atom_name) != string::npos &&
                                                                        atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case -1:
                                                            {
                                                                if(atom->GetName().find(atom_name) != string::npos &&
                                                                        atom->GetName().find(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case -2:
                                                            {
                                                                string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                                int range_selection = 0;
                                                                if(atom_name.find("-") != string::npos)
                                                                {
                                                                    range_selection = 1;

                                                                }
                                                                else
                                                                    range_selection = 0;
                                                                switch(range_selection)
                                                                {
                                                                    case 0:
                                                                    {
                                                                        if(atom_serial_number.find(atom_name) != string::npos &&
                                                                                atom_serial_number.find(atom_name) == 0)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                    case 1:
                                                                    {
                                                                        string start_index = Split(atom_name, "-").at(0);
                                                                        string end_index = Split(atom_name, "-").at(1);
                                                                        int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                        int start_index_int = ConvertString<int>(start_index);
                                                                        int end_index_int = ConvertString<int>(end_index);
                                                                        if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                }
                                                                break;
                                                            }
                                                        }

                                                    }
                                                }
                                            }
                                            else
                                            {
                                                AtomVector atoms = residue->GetAtoms();
                                                for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                {
                                                    Atom* atom = *it4;
                                                    selection.push_back(atom);
                                                }
                                            }
                                        }
                                        break;
                                    }
                                    case 1:   /// Searching the residue set by maching the end of the residue names
                                    {
                                        if(residue->GetName().find(residue_name) != string::npos &&
                                                residue->GetName().find(residue_name) == residue->GetName().size() - residue_name.size())
                                        {
                                            vector<string> atom_names = (*it1).second;
                                            if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                            {
                                                AtomVector atoms = residue->GetAtoms();
                                                for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                                {
                                                    string atom_name = *it3;
                                                    int atom_name_search_type = 0;
                                                    if(atom_name.at(0) == '^')
                                                    {
                                                        atom_name = atom_name.substr(1);
                                                        atom_name_search_type = -1;
                                                    }
                                                    else if(atom_name.at(atom_name.size() - 1) == '$')
                                                    {
                                                        atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                        atom_name_search_type = 1;
                                                    }
                                                    else if(atom_name.at(0) == '#')
                                                    {
                                                        atom_name = atom_name.substr(1);
                                                        atom_name_search_type = -2;
                                                    }
                                                    else
                                                        atom_name_search_type = 0;
                                                    for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                    {
                                                        Atom* atom = *it4;
                                                        switch(atom_name_search_type)
                                                        {
                                                            case 0:
                                                            {
                                                                if(atom->GetName().compare(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case 1:
                                                            {
                                                                if(atom->GetName().find(atom_name) != string::npos &&
                                                                        atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case -1:
                                                            {
                                                                if(atom->GetName().find(atom_name) != string::npos &&
                                                                        atom->GetName().find(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case -2:
                                                            {
                                                                string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                                int range_selection = 0;
                                                                if(atom_name.find("-") != string::npos)
                                                                {
                                                                    range_selection = 1;

                                                                }
                                                                else
                                                                    range_selection = 0;
                                                                switch(range_selection)
                                                                {
                                                                    case 0:
                                                                    {
                                                                        if(atom_serial_number.find(atom_name) != string::npos &&
                                                                                atom_serial_number.find(atom_name) == 0)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                    case 1:
                                                                    {
                                                                        string start_index = Split(atom_name, "-").at(0);
                                                                        string end_index = Split(atom_name, "-").at(1);
                                                                        int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                        int start_index_int = ConvertString<int>(start_index);
                                                                        int end_index_int = ConvertString<int>(end_index);
                                                                        if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                }
                                                                break;
                                                            }
                                                        }

                                                    }
                                                }
                                            }
                                            else
                                            {
                                                AtomVector atoms = residue->GetAtoms();
                                                for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                {
                                                    Atom* atom = *it4;
                                                    selection.push_back(atom);
                                                }
                                            }
                                        }
                                        break;
                                    }
                                    case -1:  /// Searching the residue set by matching the begining of the residue names
                                    {
                                        if(residue->GetName().find(residue_name) != string::npos &&
                                                residue->GetName().find(residue_name) == 0)
                                        {
                                            vector<string> atom_names = (*it1).second;
                                            if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                            {
                                                AtomVector atoms = residue->GetAtoms();
                                                for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                                {
                                                    string atom_name = *it3;
                                                    int atom_name_search_type = 0;
                                                    if(atom_name.at(0) == '^')
                                                    {
                                                        atom_name = atom_name.substr(1);
                                                        atom_name_search_type = -1;
                                                    }
                                                    else if(atom_name.at(atom_name.size() - 1) == '$')
                                                    {
                                                        atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                        atom_name_search_type = 1;
                                                    }
                                                    else if(atom_name.at(0) == '#')
                                                    {
                                                        atom_name = atom_name.substr(1);
                                                        atom_name_search_type = -2;
                                                    }
                                                    else
                                                        atom_name_search_type = 0;
                                                    for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                    {
                                                        Atom* atom = *it4;
                                                        switch(atom_name_search_type)
                                                        {
                                                            case 0:
                                                            {
                                                                if(atom->GetName().compare(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case 1:
                                                            {
                                                                if(atom->GetName().find(atom_name) != string::npos &&
                                                                        atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case -1:
                                                            {
                                                                if(atom->GetName().find(atom_name) != string::npos &&
                                                                        atom->GetName().find(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case -2:
                                                            {
                                                                string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                                int range_selection = 0;
                                                                if(atom_name.find("-") != string::npos)
                                                                {
                                                                    range_selection = 1;

                                                                }
                                                                else
                                                                    range_selection = 0;
                                                                switch(range_selection)
                                                                {
                                                                    case 0:
                                                                    {
                                                                        if(atom_serial_number.find(atom_name) != string::npos &&
                                                                                atom_serial_number.find(atom_name) == 0)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                    case 1:
                                                                    {
                                                                        string start_index = Split(atom_name, "-").at(0);
                                                                        string end_index = Split(atom_name, "-").at(1);
                                                                        int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                        int start_index_int = ConvertString<int>(start_index);
                                                                        int end_index_int = ConvertString<int>(end_index);
                                                                        if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                }
                                                                break;
                                                            }
                                                        }

                                                    }
                                                }
                                            }
                                            else
                                            {
                                                AtomVector atoms = residue->GetAtoms();
                                                for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                {
                                                    Atom* atom = *it4;
                                                    selection.push_back(atom);
                                                }
                                            }
                                        }
                                        break;
                                    }
                                    case -2:  /// Searching the residue set by matching the residue sequence number
                                    {
                                        string residue_sequence_number = Split(residue->GetId(), "_").at(2);
                                        int range_residue_selection = 0;
                                        if(residue_sequence_number.find("-") != string::npos)
                                        {
                                            range_residue_selection = 1;
                                        }
                                        else
                                            range_residue_selection = 0;
                                        switch(range_residue_selection)
                                        {
                                            case 0:
                                            {
                                                if(residue_sequence_number.find(residue_name) != string::npos &&
                                                        residue_sequence_number.find(residue_name) == 0)
                                                {
                                                    vector<string> atom_names = (*it1).second;
                                                    if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                                    {
                                                        AtomVector atoms = residue->GetAtoms();
                                                        for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                                        {
                                                            string atom_name = *it3;
                                                            int atom_name_search_type = 0;
                                                            if(atom_name.at(0) == '^')
                                                            {
                                                                atom_name = atom_name.substr(1);
                                                                atom_name_search_type = -1;
                                                            }
                                                            else if(atom_name.at(atom_name.size() - 1) == '$')
                                                            {
                                                                atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                                atom_name_search_type = 1;
                                                            }
                                                            else if(atom_name.at(0) == '#')
                                                            {
                                                                atom_name = atom_name.substr(1);
                                                                atom_name_search_type = -2;
                                                            }
                                                            else
                                                                atom_name_search_type = 0;
                                                            for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                            {
                                                                Atom* atom = *it4;
                                                                switch(atom_name_search_type)
                                                                {
                                                                    case 0:
                                                                    {
                                                                        if(atom->GetName().compare(atom_name) == 0)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                    case 1:
                                                                    {
                                                                        if(atom->GetName().find(atom_name) != string::npos &&
                                                                                atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                    case -1:
                                                                    {
                                                                        if(atom->GetName().find(atom_name) != string::npos &&
                                                                                atom->GetName().find(atom_name) == 0)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                    case -2:
                                                                    {
                                                                        string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                                        int range_selection = 0;
                                                                        if(atom_name.find("-") != string::npos)
                                                                        {
                                                                            range_selection = 1;

                                                                        }
                                                                        else
                                                                            range_selection = 0;
                                                                        switch(range_selection)
                                                                        {
                                                                            case 0:
                                                                            {
                                                                                if(atom_serial_number.find(atom_name) != string::npos &&
                                                                                        atom_serial_number.find(atom_name) == 0)
                                                                                    selection.push_back(atom);
                                                                                break;
                                                                            }
                                                                            case 1:
                                                                            {
                                                                                string start_index = Split(atom_name, "-").at(0);
                                                                                string end_index = Split(atom_name, "-").at(1);
                                                                                int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                                int start_index_int = ConvertString<int>(start_index);
                                                                                int end_index_int = ConvertString<int>(end_index);
                                                                                if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                                    selection.push_back(atom);
                                                                                break;
                                                                            }
                                                                        }
                                                                        break;
                                                                    }
                                                                }

                                                            }
                                                        }
                                                    }
                                                    else
                                                    {
                                                        AtomVector atoms = residue->GetAtoms();
                                                        for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                        {
                                                            Atom* atom = *it4;
                                                            selection.push_back(atom);
                                                        }
                                                    }
                                                }
                                                break;
                                            }
                                            case 1:
                                            {
                                                string residue_start_index = Split(residue_name, "-").at(0);
                                                string residue_end_index = Split(residue_name, "-").at(1);
                                                int residue_sequence_number_int = ConvertString<int>(residue_sequence_number);
                                                int residue_start_index_int = ConvertString<int>(residue_start_index);
                                                int residue_end_index_int = ConvertString<int>(residue_end_index);
                                                if(residue_sequence_number_int >= residue_start_index_int && residue_sequence_number_int <= residue_end_index_int)
                                                {
                                                    vector<string> atom_names = (*it1).second;
                                                    if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                                    {
                                                        AtomVector atoms = residue->GetAtoms();
                                                        for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                                        {
                                                            string atom_name = *it3;
                                                            int atom_name_search_type = 0;
                                                            if(atom_name.at(0) == '^')
                                                            {
                                                                atom_name = atom_name.substr(1);
                                                                atom_name_search_type = -1;
                                                            }
                                                            else if(atom_name.at(atom_name.size() - 1) == '$')
                                                            {
                                                                atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                                atom_name_search_type = 1;
                                                            }
                                                            else if(atom_name.at(0) == '#')
                                                            {
                                                                atom_name = atom_name.substr(1);
                                                                atom_name_search_type = -2;
                                                            }
                                                            else
                                                                atom_name_search_type = 0;
                                                            for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                            {
                                                                Atom* atom = *it4;
                                                                switch(atom_name_search_type)
                                                                {
                                                                    case 0:
                                                                    {
                                                                        if(atom->GetName().compare(atom_name) == 0)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                    case 1:
                                                                    {
                                                                        if(atom->GetName().find(atom_name) != string::npos &&
                                                                                atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                    case -1:
                                                                    {
                                                                        if(atom->GetName().find(atom_name) != string::npos &&
                                                                                atom->GetName().find(atom_name) == 0)
                                                                            selection.push_back(atom);
                                                                        break;
                                                                    }
                                                                    case -2:
                                                                    {
                                                                        string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                                        int range_selection = 0;
                                                                        if(atom_name.find("-") != string::npos)
                                                                        {
                                                                            range_selection = 1;

                                                                        }
                                                                        else
                                                                            range_selection = 0;
                                                                        switch(range_selection)
                                                                        {
                                                                            case 0:
                                                                            {
                                                                                if(atom_serial_number.find(atom_name) != string::npos &&
                                                                                        atom_serial_number.find(atom_name) == 0)
                                                                                    selection.push_back(atom);
                                                                                break;
                                                                            }
                                                                            case 1:
                                                                            {
                                                                                string start_index = Split(atom_name, "-").at(0);
                                                                                string end_index = Split(atom_name, "-").at(1);
                                                                                int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                                int start_index_int = ConvertString<int>(start_index);
                                                                                int end_index_int = ConvertString<int>(end_index);
                                                                                if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                                    selection.push_back(atom);
                                                                                break;
                                                                            }
                                                                        }
                                                                        break;
                                                                    }
                                                                }

                                                            }
                                                        }
                                                    }
                                                    else
                                                    {
                                                        AtomVector atoms = residue->GetAtoms();
                                                        for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                        {
                                                            Atom* atom = *it4;
                                                            selection.push_back(atom);
                                                        }
                                                    }
                                                    break;
                                                }
                                            }
                                        }
                                        break;
                                    }
                                }
                            }
                        }
                        else
                        {
                            for(ResidueVector::iterator it2 = residues_of_assembly.begin(); it2 != residues_of_assembly.end(); it2++)
                            {
                                Residue* residue = *it2;
                                vector<string> atom_names = (*it1).second;
                                if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                {
                                    AtomVector atoms = residue->GetAtoms();
                                    for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                    {
                                        string atom_name = *it3;
                                        int atom_name_search_type = 0;
                                        if(atom_name.at(0) == '^')
                                        {
                                            atom_name = atom_name.substr(1);
                                            atom_name_search_type = -1;
                                        }
                                        else if(atom_name.at(atom_name.size() - 1) == '$')
                                        {
                                            atom_name = atom_name.substr(0, atom_name.size() - 2);
                                            atom_name_search_type = 1;
                                        }
                                        else if(atom_name.at(0) == '#')
                                        {
                                            atom_name = atom_name.substr(1);
                                            atom_name_search_type = -2;
                                        }
                                        else
                                            atom_name_search_type = 0;
                                        for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                        {
                                            Atom* atom = *it4;
                                            switch(atom_name_search_type)
                                            {
                                                case 0:
                                                {
                                                    if(atom->GetName().compare(atom_name) == 0)
                                                        selection.push_back(atom);
                                                    break;
                                                }
                                                case 1:
                                                {
                                                    if(atom->GetName().find(atom_name) != string::npos &&
                                                            atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                        selection.push_back(atom);
                                                    break;
                                                }
                                                case -1:
                                                {
                                                    if(atom->GetName().find(atom_name) != string::npos &&
                                                            atom->GetName().find(atom_name) == 0)
                                                        selection.push_back(atom);
                                                    break;
                                                }
                                                case -2:
                                                {
                                                    string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                    int range_selection = 0;
                                                    if(atom_name.find("-") != string::npos)
                                                    {
                                                        range_selection = 1;

                                                    }
                                                    else
                                                        range_selection = 0;
                                                    switch(range_selection)
                                                    {
                                                        case 0:
                                                        {
                                                            if(atom_serial_number.find(atom_name) != string::npos &&
                                                                    atom_serial_number.find(atom_name) == 0)
                                                                selection.push_back(atom);
                                                            break;
                                                        }
                                                        case 1:
                                                        {
                                                            string start_index = Split(atom_name, "-").at(0);
                                                            string end_index = Split(atom_name, "-").at(1);
                                                            int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                            int start_index_int = ConvertString<int>(start_index);
                                                            int end_index_int = ConvertString<int>(end_index);
                                                            if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                selection.push_back(atom);
                                                            break;
                                                        }
                                                    }
                                                    break;
                                                }
                                            }

                                        }
                                    }
                                }
                                else
                                {
                                    AtomVector atoms = residue->GetAtoms();
                                    for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                    {
                                        Atom* atom = *it4;
                                        selection.push_back(atom);
                                    }
                                }
                            }
                        }

                    }

                }
            }
            else
            {
                string partial_key = key.substr(0, key.find(".*"));
                for(HierarchicalContainmentMap::iterator it5 = hierarchical_map.begin(); it5 != hierarchical_map.end(); it5++)
                {
                    if((*it5).first.find(partial_key) == 0)
                    {
                        ResidueVector residues_of_assembly = (*it5).second;
                        map<string, vector<string> > value = (*it).second;
                        for(map<string, vector<string> >::iterator it1 = value.begin(); it1 != value.end(); it1++)
                        {
                            string residue_name = (*it1).first;
                            int residue_name_search_type = 0;
                            if(residue_name.at(0) == '^')
                            {
                                residue_name = residue_name.substr(1);
                                residue_name_search_type = -1;
                            }
                            else if(residue_name.at(residue_name.size() - 1) == '$')
                            {
                                residue_name = residue_name.substr(0, residue_name.size() - 2);
                                residue_name_search_type = 1;
                            }
                            else if(residue_name.at(0) == '#')
                            {
                                residue_name = residue_name.substr(1);
                                residue_name_search_type = -2;
                            }
                            else
                                residue_name_search_type = 0;
                            if(residue_name.find("*") == string::npos)
                            {
                                for(ResidueVector::iterator it2 = residues_of_assembly.begin(); it2 != residues_of_assembly.end(); it2++)
                                {
                                    Residue* residue = *it2;
                                    switch(residue_name_search_type)
                                    {
                                        case 0:  /// Search in residue set by matching the whole name of the residue
                                        {
                                            if(residue->GetName().compare(residue_name) == 0)
                                            {
                                                vector<string> atom_names = (*it1).second;
                                                if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                                {
                                                    AtomVector atoms = residue->GetAtoms();
                                                    for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                                    {
                                                        string atom_name = *it3;
                                                        int atom_name_search_type = 0;
                                                        if(atom_name.at(0) == '^')
                                                        {
                                                            atom_name = atom_name.substr(1);
                                                            atom_name_search_type = -1;
                                                        }
                                                        else if(atom_name.at(atom_name.size() - 1) == '$')
                                                        {
                                                            atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                            atom_name_search_type = 1;
                                                        }
                                                        else if(atom_name.at(0) == '#')
                                                        {
                                                            atom_name = atom_name.substr(1);
                                                            atom_name_search_type = -2;
                                                        }
                                                        else
                                                            atom_name_search_type = 0;
                                                        for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                        {
                                                            Atom* atom = *it4;
                                                            switch(atom_name_search_type)
                                                            {
                                                                case 0:
                                                                {
                                                                    if(atom->GetName().compare(atom_name) == 0)
                                                                        selection.push_back(atom);
                                                                    break;
                                                                }
                                                                case 1:
                                                                {
                                                                    if(atom->GetName().find(atom_name) != string::npos &&
                                                                            atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                                        selection.push_back(atom);
                                                                    break;
                                                                }
                                                                case -1:
                                                                {
                                                                    if(atom->GetName().find(atom_name) != string::npos &&
                                                                            atom->GetName().find(atom_name) == 0)
                                                                        selection.push_back(atom);
                                                                    break;
                                                                }
                                                                case -2:
                                                                {
                                                                    string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                                    int range_selection = 0;
                                                                    if(atom_name.find("-") != string::npos)
                                                                    {
                                                                        range_selection = 1;

                                                                    }
                                                                    else
                                                                        range_selection = 0;
                                                                    switch(range_selection)
                                                                    {
                                                                        case 0:
                                                                        {
                                                                            if(atom_serial_number.find(atom_name) != string::npos &&
                                                                                    atom_serial_number.find(atom_name) == 0)
                                                                                selection.push_back(atom);
                                                                            break;
                                                                        }
                                                                        case 1:
                                                                        {
                                                                            string start_index = Split(atom_name, "-").at(0);
                                                                            string end_index = Split(atom_name, "-").at(1);
                                                                            int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                            int start_index_int = ConvertString<int>(start_index);
                                                                            int end_index_int = ConvertString<int>(end_index);
                                                                            if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                                selection.push_back(atom);
                                                                            break;
                                                                        }
                                                                    }
                                                                    break;
                                                                }
                                                            }

                                                        }
                                                    }
                                                }
                                                else
                                                {
                                                    AtomVector atoms = residue->GetAtoms();
                                                    for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                    {
                                                        Atom* atom = *it4;
                                                        selection.push_back(atom);
                                                    }
                                                }
                                            }
                                            break;
                                        }
                                        case 1:   /// Searching the residue set by maching the end of the residue names
                                        {
                                            if(residue->GetName().find(residue_name) != string::npos &&
                                                    residue->GetName().find(residue_name) == residue->GetName().size() - residue_name.size())
                                            {
                                                vector<string> atom_names = (*it1).second;
                                                if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                                {
                                                    AtomVector atoms = residue->GetAtoms();
                                                    for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                                    {
                                                        string atom_name = *it3;
                                                        int atom_name_search_type = 0;
                                                        if(atom_name.at(0) == '^')
                                                        {
                                                            atom_name = atom_name.substr(1);
                                                            atom_name_search_type = -1;
                                                        }
                                                        else if(atom_name.at(atom_name.size() - 1) == '$')
                                                        {
                                                            atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                            atom_name_search_type = 1;
                                                        }
                                                        else if(atom_name.at(0) == '#')
                                                        {
                                                            atom_name = atom_name.substr(1);
                                                            atom_name_search_type = -2;
                                                        }
                                                        else
                                                            atom_name_search_type = 0;
                                                        for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                        {
                                                            Atom* atom = *it4;
                                                            switch(atom_name_search_type)
                                                            {
                                                                case 0:
                                                                {
                                                                    if(atom->GetName().compare(atom_name) == 0)
                                                                        selection.push_back(atom);
                                                                    break;
                                                                }
                                                                case 1:
                                                                {
                                                                    if(atom->GetName().find(atom_name) != string::npos &&
                                                                            atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                                        selection.push_back(atom);
                                                                    break;
                                                                }
                                                                case -1:
                                                                {
                                                                    if(atom->GetName().find(atom_name) != string::npos &&
                                                                            atom->GetName().find(atom_name) == 0)
                                                                        selection.push_back(atom);
                                                                    break;
                                                                }
                                                                case -2:
                                                                {
                                                                    string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                                    int range_selection = 0;
                                                                    if(atom_name.find("-") != string::npos)
                                                                    {
                                                                        range_selection = 1;

                                                                    }
                                                                    else
                                                                        range_selection = 0;
                                                                    switch(range_selection)
                                                                    {
                                                                        case 0:
                                                                        {
                                                                            if(atom_serial_number.find(atom_name) != string::npos &&
                                                                                    atom_serial_number.find(atom_name) == 0)
                                                                                selection.push_back(atom);
                                                                            break;
                                                                        }
                                                                        case 1:
                                                                        {
                                                                            string start_index = Split(atom_name, "-").at(0);
                                                                            string end_index = Split(atom_name, "-").at(1);
                                                                            int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                            int start_index_int = ConvertString<int>(start_index);
                                                                            int end_index_int = ConvertString<int>(end_index);
                                                                            if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                                selection.push_back(atom);
                                                                            break;
                                                                        }
                                                                    }
                                                                    break;
                                                                }
                                                            }

                                                        }
                                                    }
                                                }
                                                else
                                                {
                                                    AtomVector atoms = residue->GetAtoms();
                                                    for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                    {
                                                        Atom* atom = *it4;
                                                        selection.push_back(atom);
                                                    }
                                                }
                                            }
                                            break;
                                        }
                                        case -1:  /// Searching the residue set by matching the begining of the residue names
                                        {
                                            if(residue->GetName().find(residue_name) != string::npos &&
                                                    residue->GetName().find(residue_name) == 0)
                                            {
                                                vector<string> atom_names = (*it1).second;
                                                if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                                {
                                                    AtomVector atoms = residue->GetAtoms();
                                                    for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                                    {
                                                        string atom_name = *it3;
                                                        int atom_name_search_type = 0;
                                                        if(atom_name.at(0) == '^')
                                                        {
                                                            atom_name = atom_name.substr(1);
                                                            atom_name_search_type = -1;
                                                        }
                                                        else if(atom_name.at(atom_name.size() - 1) == '$')
                                                        {
                                                            atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                            atom_name_search_type = 1;
                                                        }
                                                        else if(atom_name.at(0) == '#')
                                                        {
                                                            atom_name = atom_name.substr(1);
                                                            atom_name_search_type = -2;
                                                        }
                                                        else
                                                            atom_name_search_type = 0;
                                                        for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                        {
                                                            Atom* atom = *it4;
                                                            switch(atom_name_search_type)
                                                            {
                                                                case 0:
                                                                {
                                                                    if(atom->GetName().compare(atom_name) == 0)
                                                                        selection.push_back(atom);
                                                                    break;
                                                                }
                                                                case 1:
                                                                {
                                                                    if(atom->GetName().find(atom_name) != string::npos &&
                                                                            atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                                        selection.push_back(atom);
                                                                    break;
                                                                }
                                                                case -1:
                                                                {
                                                                    if(atom->GetName().find(atom_name) != string::npos &&
                                                                            atom->GetName().find(atom_name) == 0)
                                                                        selection.push_back(atom);
                                                                    break;
                                                                }
                                                                case -2:
                                                                {
                                                                    string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                                    int range_selection = 0;
                                                                    if(atom_name.find("-") != string::npos)
                                                                    {
                                                                        range_selection = 1;

                                                                    }
                                                                    else
                                                                        range_selection = 0;
                                                                    switch(range_selection)
                                                                    {
                                                                        case 0:
                                                                        {
                                                                            if(atom_serial_number.find(atom_name) != string::npos &&
                                                                                    atom_serial_number.find(atom_name) == 0)
                                                                                selection.push_back(atom);
                                                                            break;
                                                                        }
                                                                        case 1:
                                                                        {
                                                                            string start_index = Split(atom_name, "-").at(0);
                                                                            string end_index = Split(atom_name, "-").at(1);
                                                                            int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                            int start_index_int = ConvertString<int>(start_index);
                                                                            int end_index_int = ConvertString<int>(end_index);
                                                                            if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                                selection.push_back(atom);
                                                                            break;
                                                                        }
                                                                    }
                                                                    break;
                                                                }
                                                            }

                                                        }
                                                    }
                                                }
                                                else
                                                {
                                                    AtomVector atoms = residue->GetAtoms();
                                                    for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                    {
                                                        Atom* atom = *it4;
                                                        selection.push_back(atom);
                                                    }
                                                }
                                            }
                                            break;
                                        }
                                        case -2:  /// Searching the residue set by matching the residue sequence number
                                        {
                                            string residue_sequence_number = Split(residue->GetId(), "_").at(2);
                                            int range_residue_selection = 0;
                                            if(residue_sequence_number.find("-") != string::npos)
                                            {
                                                range_residue_selection = 1;
                                            }
                                            else
                                                range_residue_selection = 0;
                                            switch(range_residue_selection)
                                            {
                                                case 0:
                                                {
                                                    if(residue_sequence_number.find(residue_name) != string::npos &&
                                                            residue_sequence_number.find(residue_name) == 0)
                                                    {
                                                        vector<string> atom_names = (*it1).second;
                                                        if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                                        {
                                                            AtomVector atoms = residue->GetAtoms();
                                                            for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                                            {
                                                                string atom_name = *it3;
                                                                int atom_name_search_type = 0;
                                                                if(atom_name.at(0) == '^')
                                                                {
                                                                    atom_name = atom_name.substr(1);
                                                                    atom_name_search_type = -1;
                                                                }
                                                                else if(atom_name.at(atom_name.size() - 1) == '$')
                                                                {
                                                                    atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                                    atom_name_search_type = 1;
                                                                }
                                                                else if(atom_name.at(0) == '#')
                                                                {
                                                                    atom_name = atom_name.substr(1);
                                                                    atom_name_search_type = -2;
                                                                }
                                                                else
                                                                    atom_name_search_type = 0;
                                                                for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                                {
                                                                    Atom* atom = *it4;
                                                                    switch(atom_name_search_type)
                                                                    {
                                                                        case 0:
                                                                        {
                                                                            if(atom->GetName().compare(atom_name) == 0)
                                                                                selection.push_back(atom);
                                                                            break;
                                                                        }
                                                                        case 1:
                                                                        {
                                                                            if(atom->GetName().find(atom_name) != string::npos &&
                                                                                    atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                                                selection.push_back(atom);
                                                                            break;
                                                                        }
                                                                        case -1:
                                                                        {
                                                                            if(atom->GetName().find(atom_name) != string::npos &&
                                                                                    atom->GetName().find(atom_name) == 0)
                                                                                selection.push_back(atom);
                                                                            break;
                                                                        }
                                                                        case -2:
                                                                        {
                                                                            string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                                            int range_selection = 0;
                                                                            if(atom_name.find("-") != string::npos)
                                                                            {
                                                                                range_selection = 1;

                                                                            }
                                                                            else
                                                                                range_selection = 0;
                                                                            switch(range_selection)
                                                                            {
                                                                                case 0:
                                                                                {
                                                                                    if(atom_serial_number.find(atom_name) != string::npos &&
                                                                                            atom_serial_number.find(atom_name) == 0)
                                                                                        selection.push_back(atom);
                                                                                    break;
                                                                                }
                                                                                case 1:
                                                                                {
                                                                                    string start_index = Split(atom_name, "-").at(0);
                                                                                    string end_index = Split(atom_name, "-").at(1);
                                                                                    int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                                    int start_index_int = ConvertString<int>(start_index);
                                                                                    int end_index_int = ConvertString<int>(end_index);
                                                                                    if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                                        selection.push_back(atom);
                                                                                    break;
                                                                                }
                                                                            }
                                                                            break;
                                                                        }
                                                                    }

                                                                }
                                                            }
                                                        }
                                                        else
                                                        {
                                                            AtomVector atoms = residue->GetAtoms();
                                                            for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                            {
                                                                Atom* atom = *it4;
                                                                selection.push_back(atom);
                                                            }
                                                        }
                                                    }
                                                    break;
                                                }
                                                case 1:
                                                {
                                                    string residue_start_index = Split(residue_name, "-").at(0);
                                                    string residue_end_index = Split(residue_name, "-").at(1);
                                                    int residue_sequence_number_int = ConvertString<int>(residue_sequence_number);
                                                    int residue_start_index_int = ConvertString<int>(residue_start_index);
                                                    int residue_end_index_int = ConvertString<int>(residue_end_index);
                                                    if(residue_sequence_number_int >= residue_start_index_int && residue_sequence_number_int <= residue_end_index_int)
                                                    {
                                                        vector<string> atom_names = (*it1).second;
                                                        if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                                        {
                                                            AtomVector atoms = residue->GetAtoms();
                                                            for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                                            {
                                                                string atom_name = *it3;
                                                                int atom_name_search_type = 0;
                                                                if(atom_name.at(0) == '^')
                                                                {
                                                                    atom_name = atom_name.substr(1);
                                                                    atom_name_search_type = -1;
                                                                }
                                                                else if(atom_name.at(atom_name.size() - 1) == '$')
                                                                {
                                                                    atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                                    atom_name_search_type = 1;
                                                                }
                                                                else if(atom_name.at(0) == '#')
                                                                {
                                                                    atom_name = atom_name.substr(1);
                                                                    atom_name_search_type = -2;
                                                                }
                                                                else
                                                                    atom_name_search_type = 0;
                                                                for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                                {
                                                                    Atom* atom = *it4;
                                                                    switch(atom_name_search_type)
                                                                    {
                                                                        case 0:
                                                                        {
                                                                            if(atom->GetName().compare(atom_name) == 0)
                                                                                selection.push_back(atom);
                                                                            break;
                                                                        }
                                                                        case 1:
                                                                        {
                                                                            if(atom->GetName().find(atom_name) != string::npos &&
                                                                                    atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                                                selection.push_back(atom);
                                                                            break;
                                                                        }
                                                                        case -1:
                                                                        {
                                                                            if(atom->GetName().find(atom_name) != string::npos &&
                                                                                    atom->GetName().find(atom_name) == 0)
                                                                                selection.push_back(atom);
                                                                            break;
                                                                        }
                                                                        case -2:
                                                                        {
                                                                            string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                                            int range_selection = 0;
                                                                            if(atom_name.find("-") != string::npos)
                                                                            {
                                                                                range_selection = 1;

                                                                            }
                                                                            else
                                                                                range_selection = 0;
                                                                            switch(range_selection)
                                                                            {
                                                                                case 0:
                                                                                {
                                                                                    if(atom_serial_number.find(atom_name) != string::npos &&
                                                                                            atom_serial_number.find(atom_name) == 0)
                                                                                        selection.push_back(atom);
                                                                                    break;
                                                                                }
                                                                                case 1:
                                                                                {
                                                                                    string start_index = Split(atom_name, "-").at(0);
                                                                                    string end_index = Split(atom_name, "-").at(1);
                                                                                    int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                                    int start_index_int = ConvertString<int>(start_index);
                                                                                    int end_index_int = ConvertString<int>(end_index);
                                                                                    if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                                        selection.push_back(atom);
                                                                                    break;
                                                                                }
                                                                            }
                                                                            break;
                                                                        }
                                                                    }

                                                                }
                                                            }
                                                        }
                                                        else
                                                        {
                                                            AtomVector atoms = residue->GetAtoms();
                                                            for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                                            {
                                                                Atom* atom = *it4;
                                                                selection.push_back(atom);
                                                            }
                                                        }
                                                        break;
                                                    }
                                                }
                                            }
                                            break;
                                        }
                                    }
                                }
                            }
                            else
                            {
                                for(ResidueVector::iterator it2 = residues_of_assembly.begin(); it2 != residues_of_assembly.end(); it2++)
                                {
                                    Residue* residue = *it2;
                                    vector<string> atom_names = (*it1).second;
                                    if(find(atom_names.begin(), atom_names.end(), "*") == atom_names.end())
                                    {
                                        AtomVector atoms = residue->GetAtoms();
                                        for(vector<string>::iterator it3 = atom_names.begin(); it3 != atom_names.end(); it3++)
                                        {
                                            string atom_name = *it3;
                                            int atom_name_search_type = 0;
                                            if(atom_name.at(0) == '^')
                                            {
                                                atom_name = atom_name.substr(1);
                                                atom_name_search_type = -1;
                                            }
                                            else if(atom_name.at(atom_name.size() - 1) == '$')
                                            {
                                                atom_name = atom_name.substr(0, atom_name.size() - 2);
                                                atom_name_search_type = 1;
                                            }
                                            else if(atom_name.at(0) == '#')
                                            {
                                                atom_name = atom_name.substr(1);
                                                atom_name_search_type = -2;
                                            }
                                            else
                                                atom_name_search_type = 0;
                                            for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                            {
                                                Atom* atom = *it4;
                                                switch(atom_name_search_type)
                                                {
                                                    case 0:
                                                    {
                                                        if(atom->GetName().compare(atom_name) == 0)
                                                            selection.push_back(atom);
                                                        break;
                                                    }
                                                    case 1:
                                                    {
                                                        if(atom->GetName().find(atom_name) != string::npos &&
                                                                atom->GetName().find(atom_name) == atom->GetName().size() - atom_name.size())
                                                            selection.push_back(atom);
                                                        break;
                                                    }
                                                    case -1:
                                                    {
                                                        if(atom->GetName().find(atom_name) != string::npos &&
                                                                atom->GetName().find(atom_name) == 0)
                                                            selection.push_back(atom);
                                                        break;
                                                    }
                                                    case -2:
                                                    {
                                                        string atom_serial_number = Split(atom->GetId(), "_").at(1);
                                                        int range_selection = 0;
                                                        if(atom_name.find("-") != string::npos)
                                                        {
                                                            range_selection = 1;

                                                        }
                                                        else
                                                            range_selection = 0;
                                                        switch(range_selection)
                                                        {
                                                            case 0:
                                                            {
                                                                if(atom_serial_number.find(atom_name) != string::npos &&
                                                                        atom_serial_number.find(atom_name) == 0)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                            case 1:
                                                            {
                                                                string start_index = Split(atom_name, "-").at(0);
                                                                string end_index = Split(atom_name, "-").at(1);
                                                                int atom_serial_number_int = ConvertString<int>(atom_serial_number);
                                                                int start_index_int = ConvertString<int>(start_index);
                                                                int end_index_int = ConvertString<int>(end_index);
                                                                if(atom_serial_number_int >= start_index_int && atom_serial_number_int <= end_index_int)
                                                                    selection.push_back(atom);
                                                                break;
                                                            }
                                                        }
                                                        break;
                                                    }
                                                }

                                            }
                                        }
                                    }
                                    else
                                    {
                                        AtomVector atoms = residue->GetAtoms();
                                        for(AtomVector::iterator it4 = atoms.begin(); it4 != atoms.end(); it4++)
                                        {
                                            Atom* atom = *it4;
                                            selection.push_back(atom);
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

    return selection;
}

Assembly::SelectPatternMap Assembly::ParsePatternString(string pattern)
{
    SelectPatternMap select_pattern_map = SelectPatternMap();
    vector<string> tokens = Split(pattern, ";");
    for(int i = 0; i < tokens.size(); i++)
    {
        string token = tokens.at(i);
        vector<string> assembly_tokens = Split(token, ":");
        string assembly_id_token = assembly_tokens.at(0);
        vector<string> assembly_ids = Split(assembly_id_token, ",");
        map<string, vector<string> > residues = map<string, vector<string> >();
        for(int j = 1; j < assembly_tokens.size(); j++)
        {
            string assembly_token = assembly_tokens.at(j);
            vector<string> assembly_residue_tokens = Split(assembly_token, "@");
            string residue_token = assembly_residue_tokens.at(0);
            string atom_token = assembly_residue_tokens.at(1);
            vector<string> residue_tokens = Split(residue_token, ",");
            vector<string> atom_tokens = Split(atom_token, ",");
            for(int k = 0; k < residue_tokens.size(); k++)
            {
                for(int l = 0; l < atom_tokens.size(); l++)
                {
                    if(find(residues[residue_tokens.at(k)].begin(), residues[residue_tokens.at(k)].end(), atom_tokens.at(l)) == residues[residue_tokens.at(k)].end())
                        residues[residue_tokens.at(k)].push_back(atom_tokens.at(l));
                }
            }
        }
        for(int j = 0; j < assembly_ids.size(); j++)
            select_pattern_map[assembly_ids.at(j)] = residues;
    }
    return select_pattern_map;
}

