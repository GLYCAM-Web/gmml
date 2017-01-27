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
void Assembly::DetectShape(AtomVector cycle, Monosaccharide* mono)
{
    ///Creating a new assembly only from the ring atoms for external detect shape program
    Assembly* detect_shape_assembly = new Assembly();
    detect_shape_assembly->AddResidue(cycle.at(0)->GetResidue());
    Residue* detect_shape_residue = detect_shape_assembly->GetResidues().at(0);
    for(int i = 0; i < cycle.size(); i++)
    {
        string name = cycle.at(i)->GetName();
        string id = cycle.at(i)->GetId();

        ///Preparing atom name and id for detect shape program. It doesn't work with atoms containing special characters: C', C* etc
        //        replace( id.begin(), id.end(), '?', 'n'); // replace all '?' with 'n'
        FindReplaceString(id, "\'", "");
        FindReplaceString(id, ",", "");
        FindReplaceString(name, "*", "");
        replace( id.begin(), id.end(), '*', 's'); // replace all '*' with ''

        //        replace( name.begin(), name.end(), '?', 'n'); // replace all '?' with 'n'
        FindReplaceString(name, "\'", "");
        FindReplaceString(name, ",", "");
        FindReplaceString(name, "*", "");

        cycle.at(i)->SetName(name);
        cycle.at(i)->SetId(id);
    }
    detect_shape_residue->SetAtoms(cycle);

    ///Write a new PDB file from the new assembly
    PdbFile* pdb = detect_shape_assembly->BuildPdbFileStructureFromAssembly();
    pdb->Write("temp_gmml_pdb.pdb");

    ///Converting the written PDB file to fomrat readable by detect_shape program
    string line = "";
    ifstream gmml_pdb ("temp_gmml_pdb.pdb");
    ofstream detect_shape_pdb ("temp_detect_shape_pdb.pdb");
    int n = 0;
    if (gmml_pdb.is_open())
    {
        while (!gmml_pdb.eof()) {
            getline(gmml_pdb, line);
            if(line.find("HETATM") != string::npos)
            {
                detect_shape_pdb << line << endl;
            }
            n++;
        }
        gmml_pdb.close();
        detect_shape_pdb.close();
    }
    else cout << "Unable to open temp_gmml_pdb.pdb file" << endl;

    ///Writing a configuration file for the second argument of the detect_sugar program
    ofstream detect_shape_configuration ("temp_config");
    detect_shape_configuration << "Atom" << endl;
    for(int i = 0; i < cycle.size(); i++)
    {
        detect_shape_configuration << cycle.at(i)->GetName() << endl;
    }
    detect_shape_configuration << "Residue" << endl;
    detect_shape_configuration << "1" << endl;
    detect_shape_configuration << "Path" << endl;
    detect_shape_configuration << "apps/BFMP/canonicals.txt" << endl;
    detect_shape_configuration.close();

    ///Calling detect_shape program
    system("apps/BFMP/detect_shape temp_detect_shape_pdb.pdb temp_config > /dev/null");

    ///Adding the BFMP ring conformation infomration gained from the detect_sugar program to the monosaccharide
    ifstream shape_detection_result ("ring_conformations.txt");
    line = "";
    if (shape_detection_result.is_open())
    {
        getline (shape_detection_result,line);
        getline (shape_detection_result,line);
        vector<string> line_tokens = Split(line, "\t");
        if(line_tokens.at(1).compare("-") == 0)
            mono->bfmp_ring_conformation_ = line_tokens.at(2);
        else
            mono->bfmp_ring_conformation_ = line_tokens.at(1);
        shape_detection_result.close();
    }
    else cout << "Unable to open ring_conformations.txt file from detect shape program" << endl;

    ///Deleting temporary files
    remove("temp_detect_shape_pdb.pdb");
    remove("temp_gmml_pdb.pdb");
    remove("temp_config");
    remove("ring_conformations.txt");
}

vector<Oligosaccharide*> Assembly::ExtractSugars(vector<string> amino_lib_files, bool glyprobity_report, bool populate_ontology)
{
    ResidueNameMap dataset_residue_names = GetAllResidueNamesFromMultipleLibFilesMap(amino_lib_files);

    ///CYCLE DETECTION
    CycleMap cycles = DetectCyclesByExhaustiveRingPerception();

    //    CycleMap cycles = DetectCyclesByDFS();

    ///PRINTING ALL DETECTED CYCLES
    cout << endl << "All detected cycles" << endl;
    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
    {
        string cycle_atoms_str = (*it).first;
        cout << cycle_atoms_str << endl;
    }

    ///FILTERING OUT FUSED CYCLES. aka Cycles that are sharing an edge
    RemoveFusedCycles(cycles);

    ///FILTERING OUT OXYGENLESS CYCLES
    FilterAllCarbonCycles(cycles);

    ///ANOMERIC CARBON DETECTION and SORTING
    cout << endl << "Cycles after discarding rings that are all-carbon" << endl;
    vector<string> anomeric_carbons_status = vector<string>();
    vector<Note*> anomeric_notes = vector<Note*>();
    CycleMap sorted_cycles = CycleMap();
    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
    {
        string cycle_atoms_str = (*it).first;
        AtomVector cycle_atoms = (*it).second;

        cout << cycle_atoms_str << endl; ///e.g. C1_3810_NAG_A_1521_?_?_1-O5_3821_NAG_A_1521_?_?_1-C5_3814_NAG_A_1521_?_?_1-C4_3813_NAG_A_1521_?_?_1-C3_3812_NAG_A_1521_?_?_1-C2_3811_NAG_A_1521_?_?_1

        Note* anomeric_note = new Note();
        Atom* anomeric = FindAnomericCarbon(anomeric_note, anomeric_carbons_status, cycle_atoms, cycle_atoms_str);
        anomeric_notes.push_back(anomeric_note);
        if(anomeric != NULL) ///Sorting the cycle atoms and adding it to the sorted_cycles map if an anomeric carbon identified, otherwise the structure can't be a sugar
        {
            AtomVector sorted_cycle_atoms = AtomVector();
            stringstream sorted_cycle_stream;
            sorted_cycle_atoms = SortCycle(cycle_atoms, anomeric, sorted_cycle_stream);
            sorted_cycles[sorted_cycle_stream.str()] = sorted_cycle_atoms;
        }
    }
    cycles = sorted_cycles;

    ///CREATING MONOSACCHARIDE STRUCTURE. Ring atoms, side atoms, chemical code (Glycode), modifications/derivatives, names
    cout << endl << "Detailed information of sorted cycles after discarding fused or oxygenless rings: " << endl;
    vector<Monosaccharide*> monos = vector<Monosaccharide*>();
    int mono_id = 0;
    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
    {
        string cycle_atoms_str = (*it).first;
        AtomVector cycle = (*it).second;

        Monosaccharide* mono = new Monosaccharide();
        int status_index = distance(cycles.begin(), it);
        mono->anomeric_status_ = anomeric_carbons_status.at(status_index);

        ///ASSIGNING RING ATOMS TO MONOSACCHARIDE OBJECT
        cout << "Ring atoms: " << cycle_atoms_str << endl;
        mono->cycle_atoms_str_ = cycle_atoms_str;
        mono->cycle_atoms_ = cycle;

        ///ASSIGNING SIDE ATOMS (EXCOCYCLIC ATOMS) TO MONOSACCHARIDE OBJECT
        vector<string> orientations = GetSideGroupOrientations(mono, cycle_atoms_str);

        ///PRINTING ASSIGNED SIDE ATOMS
        cout << "Side group atoms: " << endl;
        for(vector<AtomVector>::iterator it1 = mono->side_atoms_.begin(); it1 != mono->side_atoms_.end(); it1++)
        {
            AtomVector sides = (*it1);
            if(it1 == mono->side_atoms_.begin())///side atoms of anomeric carbon
            {
                if(sides.at(0) != NULL && sides.at(1) != NULL)
                    cout << "[1] -> " << sides.at(0)->GetId() << ", " << sides.at(1)->GetId() << endl;
                else if(sides.at(1) != NULL)
                    cout << "[1] -> " << sides.at(1)->GetId() << endl;
                else if(sides.at(0) != NULL)
                    cout << "[1] -> " << sides.at(0)->GetId() << endl;
            }
            else if(it1 == mono->side_atoms_.end() - 1)//side atoms of last carbon of the ring
            {
                cout << "[" << mono->cycle_atoms_.size() - 1 << "] -> ";
                if(sides.at(0) != NULL)
                    cout << sides.at(0)->GetId() << endl;
            }
            else if(sides.at(1) != NULL)
            {
                int cycle_atom_index = distance(mono->side_atoms_.begin(), it1);
                cout << "[" << cycle_atom_index + 1 << "] -> " << sides.at(1)->GetId() << endl;
            }
        }

        ///PRINTING ANOMERIC STATUS
        cout << mono->anomeric_status_ << mono->cycle_atoms_.at(0)->GetId() << endl;

        ///CREATING CHEMICAL CODE (Glycode) OBJECT
        ChemicalCode* code = BuildChemicalCode(orientations);
        if(code != NULL)
        {
            mono->chemical_code_ = code;
        }
        cout << endl << "Stereo chemistry chemical code:"  << endl;
        code->Print(cout);
        cout << endl;

        ///DETECT SHAPE USING BFMP EXTERNAL PROGRAM. Currently, the program does not work for furanoses
        if(cycle.size() > 5)
        {
            DetectShape(cycle, mono);
            if(mono->bfmp_ring_conformation_.compare("") != 0)
                cout << "BFMP ring conformation: " << mono->bfmp_ring_conformation_ << endl << endl; ///Part of Glyprobity report
        }

        ///CHECKING FOR +2 and +3 SIDE CARBONS
        AtomVector plus_sides = ExtractAdditionalSideAtoms(mono);

        ///FINDING CHEMICAL CODE IN NAME LOOKUP TABLE
        string code_str = code->toString();
        mono->sugar_name_ = SugarStereoChemistryNameLookup(code_str);

        ///DERIVATIVE/MODIFICATION PATTERN EXTRACTION
        ExtractDerivatives(mono);
        bool minus_one = false;
        if(mono->derivatives_map_.find("-1") != mono->derivatives_map_.end())
            minus_one = true;
        for(map<string, string>::iterator it1 = mono->derivatives_map_.begin(); it1 != mono->derivatives_map_.end(); it1++)
        {
            string key = (*it1).first;
            string value = (*it1).second;
            if(minus_one)
            {
                if(key.compare("-1") == 0)
                    cout << "Carbon at position 1 is attached to " << value << endl;
                else if(key.compare("a") == 0)
                    cout << "Carbon at position 2 is attached to " << value << endl;
                else if(key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                    cout << "Carbon at position " << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) + 1 << " is attached to " << value << endl;
                else
                    cout << "Carbon at position " << (gmml::ConvertString<int>(key) + 1) << " is attached to " << value << endl;
            }
            else
            {
                if(key.compare("a") == 0)
                    cout << "Carbon at position 1 is attached to " << value << endl;
                else if(key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                    cout << "Carbon at position " << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << " is attached to " << value << endl;
                else
                    cout << "Carbon at position " << key << " is attached to " << value << endl;
            }
        }

        ///GENERATING COMPLETE NAME
        if(plus_sides.size() <= 1)
        {
            ///COMPLETE NAME GENERATION BASED ON DERIVATIVE MAP
            GenerateCompleteSugarName(mono);
        }
        else///UPDATING SIDE ATOMS
        {
            vector<string>::iterator index_it;
            if((index_it = find(mono->chemical_code_->right_up_.begin(), mono->chemical_code_->right_up_.end(), "+1")) != mono->chemical_code_->right_up_.end()){}
            else if((index_it = find(mono->chemical_code_->right_down_.begin(), mono->chemical_code_->right_down_.end(), "+1")) != mono->chemical_code_->right_down_.end()){}

            ///CHECKING R or S
            stringstream plus_one;
            string orientation = CalculateRSOrientations(mono->cycle_atoms_.at(mono->cycle_atoms_.size() - 2), plus_sides.at(0), plus_sides.at(1));
            plus_one << "+1" << orientation;
            (*index_it) = plus_one.str();
            if(plus_sides.size() == 3)
            {
                stringstream plus_two;
                orientation = CalculateRSOrientations(plus_sides.at(0), plus_sides.at(1), plus_sides.at(2));
                plus_two << "+2" << orientation;
                mono->chemical_code_->right_up_.push_back(plus_two.str());
                mono->chemical_code_->right_up_.push_back("+3");
            }

            ///UPDATING CHEMICAL CODE
            UpdateComplexSugarChemicalCode(mono);

            ///PRINTING SIDE ATOMS OF COMPLEX STRUCTURE
            cout << "Complex structure side group atoms: " << endl;
            for(vector<AtomVector>::iterator it1 = mono->side_atoms_.begin(); it1 != mono->side_atoms_.end(); it1++)
            {
                stringstream complex_structure_side;
                AtomVector sides = (*it1);
                if(it1 == mono->side_atoms_.begin())///side atoms of anomeric carbon
                {
                    if(sides.at(0) != NULL && sides.at(1) != NULL)
                        complex_structure_side << "[1] -> " << sides.at(0)->GetId() << ", " << sides.at(1)->GetId() << endl;
                    else if(sides.at(1) != NULL)
                        complex_structure_side << "[1] -> " << sides.at(1)->GetId() << endl;
                    else if(sides.at(0) != NULL)
                        complex_structure_side << "[1] -> " << sides.at(0)->GetId() << endl;
                }
                else if(it1 == mono->side_atoms_.end() - 1)///side atoms of last carbon of the ring
                {
                    complex_structure_side << "[" << mono->cycle_atoms_.size() - 1 << "]";
                    for(int i = 0; i < plus_sides.size() ; i++)
                        complex_structure_side << " -> " << sides.at(i)->GetId();
                    cout << complex_structure_side.str() << endl;
                }
                else if(sides.at(1) != NULL)
                {
                    int cycle_atom_index = distance(mono->side_atoms_.begin(), it1);
                    complex_structure_side << "[" << cycle_atom_index + 1 << "] -> " << sides.at(1)->GetId() << endl;
                }
            }

            ///PRINTING COMPLEX SUGAR CHEMICAL CODE
            cout << endl << "Complex sugar chemical code:" << endl;
            mono->chemical_code_->Print(cout);

            ///FINDING COMPLEX CHEMICAL CODE IN COMPLEX SUGAR NAME LOOKUP TABLE
            mono->sugar_name_ = ComplexSugarNameLookup(mono->chemical_code_->toString());

            if(plus_sides.size() == 2)
            {
                ///COMPLETE NAME GENERATION BASED ON DERIVATIVE MAP
                GenerateCompleteSugarName(mono);
            }
        }
        cout << endl;

        if(mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") == 0 && mono->sugar_name_.monosaccharide_name_.compare("") == 0)
        {
            ///FINDING CLOSEST MATCH FOR THE CHEMICAL CODE IN THE LOOKUP TABLE
            vector<Glycan::SugarName> closest_matches = vector<Glycan::SugarName>();
            mono->sugar_name_ = ClosestMatchSugarStereoChemistryNameLookup(mono->chemical_code_->toString(), closest_matches);
            if(mono->sugar_name_.monosaccharide_name_.compare("") == 0)
                mono->sugar_name_.monosaccharide_name_ = mono->sugar_name_.monosaccharide_stereochemistry_name_;
            if(mono->sugar_name_.monosaccharide_short_name_.compare("") == 0)
                mono->sugar_name_.monosaccharide_short_name_ = mono->sugar_name_.monosaccharide_stereochemistry_short_name_;

            ///ADDING NOTES/ISSUES OF MONOSACCHARIDE STRUCTURE
            Note* matching_note = new Note();
            matching_note->type_ = Glycan::COMMENT;
            matching_note->category_ = Glycan::MONOSACCHARIDE;
            stringstream ss;
            ss << "No exact match for " << mono->sugar_name_.monosaccharide_stereochemistry_short_name_ << ". close matches: ";
            for(vector<Glycan::SugarName>::iterator ite = closest_matches.begin(); ite != closest_matches.end(); ite++)
            {
                Glycan::SugarName sn = (*ite);
                if(ite == closest_matches.end() - 1)
                    ss << sn.monosaccharide_stereochemistry_short_name_;
                else
                    ss << sn.monosaccharide_stereochemistry_short_name_ << ", ";
            }
            matching_note->description_ = ss.str();
            this->AddNote(matching_note);


            if(mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") == 0)
            {
                mono->sugar_name_.monosaccharide_stereochemistry_name_ = "Unknown";
                mono->sugar_name_.monosaccharide_stereochemistry_short_name_ = "Unknown";
                mono->sugar_name_.monosaccharide_name_ = "Unknown";
                mono->sugar_name_.monosaccharide_short_name_ = "Unknown";
            }
            else
                cout << "No exact match found for the chemical code, the following information comes from one of the closest matches:" << endl;
        }

        ///GLYPROBITY REPORT (GEOMETRY OUTLIERS)
        if(glyprobity_report)
            CalculateGlyprobityGeometryOutliers(mono);

        ///ADDING NOTES/ISSUES OF ANOMERIC CONFIGURATION
        Note* anomeric_note = anomeric_notes.at(status_index);
        stringstream n;
        if(anomeric_note->description_.compare("") != 0)
        {
            n << mono->sugar_name_.monosaccharide_short_name_ << ": " << anomeric_note->description_;
            anomeric_note->description_ = n.str();
            this->AddNote(anomeric_note);
        }

        ///ADDING NOTES/ISSUES OF RESIDUE NAMING
        string original_residue = mono->cycle_atoms_.at(0)->GetResidue()->GetName();
        vector<string> pdb_codes = Split(mono->sugar_name_.pdb_code_, ",");
        if(pdb_codes.size() > 0)
        {
            string pdb_code = "";
            bool found_code = false;
            for(vector<string>::iterator codes_it = pdb_codes.begin(); codes_it != pdb_codes.end(); codes_it++)
            {
                pdb_code = (*codes_it);
                if(pdb_code.compare(original_residue) == 0)
                {
                    found_code = true;
                    break;
                }
            }
            if(!found_code)
            {
                Note* residue_naming_note = new Note();
                residue_naming_note->category_ = Glycan::RESIDUE_NAME;
                stringstream res_ss;
                if(mono->sugar_name_.pdb_code_.compare("") == 0)
                {
                    residue_naming_note->type_ = Glycan::WARNING;
                    res_ss << "PDB 3 letter code not found for " << mono->sugar_name_.monosaccharide_short_name_;
                }
                else
                {
                    residue_naming_note->type_ = Glycan::ERROR;
                    res_ss << "Residue name in input PDB file for " << mono->sugar_name_.monosaccharide_short_name_ << " does not match to PDB code: " << mono->sugar_name_.pdb_code_;
                }
                residue_naming_note->description_ = res_ss.str();
                this->AddNote(residue_naming_note);
            }
        }

        ///PRINTING NAMES OF MONOSACCHARIDE
        cout << "Stereochemistry name: " << mono->sugar_name_.monosaccharide_stereochemistry_name_ << endl;
        cout << "Stereochemistry short name: " << mono->sugar_name_.monosaccharide_stereochemistry_short_name_ << endl;
        cout << "Complete name: " << mono->sugar_name_.monosaccharide_name_ << endl;
        cout << "Short name: " << mono->sugar_name_.monosaccharide_short_name_ << endl;

        cout << "-------------------------------------------------------------------------------------------------------------------------------------------" << endl;

        mono_id++;
        mono->mono_id = mono_id;
        monos.push_back(mono);
    }

    ///CREATING TREE-LIKE STRUCTURE OF OLIGOSACCHARIDE
    cout << endl << "Oligosaccharides:" << endl;
    int number_of_covalent_links = 0;
    int number_of_probable_non_covalent_complexes = 0;
    vector<Oligosaccharide*> oligosaccharides = ExtractOligosaccharides(monos, dataset_residue_names, number_of_covalent_links, number_of_probable_non_covalent_complexes);

    ///BUILDING OLIGOSACCHARIDE SEQUENCE
    int number_of_oligosaccharides = 0;
    int number_of_monosaccharides = 0;
    for(vector<Oligosaccharide*>::iterator it = oligosaccharides.begin(); it != oligosaccharides.end(); it++)
    {
        Oligosaccharide* oligo = (*it);
        if(oligo->child_oligos_linkages_.size() > 0)
            number_of_oligosaccharides++;
        else
            number_of_monosaccharides++;
        oligo->Print(cout);
    }

    ///PRINTING NOTES AND ISSUES FOUND WITH THE INPUT FILE
    vector<Note*> notes = this->GetNotes();
    cout << "-------------------------------------------------------------------------------------------------------------------------------------------" << endl;
    cout << endl << "NOTES/ISSUES:" << endl;
    for(vector<Note*>::iterator note_it = notes.begin(); note_it != notes.end(); note_it++)
    {
        Note* note = (*note_it);
        cout << endl << "Category: " << note->ConvertGlycanNoteCat2String(note->category_) << endl;
        cout << "Type: " << note->ConvertGlycanNoteType2String(note->type_) << endl;
        cout << "Description: " << note->description_ << endl;
    }

    cout << "-------------------------------------------------------------------------------------------------------------------------------------------" << endl;

    ///PRINTING STATISTICAL REPORT OF GLYPROBITY
    if(glyprobity_report)
    {
        cout << endl << "GLYPROBITY REPORT" << endl;
        cout << "<-------Topology------>" << endl;
        cout << "Monosaccharide detected: " << monos.size() << endl;
        cout << "Residue Distribution " << endl;
        cout << " Monosaccharides: " << number_of_monosaccharides << endl;
        cout << " Oligosaccharides: " << number_of_oligosaccharides << endl;
        cout << "Carbohydrate Context " << endl;
        cout << " Covalently linked to protein: " << number_of_covalent_links << endl;
        cout << " Non-covalent complex: " << number_of_probable_non_covalent_complexes << endl;
        cout << "<--------------------->" << endl;
    }

    ///POPULATING GMMO ONTOLOGY
    if(populate_ontology)
    {
        if(oligosaccharides.size() > 0)
        {
            string gmmo = "gmmo.ttl";
            std::ofstream out_file;
            out_file.open(gmmo.c_str(), fstream::app);

            ifstream in("gmmo.ttl");///Checking if the file is empty
            size_t out_file_size = 0;
            in.seekg(0,ios_base::end);
            out_file_size = in.tellg();
            in.close();
            if(out_file_size == 0) ///If the file is empty add the prefixes first
                out_file << Ontology::TTL_FILE_PREFIX << endl;

            this->PopulateOntology(out_file, oligosaccharides);
            out_file.close();
        }
    }

    return oligosaccharides;
}

bool Assembly::MatchDisaccharide(queue<Oligosaccharide*> oligo_queue, double &phi_angle, double &psi_angle, string first_mono, char mono1_carbon_index, string second_mono, char mono2_carbon_index)
{
    Oligosaccharide* oligo = oligo_queue.front();
    oligo_queue.pop();
    bool found_disaccharide = false;

    Oligosaccharide* corresponding_second_oligo = NULL;
    ///second mono of the input disaccharide comes first in the tree structure of oligo (parent of first mono)
    if(oligo->root_->sugar_name_.monosaccharide_short_name_.compare(second_mono) == 0 )///found a mono with the same name as the right side mono of the disaccharide
        corresponding_second_oligo = oligo;

    OligosaccharideVector child_oligos = oligo->child_oligos_;
    for(OligosaccharideVector::iterator it1 = child_oligos.begin(); it1 != child_oligos.end(); it1++)
    {
        oligo_queue.push(*it1);

        if(corresponding_second_oligo != NULL) ///if current mono matches the disaccharide, look for a linked mono that macthes the other mono in the disaccharide
        {
            Oligosaccharide* corresponding_first_oligo = (*it1);
            oligo_queue.push(corresponding_first_oligo);

            if(corresponding_first_oligo->root_->sugar_name_.monosaccharide_short_name_.compare(first_mono) == 0)///found a mono with the same name as the left side mono of the disaccharide
            {
                vector<string> child_links = corresponding_second_oligo->child_oligos_linkages_;///links from right mono of the disaccharide to child monos

                vector<string> mono2_cycle_atom_tokens = Split(corresponding_second_oligo->root_->cycle_atoms_str_, "-");
                for(vector<string>::iterator it2 = child_links.begin(); it2 != child_links.end(); it2++)
                {
                    //                int index = distance(child_links.begin(), it2);
                    string link = (*it2);
                    vector<string> link_tokens = Split(link, "-");
                    int parent_c_index = 0;
                    int child_c_index = 0;

                    if(corresponding_second_oligo->root_->side_atoms_.at(0).at(0) != NULL)
                        parent_c_index++;
                    for(int i = 0; i < mono2_cycle_atom_tokens.size(); i++)
                    {
                        parent_c_index++;
                        if(mono2_cycle_atom_tokens.at(i).compare(link_tokens.at(0)) == 0)
                            break;
                    }
                    vector<string> mono1_cycle_atom_tokens = Split(corresponding_first_oligo->root_->cycle_atoms_str_, "-");
                    if(corresponding_first_oligo->root_->side_atoms_.at(0).at(0) != NULL)
                        child_c_index++;
                    for(int i = 0; i < mono1_cycle_atom_tokens.size(); i++)
                    {
                        child_c_index++;
                        if(mono1_cycle_atom_tokens.at(i).compare(link_tokens.at(2)) == 0)
                            break;
                    }

                    if(mono2_carbon_index == parent_c_index + '0' && mono1_carbon_index == child_c_index + '0')///indeces matched the indeces of the given disaccharide
                    {
                        found_disaccharide = true;
                        AtomVector mono1_ring_atoms = corresponding_first_oligo->root_->cycle_atoms_;

                        ///Preparing atoms for phi and psi angle
                        Atom* phi_atom1 = mono1_ring_atoms.at(mono1_ring_atoms.size() - 1); ///O5 or N5
                        Atom* phi_atom2 = new Atom(); ///C1
                        phi_atom2 = corresponding_first_oligo->root_->cycle_atoms_.at(0); ///anomeric carbon
                        Atom* phi_atom3 = NULL;
                        AtomVector phi_atom2_neighbors = phi_atom2->GetNode()->GetNodeNeighbors();
                        for(AtomVector::iterator it3 = phi_atom2_neighbors.begin(); it3 != phi_atom2_neighbors.end(); it3++)
                        {
                            Atom* atom2_neighbor= (*it3);
                            ///If the neighbor id is the same as the linkage intermediate atom id
                            if(atom2_neighbor->GetId().compare(link_tokens.at(1)) == 0)
                            {
                                phi_atom3 = atom2_neighbor; ///Ox
                                break;
                            }
                        }
                        if(phi_atom3 != NULL)
                        {
                            Atom* phi_atom4 = NULL;
                            AtomVector phi_atom3_neighbors = phi_atom3->GetNode()->GetNodeNeighbors();
                            for(AtomVector::iterator it4 = phi_atom3_neighbors.begin(); it4 != phi_atom3_neighbors.end(); it4++)
                            {
                                Atom* atom3_neighbor= (*it4);
                                ///If the neighbor id is the same as the linkage carbon at second mono side (first linkage atom from second mono side)
                                if(atom3_neighbor->GetId().compare(link_tokens.at(0)) == 0)
                                {
                                    phi_atom4 = atom3_neighbor; ///Cx
                                    phi_angle = CalculateTorsionAngleByAtoms(phi_atom1, phi_atom2, phi_atom3, phi_atom4); /// ϕ (O5′-C1′-Ox-Cx)
                                    AtomVector phi_atom4_neighbors = phi_atom4->GetNode()->GetNodeNeighbors();

                                    int atom4_index = ConvertString<int>(Split(phi_atom4->GetName(), "C*,\'").at(0));
                                    for(AtomVector::iterator it5 = phi_atom4_neighbors.begin(); it5 != phi_atom4_neighbors.end(); it5++)
                                    {
                                        Atom* atom4_neighbor = (*it5);
                                        string neighbor_name = atom4_neighbor->GetName();
                                        if(neighbor_name.find("C") != string::npos)
                                        {
                                            int neighbor_index = ConvertString<int>(Split(neighbor_name, "C*,\'").at(0));
                                            if(atom4_index > neighbor_index) ///Cx-1
                                            {
                                                psi_angle = CalculateTorsionAngleByAtoms(phi_atom2, phi_atom3, phi_atom4, atom4_neighbor); /// ψ (C1′-Ox-Cx-Cx−1)
                                                break;
                                            }
                                        }
                                    }
                                    break;
                                }
                            }
                        }
                    }
                    if(found_disaccharide)
                        break;
                }
                if(found_disaccharide)
                    break;
            }
        }
    }
    if(!found_disaccharide && !oligo_queue.empty())
        return MatchDisaccharide(oligo_queue, phi_angle, psi_angle, first_mono, mono1_carbon_index, second_mono, mono2_carbon_index);
    else
        return found_disaccharide;
}

void Assembly::ExtractRingAtomsInformation()
{
    ///CYCLE DETECTION
    CycleMap cycles = DetectCyclesByExhaustiveRingPerception();
    ///FILTERING OUT FUSED CYCLES
    RemoveFusedCycles(cycles);
    ///FILTERING OUT OXYGENLESS CYCLES
    FilterAllCarbonCycles(cycles);
    CycleMap sorted_cycles = CycleMap();
    vector<string> anomeric_carbons_status = vector<string>();
    ///ANOMERIC DETECTION and SORTING
    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
    {
        string cycle_atoms_str = (*it).first;
        AtomVector cycle_atoms = (*it).second;
        Note* anomeric_note = new Note();
        Atom* anomeric = FindAnomericCarbon(anomeric_note, anomeric_carbons_status, cycle_atoms, cycle_atoms_str);

        if(anomeric != NULL)
        {
            AtomVector sorted_cycle_atoms = AtomVector();
            stringstream sorted_cycle_stream;
            sorted_cycle_atoms = SortCycle(cycle_atoms, anomeric, sorted_cycle_stream);
            sorted_cycles[sorted_cycle_stream.str()] = sorted_cycle_atoms;
        }
    }
    cycles = sorted_cycles;
    cout << endl << "Detailed information of sorted cycles after discarding fused or oxygenless rings: " << endl;
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Detailed information of sorted cycles after discarding fused or oxygenless rings: ");

    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
    {
        Monosaccharide* mono = new Monosaccharide();
        int status_index = distance(cycles.begin(), it);
        mono->anomeric_status_ = anomeric_carbons_status.at(status_index);
        string cycle_atoms_str = (*it).first;
        AtomVector cycle = (*it).second;

        stringstream ring_atoms;
        ring_atoms << "Ring atoms: " << cycle_atoms_str;
        cout << ring_atoms.str() << endl;
        gmml::log(__LINE__, __FILE__,  gmml::INF, ring_atoms.str());

        stringstream ring_atom_names;
        ring_atom_names << "Ring atom names: " ;
        for(AtomVector::iterator it1 = cycle.begin(); it1 != cycle.end(); it1++)
        {
            Atom* atom = (*it1);
            vector<string> atom_tokens = Split(atom->GetId(), "_");
            if(it1 == cycle.end() - 1)
                ring_atom_names << atom_tokens.at(0);
            else
                ring_atom_names << atom_tokens.at(0) << ",";
        }
        cout << ring_atom_names.str() << endl;
        gmml::log(__LINE__, __FILE__,  gmml::INF, ring_atom_names.str());
    }
}

void Assembly::ReturnCycleAtoms(string src_id, Atom *current_atom, AtomIdAtomMap &atom_parent_map, AtomVector &cycle, stringstream &cycle_stream)
{
    cycle.push_back(current_atom);
    cycle_stream << current_atom->GetId() << "-";
    Atom* parent = atom_parent_map[current_atom->GetId()];
    if(src_id.compare(parent->GetId()) == 0)
    {
        cycle.push_back(parent);
        cycle_stream << parent->GetId();
        return;
    }
    ReturnCycleAtoms(src_id, parent, atom_parent_map, cycle, cycle_stream);
}

void Assembly::FilterAllCarbonCycles(CycleMap &cycles)
{
    map<string, bool> to_be_deleted_cycles = map<string, bool>();
    bool all_carbons = true;
    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
    {
        string cycle_str = (*it).first;
        AtomVector cycle_atoms = (*it).second;
        all_carbons = true;
        for(AtomVector::iterator it1 = cycle_atoms.begin(); it1 != cycle_atoms.end(); it1++)
        {
            Atom* atom = (*it1);
            atom->GetName().find("C") == string::npos;
            {
                all_carbons = false;
                break;
            }
        }
        if(all_carbons)
            to_be_deleted_cycles[cycle_str] = true;
    }
    CycleMap all_carbons_filtered_cycles = CycleMap();
    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
    {
        string cycle_str = (*it).first;
        AtomVector cycle_atoms = (*it).second;
        if(to_be_deleted_cycles.find(cycle_str) == to_be_deleted_cycles.end())
            all_carbons_filtered_cycles[cycle_str] = cycle_atoms;
    }
    cycles.clear();
    cycles = all_carbons_filtered_cycles;
}

void Assembly::RemoveFusedCycles(CycleMap &cycles)
{
    map<string, bool> to_be_deleted_cycles = map<string, bool>();
    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
    {
        string cycle_i_str = (*it).first; ///cycle i will be compared with all the other cycles
        AtomVector cycle_i_atoms = (*it).second;
        for(CycleMap::iterator it1 = cycles.begin(); it1 != cycles.end(); it1++)
        {
            if(it != it1)
            {
                string cycle_j_str = (*it1).first; ///cycle j to be compared with cycle i
                for(int i = 0; i < cycle_i_atoms.size(); i++)
                {
                    stringstream mutual_edge;
                    stringstream mutual_edge_reverse;
                    Atom* a1 = new Atom();
                    Atom* a2 = new Atom();
                    if(i == cycle_i_atoms.size() - 1)
                    {
                        a1 = cycle_i_atoms.at(i);
                        a2 = cycle_i_atoms.at(0);
                    }
                    else
                    {
                        a1 = cycle_i_atoms.at(i);
                        a2 = cycle_i_atoms.at(i + 1);
                    }
                    mutual_edge << a1->GetId() << "-" << a2->GetId();
                    mutual_edge_reverse << a2->GetId() << "-" << a1->GetId();
                    if(cycle_j_str.find(mutual_edge.str()) != string::npos || cycle_j_str.find(mutual_edge_reverse.str()) != string::npos)///mutual edge found
                    {
                        to_be_deleted_cycles[cycle_i_str] = true;
                        to_be_deleted_cycles[cycle_j_str] = true;
                        break;
                    }
                }
            }
        }
    }
    CycleMap fused_filtered_cycles = CycleMap();
    for(CycleMap::iterator it = cycles.begin(); it != cycles.end(); it++)
    {
        string cycle_str = (*it).first;
        AtomVector cycle_atoms = (*it).second;
        if(to_be_deleted_cycles.find(cycle_str) == to_be_deleted_cycles.end())
            fused_filtered_cycles[cycle_str] = cycle_atoms;
    }
    cycles.clear();
    cycles = fused_filtered_cycles;
}

Atom* Assembly::FindAnomericCarbon(Note* anomeric_note, vector<string>& anomeric_carbons_status, AtomVector cycle, string cycle_atoms_str)
{
    Atom* anomeric_carbon = new Atom();
    for(AtomVector::iterator it = cycle.begin(); it != cycle.end(); it++)
    {
        Atom* cycle_atom = (*it);
        if((cycle_atom->GetName().substr(0,1).compare("O") == 0 ))///find oxygen in ring
            //                && isdigit(ConvertString<char>(cycle_atom->GetName().substr(1,1)))))
        {
            AtomNode* node = cycle_atom->GetNode();
            AtomVector neighbors = node->GetNodeNeighbors();

            ///Check the first neighbor of oxygen
            Atom* o_neighbor1 = neighbors.at(0);
            AtomNode* o_neighbor1_node = o_neighbor1->GetNode();
            AtomVector o_neighbor1_neighbors = o_neighbor1_node->GetNodeNeighbors();
            for(AtomVector::iterator it1 = o_neighbor1_neighbors.begin(); it1 != o_neighbor1_neighbors.end(); it1++)///check if neighbor1 of oxygen has another oxygen or nitrogen neighbor
            {
                Atom* neighbor1_neighbor = (*it1);
                if(cycle_atoms_str.find(neighbor1_neighbor->GetId()) == string::npos ///if the neighbor is not one of the cycle atoms
                        && (neighbor1_neighbor->GetName().substr(0,1).compare("O") == 0 || neighbor1_neighbor->GetName().substr(0,1).compare("N") == 0)) ///if first element is "O" or "N"
                    //                        && isdigit(ConvertString<char>(neighbor1_neighbor->GetName().substr(1,1))))///if second element is a digit
                {
                    anomeric_carbon = o_neighbor1;
                    anomeric_carbons_status.push_back("Anomeric carbon: ");
                    anomeric_note->description_ = "";

                    return anomeric_carbon;
                }
            }

            ///Check the second neighbor of oxygen
            Atom* o_neighbor2 = neighbors.at(1);
            AtomNode* o_neighbor2_node = o_neighbor2->GetNode();
            AtomVector o_neighbor2_neighbors = o_neighbor2_node->GetNodeNeighbors();
            for(AtomVector::iterator it2 = o_neighbor2_neighbors.begin(); it2 != o_neighbor2_neighbors.end(); it2++)///check if neighbor2 of oxygen has another oxygen or nitrogen neighbor
            {
                Atom* neighbor2_neighbor = (*it2);
                if( cycle_atoms_str.find(neighbor2_neighbor->GetId()) == string::npos
                        && (neighbor2_neighbor->GetName().substr(0,1).compare("O") == 0 || neighbor2_neighbor->GetName().substr(0,1).compare("N") == 0))
                    //                        && isdigit(ConvertString<char>(neighbor2_neighbor->GetName().substr(1,1))))
                {
                    anomeric_carbon = o_neighbor2;
                    //                    cout << "Anomeric carbon is: " << anomeric_carbon->GetId() << endl;
                    anomeric_carbons_status.push_back("Anomeric carbon: ");
                    anomeric_note->description_ = "";
                    return anomeric_carbon;
                }
            }

            ///Check the order of the carbons based on their names to locate the anomeric
            stringstream ss1;
            for(int i = 0; i < o_neighbor1->GetName().size(); i++)
            {
                if(isdigit(o_neighbor1->GetName().at(i)) != 0)
                    ss1 << o_neighbor1->GetName().at(i);
            }
            stringstream ss2;
            for(int i = 0; i < o_neighbor2->GetName().size(); i++)
            {
                if(isdigit(o_neighbor2->GetName().at(i)) != 0)
                    ss2 << o_neighbor2->GetName().at(i);
            }
            if(ConvertString<int>(ss1.str()) < ConvertString<int>(ss2.str()))
            {
                anomeric_note->type_ = Glycan::WARNING;
                anomeric_note->category_ = Glycan::ANOMERIC;
                anomeric_note->description_ = "Anomeric oxygen is missing";
                anomeric_carbons_status.push_back("Anomeric carbon probably is: ");
                return o_neighbor1;
            }
            if(ConvertString<int>(ss2.str()) < ConvertString<int>(ss1.str()))
            {
                anomeric_note->type_ = Glycan::WARNING;
                anomeric_note->category_ = Glycan::ANOMERIC;
                anomeric_note->description_ = "Anomeric oxygen is missing";
                anomeric_carbons_status.push_back("Anomeric carbon probably is: ");
                return o_neighbor2;
            }

            ///Check non-ring neighbors of oxygen neighbors (the one without non-ring carbon is anomeric)
            bool neighbor2_is_anomeric = false;
            for(AtomVector::iterator it1 = o_neighbor1_neighbors.begin(); it1 != o_neighbor1_neighbors.end(); it1++)///check if neighbor1 of oxygen has non-ring carbon neighbor
            {
                Atom* neighbor1_neighbor = (*it1);
                if(cycle_atoms_str.find(neighbor1_neighbor->GetId()) == string::npos ///if the neighbor is not one of the cycle atoms
                        && (neighbor1_neighbor->GetName().substr(0,1).compare("C") == 0 )) ///if first element is "C"
                {
                    neighbor2_is_anomeric = true;
                    break;
                }
            }
            bool neighbor1_is_anomeric = false;
            for(AtomVector::iterator it1 = o_neighbor2_neighbors.begin(); it1 != o_neighbor2_neighbors.end(); it1++)///check if neighbor1 of oxygen has non-ring carbon neighbor
            {
                Atom* neighbor2_neighbor = (*it1);
                if(cycle_atoms_str.find(neighbor2_neighbor->GetId()) == string::npos ///if the neighbor is not one of the cycle atoms
                        && (neighbor2_neighbor->GetName().substr(0,1).compare("C") == 0 )) ///if first element is "C"
                {
                    neighbor1_is_anomeric = true;
                }
            }
            if(!neighbor1_is_anomeric)
            {
                anomeric_note->type_ = Glycan::WARNING;
                anomeric_note->category_ = Glycan::ANOMERIC;
                anomeric_note->description_ = "Anomeric oxygen is missing";
                anomeric_carbons_status.push_back("Anomeric carbon probably is: ");
                return o_neighbor2;
            }
            else if(!neighbor2_is_anomeric)
            {
                anomeric_note->type_ = Glycan::WARNING;
                anomeric_note->category_ = Glycan::ANOMERIC;
                anomeric_note->description_ = "Anomeric oxygen is missing";
                anomeric_carbons_status.push_back("Anomeric carbon probably is: ");
                return o_neighbor1;
            }
            anomeric_note->type_ = Glycan::WARNING;
            anomeric_note->category_ = Glycan::ANOMERIC;
            anomeric_note->description_ = "Anomeric oxygen is missing";
            anomeric_carbons_status.push_back("Not enough information to detect the anomeric carbon, it has been chosen randomely: ");
            return o_neighbor1;
        }
    }
    return NULL;
}

Assembly::AtomVector Assembly::SortCycle(AtomVector cycle, Atom *anomeric_atom, stringstream &sorted_cycle_stream)
{
    AtomVector sorted_cycle = AtomVector();
    for(AtomVector::iterator it = cycle.begin(); it != cycle.end(); it++)
    {
        Atom* atom = (*it);
        int index = distance(cycle.begin(), it);
        if(atom->GetId().compare(anomeric_atom->GetId()) == 0)
        {
            if(index == cycle.size() - 1)///anomeric atom is at the end of the cycle
            {
                sorted_cycle.push_back(anomeric_atom);///anomeric atom as the first atom of the cycle
                sorted_cycle_stream << anomeric_atom->GetId() << "-";
                anomeric_atom->SetIsRing(true);
                Atom* a0 = cycle.at(0);
                if(a0->GetName().substr(0,1).compare("O") == 0)///a0 is oxygen so the vector is in reverse order
                {
                    for(AtomVector::iterator it1 = it - 1; it1 != cycle.begin(); it1--)///atoms before the anomeric atom in reverse order
                    {
                        Atom* a = (*it1);
                        sorted_cycle.push_back(a);
                        sorted_cycle_stream << a->GetId() << "-";
                        a->SetIsRing(true);
                    }
                    sorted_cycle.push_back((*cycle.begin()));
                    sorted_cycle_stream << (*cycle.begin())->GetId();
                    (*cycle.begin())->SetIsRing(true);
                }
                else
                {
                    for(AtomVector::iterator it1 = cycle.begin(); it1 != it; it1++)///atoms before the anomeric atom from beginning of vector
                    {
                        Atom* a = (*it1);
                        sorted_cycle.push_back(a);
                        a->SetIsRing(true);
                        if(it1 == it - 1)
                            sorted_cycle_stream << a->GetId();
                        else
                            sorted_cycle_stream << a->GetId() << "-";
                    }
                }
            }
            else///anomeric is not at the end of the cycle
            {
                Atom* next_atom = cycle.at(index + 1);
                if(next_atom->GetName().substr(0,1).compare("O") == 0)///next atom is oxygen so the vector is in reverse order
                {
                    for(AtomVector::iterator it1 = it; it1 != cycle.begin(); it1--) ///atoms befor anomeric atom down to beginning of the vector
                    {
                        Atom* a_before = (*it1);
                        sorted_cycle.push_back(a_before);
                        sorted_cycle_stream << a_before->GetId() << "-";
                        a_before->SetIsRing(true);
                    }
                    sorted_cycle.push_back((*cycle.begin()));
                    sorted_cycle_stream << (*cycle.begin())->GetId() << "-";
                    (*cycle.begin())->SetIsRing(true);
                    for(AtomVector::iterator it2 = cycle.end() - 1; it2 != it; it2--)///atoms from end of the vector down to anomeric atom
                    {
                        Atom* atom_after = (*it2);
                        sorted_cycle.push_back(atom_after);
                        atom_after->SetIsRing(true);
                        if(it2 == it + 1)
                            sorted_cycle_stream << atom_after->GetId();
                        else
                            sorted_cycle_stream << atom_after->GetId() << "-";
                    }
                }
                else///oxygen is before the anomeric atom so the vector is in normal order
                {
                    for(AtomVector::iterator it1 = it; it1 != cycle.end(); it1++) ///atoms after anomeric atom to the end of the vector
                    {
                        Atom* atom_after = (*it1);
                        sorted_cycle.push_back(atom_after);
                        atom_after->SetIsRing(true);
                        //                        if(it1 == cycle.end() - 1)
                        //                            sorted_cycle_stream << atom_after->GetId();
                        //                        else
                        sorted_cycle_stream << atom_after->GetId() << "-";
                    }
                    for(AtomVector::iterator it2 = cycle.begin(); it2 != it; it2++)///atoms befor the anomeric atom from beginning of vector
                    {
                        Atom* atom_before = (*it2);
                        sorted_cycle.push_back(atom_before);
                        atom_before->SetIsRing(true);
                        if(it2 == it - 1)
                            sorted_cycle_stream << atom_before->GetId();
                        else
                            sorted_cycle_stream << atom_before->GetId() << "-";
                    }
                }
            }
        }
    }
    return sorted_cycle;
}

vector<string> Assembly::GetSideGroupOrientations(Monosaccharide* mono, string cycle_atoms_str)
{
    vector<string> orientations = vector<string>();
    vector<AtomVector> side_atoms = vector<AtomVector>();
    AtomVector default_atom_vector = AtomVector(3, NULL);

    for(AtomVector::iterator it = mono->cycle_atoms_.begin(); it != mono->cycle_atoms_.end() - 1; it++) ///iterate on cycle atoms except the oxygen in the ring
    {
        orientations.push_back("N");
        side_atoms.push_back(default_atom_vector);
        int index = distance(mono->cycle_atoms_.begin(), it);
        Atom* prev_atom = new Atom();
        Atom* current_atom = (*it);
        Atom* next_atom = new Atom();
        if(index == 0)///if the current atom is the anomeric atom
            prev_atom = mono->cycle_atoms_.at(mono->cycle_atoms_.size() - 1); ///previous atom is the oxygen(last atom of the sorted cycle)
        else
            prev_atom = mono->cycle_atoms_.at(index - 1);
        next_atom = mono->cycle_atoms_.at(index + 1);

        ///Calculating the plane based on the two ring neighbors of the current atom
        Coordinate prev_atom_coord = Coordinate(*prev_atom->GetCoordinates().at(model_index_));
        Coordinate current_atom_coord = Coordinate(*current_atom->GetCoordinates().at(model_index_));
        Coordinate next_atom_coord = Coordinate(*next_atom->GetCoordinates().at(model_index_));
        prev_atom_coord.operator -(current_atom_coord) ;
        next_atom_coord.operator -(current_atom_coord) ;
        Plane plane = Plane();
        plane.SetV1(prev_atom_coord);
        plane.SetV2(next_atom_coord);
        Coordinate normal_v = plane.GetUnitNormalVector();

        ///Calculating the orientation of the side atoms
        AtomNode* node = current_atom->GetNode();
        AtomVector neighbors = node->GetNodeNeighbors();
        int not_h_neighbors = 0;
        for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
        {
            Atom* neighbor = (*it1);
            string neighbor_id = neighbor->GetId();
            if(cycle_atoms_str.find(neighbor_id) == string::npos) ///if not one of the cycle atoms
            {
                if(neighbor->GetName().at(0) != 'H') ///deoxy check
                    not_h_neighbors++;
            }
        }
        if(not_h_neighbors != 0)
        {
            for(AtomVector::iterator it1 = neighbors.begin(); it1 != neighbors.end(); it1++)
            {
                Atom* neighbor = (*it1);
                string neighbor_id = neighbor->GetId();
                if(cycle_atoms_str.find(neighbor_id) == string::npos) ///if not one of the cycle atoms
                {
                    if(neighbor->GetName().at(0) != 'H') ///deoxy check
                        not_h_neighbors++;
                    Coordinate side_atom_coord = Coordinate(*neighbor->GetCoordinates().at(model_index_));
                    side_atom_coord.operator -(current_atom_coord);
                    side_atom_coord.Normalize();
                    double theta = acos(normal_v.DotProduct(side_atom_coord));

                    if(index == 0 && neighbor_id.at(0) == 'C')///if anomeric atom has a non-ring carbon neighbor
                    {
                        if(orientations.at(index).compare("N") == 0) ///if the position of non-ring oxygen or nitrogen hasn't been set yet
                        {
                            if(theta > (gmml::PI_RADIAN/2))
                            {
                                orientations.at(index) = "-1D";
                                side_atoms.at(index).at(0) = neighbor;
                            }
                            else
                            {
                                orientations.at(index) = "-1U";
                                side_atoms.at(index).at(0) = neighbor;
                            }
                            continue;
                        }
                        else
                        { ///position of non-ring oxygen or nitrogen + the non-ring carbon
                            stringstream ss;
                            if(theta > (gmml::PI_RADIAN/2))
                            {
                                ss << orientations.at(index) << "-1D";
                                orientations.at(index) = ss.str();
                                side_atoms.at(index).at(0) = neighbor;
                            }
                            else
                            {
                                ss << orientations.at(index) << "-1U";
                                orientations.at(index) = ss.str();
                                side_atoms.at(index).at(0) = neighbor;
                            }
                            break;
                        }
                    }
                    else if(neighbor_id.at(0) == 'O' || neighbor_id.at(0) == 'N')///if neighbor is a non-ring oxygen or nitrogen
                    {
                        if(index == 0)///current atom is anomeric
                        {
                            if(orientations.at(index).compare("N") == 0) ///if the position of non-ring carbon neighbor (if exist) hasn't been set yet
                            {
                                if(theta > (gmml::PI_RADIAN/2))
                                {
                                    orientations.at(index) = "D";
                                    side_atoms.at(index).at(1) = neighbor;
                                }
                                else
                                {
                                    orientations.at(index) = "U";
                                    side_atoms.at(index).at(1) = neighbor;
                                }
                                continue;
                            }
                            else
                            { ///position of non-ring oxygen or nitrogen + the non-ring carbon
                                stringstream ss;
                                if(theta > (gmml::PI_RADIAN/2))
                                {
                                    ss << "D" << orientations.at(index);
                                    orientations.at(index) = ss.str();
                                    side_atoms.at(index).at(1) = neighbor;
                                }
                                else
                                {
                                    ss << "U" << orientations.at(index);
                                    orientations.at(index) = ss.str();
                                    side_atoms.at(index).at(1) = neighbor;
                                }
                                break;
                            }
                        }
                        else
                        {
                            if(theta > (gmml::PI_RADIAN/2))
                            {
                                orientations.at(index) = "D";
                                side_atoms.at(index).at(1) = neighbor;
                            }
                            else
                            {
                                orientations.at(index) = "U";
                                side_atoms.at(index).at(1) = neighbor;
                            }
                            break;
                        }
                    }
                    else if(index == mono->cycle_atoms_.size() - 2 && neighbor_id.at(0) == 'C')///if the last ring carbon has a non-ring carbon neighbor
                    {
                        ///Check if neighbor of neighbor is oxygen or nitrogen
                        AtomNode* neighbor_node = neighbor->GetNode();
                        AtomVector neighbors_of_neighbor = neighbor_node->GetNodeNeighbors();
                        int o_neighbors = 0;
                        int not_h_neighbors = 0;
                        for(AtomVector::iterator it2 = neighbors_of_neighbor.begin(); it2 != neighbors_of_neighbor.end(); it2++)
                        {
                            Atom* neighbor_of_neighbor = (*it2);
                            string neighbor_of_neighbor_id = neighbor_of_neighbor->GetId();
                            if(cycle_atoms_str.find(neighbor_of_neighbor_id) == string::npos &&
                                    (neighbor_of_neighbor_id.at(0) == 'O' || neighbor_of_neighbor_id.at(0) == 'N'))///if neighbor of neighbor is a non-ring oxygen or nitrogen
                                o_neighbors++;
                            if(cycle_atoms_str.find(neighbor_of_neighbor_id) == string::npos &&
                                    (neighbor_of_neighbor_id.at(0) != 'H'))///if neighbor of neighbor is any non-ring atom other than hydrogen
                                not_h_neighbors++;
                        }
                        if (o_neighbors >= 1)
                        {
                            if(theta > (gmml::PI_RADIAN/2))
                            {
                                orientations.at(index) = "D";
                                side_atoms.at(index).at(0) = neighbor;
                            }
                            else
                            {
                                orientations.at(index) = "U";
                                side_atoms.at(index).at(0) = neighbor;
                            }
                        }
                        else if(not_h_neighbors == 0)///Type Deoxy
                        {
                            if(theta > (gmml::PI_RADIAN/2))
                            {
                                orientations.at(index) = "Dd";
                                side_atoms.at(index).at(0) = neighbor;
                            }
                            else
                            {
                                orientations.at(index) = "Ud";
                                side_atoms.at(index).at(0) = neighbor;
                            }
                            break;
                        }
                        if(orientations.at(index).compare("N") != 0)
                            break;
                    }
                }
            }
        }
    }
    mono->side_atoms_ = side_atoms;

    return orientations;
}

ChemicalCode* Assembly::BuildChemicalCode(vector<string> orientations)
{
    ChemicalCode* code = new ChemicalCode();
    if(orientations.size() == 5 )
        code->base_ = "P";
    else if(orientations.size() == 4 )
        code->base_ = "F";
    else
        code->base_ = "?";

    ///Side atom(s) of anomeric
    ///Has only non-ring oxygen neighbor
    if(orientations.at(0).compare("U") == 0)
        code->right_up_.push_back("a");
    else if(orientations.at(0).compare("D") == 0)
        code->right_down_.push_back("a");

    ///Has non-ring oxygen and carbon neighbors
    else if(orientations.at(0).compare("U-1U") == 0)
    {
        code->right_up_.push_back("a");
        code->right_up_.push_back("-1");
    }
    else if(orientations.at(0).compare("D-1U") == 0)
    {
        code->right_down_.push_back("a");
        code->right_up_.push_back("-1");
    }
    else if(orientations.at(0).compare("U-1D") == 0)
    {
        code->right_up_.push_back("a");
        code->right_down_.push_back("-1");
    }
    else if(orientations.at(0).compare("D-1D") == 0)
    {
        code->right_down_.push_back("a");
        code->right_down_.push_back("-1");
    }

    ///Has only non-ring carbon neighbor
    else if(orientations.at(0).compare("-1U") == 0)
        code->right_up_.push_back("-1");
    else if(orientations.at(0).compare("-1D") == 0)
        code->right_down_.push_back("-1");

    ///Side atom of other carbons of the ring
    for(vector<string>::iterator it = orientations.begin() + 1; it != orientations.end() - 1; it++)
    {
        string orientation = (*it);
        int index = distance(orientations.begin(), it);
        if(orientation.compare("U") == 0)
            code->left_up_.push_back(gmml::ConvertT(index + 1));
        else if(orientation.compare("D") == 0)
            code->left_down_.push_back(gmml::ConvertT(index + 1));
        else if(orientation.compare("N") == 0)
        {
            stringstream ss;
            ss << gmml::ConvertT(index + 1) << "d";
            code->left_middle_.push_back(ss.str() );
        }
    }

    ///Side atom(s) of last carbon
    if(orientations.at(orientations.size() - 1).compare("U") == 0)
        code->right_up_.push_back("+1");
    else if(orientations.at(orientations.size() - 1).compare("D") == 0)
        code->right_down_.push_back("+1");
    ///Type Deoxy
    else if(orientations.at(orientations.size() - 1).compare("Ud") == 0)
        code->right_up_.push_back("+1d");
    else if(orientations.at(orientations.size() - 1).compare("Dd") == 0)
        code->right_down_.push_back("+1d");

    return code;
}

Assembly::AtomVector Assembly::ExtractAdditionalSideAtoms(Monosaccharide *mono)
{
    AtomVector plus_sides = AtomVector();
    if(mono->side_atoms_.at(mono->side_atoms_.size() - 1).at(0) != NULL)///if there exist a +1 carbon atom. in side_atoms_ structure (vector<AtomVector>) the first index of the last element is dedicated to +1 atom
    {
        plus_sides.push_back(mono->side_atoms_.at(mono->side_atoms_.size() - 1).at(0));
        AtomVector plus_one_atom_neighbors = mono->side_atoms_.at(mono->side_atoms_.size() - 1).at(0)->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it1 = plus_one_atom_neighbors.begin(); it1 != plus_one_atom_neighbors.end(); it1++)
        {
            if((*it1)->GetName().at(0) == 'C' && mono->cycle_atoms_str_.find((*it1)->GetId()) == string::npos)///+2 carbon atom found
            {
                Atom* plus_two = (*it1);
                plus_sides.push_back(plus_two);
                mono->side_atoms_.at(mono->side_atoms_.size() - 1).at(1) = plus_two;///in side_atoms_ structure (vector<AtomVector>) the second index of the last element is dedicated to +2 atom

                AtomVector plus_two_atom_neighbors = plus_two->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it2 = plus_two_atom_neighbors.begin(); it2 != plus_two_atom_neighbors.end(); it2++)
                {
                    Atom* plus_three = (*it2);
                    if(plus_three->GetName().at(0) == 'C' && plus_three->GetId().compare(mono->side_atoms_.at(mono->side_atoms_.size() - 1).at(0)->GetId()) != 0)///+3 carbon atom found
                    {
                        plus_sides.push_back(plus_three);
                        mono->side_atoms_.at(mono->side_atoms_.size() - 1).at(2) = plus_three;///in side_atoms_ structure (vector<AtomVector>) the third index of the last element is dedicated to +3 atom
                        break;
                    }
                }
                break;
            }
        }
    }
    return plus_sides;
}

void Assembly::ExtractDerivatives(Monosaccharide * mono)
{
    for(AtomVector::iterator it = mono->cycle_atoms_.begin(); it != mono->cycle_atoms_.end() - 1; it++) ///iterate on cycle atoms except the oxygen of the ring
    {
        int index = distance(mono->cycle_atoms_.begin(), it);
        Atom* target = (*it);
        string key = "";
        string value = "";
        AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it1 = t_neighbors.begin(); it1 != t_neighbors.end(); it1++)
        {
            AtomVector pattern_atoms = AtomVector();
            Atom* t_neighbor = (*it1);
            if(t_neighbor->GetName().at(0) == 'N' && mono->cycle_atoms_str_.find(t_neighbor->GetId()) == string::npos)///check formulas with nitrogen
            {
                if((value = CheckxC_N(target, mono->cycle_atoms_str_, pattern_atoms)).compare("") != 0)///xCH-N
                    break;
                if((value = CheckxC_NxO_CO_C(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-C=OCH3
                    break;
                if((value = CheckxC_NxO_CO_CO(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-C=OCH2OH
                    break;
                if((value = CheckxC_NxO_SO3(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-SO3
                    break;
                if((value = CheckxC_NxO_PO3(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-PO3
                    break;
                if((value = CheckxC_NxO_C(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-CH3
                    break;
            }
            if(t_neighbor->GetName().at(0) == 'O' && mono->cycle_atoms_str_.find(t_neighbor->GetId()) == string::npos)///check formulas with oxygen
            {
                if((value = CheckxC_NxO_CO_C(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-C=OCH3
                    break;
                if((value = CheckxC_NxO_CO_CO(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-C=OCH2OH
                    break;
                if((value = CheckxC_NxO_SO3(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-SO3
                    break;
                if((value = CheckxC_NxO_PO3(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-PO3
                    break;
                if((value = CheckxC_NxO_C(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-CH3
                    break;
                if((value = CheckxCOO(target, mono->cycle_atoms_str_, pattern_atoms)).compare("") != 0)///xC-(O,O) and xC-(O,OH)
                    break;
            }
        }
        if(value.compare("") != 0)///if any pattern matched add it to the index-derivative map
        {
            if(index == 0)
                key = "a";
            else
                key = gmml::ConvertT(index + 1);
            mono->derivatives_map_[key] = value;
        }
    }
    for(vector<AtomVector>::iterator it = mono->side_atoms_.begin(); it != mono->side_atoms_.end(); it++) ///iterate on side atoms
    {
        int index = distance(mono->side_atoms_.begin(), it);
        int side_branch_last_carbon_index = 0;
        AtomVector sides = (*it);
        Atom* target = NULL;
        string key = "";
        string value = "";
        if(it == mono->side_atoms_.begin())///side atoms of anomeric carbon
        {
            if(sides.at(0) != NULL)
                target = sides.at(0);///first index of each side is for carbon atoms in the vector<AtomVector> structure
        }
        if(it == mono->side_atoms_.end() - 1)//side atoms of last carbon of the ring
        {
            if(sides.at(0) != NULL)
            {
                for(side_branch_last_carbon_index = sides.size() - 1; sides.at(side_branch_last_carbon_index) == NULL; side_branch_last_carbon_index-- ){}
                target = sides.at(side_branch_last_carbon_index);
            }
        }
        if(target != NULL)
        {
            AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it1 = t_neighbors.begin(); it1 != t_neighbors.end(); it1++)
            {
                AtomVector pattern_atoms = AtomVector();
                Atom* t_neighbor = (*it1);
                if(t_neighbor->GetName().at(0) == 'N' && mono->cycle_atoms_str_.find(t_neighbor->GetId()) == string::npos)///check formulas with nitrogen
                {
                    if((value = CheckxC_N(target, mono->cycle_atoms_str_, pattern_atoms)).compare("") != 0)///xCH-N
                        break;
                    if((value = CheckxC_NxO_CO_C(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-C=OCH3
                        break;
                    if((value = CheckxC_NxO_CO_CO(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-C=OCH2OH
                        break;
                    if((value = CheckxC_NxO_SO3(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-SO3
                        break;
                    if((value = CheckxC_NxO_PO3(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-PO3
                        break;
                    if((value = CheckxC_NxO_C(target, mono->cycle_atoms_str_, 'N', pattern_atoms)).compare("") != 0)///xC-N-CH3
                        break;
                }
                if(t_neighbor->GetName().at(0) == 'O' && mono->cycle_atoms_str_.find(t_neighbor->GetId()) == string::npos)///check formulas with oxygen
                {
                    if((value = CheckxC_NxO_CO_C(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-C=OCH3
                        break;
                    if((value = CheckxC_NxO_CO_CO(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-C=OCH2OH
                        break;
                    if((value = CheckxC_NxO_SO3(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-SO3
                        break;
                    if((value = CheckxC_NxO_PO3(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-PO3
                        break;
                    if((value = CheckxC_NxO_C(target, mono->cycle_atoms_str_, 'O', pattern_atoms)).compare("") != 0)///xC-O-CH3
                        break;
                    if((value = CheckxCOO(target, mono->cycle_atoms_str_, pattern_atoms)).compare("") != 0)///xC-(O,O) and xC-(O,OH)
                        break;
                }
            }
            if(value.compare("") != 0)///if any pattern matched add it to the index-derivative map
            {
                if(index == 0)
                    key = "-1";
                else
                {
                    switch (side_branch_last_carbon_index)
                    {
                        case 0:
                            key = "+1";
                            break;
                        case 1:
                            key = "+2";
                            break;
                        case 2:
                            key = "+3";
                            break;
                    }

                }
                mono->derivatives_map_[key] = value;
            }
        }
    }
}

string Assembly::CheckxC_N(Atom* target, string cycle_atoms_str, AtomVector& pattern_atoms)
{
    stringstream pattern;
    pattern << "xC";
    Atom* N = NULL;
    AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
    for(AtomVector::iterator it1 = t_neighbors.begin(); it1 != t_neighbors.end(); it1++)
    {
        Atom* t_neighbor = (*it1);
        if(cycle_atoms_str.find(t_neighbor->GetId()) == string::npos)
        {
            if(t_neighbor->GetName().at(0) != 'N')
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == 'N' && N != NULL)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == 'N' && N == NULL)
                N = t_neighbor;
        }
    }
    if(N != NULL)
    {
        pattern << "-N";
        AtomVector n_neighbors = N->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = n_neighbors.begin(); it != n_neighbors.end(); it++)
        {
            Atom* n_neighbor = (*it);
            if(n_neighbor->GetId().compare(target->GetId()) != 0)
                pattern << n_neighbor->GetName().at(0);
        }
    }
    //        cout << "CheckxC_N:" << pattern.str() << endl;
    if(pattern.str().compare("xCH-NHH") == 0 || pattern.str().compare("xC-N") == 0 || pattern.str().compare("xCHH-NHH") == 0 || pattern.str().compare("xCH-N") == 0  ||
            pattern.str().compare("xCHH-N") == 0 || pattern.str().compare("xC-NHH") == 0)
        return "xCH-N";
    else
        return "";
}

string Assembly::CheckxC_NxO_CO_C(Atom *target, string cycle_atoms_str, char NxO, AtomVector& pattern_atoms)
{
    stringstream pattern;
    pattern << "xC";
    Atom* N_or_O = NULL;
    AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
    for(AtomVector::iterator it = t_neighbors.begin(); it != t_neighbors.end(); it++)
    {
        Atom* t_neighbor = (*it);
        if(cycle_atoms_str.find(t_neighbor->GetId()) == string::npos)
        {
            if(t_neighbor->GetName().at(0) != NxO)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O != NULL)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O == NULL)
                N_or_O = t_neighbor;
        }
    }
    if(N_or_O != NULL)
    {
        Atom* C = NULL;
        pattern << "-" << NxO;
        AtomVector n_neighbors = N_or_O->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = n_neighbors.begin(); it != n_neighbors.end(); it++)
        {
            Atom* n_neighbor = (*it);
            if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) != 'C')
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'C' && C != NULL)
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'C' && C == NULL)
                C = n_neighbor;
        }
        if(C != NULL)
        {
            pattern_atoms.push_back(C);
            Atom* CC = NULL;
            Atom* CO = NULL;
            pattern << "-C";
            AtomVector c_neighbors = C->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = c_neighbors.begin(); it != c_neighbors.end(); it++)
            {
                Atom* c_neighbor = (*it);
                if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 && c_neighbor->GetName().at(0) != 'C' && c_neighbor->GetName().at(0) != 'O')
                    pattern << c_neighbor->GetName().at(0);
                else if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  c_neighbor->GetName().at(0) == 'C' && CC != NULL)
                    pattern << c_neighbor->GetName().at(0);
                else if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  c_neighbor->GetName().at(0) == 'O' && CO != NULL)
                    pattern << c_neighbor->GetName().at(0);
                else if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  c_neighbor->GetName().at(0) == 'O' && CO == NULL)
                    CO = c_neighbor;
                else if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  c_neighbor->GetName().at(0) == 'C' && CC == NULL)
                    CC = c_neighbor;
            }
            if(CO != NULL)
            {
                pattern_atoms.push_back(CO);
                pattern << "O";
                AtomVector co_neighbors = CO->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = co_neighbors.begin(); it != co_neighbors.end(); it++)
                {
                    Atom* co_neighbor = (*it);
                    if(co_neighbor->GetId().compare(C->GetId()) != 0)
                        pattern << co_neighbor->GetName().at(0);
                }
            }
            if(CC != NULL)
            {
                pattern_atoms.push_back(CC);
                pattern << "-C";
                AtomVector cc_neighbors = CC->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = cc_neighbors.begin(); it != cc_neighbors.end(); it++)
                {
                    Atom* cc_neighbor = (*it);
                    if(cc_neighbor->GetId().compare(C->GetId()) != 0)
                        pattern << cc_neighbor->GetName().at(0);
                }
            }
        }
    }
    //            cout << "CheckxC_NxO_CO_C:" << pattern.str() << endl;
    if(NxO == 'N')
    {
        if(pattern.str().compare("xCH-NH-CO-CHHH") == 0 || pattern.str().compare("xC-N-CO-C") == 0 || pattern.str().compare("xCHH-NH-CO-CHHH") == 0 || pattern.str().compare("xC-NH-CO-CHHH") == 0 ||
                pattern.str().compare("xC-N-CO-CHHH") == 0 || pattern.str().compare("xC-NH-CO-C") == 0 || pattern.str().compare("xCHH-N-CO-CHHH") == 0 || pattern.str().compare("xCHH-NH-CO-C") == 0 ||
                pattern.str().compare("xCHH-N-CO-C") == 0 || pattern.str().compare("xCH-N-CO-CHHH") == 0 || pattern.str().compare("xCH-NH-CO-C") == 0 || pattern.str().compare("xCH-N-CO-C") == 0 )
            return "xC-N-C=OCH3";
        else
            return "";
    }
    else if(NxO == 'O')
    {
        if(pattern.str().compare("xCH-OH-CO-CHHH") == 0 || pattern.str().compare("xC-O-CO-C") == 0 || pattern.str().compare("xCHH-OH-CO-CHHH") == 0 || pattern.str().compare("xC-OH-CO-CHHH") == 0 ||
                pattern.str().compare("xC-O-CO-CHHH") == 0 || pattern.str().compare("xC-OH-CO-C") == 0 || pattern.str().compare("xCHH-O-CO-CHHH") == 0 || pattern.str().compare("xCHH-OH-CO-C") == 0 ||
                pattern.str().compare("xCHH-O-CO-C") == 0 || pattern.str().compare("xCH-O-CO-CHHH") == 0 || pattern.str().compare("xCH-OH-CO-C") == 0 || pattern.str().compare("xCH-O-CO-C") == 0)
            return "xC-O-C=OCH3";
        else
            return "";
    }
    else
        return "";
}

string Assembly::CheckxC_NxO_CO_CO(Atom *target, string cycle_atoms_str, char NxO, AtomVector& pattern_atoms)
{
    stringstream pattern;
    pattern << "xC";
    Atom* N_or_O = NULL;
    AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
    for(AtomVector::iterator it = t_neighbors.begin(); it != t_neighbors.end(); it++)
    {
        Atom* t_neighbor = (*it);
        if(cycle_atoms_str.find(t_neighbor->GetId()) == string::npos)
        {
            if(t_neighbor->GetName().at(0) != NxO)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O != NULL)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O == NULL)
                N_or_O = t_neighbor;
        }
    }
    if(N_or_O != NULL)
    {
        Atom* C = NULL;
        pattern << "-" << NxO;
        AtomVector n_neighbors = N_or_O->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = n_neighbors.begin(); it != n_neighbors.end(); it++)
        {
            Atom* n_neighbor = (*it);
            if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) != 'C')
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'C' && C != NULL)
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'C' && C == NULL)
                C = n_neighbor;
        }
        if(C != NULL)
        {
            pattern_atoms.push_back(C);
            Atom* CC = NULL;
            Atom* CO = NULL;
            pattern << "-C";
            AtomVector c_neighbors = C->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = c_neighbors.begin(); it != c_neighbors.end(); it++)
            {
                Atom* c_neighbor = (*it);
                if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 && c_neighbor->GetName().at(0) != 'C' && c_neighbor->GetName().at(0) != 'O')
                    pattern << c_neighbor->GetName().at(0);
                else if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  c_neighbor->GetName().at(0) == 'C' && CC != NULL)
                    pattern << c_neighbor->GetName().at(0);
                else if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  c_neighbor->GetName().at(0) == 'O' && CO != NULL)
                    pattern << c_neighbor->GetName().at(0);
                else if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  c_neighbor->GetName().at(0) == 'O' && CO == NULL)
                    CO = c_neighbor;
                else if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  c_neighbor->GetName().at(0) == 'C' && CC == NULL)
                    CC = c_neighbor;
            }
            if(CO != NULL)
            {
                pattern_atoms.push_back(CO);
                pattern << "O";
                AtomVector co_neighbors = CO->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = co_neighbors.begin(); it != co_neighbors.end(); it++)
                {
                    Atom* co_neighbor = (*it);
                    if(co_neighbor->GetId().compare(C->GetId()) != 0)
                        pattern << co_neighbor->GetName().at(0);
                }
            }
            if(CC != NULL)
            {
                pattern_atoms.push_back(CC);
                Atom* O = NULL;
                pattern << "-C";
                AtomVector cc_neighbors = CC->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = cc_neighbors.begin(); it != cc_neighbors.end(); it++)
                {
                    Atom* cc_neighbor = (*it);
                    if(cc_neighbor->GetId().compare(C->GetId()) != 0 && cc_neighbor->GetName().at(0) != 'O')
                        pattern << cc_neighbor->GetName().at(0);
                    else if(cc_neighbor->GetId().compare(C->GetId()) != 0 && cc_neighbor->GetName().at(0) == 'O' && O != NULL)
                        pattern << cc_neighbor->GetName().at(0);
                    else if(cc_neighbor->GetId().compare(C->GetId()) != 0 && cc_neighbor->GetName().at(0) == 'O' && O == NULL)
                        O = cc_neighbor;
                }
                if(O != NULL)
                {
                    pattern_atoms.push_back(O);
                    pattern << "-O";
                    AtomVector o_neighbors = O->GetNode()->GetNodeNeighbors();
                    for(AtomVector::iterator it = o_neighbors.begin(); it != o_neighbors.end(); it++)
                    {
                        Atom* o_neighbor = (*it);
                        if(o_neighbor->GetId().compare(CC->GetId()) != 0)
                            pattern << o_neighbor->GetName().at(0);
                    }
                }
            }
        }
    }
    //        cout << "CheckxC_NxO_CO_CO: " << pattern.str() << endl;
    if(NxO == 'N')
    {
        if(pattern.str().compare("xCH-NH-CO-CHH-OH") == 0 || pattern.str().compare("xCH-N-CO-CHH-OH") == 0 || pattern.str().compare("xCH-NH-CO-C-OH") == 0 ||
                pattern.str().compare("xCH-NH-CO-CHH-O") == 0 || pattern.str().compare("xCH-N-CO-C-OH") == 0 || pattern.str().compare("xCH-N-CO-CHH-O") == 0 ||
                pattern.str().compare("xCH-NH-CO-C-OH") == 0 || pattern.str().compare("xCH-N-CO-C-O") == 0 || pattern.str().compare("xC-N-CO-C-O") == 0 ||
                pattern.str().compare("xC-NH-CO-C-O") == 0 || pattern.str().compare("xC-NH-CO-C-OH") == 0 || pattern.str().compare("xC-NH-CO-CHH-O") == 0 ||
                pattern.str().compare("xC-N-CO-CHH-O") == 0 || pattern.str().compare("xC-N-CO-C-OH") == 0 || pattern.str().compare("xC-NH-CO-CHH-OH") == 0 ||
                pattern.str().compare("xC-N-CO-CHH-OH") == 0 || pattern.str().compare("xCHH-NH-CO-CHH-OH") == 0 || pattern.str().compare("xCHH-N-CO-C-O") == 0 ||
                pattern.str().compare("xCHH-N-CO-CHH-OH") == 0 || pattern.str().compare("xCHH-N-CO-CHH-O") == 0 || pattern.str().compare("xCHH-N-CO-C-OH") == 0 ||
                pattern.str().compare("xCHH-NH-CO-C-OH") == 0 || pattern.str().compare("xCHH-NH-CO-CHH-O") == 0 || pattern.str().compare("xCHH-NH-CO-C-O") == 0 )
            return "xC-N-C=OCH2OH";
        else
            return "";
    }
    else if(NxO == 'O')
    {
        if(pattern.str().compare("xCH-OH-CO-CHH-OH") == 0 || pattern.str().compare("xCH-O-CO-CHH-OH") == 0 || pattern.str().compare("xCH-OH-CO-C-OH") == 0 ||
                pattern.str().compare("xCH-OH-CO-CHH-O") == 0 || pattern.str().compare("xCH-O-CO-C-OH") == 0 || pattern.str().compare("xCH-O-CO-CHH-O") == 0 ||
                pattern.str().compare("xCH-OH-CO-C-OH") == 0 || pattern.str().compare("xCH-O-CO-C-O") == 0 || pattern.str().compare("xC-O-CO-C-O") == 0 ||
                pattern.str().compare("xC-OH-CO-C-O") == 0 || pattern.str().compare("xC-OH-CO-C-OH") == 0 || pattern.str().compare("xC-OH-CO-CHH-O") == 0 ||
                pattern.str().compare("xC-O-CO-CHH-O") == 0 || pattern.str().compare("xC-O-CO-C-OH") == 0 || pattern.str().compare("xC-OH-CO-CHH-OH") == 0 ||
                pattern.str().compare("xC-O-CO-CHH-OH") == 0 || pattern.str().compare("xCHH-OH-CO-CHH-OH") == 0 || pattern.str().compare("xCHH-O-CO-C-O") == 0 ||
                pattern.str().compare("xCHH-O-CO-CHH-OH") == 0 || pattern.str().compare("xCHH-O-CO-CHH-O") == 0 || pattern.str().compare("xCHH-O-CO-C-OH") == 0 ||
                pattern.str().compare("xCHH-OH-CO-C-OH") == 0 || pattern.str().compare("xCHH-OH-CO-CHH-O") == 0 || pattern.str().compare("xCHH-OH-CO-C-O") == 0 )

            return "xC-O-C=OCH2OH";
        else
            return "";
    }
    else
        return "";
}

string Assembly::CheckxC_NxO_SO3(Atom *target, string cycle_atoms_str, char NxO, AtomVector& pattern_atoms)
{
    stringstream pattern;
    pattern << "xC";
    Atom* N_or_O = NULL;
    AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
    for(AtomVector::iterator it = t_neighbors.begin(); it != t_neighbors.end(); it++)
    {
        Atom* t_neighbor = (*it);
        if(cycle_atoms_str.find(t_neighbor->GetId()) == string::npos)
        {
            if(t_neighbor->GetName().at(0) != NxO)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O != NULL)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O == NULL)
                N_or_O = t_neighbor;
        }
    }
    if(N_or_O != NULL)
    {
        Atom* S = NULL;
        pattern << "-" << NxO;
        AtomVector n_neighbors = N_or_O->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = n_neighbors.begin(); it != n_neighbors.end(); it++)
        {
            Atom* n_neighbor = (*it);
            if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) != 'S')
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'S' && S != NULL)
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'S' && S == NULL)
                S = n_neighbor;
        }
        if(S != NULL)
        {
            pattern_atoms.push_back(S);
            Atom* O1 = NULL;
            Atom* O2 = NULL;
            Atom* O3 = NULL;
            pattern << "-S";
            AtomVector s_neighbors = S->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = s_neighbors.begin(); it != s_neighbors.end(); it++)
            {
                Atom* s_neighbor = (*it);
                if(s_neighbor->GetId().compare(N_or_O->GetId()) != 0 && s_neighbor->GetName().at(0) != 'O')
                    pattern << s_neighbor->GetName().at(0);
                else if(s_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  s_neighbor->GetName().at(0) == 'O' && O1 == NULL)
                    O1 = s_neighbor;
                else if(s_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  s_neighbor->GetName().at(0) == 'O' && O2 == NULL)
                    O2 = s_neighbor;
                else if(s_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  s_neighbor->GetName().at(0) == 'O' && O3 == NULL)
                    O3 = s_neighbor;
            }
            if(O1 != NULL && O2 != NULL && O3 != NULL)
            {
                pattern_atoms.push_back(O1);
                pattern_atoms.push_back(O2);
                pattern_atoms.push_back(O3);
                pattern << "OOO";
                AtomVector o1_neighbors = O1->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = o1_neighbors.begin(); it != o1_neighbors.end(); it++)
                {
                    Atom* o1_neighbor = (*it);
                    if(o1_neighbor->GetId().compare(S->GetId()) != 0)
                        pattern << o1_neighbor->GetName().at(0);
                }
                AtomVector o2_neighbors = O2->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = o2_neighbors.begin(); it != o2_neighbors.end(); it++)
                {
                    Atom* o2_neighbor = (*it);
                    if(o2_neighbor->GetId().compare(S->GetId()) != 0)
                        pattern << o2_neighbor->GetName().at(0);
                }
                AtomVector o3_neighbors = O3->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = o3_neighbors.begin(); it != o3_neighbors.end(); it++)
                {
                    Atom* o3_neighbor = (*it);
                    if(o3_neighbor->GetId().compare(S->GetId()) != 0)
                        pattern << o3_neighbor->GetName().at(0);
                }
            }
        }
    }
    //        cout << "CheckxC_NxO_SO3: " << pattern.str() << endl;
    if(NxO == 'N')
    {
        if(pattern.str().compare("xCH-NH-SOOOH") == 0 || pattern.str().compare("xCH-N-SOOOH") == 0 || pattern.str().compare("xCH-NH-SOOO") == 0 || pattern.str().compare("xCH-N-SOOO") == 0 ||
                pattern.str().compare("xCHH-NH-SOOOH") == 0 || pattern.str().compare("xCHH-NH-SOOO") == 0 || pattern.str().compare("xCHH-N-SOOOH") == 0 || pattern.str().compare("xCHH-N-SOOO") == 0 ||
                pattern.str().compare("xC-N-SOOO") == 0 || pattern.str().compare("xC-NH-SOOOH") == 0 || pattern.str().compare("xC-N-SOOOH") == 0 || pattern.str().compare("xC-NH-SOOO") == 0 )
            return "xC-N-SO3";
        else
            return "";
    }
    else if(NxO == 'O')
    {
        if(pattern.str().compare("xCH-OH-SOOOH") == 0 || pattern.str().compare("xCH-O-SOOOH") == 0 || pattern.str().compare("xCH-OH-SOOO") == 0 || pattern.str().compare("xCH-O-SOOO") == 0 ||
                pattern.str().compare("xCHH-OH-SOOOH") == 0 || pattern.str().compare("xCHH-OH-SOOO") == 0 || pattern.str().compare("xCHH-O-SOOOH") == 0 || pattern.str().compare("xCHH-O-SOOO") == 0 ||
                pattern.str().compare("xC-O-SOOO") == 0 || pattern.str().compare("xC-OH-SOOOH") == 0 || pattern.str().compare("xC-O-SOOOH") == 0 || pattern.str().compare("xC-OH-SOOO") == 0 )
            return "xC-O-SO3";
        else
            return "";
    }
    else
        return "";
}

string Assembly::CheckxC_NxO_PO3(Atom *target, string cycle_atoms_str, char NxO, AtomVector& pattern_atoms)
{
    stringstream pattern;
    pattern << "xC";
    Atom* N_or_O = NULL;
    AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
    for(AtomVector::iterator it = t_neighbors.begin(); it != t_neighbors.end(); it++)
    {
        Atom* t_neighbor = (*it);
        if(cycle_atoms_str.find(t_neighbor->GetId()) == string::npos)
        {
            if(t_neighbor->GetName().at(0) != NxO)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O != NULL)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O == NULL)
                N_or_O = t_neighbor;
        }
    }
    if(N_or_O != NULL)
    {
        Atom* P = NULL;
        pattern << "-" << NxO;
        AtomVector n_neighbors = N_or_O->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = n_neighbors.begin(); it != n_neighbors.end(); it++)
        {
            Atom* n_neighbor = (*it);
            if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) != 'P')
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'P' && P != NULL)
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'P' && P == NULL)
                P = n_neighbor;
        }
        if(P != NULL)
        {
            pattern_atoms.push_back(P);
            Atom* O1 = NULL;
            Atom* O2 = NULL;
            Atom* O3 = NULL;
            pattern << "-P";
            AtomVector p_neighbors = P->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = p_neighbors.begin(); it != p_neighbors.end(); it++)
            {
                Atom* p_neighbor = (*it);
                if(p_neighbor->GetId().compare(N_or_O->GetId()) != 0 && p_neighbor->GetName().at(0) != 'O')
                    pattern << p_neighbor->GetName().at(0);
                else if(p_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  p_neighbor->GetName().at(0) == 'O' && O1 == NULL)
                    O1 = p_neighbor;
                else if(p_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  p_neighbor->GetName().at(0) == 'O' && O2 == NULL)
                    O2 = p_neighbor;
                else if(p_neighbor->GetId().compare(N_or_O->GetId()) != 0 &&  p_neighbor->GetName().at(0) == 'O' && O3 == NULL)
                    O3 = p_neighbor;
            }
            if(O1 != NULL && O2 != NULL && O3 != NULL)
            {
                pattern_atoms.push_back(O1);
                pattern_atoms.push_back(O2);
                pattern_atoms.push_back(O3);
                pattern << "OOO";
                AtomVector o1_neighbors = O1->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = o1_neighbors.begin(); it != o1_neighbors.end(); it++)
                {
                    Atom* o1_neighbor = (*it);
                    if(o1_neighbor->GetId().compare(P->GetId()) != 0)
                        pattern << o1_neighbor->GetName().at(0);
                }
                AtomVector o2_neighbors = O2->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = o2_neighbors.begin(); it != o2_neighbors.end(); it++)
                {
                    Atom* o2_neighbor = (*it);
                    if(o2_neighbor->GetId().compare(P->GetId()) != 0)
                        pattern << o2_neighbor->GetName().at(0);
                }
                AtomVector o3_neighbors = O3->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it = o3_neighbors.begin(); it != o3_neighbors.end(); it++)
                {
                    Atom* o3_neighbor = (*it);
                    if(o3_neighbor->GetId().compare(P->GetId()) != 0)
                        pattern << o3_neighbor->GetName().at(0);
                }
            }
        }
    }
    //        cout << "CheckxC_NxO_PO3: " << pattern.str() << endl;
    if(NxO == 'N')
    {
        if(pattern.str().compare("xCH-NH-POOOH") == 0 || pattern.str().compare("xCH-N-POOOH") == 0 || pattern.str().compare("xCH-NH-POOO") == 0 || pattern.str().compare("xCH-N-POOO") == 0 ||
                pattern.str().compare("xCHH-NH-POOOH") == 0 || pattern.str().compare("xCHH-NH-POOO") == 0 || pattern.str().compare("xCHH-N-POOOH") == 0 || pattern.str().compare("xCHH-N-POOO") == 0 ||
                pattern.str().compare("xC-N-POOO") == 0 || pattern.str().compare("xC-NH-POOOH") == 0 || pattern.str().compare("xC-N-POOOH") == 0 || pattern.str().compare("xC-NH-POOO") == 0 )
            return "xC-N-PO3";
        else
            return "";
    }
    else if(NxO = 'O')
    {
        if(pattern.str().compare("xCH-OH-POOOH") == 0 || pattern.str().compare("xCH-O-POOOH") == 0 || pattern.str().compare("xCH-OH-POOO") == 0 || pattern.str().compare("xCH-O-POOO") == 0 ||
                pattern.str().compare("xCHH-OH-POOOH") == 0 || pattern.str().compare("xCHH-OH-POOO") == 0 || pattern.str().compare("xCHH-O-POOOH") == 0 || pattern.str().compare("xCHH-O-POOO") == 0 ||
                pattern.str().compare("xC-O-POOO") == 0 || pattern.str().compare("xC-OH-POOOH") == 0 || pattern.str().compare("xC-O-POOOH") == 0 || pattern.str().compare("xC-OH-POOO") == 0 )            return "xC-O-PO3";
        else
            return "";
    }
    else
        return "";
}

string Assembly::CheckxC_NxO_C(Atom *target, string cycle_atoms_str, char NxO, AtomVector& pattern_atoms)
{
    stringstream pattern;
    pattern << "xC";
    Atom* N_or_O = NULL;
    AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
    for(AtomVector::iterator it = t_neighbors.begin(); it != t_neighbors.end(); it++)
    {
        Atom* t_neighbor = (*it);
        if(cycle_atoms_str.find(t_neighbor->GetId()) == string::npos)
        {
            if(t_neighbor->GetName().at(0) != NxO)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O != NULL)
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == NxO && N_or_O == NULL)
                N_or_O = t_neighbor;
        }
    }
    if(N_or_O != NULL)
    {
        Atom* C = NULL;
        pattern << "-" << NxO;
        AtomVector n_neighbors = N_or_O->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = n_neighbors.begin(); it != n_neighbors.end(); it++)
        {
            Atom* n_neighbor = (*it);
            if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) != 'C')
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'C' && C != NULL)
                pattern << n_neighbor->GetName().at(0);
            else if(n_neighbor->GetId().compare(target->GetId()) != 0 && n_neighbor->GetName().at(0) == 'C' && C == NULL)
                C = n_neighbor;
        }
        if(C != NULL)
        {
            pattern_atoms.push_back(C);
            pattern << "-C";
            AtomVector c_neighbors = C->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = c_neighbors.begin(); it != c_neighbors.end(); it++)
            {
                Atom* c_neighbor = (*it);
                if(c_neighbor->GetId().compare(N_or_O->GetId()) != 0)
                    pattern << c_neighbor->GetName().at(0);
            }
        }
    }
    //            cout << "CheckxC_NxO_C: " << pattern.str() << endl;
    if(NxO == 'N')
    {
        if(pattern.str().compare("xCH-N-CHHH") == 0 ||  pattern.str().compare("xCH-NH-CHHH") == 0 || pattern.str().compare("xCH-NH-C") == 0 || pattern.str().compare("xCH-N-C") == 0 ||
                pattern.str().compare("xCHH-N-CHHH") == 0 || pattern.str().compare("xCHH-NH-CHHH") == 0 || pattern.str().compare("xCHH-NH-C") == 0 || pattern.str().compare("xCHH-N-C") == 0 ||
                pattern.str().compare("xC-NH-CHHH") == 0 || pattern.str().compare("xC-NH-C") == 0 || pattern.str().compare("xC-N-CHHH") == 0 || pattern.str().compare("xC-N-C") == 0)
            return "xC-N-CH3";
        else
            return "";
    }
    else if(NxO == 'O')
    {
        if(pattern.str().compare("xCH-O-CHHH") == 0 ||  pattern.str().compare("xCH-OH-CHHH") == 0 || pattern.str().compare("xCH-OH-C") == 0 || pattern.str().compare("xCH-O-C") == 0 ||
                pattern.str().compare("xCHH-O-CHHH") == 0 || pattern.str().compare("xCHH-OH-CHHH") == 0 || pattern.str().compare("xCHH-OH-C") == 0 || pattern.str().compare("xCHH-O-C") == 0 ||
                pattern.str().compare("xC-OH-CHHH") == 0 || pattern.str().compare("xC-OH-C") == 0 || pattern.str().compare("xC-O-CHHH") == 0 || pattern.str().compare("xC-O-C") == 0)
            return "xC-O-CH3";
        else
            return "";
    }
    else
        return "";
}

string Assembly::CheckxCOO(Atom *target, string cycle_atoms_str, AtomVector& pattern_atoms)
{
    stringstream pattern;
    pattern << "xC";
    Atom* O1 = NULL;
    Atom* O2 = NULL;
    AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
    for(AtomVector::iterator it = t_neighbors.begin(); it != t_neighbors.end(); it++)
    {
        Atom* t_neighbor = (*it);
        if(cycle_atoms_str.find(t_neighbor->GetId()) == string::npos)
        {
            if(t_neighbor->GetName().at(0) != 'O')
                pattern << t_neighbor->GetName().at(0);
            else if(t_neighbor->GetName().at(0) == 'O' && O1 == NULL)
                O1 = t_neighbor;
            else if(t_neighbor->GetName().at(0) == 'O' && O2 == NULL)
                O2 = t_neighbor;
        }
    }
    if(O1 != NULL && O2 != NULL)
    {
        pattern << "OO";
        AtomVector o1_neighbors = O1->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = o1_neighbors.begin(); it != o1_neighbors.end(); it++)
        {
            Atom* o1_neighbor = (*it);
            if(o1_neighbor->GetId().compare(target->GetId()) != 0)
                pattern << o1_neighbor->GetName().at(0);
        }
        AtomVector o2_neighbors = O2->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = o2_neighbors.begin(); it != o2_neighbors.end(); it++)
        {
            Atom* o2_neighbor = (*it);
            if(o2_neighbor->GetId().compare(target->GetId()) != 0)
                pattern << o2_neighbor->GetName().at(0);
        }
    }
    //        cout << "CheckxCOO: " << pattern.str() << endl;
    if(pattern.str().compare("xCOO") == 0 || pattern.str().compare("xCHOO") == 0)
        return "xC-(O,O)";
    else if(pattern.str().compare("xCOOH") == 0 || pattern.str().compare("xCHOOH") == 0)
        return "xC-(O,OH)";
    else
        return "";
}

void Assembly::GenerateCompleteSugarName(Monosaccharide *mono)
{
    stringstream in_bracket;
    stringstream head;
    stringstream tail;
    bool minus_one = false;
    if(mono->derivatives_map_.find("-1") != mono->derivatives_map_.end())
        minus_one = true;
    for(map<string, string>::iterator it1 = mono->derivatives_map_.begin(); it1 != mono->derivatives_map_.end(); it1++)
    {
        string key = (*it1).first;
        string value = (*it1).second;
        string long_name_pattern = "";
        string cond_name_pattern = "";
        string long_name_pattern_at_minus_one = "";
        string long_name_pattern_at_plus_one = "";
        string pattern = "";
        if(value.compare("xCH-N") == 0)
        {
            long_name_pattern = "-osamine";
            cond_name_pattern = "N";
            pattern = "CH-N";
            AddModificationRuleOneInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, tail, minus_one, in_bracket);
        }
        if(value.compare("xC-N-C=OCH3") == 0)
        {
            long_name_pattern = "N-acetyl-";
            cond_name_pattern = "NAc";
            pattern = "C-N-C=OCH3";
            AddModificationRuleOneInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, tail, minus_one, in_bracket);
        }
        if(value.compare("xC-N-C=OCH2OH") == 0)
        {
            long_name_pattern = "N-glycolyl-";
            cond_name_pattern = "NGc";
            pattern = "C-N-C=OCH2OH";
            AddModificationRuleOneInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, tail, minus_one, in_bracket);
        }
        if(value.compare("xC-N-SO3") == 0)
        {
            long_name_pattern = "N-sulfo-";
            cond_name_pattern = "NS";
            pattern = "C-N-SO3";
            AddModificationRuleOneInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, tail, minus_one, in_bracket);
        }
        if(value.compare("xC-N-PO3") == 0)
        {
            long_name_pattern = "N-phospho-";
            cond_name_pattern = "NP";
            pattern = "C-N-PO3";
            AddModificationRuleOneInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, tail, minus_one, in_bracket);
        }
        if(value.compare("xC-N-CH3") == 0)
        {
            long_name_pattern = "N-methyl-";
            cond_name_pattern = "NMe";
            pattern = "C-N-CH3";
            AddModificationRuleOneInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, tail, minus_one, in_bracket);
        }
        if(value.compare("xC-O-C=OCH3") == 0)
        {
            long_name_pattern = "-acetyl-";
            cond_name_pattern = "Ac";
            pattern = "C-O-C=OCH3";
            AddDerivativeRuleInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, minus_one, in_bracket);
        }
        if(value.compare("xC-O-C=OCH2OH") == 0)
        {
            long_name_pattern = "-glycolyl-";
            cond_name_pattern = "Gc";
            pattern = "C-O-C=OCH2OH";
            AddDerivativeRuleInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, minus_one, in_bracket);
        }
        if(value.compare("xC-O-SO3") == 0)
        {
            long_name_pattern = "-sulfo-";
            cond_name_pattern = "S";
            pattern = "C-O-SO3";
            AddDerivativeRuleInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, minus_one, in_bracket);
        }
        if(value.compare("xC-O-PO3") == 0)
        {
            long_name_pattern = "-phospho-";
            cond_name_pattern = "P";
            pattern = "C-O-PO3";
            AddDerivativeRuleInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, minus_one, in_bracket);
        }
        if(value.compare("xC-O-CH3") == 0)
        {
            long_name_pattern = "-methyl-";
            cond_name_pattern = "Me";
            pattern = "C-O-CH3";
            AddDerivativeRuleInfo(key, pattern, mono, long_name_pattern, cond_name_pattern, head, minus_one, in_bracket);
        }
        if(value.compare("xC-(O,OH)") == 0)
        {
            long_name_pattern_at_minus_one = "-ulosonic acid";
            long_name_pattern_at_plus_one = "-uronic acid";
            cond_name_pattern = "AH";
            pattern = "C-(O,OH)";
            AddModificationRuleTwoInfo(key, pattern, mono, long_name_pattern_at_minus_one, long_name_pattern_at_plus_one, cond_name_pattern, tail, minus_one, in_bracket);
        }
        if(value.compare("xC-(O,O)") == 0)
        {
            long_name_pattern_at_minus_one = "-ulosonate";
            long_name_pattern_at_plus_one = "-uronate";
            cond_name_pattern = "A";
            pattern = "C-(O,O)";
            AddModificationRuleTwoInfo(key, pattern, mono, long_name_pattern_at_minus_one, long_name_pattern_at_plus_one, cond_name_pattern, tail, minus_one, in_bracket);
        }
    }
    if(in_bracket.str().size() != 0)
    {
        stringstream short_name;
        string sn;
        if(mono->sugar_name_.monosaccharide_short_name_.compare("") != 0)
            sn = mono->sugar_name_.monosaccharide_short_name_;
        else
            sn = mono->sugar_name_.monosaccharide_stereochemistry_short_name_;

        ///moving a, b or x to after the bracket: short-name + [...] + a/b/x and removing ", " from the end of bracket stream
        int condensed_name_size = sn.size();
        string condensed_name = sn;
        string new_name_part1 = condensed_name.substr(0, (condensed_name_size - 1));///short_name
        char new_name_part2 = condensed_name.at(condensed_name_size - 1);///a/b/x
        short_name << new_name_part1 << "[" << in_bracket.str().substr(0, in_bracket.str().size() - 1) << "]" << new_name_part2;

        mono->sugar_name_.monosaccharide_short_name_ = short_name.str();
    }
    else if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0 && mono->sugar_name_.monosaccharide_short_name_.compare("") == 0)
    {
        mono->sugar_name_.monosaccharide_short_name_ = mono->sugar_name_.monosaccharide_stereochemistry_short_name_;
    }
    stringstream long_name;
    if(mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
    {
        long_name << head.str() << mono->sugar_name_.monosaccharide_stereochemistry_name_ << tail.str();
        mono->sugar_name_.monosaccharide_name_ = long_name.str();
    }
}

void Assembly::AddModificationRuleOneInfo(string key, string pattern, Monosaccharide* mono, string long_name_pattern, string cond_name_pattern, stringstream& head,
                                          stringstream& tail, bool minus_one, stringstream& in_bracket)
{
    stringstream ss;
    ss << pattern;
    if(key.compare("a") == 0)
    {
        ss << " is at warning position: anomeric";
        gmml::log(__LINE__, __FILE__,  gmml::WAR, ss.str());
        Note* der_mod_note = new Note();
        der_mod_note->type_ = Glycan::WARNING;
        der_mod_note->category_ = Glycan::DER_MOD;
        stringstream note;
        note << mono->sugar_name_.monosaccharide_short_name_ << ": " << ss.str();
        der_mod_note->description_ = note.str();
        this->AddNote(der_mod_note);
        cout << ss.str() << endl;
    }
    else if(key.compare("2") == 0 && mono->sugar_name_.ring_type_.compare("P") == 0 &&
            find(mono->chemical_code_->right_down_.begin(), mono->chemical_code_->right_down_.end(), "-1") == mono->chemical_code_->right_down_.end() &&
            find(mono->chemical_code_->right_up_.begin(), mono->chemical_code_->right_up_.end(), "-1") == mono->chemical_code_->right_up_.end() &&
            mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
    {
        if(long_name_pattern.compare("-osamine") == 0)
            tail << long_name_pattern;
        else
            head << long_name_pattern;
        stringstream short_name;
        if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
        {
            ///moving a, b or x to after the N expression: short-name + Condensed name pattern + a/b/x
            int stereo_condensed_name_size = mono->sugar_name_.monosaccharide_stereochemistry_short_name_.size();
            string stereo_condensed_name = mono->sugar_name_.monosaccharide_stereochemistry_short_name_;
            string new_name_part1 = stereo_condensed_name.substr(0, (stereo_condensed_name_size - 1));///short_name
            char new_name_part2 = stereo_condensed_name.at(stereo_condensed_name_size - 1);///a/b/x
            short_name << new_name_part1 << cond_name_pattern << new_name_part2;

            mono->sugar_name_.monosaccharide_short_name_ = short_name.str();
        }
    }
    else if(mono->sugar_name_.ring_type_.compare("F") == 0 && key.compare("4") == 0)
    {
        if(!minus_one)
            ss << " is at error position: 4";
        else
            ss << " is at error position: 5";
        gmml::log(__LINE__, __FILE__,  gmml::ERR, ss.str());
        Note* der_mod_note = new Note();
        der_mod_note->type_ = Glycan::ERROR;
        der_mod_note->category_ = Glycan::DER_MOD;
        stringstream note;
        note << mono->sugar_name_.monosaccharide_short_name_ << ": " << ss.str();
        der_mod_note->description_ = note.str();
        this->AddNote(der_mod_note);
        cout << ss.str() << endl;
    }
    else if(mono->sugar_name_.ring_type_.compare("P") == 0 && key.compare("5") == 0)
    {
        if(!minus_one)
            ss << " is at error position: 5";
        else
            ss << " is at error position: 6";
        gmml::log(__LINE__, __FILE__,  gmml::ERR, ss.str());
        Note* der_mod_note = new Note();
        der_mod_note->type_ = Glycan::ERROR;
        der_mod_note->category_ = Glycan::DER_MOD;
        stringstream note;
        note << mono->sugar_name_.monosaccharide_short_name_ << ": " << ss.str();
        der_mod_note->description_ = note.str();
        this->AddNote(der_mod_note);
        cout << ss.str() << endl;
    }
    else
    {
        if(!minus_one)
        {
            if( key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                in_bracket << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << cond_name_pattern << ",";
            else
                in_bracket << key << cond_name_pattern << ",";
        }
        else
        {
            if(key.compare("-1") == 0)
                in_bracket << "1" << cond_name_pattern << ",";
            else if( key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                in_bracket << mono->cycle_atoms_.size() + ConvertString<int>(key) << cond_name_pattern << ",";
            else
                in_bracket << ConvertString<int>(key) + 1 << cond_name_pattern << ",";
        }
    }
}
void Assembly::AddDerivativeRuleInfo(string key, string pattern, Monosaccharide *mono, string long_name_pattern, string cond_name_pattern, stringstream &head,
                                     bool minus_one, stringstream &in_bracket)
{
    stringstream ss;
    ss << pattern;
    if(mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
    {
        if(!minus_one)
        {
            if(key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                head << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << long_name_pattern;
            //            else if(key.compare("a") == 0)
            //                head << "2" << long_name_pattern;
            else if(key.compare("a") != 0)
                head << ConvertString<int>(key) << long_name_pattern;
        }
        else
        {
            if(key.compare("-1") == 0)
                head << "1" << long_name_pattern;
            else if(key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                head << mono->cycle_atoms_.size() + ConvertString<int>(key) << long_name_pattern;
            //            else if(key.compare("a") == 0)
            //                head << "2" << long_name_pattern;
            else if(key.compare("a") != 0)
                head << ConvertString<int>(key) + 1 << long_name_pattern;
        }
    }
    if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
    {
        if(mono->sugar_name_.ring_type_.compare("F") == 0 && key.compare("4") == 0)
        {
            if(!minus_one)
                ss << " is at error position: 4";
            else
                ss << " is at error position: 5";
            gmml::log(__LINE__, __FILE__,  gmml::ERR, ss.str());
            Note* der_mod_note = new Note();
            der_mod_note->type_ = Glycan::ERROR;
            der_mod_note->category_ = Glycan::DER_MOD;
            stringstream note;
            note << mono->sugar_name_.monosaccharide_short_name_ << ": " << ss.str();
            der_mod_note->description_ = note.str();
            this->AddNote(der_mod_note);
            cout << ss.str() << endl;
        }
        else if(mono->sugar_name_.ring_type_.compare("P") == 0 && key.compare("5") == 0)
        {
            if(!minus_one)
                ss << " is at error position: 5";
            else
                ss << " is at error position: 6";
            gmml::log(__LINE__, __FILE__,  gmml::ERR, ss.str());
            Note* der_mod_note = new Note();
            der_mod_note->type_ = Glycan::ERROR;
            der_mod_note->category_ = Glycan::DER_MOD;
            stringstream note;
            note << mono->sugar_name_.monosaccharide_short_name_ << ": " << ss.str();
            der_mod_note->description_ = note.str();
            this->AddNote(der_mod_note);
            cout << ss.str() << endl;
        }
        else
        {
            if(!minus_one)
            {
                if(key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                    in_bracket << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << cond_name_pattern << ",";
                //                else if(key.compare("a") == 0)
                //                    in_bracket << "2" << cond_name_pattern << ",";
                else if(key.compare("a") != 0)
                    in_bracket << ConvertString<int>(key) << cond_name_pattern << ",";
            }
            else
            {
                if(key.compare("-1") == 0)
                    in_bracket << "1" << cond_name_pattern << ",";
                else if( key.compare("+1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0)
                    in_bracket << mono->cycle_atoms_.size() + ConvertString<int>(key) << cond_name_pattern << ",";
                //                else if(key.compare("a") == 0)
                //                    in_bracket << "2" << cond_name_pattern << ",";
                else if(key.compare("a") != 0)
                    in_bracket << ConvertString<int>(key) + 1 << cond_name_pattern << ",";
            }
        }
    }
}
void Assembly::AddModificationRuleTwoInfo(string key, string pattern, Monosaccharide *mono, string long_name_pattern_at_minus_one, string long_name_pattern_at_plus_one,
                                          string cond_name_pattern, stringstream &tail, bool minus_one, stringstream &in_bracket)
{
    Note* der_mod_note = new Note();
    stringstream ss;
    ss << pattern;
    if((key.compare("-1") == 0 || key.compare("+2") == 0 || key.compare("+3") == 0) && mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
    {
        if(key.compare("-1") == 0)
        {
            ss << " is at warning position: " << ConvertString<int>(key) + 1;
            der_mod_note->type_ = Glycan::WARNING;
        }
        if(!minus_one)
        {
            ss << " is at error position: " << key;
            der_mod_note->type_ = Glycan::ERROR;
        }
        else
        {
            ss << " is at error position: " << ConvertString<int>(key) + 1;
            der_mod_note->type_ = Glycan::ERROR;
        }
        cout << ss.str() << endl;
        gmml::log(__LINE__, __FILE__,  gmml::ERR, ss.str());
        tail << long_name_pattern_at_minus_one;

        if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
        {
            if(!minus_one)
            {
                if( key.compare("+2") == 0 || key.compare("+3") == 0)
                    in_bracket << mono->cycle_atoms_.size() - 1 + ConvertString<int>(key) << cond_name_pattern << ",";
            }
            else
            {
                if(key.compare("-1") == 0)
                    in_bracket << "1" << cond_name_pattern << ",";
                else if( key.compare("+2") == 0 || key.compare("+3") == 0)
                    in_bracket << mono->cycle_atoms_.size() + ConvertString<int>(key) << cond_name_pattern << ",";
            }
        }
    }
    else if(key.compare("+1") == 0 && mono->sugar_name_.monosaccharide_stereochemistry_name_.compare("") != 0)
    {
        tail << long_name_pattern_at_plus_one;
        stringstream short_name;
        if(mono->sugar_name_.monosaccharide_stereochemistry_short_name_.compare("") != 0)
        {
            ///moving a, b or x to after the AH expression: short-name + AH + a/b/x
            int stereo_condensed_name_size = mono->sugar_name_.monosaccharide_stereochemistry_short_name_.size();
            string stereo_condensed_name = mono->sugar_name_.monosaccharide_stereochemistry_short_name_;
            string new_name_part1 = stereo_condensed_name.substr(0, (stereo_condensed_name_size - 1));///short_name
            char new_name_part2 = stereo_condensed_name.at(stereo_condensed_name_size - 1);///a/b/x
            short_name << new_name_part1 << cond_name_pattern << new_name_part2;

            mono->sugar_name_.monosaccharide_short_name_ = short_name.str();
        }
    }
    else
    {
        if(!minus_one)
        {
            if(key.compare("a") != 0)
                ss << " is at warning position: " << key;
            else
                ss << " is at warning position: 1";
        }
        else
        {
            if(key.compare("a") != 0)
                ss << " is at warning position: " << ConvertString<int>(key) + 1;
            else
                ss << " is at warning position: 2";
        }
        der_mod_note->type_ = Glycan::WARNING;
    }
    der_mod_note->category_ = Glycan::DER_MOD;
    stringstream note;
    note << mono->sugar_name_.monosaccharide_short_name_ << ": " << ss.str();
    der_mod_note->description_ = note.str();
    this->AddNote(der_mod_note);
    cout << ss.str() << endl;
    gmml::log(__LINE__, __FILE__,  gmml::WAR, ss.str());
}

void Assembly::UpdateComplexSugarChemicalCode(Monosaccharide *mono)
{
    if(mono->side_atoms_.at(0).at(0) != NULL)
    {
        if(mono->derivatives_map_["-1"].compare("xC-(O,O)") == 0)
        {
            vector<string>::iterator index_it;
            if((index_it = find(mono->chemical_code_->right_up_.begin(), mono->chemical_code_->right_up_.end(), "-1")) != mono->chemical_code_->right_up_.end())
                (*index_it) = "-1A";
            else if((index_it = find(mono->chemical_code_->right_down_.begin(), mono->chemical_code_->right_down_.end(), "-1")) != mono->chemical_code_->right_down_.end())
                (*index_it) = "-1A";
        }
    }
    if(mono->cycle_atoms_.at(2) != NULL)
    {
        if(mono->derivatives_map_["-1"].compare("xC-O-C=OCH3") == 0)
        {
            vector<string>::iterator index_it;
            if((index_it = find(mono->chemical_code_->left_up_.begin(), mono->chemical_code_->left_up_.end(), "3")) != mono->chemical_code_->left_up_.end())
                (*index_it) = "3Ac";
            else if((index_it = find(mono->chemical_code_->left_down_.begin(), mono->chemical_code_->left_down_.end(), "3")) != mono->chemical_code_->left_down_.end())
                (*index_it) = "3Ac";
        }
    }
    if(mono->cycle_atoms_.at(3) != NULL)
    {
        if(mono->derivatives_map_["4"].compare("xCH-N") == 0)
        {
            vector<string>::iterator index_it;
            if((index_it = find(mono->chemical_code_->left_up_.begin(), mono->chemical_code_->left_up_.end(), "4")) != mono->chemical_code_->left_up_.end())
                (*index_it) = "4N";
            else if((index_it = find(mono->chemical_code_->left_down_.begin(), mono->chemical_code_->left_down_.end(), "4")) != mono->chemical_code_->left_down_.end())
                (*index_it) = "4N";
        }
        if(mono->derivatives_map_["4"].compare("xC-N-C=OCH3") == 0)
        {
            vector<string>::iterator index_it;
            if((index_it = find(mono->chemical_code_->left_up_.begin(), mono->chemical_code_->left_up_.end(), "4")) != mono->chemical_code_->left_up_.end())
                (*index_it) = "4NAc";
            else if((index_it = find(mono->chemical_code_->left_down_.begin(), mono->chemical_code_->left_down_.end(), "4")) != mono->chemical_code_->left_down_.end())
                (*index_it) = "4NAc";
        }
        if(mono->derivatives_map_["4"].compare("xC-N-C=OCH2OH") == 0)
        {
            vector<string>::iterator index_it;
            if((index_it = find(mono->chemical_code_->left_up_.begin(), mono->chemical_code_->left_up_.end(), "4")) != mono->chemical_code_->left_up_.end())
                (*index_it) = "4NGc";
            else if((index_it = find(mono->chemical_code_->left_down_.begin(), mono->chemical_code_->left_down_.end(), "4")) != mono->chemical_code_->left_down_.end())
                (*index_it) = "4NGc";
        }
    }
}

vector<Oligosaccharide*> Assembly::ExtractOligosaccharides(vector<Monosaccharide*> monos, ResidueNameMap dataset_residue_names,
                                                           int& number_of_covalent_links, int& number_of_probable_non_covalent_complexes)
{
    string terminal_residue_name = "";
    ResidueNameMap common_terminal_residues = gmml::InitializeCommonTerminalResidueMap();
    map<Monosaccharide*, vector<Monosaccharide*> > monos_table = map<Monosaccharide*, vector<Monosaccharide*> >();
    map<Monosaccharide*, vector<string> > monos_table_linkages = map<Monosaccharide*, vector<string> >();

    ///Iterating on list of monos to check if there is a connection to another mono in the list
    for(vector<Monosaccharide*>::iterator it = monos.begin(); it != monos.end(); it++)
    {
        Monosaccharide* mono1 = (*it);

        monos_table[mono1] = vector<Monosaccharide*>();
        monos_table_linkages[mono1] = vector<string>();

        for(vector<AtomVector>::iterator it1 = mono1->side_atoms_.begin(); it1 != mono1->side_atoms_.end(); it1++) ///iterate on side atoms
        {
            int index = distance(mono1->side_atoms_.begin(), it1);
            AtomVector sides = (*it1);
            map<Atom*, Atom*> target_parent_map = map<Atom*, Atom*>(); /// A map of target atom to it's parent atom. Target atom is a non ring oxygen or nitrogen

            if(it1 == mono1->side_atoms_.begin())///side atoms of anomeric
            {
                if(sides.at(1) != NULL)
                    target_parent_map[sides.at(1)] = mono1->cycle_atoms_.at(0);
            }
            else if(it1 == mono1->side_atoms_.end() - 1) ///side atoms of last carbon of the ring
            {
                for(AtomVector::iterator last_c_side_it = sides.begin(); last_c_side_it != sides.end(); last_c_side_it++)
                {
                    Atom* side_of_last_carbon = (*last_c_side_it);
                    if(side_of_last_carbon != NULL)
                    {
                        AtomVector last_c_side_neighbors = side_of_last_carbon->GetNode()->GetNodeNeighbors();
                        for(AtomVector::iterator it2 = last_c_side_neighbors.begin(); it2 != last_c_side_neighbors.end(); it2++)
                        {
                            if((*it2)->GetId().at(0) == 'O' || (*it2)->GetId().at(0) == 'N')
                            {
                                target_parent_map[(*it2)] = side_of_last_carbon;
                                break;
                            }
                        }
                    }
                }
            }
            else
            {
                if(sides.at(1) != NULL)
                    target_parent_map[sides.at(1)] = mono1->cycle_atoms_.at(index);///index 1 of each side is for non-carbon side atoms in the vector<AtomVector> structure
            }
            ///Examine neighbors of each target atom to check if they can be found in other monos side/ring atoms
            for(map<Atom*, Atom*>::iterator map_it = target_parent_map.begin(); map_it != target_parent_map.end(); map_it++)
            {
                bool found_in_other_mono = false;
                Atom* target = (*map_it).first;
                Atom* target_parent = (*map_it).second;
                AtomVector t_neighbors = target->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it2 = t_neighbors.begin(); it2 != t_neighbors.end(); it2++)
                {
                    Atom* t_neighbor = (*it2);
                    if(t_neighbor->GetId().compare(target_parent->GetId()) != 0)///making sure neighbor is not the parent of target atom
                    {
                        for(vector<Monosaccharide*>::iterator it3 = monos.begin(); it3 != monos.end(); it3++)
                        {
                            if(it3 != it)///Cheking monos other than the current mono
                            {
                                Monosaccharide* mono2 = (*it3);
                                AtomVector mono2_sides = mono2->side_atoms_.at(mono2->side_atoms_.size() - 1); ///side of last ring carbon

                                bool found_in_side = false;
                                for(AtomVector::iterator mono2_last_c_side_it = mono2_sides.begin(); mono2_last_c_side_it != mono2_sides.end(); mono2_last_c_side_it++)
                                {
                                    Atom* mono2_last_c_side = (*mono2_last_c_side_it);
                                    if(mono2_last_c_side != NULL)
                                    {
                                        if(t_neighbor->GetId().compare(mono2_last_c_side->GetId()) == 0) ///target atom has been attached to another cycle's side atom
                                            found_in_side = true;
                                    }
                                }
                                if(found_in_side || mono2->cycle_atoms_str_.find(t_neighbor->GetId()) != string::npos) //if target's neighbor found in another mono's side or ring atoms
                                {
                                    found_in_other_mono = true;
                                    monos_table[mono1].push_back(mono2);

                                    string mono1_carbon = target_parent->GetId();
                                    string mono1_name = "";
                                    string mono2_carbon = t_neighbor->GetId();
                                    string mono2_name = "";
                                    if(mono1->sugar_name_.monosaccharide_short_name_.compare("") != 0)
                                        mono1_name = mono1->sugar_name_.monosaccharide_short_name_;
                                    else
                                        mono1_name = mono1->sugar_name_.monosaccharide_stereochemistry_short_name_;

                                    if(mono2->sugar_name_.monosaccharide_short_name_.compare("") != 0)
                                        mono2_name = mono2->sugar_name_.monosaccharide_short_name_;
                                    else
                                        mono2_name = mono2->sugar_name_.monosaccharide_stereochemistry_short_name_;

                                    stringstream linkage;
                                    linkage << mono1_carbon << "-" << target->GetId() << "-" << mono2_carbon;
                                    monos_table_linkages[mono1].push_back(linkage.str());
                                    break;
                                }
                            }
                        }
                    }
                }
                if(found_in_other_mono)
                    break;
            }
        }
    }
    vector<int> visited_monos = vector<int>();
    vector<Oligosaccharide*> oligosaccharides = vector<Oligosaccharide*>();

    vector<string> checked_linkages = vector<string>();
    for(map<Monosaccharide*, vector<Monosaccharide*> >::iterator it = monos_table.begin(); it != monos_table.end(); it++)
    {
        Monosaccharide* key = (*it).first;
        vector<Monosaccharide*> values = (*it).second;

        vector<string> visited_linkages = vector<string>();
        if(find(visited_monos.begin(), visited_monos.end(), key->mono_id) == visited_monos.end())///if the mono is not visited
        {
            bool isRoot = false;
            stringstream anomeric_linkage;
            anomeric_linkage << key->cycle_atoms_.at(0)->GetId() << "-";
            vector<string> mono_linkages = monos_table_linkages[key];
            AtomVector terminal_atoms = AtomVector();
            if(values.size() == 0) ///mono is not attached to any other mono
            {
                Atom* anomeric_o = NULL;
                if(key->side_atoms_.at(0).at(1) != NULL)
                    anomeric_o = key->side_atoms_.at(0).at(1);
                if(anomeric_o != NULL)
                {
                    if(dataset_residue_names.find(anomeric_o->GetResidue()->GetName()) != dataset_residue_names.end() ||
                            common_terminal_residues.find(anomeric_o->GetResidue()->GetName()) != common_terminal_residues.end())///check if there is any terminal
                    {
                        terminal_residue_name = anomeric_o->GetResidue()->GetName();
                    }
                    else
                    {
                        terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
                    }
                }
                isRoot = true;
            }
            else if (values.size() == 1 && ///mono is attached to one other mono
                     (monos_table_linkages[values.at(0)].size() == 1))///the other mono is only attached to this mono
            {
                ///CHECKING LINKAGE ISSUES, e.g. C1-O3-C4 is an issue
                CheckLinkageNote(key, values.at(0), mono_linkages.at(0), checked_linkages);
                stringstream other_mono_anomeric_linkage_as_right_side;
                other_mono_anomeric_linkage_as_right_side << "-" << values.at(0)->cycle_atoms_.at(0)->GetId();///atom id on the right side of the linkage c-o-c

                Atom* anomeric_o = NULL;
                Atom* o_neighbor_1 = NULL;
                Atom* o_neighbor_2 = NULL;
                AtomVector o_neighbors = AtomVector();
                if(key->side_atoms_.at(0).at(1) != NULL)///Getting the information of anomeric oxygen's neighbors is needed for choosing the root
                {
                    anomeric_o = key->side_atoms_.at(0).at(1);
                    o_neighbors = anomeric_o->GetNode()->GetNodeNeighbors();
                    if(o_neighbors.size() > 1)
                    {
                        o_neighbor_1 = o_neighbors.at(0);
                        o_neighbor_2 = o_neighbors.at(1);
                    }
                }
                if(anomeric_o != NULL)
                {
                    ///RULE1: anomeric to anomeric linkage
                    if(((mono_linkages.at(0)).find(anomeric_linkage.str()) != string::npos) && ///this mono is attached to other mono through anomeric
                            (mono_linkages.at(0).find(other_mono_anomeric_linkage_as_right_side.str()) != string::npos))///the other mono is only attached to this mono through anomeric
                        isRoot = true;
                    ///RULE2: Directed graph
                    else if(((mono_linkages.at(0)).find(anomeric_linkage.str()) == string::npos) && ///this mono is not attached to other mono through anomeric
                            (mono_linkages.at(0).find(other_mono_anomeric_linkage_as_right_side.str()) != string::npos)) ///the other mono is attached to this mono through anomeric
                    {
                        isRoot = true;
                        terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
                    }
                    ///RULE3: Terminal
                    else if(o_neighbors.size() == 1) ///anomeric oxygen is not attached to anything else except the carbon of the ring
                    {
                        isRoot = true;
                        terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
                    }
                    else if(o_neighbors.size() == 2 && (((o_neighbor_1->GetDescription().find("Het;") != string::npos) && (o_neighbor_2->GetDescription().find("Het;") == string::npos)) ||
                                                        ((o_neighbor_2->GetDescription().find("Het;") != string::npos) && (o_neighbor_1->GetDescription().find("Het;") == string::npos))) )
                    {
                        ///anomeric oxygen is attached to protein
                        isRoot = true;
                        terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
                        number_of_covalent_links++;
                        if(terminal_residue_name.compare("NLN") != 0 && terminal_residue_name.compare("OLS") != 0 && terminal_residue_name.compare("OLT") != 0)
                        {
                            stringstream ss;
                            ss << "Root anomeric atom is attached to a non-standard " << terminal_residue_name << " protein residue!";
                            gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
                            cout << ss.str() << endl;
                            terminal_residue_name = "";
                        }
                    }
                    else if(dataset_residue_names.find(anomeric_o->GetResidue()->GetName()) != dataset_residue_names.end() ||
                            common_terminal_residues.find(anomeric_o->GetResidue()->GetName()) != common_terminal_residues.end())///mono is attached to a terminal through anomeric oxygen
                    {
                        terminal_residue_name = anomeric_o->GetResidue()->GetName();
                        isRoot = true;
                    }
                    else if((terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms)).compare("") != 0)
                        isRoot = true;
                }
                ///RULE2: Directed graph
                else if((mono_linkages.at(0).find(other_mono_anomeric_linkage_as_right_side.str()) != string::npos)) ///this mono doesn't have anomeric oxygen and the other mono is attached to this mono through anomeric
                {
                    isRoot = true;
                    terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
                    number_of_probable_non_covalent_complexes++;
                }
            }
            else
            {
                Atom* anomeric_o = NULL;
                Atom* o_neighbor_1 = NULL;
                Atom* o_neighbor_2 = NULL;
                AtomVector o_neighbors = AtomVector();
                if(key->side_atoms_.at(0).at(1) != NULL)///Getting the information of anomeric oxygen's neighbors is needed for choosing the root
                {
                    anomeric_o = key->side_atoms_.at(0).at(1);
                    o_neighbors = anomeric_o->GetNode()->GetNodeNeighbors();
                    if(o_neighbors.size() > 1)
                    {
                        o_neighbor_1 = o_neighbors.at(0);
                        o_neighbor_2 = o_neighbors.at(1);
                    }
                }
                if(anomeric_o != NULL)
                {
                    ///RULE1: anomeric to anomeric linkage
                    for(int i = 0; i < values.size(); i++)
                    {
                        CheckLinkageNote(key, values.at(i), mono_linkages.at(i), checked_linkages);
                        stringstream other_mono_anomeric_linkage_as_right_side;
                        other_mono_anomeric_linkage_as_right_side << "-" << values.at(i)->cycle_atoms_.at(0)->GetId();///atom id on the right side of the linkage c-o-c

                        if(((mono_linkages.at(i)).find(anomeric_linkage.str()) != string::npos) && ///this mono is attached to another mono through anomeric
                                (mono_linkages.at(i).find(other_mono_anomeric_linkage_as_right_side.str()) != string::npos))///the other mono is attached to this mono through anomeric
                        {
                            isRoot = true;
                            break;
                        }
                    }
                    if(!isRoot) ///RULE2: Directed graph
                    {
                        for(int i = 0; i < values.size(); i++)
                        {
                            CheckLinkageNote(key, values.at(i), mono_linkages.at(i), checked_linkages);
                            stringstream other_mono_anomeric_linkage_as_right_side;
                            other_mono_anomeric_linkage_as_right_side << "-" << values.at(i)->cycle_atoms_.at(0)->GetId();///atom id on the right side of the linkage c-o-c
                            if(((mono_linkages.at(i)).find(anomeric_linkage.str()) != string::npos)) ///this mono is attached to other mono through anomeric
                            {
                                isRoot = false;
                                break;
                            }
                            else if((mono_linkages.at(i).find(other_mono_anomeric_linkage_as_right_side.str()) != string::npos)) ///the other mono is attached to this mono through anomeric
                            {
                                isRoot = true;
                                terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
                            }
                        }
                    }
                    else if(!isRoot)///RULE3: Terminal
                    {
                        if(o_neighbors.size() == 1) ///anomeric oxygen is not attached to anything else, except the carbon of the ring
                        {
                            isRoot = true;
                            terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
                        }
                        else if(o_neighbors.size() == 2 && (((o_neighbor_1->GetDescription().find("Het;") != string::npos) && (o_neighbor_2->GetDescription().find("Het;") == string::npos)) ||
                                                            ((o_neighbor_2->GetDescription().find("Het;") != string::npos) && (o_neighbor_1->GetDescription().find("Het;") == string::npos))) )
                        {
                            ///anomeric oxygen is attached to protein
                            isRoot = true;
                            terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
                            number_of_covalent_links++;
                            if(terminal_residue_name.compare("NLN") != 0 && terminal_residue_name.compare("OLS") != 0 && terminal_residue_name.compare("OLT") != 0)
                            {
                                stringstream ss;
                                ss << "Root anomeric atom is attached to a non-standard " << terminal_residue_name << " protein residue!";
                                gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
                                cout << ss.str() << endl;
                                terminal_residue_name = "";
                            }
                        }
                        else if(dataset_residue_names.find(anomeric_o->GetResidue()->GetName()) != dataset_residue_names.end() ||
                                common_terminal_residues.find(anomeric_o->GetResidue()->GetName()) != common_terminal_residues.end())///mono is attached to a terminal through anomeric oxygen
                        {
                            terminal_residue_name = anomeric_o->GetResidue()->GetName();
                            isRoot = true;
                        }
                        else if((terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms)).compare("") != 0)
                            isRoot = true;
                    }
                }
                //this mono doesn't have anomeric oxygen
                else ///RULE2: Directed graph
                {
                    for(int i = 0; i < values.size(); i++)
                    {
                        CheckLinkageNote(key, values.at(i), mono_linkages.at(i), checked_linkages);
                        vector<string> other_mono_linkage = monos_table_linkages[values.at(i)];
                        stringstream other_mono_anomeric_linkage;
                        other_mono_anomeric_linkage << values.at(i)->cycle_atoms_.at(0)->GetId() << "-";///atom id on the left side of the linkage c-o-c
                        if((other_mono_linkage.at(0).find(other_mono_anomeric_linkage.str()) != string::npos)) ///the other mono is attached to this mono through anomeric
                        {
                            isRoot = true;
                            terminal_residue_name = CheckTerminals(anomeric_o, terminal_atoms);
                            number_of_probable_non_covalent_complexes++;
                            break;
                        }
                    }
                }
            }
            if(isRoot)
            {
                Oligosaccharide* oligo = new Oligosaccharide();
                BuildOligosaccharideTreeStructure(key, values, oligo, visited_monos, monos_table, monos_table_linkages, visited_linkages);
                oligo->terminal_ = terminal_residue_name;
                oligosaccharides.push_back(oligo);
            }
        }
    }

    for(map<Monosaccharide*, vector<Monosaccharide*> >::iterator it = monos_table.begin(); it != monos_table.end(); it++)
    {
        Monosaccharide* key = (*it).first;
        vector<Monosaccharide*> values = (*it).second;
        if(values.size() > 1)
        {
            vector<string> visited_linkages = vector<string>();
            if(find(visited_monos.begin(), visited_monos.end(), key->mono_id) == visited_monos.end())///if the mono is not visited
            {
                vector<string> mono_linkages = monos_table_linkages[key];
                stringstream anomeric_linkage;
                anomeric_linkage << key->cycle_atoms_.at(0)->GetId() << "-";

                for(vector<string>::iterator it1 = mono_linkages.begin(); it1 != mono_linkages.end(); it1++)
                {
                    if((*it1).find(anomeric_linkage.str()) != string::npos)///mono is attached to another mono through anomeric
                    {
                        Oligosaccharide* oligo = new Oligosaccharide();
                        BuildOligosaccharideTreeStructure(key, values, oligo, visited_monos, monos_table, monos_table_linkages, visited_linkages);
                        oligosaccharides.push_back(oligo);
                        break;
                    }
                }
            }
        }
    }
    return oligosaccharides;
}

string Assembly::CheckOMETerminal(Atom* target, AtomVector& terminal_atoms)
{
    terminal_atoms = AtomVector();
    AtomVector atoms_1 = AtomVector();
    AtomVector atoms_2 = AtomVector();
    stringstream pattern;
    pattern << "O";
    atoms_1.push_back(target);
    atoms_2.push_back(target);
    Atom* C1 = NULL;
    Atom* C2 = NULL;
    AtomVector o_neighbors = target->GetNode()->GetNodeNeighbors();
    for(AtomVector::iterator it = o_neighbors.begin(); it != o_neighbors.end(); it++)
    {
        Atom* o_neighbor = (*it);

        if(o_neighbor->GetName().at(0) == 'C' && C1 == NULL && C2 == NULL)
            C1 = o_neighbor;
        else if(o_neighbor->GetName().at(0) == 'C' && C1 != NULL && C2 == NULL)
            C2 = o_neighbor;
    }
    stringstream temp;
    if(C1 != NULL)
    {
        temp << pattern.str() << "-C";
        atoms_1.push_back(C1);
        AtomVector c1_neighbors = C1->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = c1_neighbors.begin(); it != c1_neighbors.end(); it++)
        {
            Atom* c1_neighbor = (*it);
            if(c1_neighbor->GetId().compare(target->GetId()) != 0)
            {
                temp << c1_neighbor->GetName().at(0);
                atoms_1.push_back(c1_neighbor);
            }
        }
    }
    if(temp.str().compare("O-C") == 0 || temp.str().compare("O-CHHH") == 0)
    {
        terminal_atoms = atoms_1;
        return "OME";
    }
    if(C2 != NULL)
    {
        pattern << "-C";
        atoms_2.push_back(C2);
        AtomVector c2_neighbors = C2->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = c2_neighbors.begin(); it != c2_neighbors.end(); it++)
        {
            Atom* c2_neighbor = (*it);
            if(c2_neighbor->GetId().compare(target->GetId()) != 0)
            {
                pattern << c2_neighbor->GetName().at(0);
                atoms_2.push_back(c2_neighbor);
            }
        }
        if(pattern.str().compare("O-C") == 0 || pattern.str().compare("O-CHHH") == 0)
        {
            terminal_atoms = atoms_2;
            return "OME";
        }
    }
    return "";
}

string Assembly::CheckROHTerminal(Atom* target, AtomVector& terminal_atoms)
{
    terminal_atoms = AtomVector();
    AtomVector o_neighbors = target->GetNode()->GetNodeNeighbors();
    if(o_neighbors.size() == 1)
    {
        terminal_atoms.push_back(target);
        return "ROH";
    }
    else if (o_neighbors.size() > 1)
    {
        if((o_neighbors.at(0)->GetName().at(0) == 'H' && o_neighbors.at(1)->GetName().at(0) != 'H'))
        {
            terminal_atoms.push_back(target);
            terminal_atoms.push_back(o_neighbors.at(0));
            return "ROH";
        }
        else if(o_neighbors.at(1)->GetName().at(0) == 'H' && o_neighbors.at(0)->GetName().at(0) != 'H')
        {
            terminal_atoms.push_back(target);
            terminal_atoms.push_back(o_neighbors.at(1));
            return "ROH";
        }
    }
    return "";
}

string Assembly::CheckTBTTerminal(Atom *target, AtomVector& terminal_atoms)
{
    terminal_atoms = AtomVector();
    AtomVector atoms_1 = AtomVector();
    AtomVector atoms_2 = AtomVector();
    stringstream pattern;
    pattern << "O";
    atoms_1.push_back(target);
    atoms_2.push_back(target);
    Atom* C1 = NULL;
    Atom* C2 = NULL;
    AtomVector o_neighbors = target->GetNode()->GetNodeNeighbors();
    for(AtomVector::iterator it = o_neighbors.begin(); it != o_neighbors.end(); it++)
    {
        Atom* o_neighbor = (*it);
        if(o_neighbor->GetName().at(0) == 'C' && C1 == NULL && C2 == NULL)
            C1 = o_neighbor;
        else if(o_neighbor->GetName().at(0) == 'C' && C1 != NULL && C2 == NULL)
            C2 = o_neighbor;
    }
    stringstream temp;
    if(C1 != NULL)
    {
        Atom* C1C1 = NULL;
        Atom* C1C2 = NULL;
        Atom* C1C3 = NULL;
        temp << pattern.str() << "-" << "C";
        atoms_1.push_back(C1);
        AtomVector c1_neighbors = C1->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = c1_neighbors.begin(); it != c1_neighbors.end(); it++)
        {
            Atom* c1_neighbor = (*it);
            if(c1_neighbor->GetId().compare(target->GetId()) != 0 && c1_neighbor->GetName().at(0) == 'C' && C1C1 == NULL)
                C1C1 = c1_neighbor;
            else if(c1_neighbor->GetId().compare(target->GetId()) != 0 && c1_neighbor->GetName().at(0) == 'C' && C1C1 != NULL && C1C2 == NULL)
                C1C2 = c1_neighbor;
            else if(c1_neighbor->GetId().compare(target->GetId()) != 0 && c1_neighbor->GetName().at(0) == 'C' && C1C1 != NULL && C1C2 != NULL && C1C3 == NULL)
                C1C3 = c1_neighbor;
        }
        if(C1C1 != NULL && C1C2 != NULL && C1C3 != NULL)
        {
            temp << "C";
            atoms_1.push_back(C1C1);
            AtomVector c1c1_neighbors = C1C1->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = c1c1_neighbors.begin(); it != c1c1_neighbors.end(); it++)
            {
                Atom* c1c1_neighbor = (*it);
                if(c1c1_neighbor->GetId().compare(C1->GetId()) != 0)
                {
                    temp << c1c1_neighbor->GetName().at(0);
                    atoms_1.push_back(c1c1_neighbor);
                }
            }
            temp << "C";
            atoms_1.push_back(C1C2);
            AtomVector c1c2_neighbors = C1C2->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = c1c2_neighbors.begin(); it != c1c2_neighbors.end(); it++)
            {
                Atom* c1c2_neighbor = (*it);
                if(c1c2_neighbor->GetId().compare(C1->GetId()) != 0)
                {
                    temp << c1c2_neighbor->GetName().at(0);
                    atoms_1.push_back(c1c2_neighbor);
                }
            }
            temp << "C";
            atoms_1.push_back(C1C3);
            AtomVector c1c3_neighbors = C1C3->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = c1c3_neighbors.begin(); it != c1c3_neighbors.end(); it++)
            {
                Atom* c1c3_neighbor = (*it);
                if(c1c3_neighbor->GetId().compare(C1->GetId()) != 0)
                {
                    temp << c1c3_neighbor->GetName().at(0);
                    atoms_1.push_back(c1c3_neighbor);
                }
            }
        }
        if(temp.str().compare("O-CCCC") == 0 || temp.str().compare("O-CCHHHCHHHCHHH") == 0 )
        {
            terminal_atoms = atoms_1;
            return "TBT";
        }
    }
    if(C2 != NULL)
    {
        Atom* C2C1 = NULL;
        Atom* C2C2 = NULL;
        Atom* C2C3 = NULL;
        pattern << "-" << "C";
        atoms_2.push_back(C2);
        AtomVector c2_neighbors = C2->GetNode()->GetNodeNeighbors();
        for(AtomVector::iterator it = c2_neighbors.begin(); it != c2_neighbors.end(); it++)
        {
            Atom* c2_neighbor = (*it);
            if(c2_neighbor->GetId().compare(target->GetId()) != 0 && c2_neighbor->GetName().at(0) == 'C' && C2C1 == NULL)
                C2C1 = c2_neighbor;
            else if(c2_neighbor->GetId().compare(target->GetId()) != 0 && c2_neighbor->GetName().at(0) == 'C' && C2C1 != NULL && C2C2 == NULL)
                C2C2 = c2_neighbor;
            else if(c2_neighbor->GetId().compare(target->GetId()) != 0 && c2_neighbor->GetName().at(0) == 'C' && C2C1 != NULL && C2C2 != NULL && C2C3 == NULL)
                C2C3 = c2_neighbor;
        }
        if(C2C1 != NULL && C2C2 != NULL && C2C3 != NULL)
        {
            pattern << "C";
            atoms_2.push_back(C2C1);
            AtomVector c2c1_neighbors = C2C1->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = c2c1_neighbors.begin(); it != c2c1_neighbors.end(); it++)
            {
                Atom* c2c1_neighbor = (*it);
                if(c2c1_neighbor->GetId().compare(C2->GetId()) != 0)
                {
                    pattern << c2c1_neighbor->GetName().at(0);
                    atoms_2.push_back(c2c1_neighbor);
                }
            }
            pattern << "C";
            atoms_2.push_back(C2C2);
            AtomVector c2c2_neighbors = C2C2->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = c2c2_neighbors.begin(); it != c2c2_neighbors.end(); it++)
            {
                Atom* c2c2_neighbor = (*it);
                if(c2c2_neighbor->GetId().compare(C2->GetId()) != 0)
                {
                    pattern << c2c2_neighbor->GetName().at(0);
                    atoms_2.push_back(c2c2_neighbor);
                }
            }
            pattern << "C";
            atoms_2.push_back(C2C3);
            AtomVector c2c3_neighbors = C2C3->GetNode()->GetNodeNeighbors();
            for(AtomVector::iterator it = c2c3_neighbors.begin(); it != c2c3_neighbors.end(); it++)
            {
                Atom* c2c3_neighbor = (*it);
                if(c2c3_neighbor->GetId().compare(C2->GetId()) != 0)
                {
                    pattern << c2c3_neighbor->GetName().at(0);
                    atoms_2.push_back(c2c3_neighbor);
                }
            }
        }
        if(pattern.str().compare("O-CCCC") == 0 || pattern.str().compare("O-CCHHHCHHHCHHH") == 0 )
        {
            terminal_atoms = atoms_2;
            return "TBT";
        }
    }
    return "";
}

string Assembly::CheckTerminals(Atom* target, AtomVector& terminal_atoms)
{
    if(target != NULL)
    {
        AtomVector o_neighbors = target->GetNode()->GetNodeNeighbors();
        if(CheckROHTerminal(target, terminal_atoms).compare("") != 0)
            return "ROH";
        else if(CheckOMETerminal(target, terminal_atoms).compare("") != 0)
            return "OME";
        else if(CheckTBTTerminal(target, terminal_atoms).compare("") != 0)
            return "TBT";
        else if(o_neighbors.size() == 2)
        {
            Atom* target_o_neighbor = NULL;
            if(o_neighbors.at(0)->GetDescription().find("Het;") != string::npos && o_neighbors.at(1)->GetDescription().find("Het;") == string::npos)
                target_o_neighbor = o_neighbors.at(1);
            else if(o_neighbors.at(0)->GetDescription().find("Het;") == string::npos && o_neighbors.at(1)->GetDescription().find("Het;") != string::npos)
                target_o_neighbor = o_neighbors.at(0);

            if(target_o_neighbor != NULL)
            {
                ResidueVector residues = this->GetAllResiduesOfAssembly();
                Residue* target_residue = NULL;
                for(ResidueVector::iterator it = residues.begin(); it != residues.end(); it++)
                {
                    Residue* residue = *it;
                    if(residue->GetId().compare(target_o_neighbor->GetResidue()->GetId()) == 0)
                    {
                        target_residue = residue;
                        break;
                    }
                }
                if(target_residue != NULL)
                    terminal_atoms = target_residue->GetAtoms();


                AmberGlycamMap amber_glycam = AmberGlycamLookup(target_o_neighbor->GetResidue()->GetName());
                AmberGlycamMap glycam_amber = GlycamAmberLookup(target_o_neighbor->GetResidue()->GetName());

                if(amber_glycam.amber_name_.compare("") != 0)
                    return amber_glycam.amber_name_;
                else if(glycam_amber.glycam_name_.compare("") != 0)
                    return glycam_amber.amber_name_;
                else
                    return target_o_neighbor->GetResidue()->GetName();
            }
            else
                return "";
        }
        else
            return "";
    }
    else
        return "";
}

void Assembly::CheckLinkageNote(Monosaccharide* mono1, Monosaccharide* mono2, string linkage, vector<string>& checked_linkages)
{
    if(find(checked_linkages.begin(), checked_linkages.end(), linkage) == checked_linkages.end())///If this linkage hasn't been checked before by calling the function on other side of the linkage
    {
        vector<string> linkage_tokens = Split(linkage, "-");
        stringstream reverse_linkage;
        reverse_linkage << linkage_tokens.at(2) << "-" << linkage_tokens.at(1) << "-" << linkage_tokens.at(0);
        checked_linkages.push_back(linkage);
        checked_linkages.push_back(reverse_linkage.str());

        int left_c_index = ConvertString<int>(Split(Split(linkage_tokens.at(0), "_").at(0), "C*,\'").at(0));
        int right_c_index = ConvertString<int>(Split(Split(linkage_tokens.at(2), "_").at(0), "C*,\'").at(0));
        int glycosidic_o_index = ConvertString<int>(Split(Split(linkage_tokens.at(1), "_").at(0), "ON*,\'").at(0));
        if(left_c_index != glycosidic_o_index && right_c_index != glycosidic_o_index)
        {
            Note* linkage_note = new Note();
            linkage_note->type_ = Glycan::ERROR;
            linkage_note->category_ = Glycan::GLYCOSIDIC;
            stringstream n;
            n << mono1->sugar_name_.monosaccharide_short_name_ << ": Glycosidic oxygen/nitrogen index does not conform to carbon index in the linkage to "
              << mono2->sugar_name_.monosaccharide_short_name_ << ". " << Split(linkage_tokens.at(0), "_").at(0) << "-" << Split(linkage_tokens.at(1), "_").at(0)
              << "-" << Split(linkage_tokens.at(2), "_").at(0);
            linkage_note->description_ = n.str();
            this->AddNote(linkage_note);
        }
    }
}

void Assembly::BuildOligosaccharideTreeStructure(Monosaccharide *key, vector<Monosaccharide*> values, Oligosaccharide *oligo,
                                                 vector<int>& visited_monos, map<Monosaccharide*, vector<Monosaccharide*> > monos_table,
                                                 map<Monosaccharide*, vector<string> > monos_table_linkages, vector<string>& visited_linkages)
{
    oligo->root_ = key;
    if(values.size() == 0)
    {
        oligo->child_oligos_ = vector<Oligosaccharide*>();
        oligo->child_oligos_linkages_ = vector<string>();
        visited_monos.push_back(key->mono_id);
        return;
    }
    else
    {
        for(vector<Monosaccharide*>::iterator it = values.begin(); it != values.end(); it++)
        {
            Monosaccharide* value_mono = (*it);
            if(find(visited_monos.begin(), visited_monos.end(), value_mono->mono_id) == visited_monos.end())
            {
                int it_index = distance(values.begin(), it);
                vector<string> key_mono_linkages = monos_table_linkages[key];
                string link = key_mono_linkages.at(it_index);
                stringstream reverse_link;
                reverse_link << Split(link, "-").at(2) << "-" << Split(link, "-").at(1) << "-" << Split(link, "-").at(0);
                if(find(visited_linkages.begin(), visited_linkages.end(), link) == visited_linkages.end() &&
                        find(visited_linkages.begin(), visited_linkages.end(), reverse_link.str()) == visited_linkages.end())
                {
                    //                    cout << "key id " << key->mono_id  << ", value id " << value_mono->mono_id << endl;
                    Oligosaccharide* child_oligo = new Oligosaccharide();
                    vector<Monosaccharide*> value_mono_values = monos_table[value_mono];
                    visited_linkages.push_back(link);
                    //                    cout << "call " << value_mono->mono_id << endl;
                    BuildOligosaccharideTreeStructure(value_mono, value_mono_values, child_oligo, visited_monos, monos_table, monos_table_linkages, visited_linkages);
                    oligo->child_oligos_.push_back(child_oligo);
                    oligo->child_oligos_linkages_.push_back(link);
                }
            }
        }
        visited_monos.push_back(key->mono_id);
        return;
    }
}

