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
void Assembly::BuildStructure(gmml::BuildingStructureOption building_option, vector<string> options, vector<string> file_paths)
{
    switch(building_option)
    {
        case gmml::DISTANCE:
            if(options.size() == 1)
            {
                stringstream ss(description_);
                ss << "Building option: Distance;";
                description_ = ss.str();
                vector<string> tokens = gmml::Split(options.at(0), ":");
                if(tokens.at(0).compare("cutoff") == 0)
                {
                    double cutoff = gmml::ConvertString<double>(tokens.at(1));
                    this->BuildStructureByDistance(cutoff);
                }
                if(tokens.at(0).compare("model_index") == 0)
                {
                    double cutoff = gmml::dCutOff;
                    int model_index = gmml::ConvertString<int>(tokens.at(1));
                    this->BuildStructureByDistance(cutoff, model_index);
                }
            }
            if(options.size() == 2)
            {
                stringstream ss(description_);
                ss << "Building option: Distance;";
                description_ = ss.str();
                vector<string> cutoff_tokens = gmml::Split(options.at(0), ":");
                vector<string> model_tokens = gmml::Split(options.at(1), ":");
                double cutoff = gmml::dCutOff;
                int model_index = 0;
                if(cutoff_tokens.at(0).compare("cutoff") == 0)
                {
                    cutoff = gmml::ConvertString<double>(cutoff_tokens.at(1));
                }
                if(model_tokens.at(0).compare("model_index") == 0)
                {
                    model_index = gmml::ConvertString<int>(model_tokens.at(1));
                }
                BuildStructureByDistance(cutoff, model_index);
            }
            else
            {
                BuildStructureByDistance();
            }
            break;
        case gmml::ORIGINAL:
        {
            stringstream ss(description_);
            ss << "Building option: Original;";
            ss << "File type: " << gmml::ConvertAssemblyInputFileType2String(this->GetSourceFileType()) << ";"
               << "File path: " << this->GetSourceFile() << ";";
            description_ = ss.str();
            this->BuildStructureByOriginalFileBondingInformation();
            break;
        }
        case gmml::DATABASE:
            vector<gmml::InputFileType> types = vector<gmml::InputFileType>();
            for(unsigned int i = 0; i < options.size(); i++)
            {
                vector<string> tokens = gmml::Split(options.at(i), ":");
                if(tokens.at(0).compare("type") == 0)
                {
                    gmml::InputFileType type = gmml::ConvertString2AssemblyInputFileType(tokens.at(1));
                    types.push_back(type);
                }
            }
            if(types.size() == file_paths.size())
            {
                stringstream ss(description_);
                ss << "Building option: Distance;";
                for(unsigned int i = 0; i < types.size(); i++)
                {
                    ss << "File type: " << types.at(i) << ";" << "File path: " << file_paths.at(i) << ";";
                }
                description_ = ss.str();
                this->BuildStructureByDatabaseFilesBondingInformation(types, file_paths);
            }
            break;
    }
}

pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex2 = PTHREAD_MUTEX_INITIALIZER;
///7*7 half matrix 6^2/2 + 4, 2 types of thread complete tile, half tile _ 6*6 complete tiles, 22 threads
///or first chunk small next one bigger and so on

void* BuildStructureByDistanceThread(void* args){
    DistanceCalculationThreadArgument* arg = (DistanceCalculationThreadArgument*)args;
    double cutoff = arg->cutoff;
    int model_index = arg->model_index;
    int ti = arg->thread_index;
    int t = arg->number_of_threads;

    //    cout << "Thread" << ti << " start" << endl;
    Assembly::AtomVector all_atoms_of_assembly = arg->a->GetAllAtomsOfAssembly();
    int atoms_size = all_atoms_of_assembly.size();
    int i = ti * (atoms_size/t);

    for(Assembly::AtomVector::iterator it = all_atoms_of_assembly.begin() + ti * (atoms_size/t); ; it++)
    {
        if(ti + 1 < t)
        {
            if(it == all_atoms_of_assembly.begin() + (ti+1) * (atoms_size/t))
                break;
        }
        else
        {
            if(it == all_atoms_of_assembly.end() - 1)
            {
                if((*it)->GetNode() == NULL)
                {
                    Atom* atom = (*it);
                    AtomNode* atom_node = new AtomNode();
                    atom_node->SetAtom(atom);
                    atom->SetNode(atom_node);
                }
                break;
            }
        }
        Atom* atom = (*it);
        AtomNode* atom_node;
        pthread_mutex_lock(&mutex1);
        if(atom->GetNode() == NULL)
        {
            atom_node = new AtomNode();
            atom_node->SetAtom(atom);
            atom->SetNode(atom_node);
        }
        else
            atom_node = atom->GetNode();
        atom_node->SetId(i);
        //        cout << "Thread" << ti << " atom id " << i << endl;
        i++;
        pthread_mutex_unlock(&mutex1);
        for(Assembly::AtomVector::iterator it1 = it + 1; it1 != all_atoms_of_assembly.end(); it1++)
        {
            Atom* neighbor_atom = (*it1);
            // X distance
            if(atom->GetCoordinates().at(model_index)->GetX() - neighbor_atom->GetCoordinates().at(model_index)->GetX() < cutoff)
            {
                // Y distance
                if(atom->GetCoordinates().at(model_index)->GetY() - neighbor_atom->GetCoordinates().at(model_index)->GetY() < cutoff)
                {
                    // Z distance
                    if(atom->GetCoordinates().at(model_index)->GetZ() - neighbor_atom->GetCoordinates().at(model_index)->GetZ() < cutoff)
                    {
                        if((atom->GetCoordinates().at(model_index)->Distance(*(neighbor_atom->GetCoordinates().at(model_index)))) < cutoff)
                        {
                            AtomNode* neighbor_node;
                            pthread_mutex_lock(&mutex1);
                            if (neighbor_atom->GetNode() == NULL)
                            {
                                neighbor_node = new AtomNode();
                                neighbor_node->SetAtom(neighbor_atom);
                            }
                            else
                                neighbor_node = neighbor_atom->GetNode();
                            atom_node->AddNodeNeighbor(neighbor_atom);
                            neighbor_node->AddNodeNeighbor(atom);
                            neighbor_atom->SetNode(neighbor_node);
                            pthread_mutex_unlock(&mutex1);
                        }
                    }
                }
            }
        }
        atom->SetNode(atom_node);
    }

    //    cout << "Thread" << ti << " END" << endl;
    pthread_exit((void*) ti);
}

void* BuildStructureByDistanceByOptimizedThread(void* args){

    DistanceCalculationThreadArgument* arg = (DistanceCalculationThreadArgument*)args;
    double cutoff = arg->cutoff;
    int model_index = arg->model_index;
    int ti = arg->thread_index;
    int t = arg->number_of_threads;


    //    cout << "Thread" << ti << " start" << endl;
    Assembly::AtomVector all_atoms_of_assembly = arg->a->GetAllAtomsOfAssembly();
    int atoms_size = all_atoms_of_assembly.size();

    int increase_factor = 0;
    increase_factor = atoms_size / ((t * t+1) / 2);
    int end_index = 0;
    int begin_index = 0;
    int i = 1;
    for(i = 1; i < ti+1; i++)
    {
        begin_index = begin_index + (i*increase_factor);
    }
    end_index = begin_index + (i*increase_factor);

    for(Assembly::AtomVector::iterator it = all_atoms_of_assembly.begin() + begin_index; ; it++)
    {
        if(ti + 1 < t)///if it's not the last thread
        {
            if(it == all_atoms_of_assembly.begin() + end_index)
                break;
        }
        else
        {
            if(it == all_atoms_of_assembly.end() - 1)
            {
                if((*it)->GetNode() == NULL)
                {
                    Atom* atom = (*it);
                    AtomNode* atom_node = new AtomNode();
                    atom_node->SetAtom(atom);
                    atom->SetNode(atom_node);
                }
                break;
            }
        }
        int index = distance(all_atoms_of_assembly.begin(), it);
        Atom* atom = (*it);
        AtomNode* atom_node;
        pthread_mutex_lock(&mutex1);
        if(atom->GetNode() == NULL)
        {
            atom_node = new AtomNode();
            atom_node->SetAtom(atom);
            atom->SetNode(atom_node);
        }
        else
            atom_node = atom->GetNode();
        atom_node->SetId(index);
        //        cout << "Thread" << ti << " atom id " << i << endl;
        pthread_mutex_unlock(&mutex1);
        for(Assembly::AtomVector::iterator it1 = it + 1; it1 != all_atoms_of_assembly.end(); it1++)
        {
            Atom* neighbor_atom = (*it1);
            // X distance
            if(atom->GetCoordinates().at(model_index)->GetX() - neighbor_atom->GetCoordinates().at(model_index)->GetX() < cutoff)
            {
                // Y distance
                if(atom->GetCoordinates().at(model_index)->GetY() - neighbor_atom->GetCoordinates().at(model_index)->GetY() < cutoff)
                {
                    // Z distance
                    if(atom->GetCoordinates().at(model_index)->GetZ() - neighbor_atom->GetCoordinates().at(model_index)->GetZ() < cutoff)
                    {
                        if((atom->GetCoordinates().at(model_index)->Distance(*(neighbor_atom->GetCoordinates().at(model_index)))) < cutoff)
                        {
                            AtomNode* neighbor_node;
                            pthread_mutex_lock(&mutex1);
                            if (neighbor_atom->GetNode() == NULL)
                            {
                                neighbor_node = new AtomNode();
                                neighbor_node->SetAtom(neighbor_atom);
                            }
                            else
                                neighbor_node = neighbor_atom->GetNode();
                            atom_node->AddNodeNeighbor(neighbor_atom);
                            neighbor_node->AddNodeNeighbor(atom);
                            neighbor_atom->SetNode(neighbor_node);
                            pthread_mutex_unlock(&mutex1);
                        }
                    }
                }
            }
        }
        atom->SetNode(atom_node);
    }

    //    cout << "Thread" << ti << " END" << endl;
    pthread_exit((void*) ti);
}

void* BuildStructureByDistanceByMatrixThread(void* args){

    DistanceCalculationByMatrixThreadArgument* arg = (DistanceCalculationByMatrixThreadArgument*)args;
    int ti = arg->thread_index;
    double cutoff = arg->cutoff;
    int model_index = arg->model_index;
    Assembly::AtomVector* first_chunk = arg->first_chunk;
    Assembly::AtomVector* second_chunk = arg->second_chunk;

    for(Assembly::AtomVector::iterator it = first_chunk->begin(); it != first_chunk->end(); it++)
    {
        int index = distance(first_chunk->begin(), it);
        Atom* atom = (*it);
        AtomNode* atom_node;
        pthread_mutex_lock(&mutex1);
        if(atom->GetNode() == NULL)
        {
            atom_node = new AtomNode();
            atom_node->SetAtom(atom);
            atom->SetNode(atom_node);
        }
        else
            atom_node = atom->GetNode();
        atom_node->SetId(index);
        pthread_mutex_unlock(&mutex1);
        for(Assembly::AtomVector::iterator it1 = second_chunk->begin(); it1 != second_chunk->end(); it1++)
        {
            Atom* neighbor_atom = (*it1);
            // X distance
            if(atom->GetCoordinates().at(model_index)->GetX() - neighbor_atom->GetCoordinates().at(model_index)->GetX() < cutoff)
            {
                // Y distance
                if(atom->GetCoordinates().at(model_index)->GetY() - neighbor_atom->GetCoordinates().at(model_index)->GetY() < cutoff)
                {
                    // Z distance
                    if(atom->GetCoordinates().at(model_index)->GetZ() - neighbor_atom->GetCoordinates().at(model_index)->GetZ() < cutoff)
                    {
                        if((atom->GetCoordinates().at(model_index)->Distance(*(neighbor_atom->GetCoordinates().at(model_index)))) < cutoff)
                        {
                            AtomNode* neighbor_node;
                            pthread_mutex_lock(&mutex1);
                            if (neighbor_atom->GetNode() == NULL)
                            {
                                neighbor_node = new AtomNode();
                                neighbor_node->SetAtom(neighbor_atom);
                            }
                            else
                                neighbor_node = neighbor_atom->GetNode();
                            atom_node->AddNodeNeighbor(neighbor_atom);
                            neighbor_node->AddNodeNeighbor(atom);
                            neighbor_atom->SetNode(neighbor_node);
                            pthread_mutex_unlock(&mutex1);
                        }
                    }
                }
            }
        }
        atom->SetNode(atom_node);
    }
    pthread_exit((void*) ti);
}

void* BuildStructureByDistanceByMatrixDiameterThread(void* args){

    DistanceCalculationByMatrixThreadArgument* arg = (DistanceCalculationByMatrixThreadArgument*)args;
    int ti = arg->thread_index;
    double cutoff = arg->cutoff;
    int model_index = arg->model_index;
    Assembly::AtomVector* first_chunk = arg->first_chunk;
    Assembly::AtomVector* second_chunk = arg->second_chunk;
    int j = ti * 10;

    vector<Assembly::AtomVector*> chunks = vector<Assembly::AtomVector*>();
    chunks.push_back(first_chunk);
    chunks.push_back(second_chunk);
    for(int i = 0; i < 2; i++)
    {
        Assembly::AtomVector* chunk = chunks.at(i);
        if(chunk->size() != 0)
        {
            for(Assembly::AtomVector::iterator it = chunk->begin(); it != chunk->end(); it++)
            {
                Atom* atom = (*it);
                //                    cout << "chunk" << i << " start " << atom->GetId() << endl;
                AtomNode* atom_node;
                pthread_mutex_lock(&mutex1);
                if(atom->GetNode() == NULL)
                {
                    atom_node = new AtomNode();
                    atom_node->SetAtom(atom);
                    atom->SetNode(atom_node);
                }
                else
                    atom_node = atom->GetNode();
                atom_node->SetId(j);
                j++;
                if(it != chunk->end())
                {
                    pthread_mutex_unlock(&mutex1);
                    for(Assembly::AtomVector::iterator it1 = it+1; it1 != chunk->end(); it1++)
                    {
                        Atom* neighbor_atom = (*it1);
                        // X distance
                        if(atom->GetCoordinates().at(model_index)->GetX() - neighbor_atom->GetCoordinates().at(model_index)->GetX() < cutoff)
                        {
                            // Y distance
                            if(atom->GetCoordinates().at(model_index)->GetY() - neighbor_atom->GetCoordinates().at(model_index)->GetY() < cutoff)
                            {
                                // Z distance
                                if(atom->GetCoordinates().at(model_index)->GetZ() - neighbor_atom->GetCoordinates().at(model_index)->GetZ() < cutoff)
                                {
                                    if((atom->GetCoordinates().at(model_index)->Distance(*(neighbor_atom->GetCoordinates().at(model_index)))) < cutoff)
                                    {
                                        AtomNode* neighbor_node;
                                        pthread_mutex_lock(&mutex1);
                                        if (neighbor_atom->GetNode() == NULL)
                                        {
                                            neighbor_node = new AtomNode();
                                            neighbor_node->SetAtom(neighbor_atom);
                                        }
                                        else
                                            neighbor_node = neighbor_atom->GetNode();
                                        atom_node->AddNodeNeighbor(neighbor_atom);
                                        neighbor_node->AddNodeNeighbor(atom);
                                        neighbor_atom->SetNode(neighbor_node);
                                        pthread_mutex_unlock(&mutex1);
                                    }
                                }
                            }
                        }
                    }
                }
                atom->SetNode(atom_node);
            }
        }
    }
    pthread_exit((void*) ti);
}

// MATRIX VERSION
/*
void Assembly::BuildStructureByDistance(int number_of_threads, double cutoff, int model_index)
{
    cout << "Building structure by distance ..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building structure by distance ...");
    model_index_ = model_index;

    Assembly::AtomVector all_atoms_of_assembly = this->GetAllAtomsOfAssembly();
    int atoms_size = all_atoms_of_assembly.size();

    number_of_threads = 5;
    int matrix_size = 3;
//    number_of_threads = 2;
//    int matrix_size = 2;
    pthread_t threads[number_of_threads];
    DistanceCalculationByMatrixThreadArgument arg[number_of_threads];
    int k = 0;
    for(int i = 0; i < matrix_size; i++)
    {
        AtomVector* first_chunk = new AtomVector();
        for(AtomVector::iterator it = all_atoms_of_assembly.begin() + i * (atoms_size/matrix_size); ; it++)
        {
            if(i + 1 < matrix_size)
            {
                if(it == all_atoms_of_assembly.begin() + (i+1) * (atoms_size/matrix_size))
                    break;
            }
            else
            {
                if(it == all_atoms_of_assembly.end())
                    break;
            }
            Atom* atom = (*it);
            first_chunk->push_back(atom);
        }
        for(int j = i+1; j < matrix_size; j++)
        {
            AtomVector* second_chunk = new AtomVector();
            for(AtomVector::iterator it = all_atoms_of_assembly.begin() + (j * (atoms_size/matrix_size)); ; it++)
            {
                if(j + 1 < matrix_size)
                {
                    if(it == all_atoms_of_assembly.begin() + ((j+1)*(atoms_size/matrix_size)) )
                        break;
                }
                else
                {
                    if(it == all_atoms_of_assembly.end())
                        break;
                }
                Atom* atom = (*it);
                second_chunk->push_back(atom);
            }
//            cout << "thread=" << k << ", first chunk size = " << first_chunk->size() << ", starting from " << i * (atoms_size/matrix_size);
//            cout << ", 2nd chunk size = " << second_chunk->size() << ", starting from " << (j * ((atoms_size/matrix_size) )) << endl;
            arg[k] = DistanceCalculationByMatrixThreadArgument(k, model_index, cutoff, first_chunk, second_chunk);
            pthread_create(&threads[k], NULL, &BuildStructureByDistanceByMatrixThread, &arg[k]);
            k++;
        }
    }

    bool no_atoms_for_second_chunk = false;
    for(int i = 0; i < matrix_size; i++)///DIAMETER threads
    {
        AtomVector* first_chunk = new AtomVector();
        AtomVector* second_chunk = new AtomVector();
        for(AtomVector::iterator it = all_atoms_of_assembly.begin() + (i * (atoms_size/matrix_size)); ; it++)
        {
            if(i + 1 < 4)
            {
                if(it == all_atoms_of_assembly.begin() + ((i+1) * (atoms_size/matrix_size)))
                    break;
            }
            else
            {
                if(it == all_atoms_of_assembly.end())
                {
                    no_atoms_for_second_chunk = true;
                    break;
                }
            }
            Atom* atom = (*it);
            first_chunk->push_back(atom);
        }
//        cout << "thread=" << k << ", first chunk size = " << first_chunk->size() << ", starting from " << i * (atoms_size/matrix_size);
        i++;
        if(!no_atoms_for_second_chunk)
        {
            for(AtomVector::iterator it = all_atoms_of_assembly.begin() + (i * (atoms_size/matrix_size)); ; it++)
            {
                if(i + 1 < matrix_size)
                {
                    if(it == all_atoms_of_assembly.begin() + (i+1) * (atoms_size/matrix_size))
                        break;
                }
                else
                    break;

                Atom* atom = (*it);
                second_chunk->push_back(atom);
            }
        }
        //cout << ", 2nd chunk size = " << second_chunk->size() << ", starting from " << (i * (atoms_size/matrix_size)) << endl;
        arg[k] = DistanceCalculationByMatrixThreadArgument(k, model_index, cutoff, first_chunk, second_chunk);
            pthread_create(&threads[k], NULL, &BuildStructureByDistanceByMatrixDiameterThread, &arg[k]);
        k++;
    }
    for(int i = 0; i < number_of_threads; i++)
    {
        pthread_join(threads[i], NULL);
    }
}
*/

///First and second version
void Assembly::BuildStructureByDistance(int number_of_threads, double cutoff, int model_index)
{
    cout << "Building structure by distance ..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building structure by distance ...");
    model_index_ = model_index;

    pthread_t threads[number_of_threads];
    DistanceCalculationThreadArgument arg[number_of_threads];
    for(int i = 0; i < number_of_threads; i++)
    {
        arg[i] = DistanceCalculationThreadArgument(i, number_of_threads, model_index, cutoff, this);
        //        pthread_create(&threads[i], NULL, &BuildStructureByDistanceThread, &arg[i]); ///First version. Workload of threads are not equal
        pthread_create(&threads[i], NULL, &BuildStructureByDistanceByOptimizedThread, &arg[i]); ///Second version. Workload of threads are roughly equal.
    }
    for(int i = 0; i < number_of_threads; i++)
    {
        pthread_join(threads[i], NULL);
    }
}

void Assembly::BuildStructureByOriginalFileBondingInformation()
{
    gmml::InputFileType type = this->GetSourceFileType();
    switch(type)
    {
        case gmml::PDB:
            this->BuildStructureByPDBFileInformation();
            break;
        case gmml::PDBQT:
            break;
        case gmml::TOP:
            this->BuildStructureByTOPFileInformation();
            break;
        case gmml::LIB:
            this->BuildStructureByLIBFileInformation();
            break;
        case gmml::PREP:
            this->BuildStructureByPrepFileInformation();
            break;
        case gmml::TOP_CRD:
            break;
        case gmml::MULTIPLE:
            break;
        case gmml::UNKNOWN:
            break;
    }
}

void Assembly::BuildStructureByPDBFileInformation()
{
    try{
        cout << "Building structure by pdb file information ..." << endl;
        gmml::log(__LINE__, __FILE__, gmml::INF, "Building structure by pdb file information ...");
        PdbFile* pdb_file = new PdbFile(this->GetSourceFile());
        AtomVector all_atoms_of_assembly = this->GetAllAtomsOfAssembly();
        int i = 0;
        for(AtomVector::iterator it = all_atoms_of_assembly.begin(); it != all_atoms_of_assembly.end(); it++)
        {
            Atom* atom = (*it);
            AtomNode* atom_node = new AtomNode();
            atom_node->SetAtom(atom);
            atom_node->SetId(i);
            i++;
            PdbAtom* pdb_atom = pdb_file->GetAtomOfResidueByAtomKey(atom->GetId());
            if(pdb_atom != NULL)
            {
                int atom_serial_number = pdb_atom->GetAtomSerialNumber();
                PdbConnectCard* connectivities = pdb_file->GetConnectivities();
                PdbConnectCard::BondedAtomsSerialNumbersMap bonded_atoms_map = connectivities->GetBondedAtomsSerialNumbers();
                vector<int> bonded_atoms_serial_number = bonded_atoms_map[atom_serial_number];
                for(vector<int>::iterator it1 = bonded_atoms_serial_number.begin(); it1 != bonded_atoms_serial_number.end(); it1++)
                {
                    int bonded_atom_serial_number = *it1;
                    PdbAtom* pdb_bonded_atom = pdb_file->GetAtomBySerialNumber(bonded_atom_serial_number);
                    stringstream sss;
                    sss << pdb_bonded_atom->GetAtomName() << "_" << pdb_bonded_atom->GetAtomSerialNumber() << "_" << pdb_bonded_atom->GetAtomResidueName()
                        << "_" << pdb_bonded_atom->GetAtomChainId() << "_" << pdb_bonded_atom->GetAtomResidueSequenceNumber()
                        << "_" << pdb_bonded_atom->GetAtomInsertionCode() << "_" << pdb_bonded_atom->GetAtomAlternateLocation();
                    string pdb_bonded_atom_key = sss.str();
                    for(AtomVector::iterator it2 = all_atoms_of_assembly.begin(); it2 != all_atoms_of_assembly.end(); it2++)
                    {
                        if(it != it2)
                        {
                            Atom* assembly_atom = (*it2);
                            string assembly_atom_key = assembly_atom->GetId();
                            if(assembly_atom_key.compare(pdb_bonded_atom_key) == 0)
                            {
                                atom_node->AddNodeNeighbor(assembly_atom);
                                break;
                            }
                        }
                    }
                }
            }
            atom->SetNode(atom_node);
        }
    }
    catch(PdbFileSpace::PdbFileProcessingException &ex)
    {}
}

void Assembly::BuildStructureByTOPFileInformation()
{
    cout << "Building structure by topology file information ..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building structure by topology file information ...");
    TopologyFile* topology_file = new TopologyFile(gmml::Split(this->GetSourceFile(), ";")[0]);
    AtomVector all_atoms_of_assembly = this->GetAllAtomsOfAssembly();
    int i = 0;
    for (AtomVector::iterator it = all_atoms_of_assembly.begin(); it != all_atoms_of_assembly.end(); it++)
    {
        Atom* atom_1 = (*it);
        AtomNode* atom_node = new AtomNode();
        atom_node->SetAtom(atom_1);
        atom_node->SetId(i);
        i++;
        for(AtomVector::iterator it1 = all_atoms_of_assembly.begin(); it1 != all_atoms_of_assembly.end(); it1++)
        {
            if(it != it1)
            {
                Atom* atom_2 = (*it1);
                stringstream ss;
                ss << gmml::Split(atom_1->GetId(), "_").at(2) << "(" << gmml::Split(atom_1->GetId(), "_").at(4) << ")"
                   << ":" << gmml::Split(atom_1->GetId(), "_").at(0) << "(" <<  gmml::Split(atom_1->GetId(), "_").at(1) << ")" << "-"
                   << gmml::Split(atom_2->GetId(), "_").at(2) << "(" << gmml::Split(atom_2->GetId(), "_").at(4) << ")"
                   << ":" << gmml::Split(atom_2->GetId(), "_").at(0) << "(" << gmml::Split(atom_2->GetId(), "_").at(1) << ")";
                string key = ss.str();
                TopologyFile::TopologyBondMap topology_bond = topology_file->GetBonds();
                for(TopologyFile::TopologyBondMap::iterator it2 = topology_bond.begin(); it2 != topology_bond.end(); it2++)
                {
                    TopologyBond* bond = (*it2).second;
                    stringstream sss;
                    sss << bond->GetResidueNames().at(0) << ":" << bond->GetBonds().at(0) << "-" << bond->GetResidueNames().at(1) << ":" << bond->GetBonds().at(1);
                    string topology_bond_key = sss.str();
                    if(key.compare(topology_bond_key) == 0)
                    {
                        atom_node->AddNodeNeighbor(atom_2);
                        break;
                    }
                    stringstream ssss;
                    ssss << bond->GetResidueNames().at(1) << ":" << bond->GetBonds().at(1) << "-" << bond->GetResidueNames().at(0) << ":" << bond->GetBonds().at(0);
                    topology_bond_key = ssss.str();
                    if(key.compare(topology_bond_key) == 0)
                    {
                        atom_node->AddNodeNeighbor(atom_2);
                        break;
                    }
                }
            }
        }
        atom_1->SetNode(atom_node);
    }
}

void Assembly::BuildStructureByLIBFileInformation()
{
    cout << "Building structure by library file information..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building structure by library file information ...");
    LibraryFile* library_file = new LibraryFile(this->GetSourceFile());
    AtomVector all_atoms_of_assembly = this->GetAllAtomsOfAssembly();
    int i = 0;
    for(AtomVector::iterator it = all_atoms_of_assembly.begin(); it != all_atoms_of_assembly.end(); it++)
    {
        Atom* atom = (*it);
        AtomNode* atom_node = new AtomNode();
        atom_node->SetAtom(atom);
        atom_node->SetId(i);
        i++;
        Residue* assembly_residue = atom->GetResidue();
        vector<string> atom_id_tokens = Split(atom->GetId(), "_");
        LibraryFileResidue* library_residue = library_file->GetLibraryResidueByResidueName(assembly_residue->GetName());
        if(library_residue != NULL)
        {
            LibraryFileAtom* library_atom = library_residue->GetAtomByOrder(ConvertString<int>(atom_id_tokens.at(1)));
            if(library_atom != NULL)
            {
                vector<int> library_bonded_atom_indices = library_atom->GetBondedAtomsIndices();
                for(vector<int>::iterator it1 = library_bonded_atom_indices.begin(); it1 != library_bonded_atom_indices.end(); it1++)
                {
                    int library_bonded_atom_index = (*it1);
                    LibraryFileAtom* library_atom = library_residue->GetAtomByOrder(library_bonded_atom_index);
                    for(AtomVector::iterator it2 = all_atoms_of_assembly.begin(); it2 != all_atoms_of_assembly.end(); it2++)
                    {
                        Atom* assembly_atom = (*it2);
                        string assembly_atom_id = assembly_atom->GetId();
                        stringstream ss;
                        ss << library_residue->GetName() << ":" << library_atom->GetName() << "(" << library_atom->GetAtomOrder() << ")";
                        string library_atom_id = ss.str();
                        vector<string> assembly_atom_id_tokens = gmml::Split(assembly_atom_id, "_");
                        stringstream sss;
                        sss << assembly_atom_id_tokens.at(2) << ":" << assembly_atom_id_tokens.at(0) << "(" << assembly_atom_id_tokens.at(1) << ")";
                        if(sss.str().compare(library_atom_id) == 0)
                        {
                            atom_node->AddNodeNeighbor(assembly_atom);
                            break;
                        }
                    }
                }
            }
        }
        atom->SetNode(atom_node);
    }
}

void Assembly::BuildStructureByPrepFileInformation()
{
    cout << "Building structure by prep file information ..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building structure by prep file information ...");
    PrepFile* prep_file = new PrepFile(this->GetSourceFile());
    AtomVector all_atoms_of_assembly = this->GetAllAtomsOfAssembly();
    int i = 0;
    for(AtomVector::iterator it = all_atoms_of_assembly.begin(); it != all_atoms_of_assembly.end(); it++)
    {
        Atom* atom = (*it);
        AtomNode* atom_node = new AtomNode();
        atom_node->SetAtom(atom);
        atom_node->SetId(i);
        i++;
        Residue* assembly_residue = atom->GetResidue();
        //        cout << assembly_residue->GetName() << endl;
        PrepFileResidue* prep_residue = prep_file->GetResidues()[assembly_residue->GetName()];
        if(prep_residue != NULL)
        {
            PrepFileAtom* prep_atom = prep_residue->GetPrepAtomByName(atom->GetName());
            if(prep_atom != NULL)
            {
                vector<int> bonded_atoms_index = prep_residue->GetBondingsOfResidue()[prep_residue->GetAtomIndexByName(atom->GetName())];
                for(vector<int>::iterator it1 = bonded_atoms_index.begin(); it1 != bonded_atoms_index.end(); it1++)
                {
                    int bonded_atom_index = (*it1);
                    PrepFileAtom* bonded_atom = prep_residue->GetPrepAtomByName(prep_residue->GetAtomNameByIndex(bonded_atom_index));
                    stringstream ss;
                    ss << prep_residue->GetName() << ":" << bonded_atom->GetName();
                    for(AtomVector::iterator it2 = all_atoms_of_assembly.begin(); it2 != all_atoms_of_assembly.end(); it2++)
                    {
                        Atom* assembly_atom = (*it2);
                        vector<string> atom_id_tokens = gmml::Split(assembly_atom->GetId(), "_");
                        stringstream sss;
                        sss << atom_id_tokens.at(2) << ":" << atom_id_tokens.at(0);
                        if(sss.str().compare(ss.str()) == 0)
                        {
                            atom_node->AddNodeNeighbor(assembly_atom);
                            break;
                        }
                    }
                }
            }
        }
        atom->SetNode(atom_node);
    }
}

void Assembly::BuildStructureByDatabaseFilesBondingInformation(vector<gmml::InputFileType> types, vector<string> file_paths)
{
    cout << "Building structure by dataset files information ..." << endl;
    gmml::log(__LINE__, __FILE__, gmml::INF, "Building structure by dataset files information ...");
    AtomVector all_atoms_of_assembly = this->GetAllAtomsOfAssembly();
    int i = 0;
    for(AtomVector::iterator it = all_atoms_of_assembly.begin(); it != all_atoms_of_assembly.end(); it++)
    {
        Atom* atom = (*it);
        AtomNode* atom_node = new AtomNode();
        atom_node->SetAtom(atom);
        atom_node->SetId(i);
        i++;
        Residue* assembly_residue = atom->GetResidue();
        for(unsigned int i = 0; i < types.size(); i++)
        {
            if(types.at(i) == gmml::LIB)
            {
                string lib_path = file_paths.at(i);
                LibraryFile* library_file = new LibraryFile(lib_path);
                LibraryFileResidue* library_residue = library_file->GetLibraryResidueByResidueName(assembly_residue->GetName());
                if(library_residue != NULL)
                {
                    LibraryFileAtom* library_atom = library_residue->GetLibraryAtomByAtomName(atom->GetName());
                    if(library_atom != NULL)
                    {
                        vector<int> library_bonded_atom_indices = library_atom->GetBondedAtomsIndices();
                        for(vector<int>::iterator it1 = library_bonded_atom_indices.begin(); it1 != library_bonded_atom_indices.end(); it1++)
                        {
                            int library_bonded_atom_index = (*it1);
                            LibraryFileAtom* library_atom = library_residue->GetAtomByIndex(library_bonded_atom_index);
                            for(AtomVector::iterator it2 = all_atoms_of_assembly.begin(); it2 != all_atoms_of_assembly.end(); it2++)
                            {
                                Atom* assembly_atom = (*it2);
                                string assembly_atom_id = assembly_atom->GetId();
                                stringstream ss;
                                ss << library_residue->GetName() << ":" << library_atom->GetName();
                                string library_atom_id = ss.str();
                                if(assembly_atom_id.compare(library_atom_id) == 0)
                                {
                                    atom_node->AddNodeNeighbor(assembly_atom);
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            if(types.at(i) == gmml::PREP)
            {
                string prep_path = file_paths.at(i);
                PrepFile* prep_file = new PrepFile(prep_path);
                PrepFileResidue* prep_residue = prep_file->GetResidues()[assembly_residue->GetName()];
                if(prep_residue != NULL)
                {
                    PrepFileAtom* prep_atom = prep_residue->GetPrepAtomByName(atom->GetName());
                    if(prep_atom != NULL)
                    {
                        vector<int> bonded_atoms_index = prep_residue->GetBondingsOfResidue()[prep_residue->GetAtomIndexByName(atom->GetName())];
                        for(vector<int>::iterator it1 = bonded_atoms_index.begin(); it1 != bonded_atoms_index.end(); it1++)
                        {
                            int bonded_atom_index = (*it1);
                            PrepFileAtom* bonded_atom = prep_residue->GetPrepAtomByName(prep_residue->GetAtomNameByIndex(bonded_atom_index));
                            stringstream ss;
                            ss << prep_residue->GetName() << ":" << bonded_atom->GetName();
                            for(AtomVector::iterator it2 = all_atoms_of_assembly.begin(); it2 != all_atoms_of_assembly.end(); it2++)
                            {
                                Atom* assembly_atom = (*it2);
                                if(assembly_atom->GetId().compare(ss.str()) == 0)
                                {
                                    atom_node->AddNodeNeighbor(assembly_atom);
                                    break;
                                }
                            }
                        }
                    }
                }

            }
            atom->SetNode(atom_node);
        }
    }
}
