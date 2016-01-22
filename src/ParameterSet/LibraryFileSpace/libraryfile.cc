#include <fstream>
#include <iostream>
#include <iomanip>

#include "utils.hpp"
#include "common.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfile.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileprocessingexception.hpp"

using namespace std;
using namespace gmml;
using namespace LibraryFileSpace;
using namespace GeometryTopology;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
LibraryFile::LibraryFile() : path_("GMML-Generated"){}

LibraryFile::LibraryFile(const std::string &lib_file)
{
    path_ = lib_file;
    std::ifstream in_file;
    if(std::ifstream(lib_file.c_str()))
        in_file.open(lib_file.c_str());
    else
    {
        throw LibraryFileProcessingException(__LINE__, "Library file not found");
    }
    Read(in_file);
    in_file.close();            /// Close the parameter files
}

LibraryFile::~LibraryFile(){}
//////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////
/// Return the library file path
const std::string& LibraryFile::GetFilePath() const
{
    return path_;
}

/// Return a map of residues including in the file mapped to the names of the residues
const LibraryFile::ResidueMap& LibraryFile::GetResidues() const
{
    return residues_;
}
vector<string> LibraryFile::GetAllResidueNames()
{
    vector<string> residue_names;
    for(LibraryFile::ResidueMap::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        string residue_name = (*it).first;
        residue_names.push_back(residue_name);
    }
    return residue_names;
}

ResidueNameMap LibraryFile::GetAllResidueNamesMap()
{
    ResidueNameMap residue_names = ResidueNameMap();
    for(LibraryFile::ResidueMap::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        string residue_name = (*it).first;
        residue_names[residue_name] = residue_name;
    }
    return residue_names;
}

vector<string> LibraryFile::GetAllAtomNamesOfResidue(string residue_name)
{
    vector<string> atom_names_of_residue;
    ResidueMap residue_map = GetResidues();
    LibraryFileResidue* library_file_residue = residue_map[residue_name];
    LibraryFileResidue::AtomMap atoms = library_file_residue->GetAtoms();
    for(LibraryFileResidue::AtomMap::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        LibraryFileAtom* atom = (*it).second;
        atom_names_of_residue.push_back(atom->GetName());
    }
    return atom_names_of_residue;
}

//////////////////////////////////////////////////////////
//                         MUTATORS                     //
//////////////////////////////////////////////////////////
void LibraryFile::SetPath(string path)
{
    path_ = path;
}

void LibraryFile::SetResidues(ResidueMap residues)
{
    residues_.clear();
    for(ResidueMap::iterator it = residues.begin(); it != residues.end(); it++)
    {
        LibraryFileResidue* residue = (*it).second;
        string residue_name = (*it).first;
        residues_[residue_name] = residue;
    }
}

//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
/// Read the given file and extract all required fileds into the library file data structure
void LibraryFile::Read(std::ifstream& in_file)
{
    string line;

    /// Unable to read file
    if (!getline(in_file, line))
    {
        throw LibraryFileProcessingException("Error reading file");
    }
    /// Skip blank lines at the begining of the file
    while(line[0] != '!')
    {
        getline(in_file, line);
    }

    if(line.find("index") != string::npos)
    {
        getline(in_file, line);
        while(line[0] != '!')
        {
            try
            {
                /// Process index section
                RemoveQuotes(line);
                RemoveSpaces(line);
                residues_[line] = new LibraryFileResidue(line);
                getline(in_file,line);      /// Read the next line
            } catch(...)
            {
                throw LibraryFileProcessingException(__LINE__, "Error processing index section");
            }
        }
    }

    /// Iterate on all residues indicated in the index section of the file
    for(LibraryFile::ResidueMap::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        /// Process the atom section of the file for the corresponding residue
        if(line.find("atoms") != string::npos)
        {
            int order = 0;
            getline(in_file, line);                 /// Get the first line of the section
            while(line[0] != '!')                   /// Iterate until to the next section that indicates by ! at the begining of the read line
            {
                try
                {
                    /// Process atom section
                    order++;
                    LibraryFileAtom* newAtom = ProcessAtom(line);
                    newAtom->SetAtomOrder(order);;
                    it->second->AddAtom(newAtom);               /// Add the atom into the list of atoms of the corresponding residue

                    getline(in_file,line);      /// Read the next line
                } catch(...)
                {
                    throw LibraryFileProcessingException(__LINE__, "Error processing atom section");
                }
            }
        }
//        else
//        {
//            throw LibraryFileProcessingException(__LINE__, "Unknown section or missing section");
//        }

        /// Process the atompertinfo section of the file for the corresponding residue
        if(line.find("atomspertinfo") != string::npos)
        {
            getline(in_file, line);                 /// Get the first line of the section
            while(line[0] != '!')                   /// Iterate until to the next section that indicates by ! at the begining of the read line
            {
                try
                {
                    /// Process atompertinfo section -> This section doesn't have any useful information -> Ignore the section by reading the file until the next section
                    getline(in_file,line);      /// Read the next line
                } catch(...)
                {
                    throw LibraryFileProcessingException(__LINE__, "Error processing atom pert info section");
                }
            }
        }
//        else
//        {
//            throw LibraryFileProcessingException(__LINE__, "Unknown section or missing section");
//        }
        /// Process the boundbox section of the file for the corresponding residue

        if(line.find("boundbox") != string::npos)
        {
            getline(in_file, line);                 /// Get the first line of the section
            while(line[0] != '!')                   /// Iterate until to the next section that indicates by ! at the begining of the read line
            {
                try
                {
                    /// Process boundbox section
                    double is_box_set;
                    istringstream ss(line);
                    ss >> is_box_set;
                    if(is_box_set < 0)              /// If the number written in the first line of the section is negative then boundbox attributes have not been defined in the file
                    {
                        getline(in_file, line);
                        it->second->SetBoxAngle(dNotSet);
                        getline(in_file, line);
                        it->second->SetBoxLength(dNotSet);
                        getline(in_file, line);
                        it->second->SetBoxWidth(dNotSet);
                        getline(in_file, line);
                        it->second->SetBoxHeight(dNotSet);
                    }
                    else                        /// If the number written in the first line of the section is positive then set the attributes of the boundbox section
                    {
                        getline(in_file, line);
                        double val;

                        ss.str(line);
                        ss >> val;
                        it->second->SetBoxAngle(val);
                        getline(in_file, line);

                        ss.str(line);
                        ss >> val;
                        it->second->SetBoxLength(val);
                        getline(in_file, line);

                        ss.str(line);
                        ss >> val;
                        it->second->SetBoxWidth(val);
                        getline(in_file, line);

                        ss.str(line);
                        ss >> val;
                        it->second->SetBoxHeight(val);
                    }
                    getline(in_file,line);      /// Read the next line
                } catch(...)
                {
                    throw LibraryFileProcessingException(__LINE__, "Error processing bound box section");
                }
            }
        }
//        else
//        {
//            throw LibraryFileProcessingException(__LINE__, "Unknown section or missing section");
//        }

        /// Process the childsequence section of the file for the corresponding residue
        if(line.find("childsequence") != string::npos)
        {
            getline(in_file, line);                 /// Get the first line of the section
            while(line[0] != '!')                   /// Iterate until to the next section that indicates by ! at the begining of the read line
            {
                try
                {
                    /// Process childsequence section -> This section doesn't have any useful information -> Ignore the section by reading the file until the next section
                    getline(in_file,line);
                } catch(...)
                {
                    throw LibraryFileProcessingException(__LINE__, "Error processing child sequence section");
                }
            }
        }
//        else
//        {
//            throw LibraryFileProcessingException(__LINE__, "Unknown section or missing section");
//        }

        /// Process the connect section of the file for the corresponding residue
        if(line.find("connect") != string::npos)
        {
            getline(in_file, line);                 /// Get the first line of the section
            while(line[0] != '!')                   /// Iterate until to the next section that indicates by ! at the begining of the read line
            {
                try
                {
                    /// Process connect section
                    int head_index;
                    istringstream ss(line);
                    ss >> head_index;
                    getline(in_file,line);

                    int tail_index;
                    istringstream sss(line);
                    sss >> tail_index;
                    getline(in_file,line);

                    it->second->SetHeadAtomIndex(head_index);
                    it->second->SetTailAtomIndex(tail_index);

                } catch(...)
                {
                    throw LibraryFileProcessingException(__LINE__, "Error processing connect section");
                }
            }
        }
//        else
//        {
//            throw LibraryFileProcessingException(__LINE__, "Unknown section or missing section");
//        }

        /// Process the connectivity section of the file for the corresponding residue
        if(line.find("connectivity") != string::npos)
        {
            getline(in_file, line);                 /// Get the first line of the section
            while(line[0] != '!')                   /// Iterate until to the next section that indicates by ! at the begining of the read line
            {
                try
                {
                    /// Process connectivity section
                    istringstream ss(line);
                    int from;
                    int to;
                    int t_int;
                    ss >> from >> to >> t_int;
                    it -> second -> GetAtomByIndex(from)->AddBondedAtomIndex(to);
                    it -> second -> GetAtomByIndex(to)->AddBondedAtomIndex(from);
                    getline(in_file,line);      /// Read the next line
                } catch(...)
                {
                    throw LibraryFileProcessingException(__LINE__, "Error processing connectivity section");
                }
            }
        }
//        else
//        {
//            throw LibraryFileProcessingException(__LINE__, "Unknown section or missing section");
//        }

        /// Process the hierarchy section of the file for the corresponding residue
        if(line.find("hierarchy") != string::npos)
        {
            getline(in_file, line);                 /// Get the first line of the section
            while(line[0] != '!')                   /// Iterate until to the next section that indicates by ! at the begining of the read line
            {
                try
                {
                    /// Process hierarchy section -> This section doesn't have any useful information -> Ignore the section by reading the file until the next section
                    getline(in_file,line);      /// Read the next line
                } catch(...)
                {
                    throw LibraryFileProcessingException(__LINE__, "Error processing hierarchy section");
                }
            }
        }
//        else
//        {
//            throw LibraryFileProcessingException(__LINE__, "Unknown section or missing section");
//        }

        /// Process the name section of the file for the corresponding residue
        if(line.find("name") != string::npos)
        {
            getline(in_file, line);                 /// Get the first line of the section
            while(line[0] != '!')                   /// Iterate until to the next section that indicates by ! at the begining of the read line
            {
                try
                {
                    /// Process name section -> This section doesn't have any useful information -> Ignore the section by reading the file until the next section
                    getline(in_file,line);      /// Read the next line
                } catch(...)
                {
                    throw LibraryFileProcessingException(__LINE__, "Error processing name section");
                }
            }
        }
//        else
//        {
//            throw LibraryFileProcessingException(__LINE__, "Unknown section or missing section");
//        }

        /// Process the positions section of the file for the corresponding residue
        if(line.find("positions") != string::npos)
        {
            int order = 0;
            getline(in_file, line);                 /// Get the first line of the section
            while(line[0] != '!')                   /// Iterate until to the next section that indicates by ! at the begining of the read line
            {
                try
                {
                    /// Process positions section
                    order ++;
                    istringstream ss(line);
                    double x, y, z;
                    ss >> x >> y >> z;
                    Coordinate crd(x, y, z);
                    it->second->GetAtomByOrder(order)->SetCoordinate(crd);
                    getline(in_file,line);      /// Read the next line
                } catch(...)
                {
                    throw LibraryFileProcessingException(__LINE__, "Error processing positions section");
                }
            }
        }
//        else
//        {
//            throw LibraryFileProcessingException(__LINE__, "Unknown section or missing section");
//        }

        /// Process the residueconnect section of the file for the corresponding residue
        if(line.find("residueconnect") != string::npos)
        {
            getline(in_file, line);                 /// Get the first line of the section
            while(line[0] != '!')                   /// Iterate until to the next section that indicates by ! at the begining of the read line
            {
                try
                {
                    /// Process residueconnect section -> This section doesn't have any useful information -> Ignore the section by reading the file until the next section
                    getline(in_file,line);      /// Read the next line
                } catch(...)
                {
                    throw LibraryFileProcessingException(__LINE__, "Error processing residue connect section");
                }
            }
        }
//        else
//        {
//            throw LibraryFileProcessingException(__LINE__, "Unknown section or missing section");
//        }

        /// Process the residues section of the file for the corresponding residue
        if(line.find("residues") != string::npos)
        {
            getline(in_file, line);                 /// Get the first line of the section
            while(line[0] != '!')                   /// Iterate until to the next section that indicates by ! at the begining of the read line
            {
                try
                {
                    /// Process residues section -> This section doesn't have any useful information -> Ignore the section by reading the file until the next section
                    getline(in_file,line);      /// Read the next line
                } catch(...)
                {
                    throw LibraryFileProcessingException(__LINE__, "Error processing residues section");
                }
            }
        }
//        else
//        {
//            throw LibraryFileProcessingException(__LINE__, "Unknown section or missing section");
//        }

        /// Process the residuesPdbSequenceNumber section of the file for the corresponding residue
        if(line.find("residuesPdbSequenceNumber") != string::npos)
        {
            getline(in_file, line);                 /// Get the first line of the section
            while(line[0] != '!')                   /// Iterate until to the next section that indicates by ! at the begining of the read line
            {
                try
                {
                    /// Process residuesPdbSequenceNumber section -> This section doesn't have any useful information -> Ignore the section by reading the file until the next section
                    getline(in_file,line);      /// Read the next line
                } catch(...)
                {
                    throw LibraryFileProcessingException(__LINE__, "Error processing residuePdbSequenceNumber section");
                }
            }
        }
//        else
//        {
//            throw LibraryFileProcessingException(__LINE__, "Unknown section or missing section");
//        }

        /// Process the solventcap section of the file for the corresponding residue
        if(line.find("solventcap") != string::npos)
        {
            getline(in_file, line);                 /// Get the first line of the section
            while(line[0] != '!')                   /// Iterate until to the next section that indicates by ! at the begining of the read line
            {
                try
                {
                    /// Process solventcap section -> This section doesn't have any useful information -> Ignore the section by reading the file until the next section
                    getline(in_file,line);      /// Read the next line
                } catch(...)
                {
                    throw LibraryFileProcessingException(__LINE__, "Error processing solventcap section");
                }
            }
        }
//        else
//        {
//            throw LibraryFileProcessingException(__LINE__, "Unknown section or missing section");
//        }

        /// Process the velocities section of the file for the corresponding residue
        if(line.find("velocities") != string::npos)
        {
            getline(in_file, line);                 /// Get the first line of the section
            while(line[0] != '!')                   /// Iterate until to the next section that indicates by ! at the begining of the read line
            {
                try
                {
                    /// Process velocities section -> This section doesn't have any useful information -> Ignore the section by reading the file until the next section
                    if(!getline(in_file,line))
                        return;      /// Read the next line
                } catch(...)
                {
                    throw LibraryFileProcessingException(__LINE__, "Error processing velocities section");
                }
            }
        }
//        else
//        {
//            throw LibraryFileProcessingException(__LINE__, "Unknown section or missing section");
//        }
    }
}

/// Process a line from the atom section of the file and return a new atom object
LibraryFileAtom* LibraryFile::ProcessAtom(string &line)
{
    string name;
    string type;
    int residue_index;
    int atom_index;
    int atomic_number;
    double charge;
    int int_t;

    istringstream ss(line);             /// Create a stream from the given line
    ss >> name >> type >> int_t >> residue_index >> int_t >> atom_index >> atomic_number >> charge;     /// Split the line by space to extract attributes of the atom

    RemoveQuotes(name);
    RemoveQuotes(type);

    LibraryFileAtom* atom = new LibraryFileAtom(type, name, residue_index, atom_index, atomic_number, charge);
    return atom;
}

LibraryFileResidue* LibraryFile::GetLibraryResidueByResidueName(string residue_name)
{
    ResidueMap residue_map = this->GetResidues();
    for(ResidueMap::iterator it = residue_map.begin(); it != residue_map.end(); it++)
    {
        string name = (*it).first;
        LibraryFileResidue* residue = (*it).second;
        if(name.compare(residue_name) == 0)
            return residue;
    }
    return NULL;
}
void LibraryFile::Write(const string& library_file)
{
    std::ofstream out_file;
    try
    {
        out_file.open(library_file.c_str());
    }
    catch(...)
    {
        throw LibraryFileProcessingException(__LINE__,"File could not be created");
    }
    try
    {
        this->BuildLibraryFile(out_file);
    }
    catch(...)
    {
        out_file.close();
    }
}
void LibraryFile::BuildLibraryFile(ofstream& out_stream)
{
    out_stream << "!!index array str" << endl;
    for(ResidueMap::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        LibraryFileResidue* residue = (*it).second;
        out_stream << "\"" << residue->GetName() << "\"" << endl;
    }
    for(ResidueMap::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        LibraryFileResidue* residue = (*it).second;
        ResolveAtomSection(out_stream, residue);
        ResolveAtomPertInfoSection(out_stream, residue);
        ResolveBoundBoxSection(out_stream, residue);
        ResolveChildSequenceSection(out_stream, residue);
        ResolveConnectSection(out_stream, residue);
        ResolveConnectivitySection(out_stream, residue);
        ResolveHierarchySection(out_stream, residue);
        ResolveNameSection(out_stream, residue);
        ResolvePositionSection(out_stream, residue);
        ResolveResidueConnectSection(out_stream, residue);
        ResolveResiduesSection(out_stream, residue);
        ResolveResiduePdbSequenceNumberSection(out_stream, residue);
        ResolveSolventCapSection(out_stream, residue);
        ResolveVelocitiesSection(out_stream, residue);
    }
}
void LibraryFile::ResolveAtomSection(ofstream &stream, LibraryFileResidue *residue)
{
    LibraryFileResidue::AtomMap atoms = residue->GetAtoms();
    const string FLAG = "131072";
    stream << "!entry." << residue->GetName() << ".unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg" << endl;
    for(unsigned int i = 0; i < atoms.size(); i++)
    {
        LibraryFileAtom* atom = residue->GetAtomByOrder(i+1);
        stream << "\"" << atom->GetName() << "\" " << "\"" << atom->GetType() << "\" " << "0" << " " << atom->GetResidueIndex() << " " << FLAG << " "
               << atom->GetAtomIndex() << " " << atom->GetAtomicNumber() << " " << fixed << atom->GetCharge() << endl;
    }
}
void LibraryFile::ResolveAtomPertInfoSection(ofstream& stream, LibraryFileResidue* residue)
{
    stream << "!entry." << residue->GetName() << ".unit.atomspertinfo table  str pname  str ptype  int ptypex  int pelmnt  dbl pchg" << endl;
    stream << endl;

}
void LibraryFile::ResolveBoundBoxSection(ofstream& stream, LibraryFileResidue* residue)
{
    stream << "!entry." << residue->GetName() << ".unit.boundbox array dbl" << endl;
    stream << "-1.000000" << endl;
    if(residue->GetBoxAngle() == dNotSet)
        stream << "0.0" << endl;
    else
        stream << residue->GetBoxAngle() << endl;
    if(residue->GetBoxLength() == dNotSet)
        stream << "0.0" << endl;
    else
        stream << fixed << residue->GetBoxLength() << endl;
    if(residue->GetBoxWidth() == dNotSet)
        stream << "0.0" << endl;
    else
        stream << fixed << residue->GetBoxWidth() << endl;
    if(residue->GetBoxHeight() == dNotSet)
        stream << "0.0" << endl;
    else
        stream << fixed << residue->GetBoxHeight() << endl;
}
void LibraryFile::ResolveChildSequenceSection(ofstream& stream, LibraryFileResidue* residue)
{
    stream << "!entry." << residue->GetName() << ".unit.childsequence single int" << endl;
    stream << endl;
}
void LibraryFile::ResolveConnectSection(ofstream& stream, LibraryFileResidue* residue)
{
    stream << "!entry." << residue->GetName() << ".unit.connect array int" << endl;
    if(residue->GetHeadAtomIndex() != iNotSet)
        stream << residue->GetHeadAtomIndex() << endl;
    if(residue->GetTailAtomIndex() != iNotSet)
        stream << residue->GetTailAtomIndex() << endl;
}
void LibraryFile::ResolveConnectivitySection(ofstream& stream, LibraryFileResidue* residue)
{
    stream << "!entry." << residue->GetName() << ".unit.connectivity table  int atom1x  int atom2x  int flags" << endl;
    LibraryFileResidue::AtomMap atoms = residue->GetAtoms();
    for(unsigned int i = 0; i < atoms.size(); i++)
    {
        LibraryFileAtom* atom = residue->GetAtomByOrder(i+1);
        vector<int> bonded_atoms_indices = atom->GetBondedAtomsIndices();
//        cout << bonded_atoms_indices.size() << endl;
        for(vector<int>::iterator it = bonded_atoms_indices.begin(); it != bonded_atoms_indices.end(); it++)
        {
            int bonded_atom_index = (*it);
//            cout << bonded_atom_index << endl;
            if(bonded_atom_index > atom->GetAtomIndex())
            {
                stream << atom->GetAtomIndex() << " " << bonded_atom_index << " " << "1" << endl;
            }
        }
    }
}
void LibraryFile::ResolveHierarchySection(ofstream& stream, LibraryFileResidue* residue)
{
    stream << "!entry." << residue->GetName() << ".unit.hierarchy table  str abovetype  int abovex  str belowtype  int belowx" << endl;
    stream << endl;
}
void LibraryFile::ResolveNameSection(ofstream& stream, LibraryFileResidue* residue)
{
    stream << "!entry." << residue->GetName() << ".unit.name single str" << endl;
    stream << "\"" << residue->GetName() << "\"" << endl;
}
void LibraryFile::ResolvePositionSection(ofstream& stream, LibraryFileResidue* residue)
{
    stream << "!entry." << residue->GetName() << ".unit.positions table  dbl x  dbl y  dbl z" << endl;
    LibraryFileResidue::AtomMap atoms = residue->GetAtoms();
    for(unsigned int i = 0; i < atoms.size(); i++)
    {
        LibraryFileAtom* atom = residue->GetAtomByOrder(i+1);
        GeometryTopology::Coordinate coordinate = atom->GetCoordinate();
        stream << fixed << coordinate.GetX() << " " << fixed << coordinate.GetY() << " " << fixed << coordinate.GetZ() << endl;
    }
}
void LibraryFile::ResolveResidueConnectSection(ofstream& stream, LibraryFileResidue* residue)
{
    stream << "!entry." << residue->GetName() << ".unit.residueconnect table  int c1x  int c2x  int c3x  int c4x  int c5x  int c6x" << endl;
    stream << endl;
}
void LibraryFile::ResolveResiduesSection(ofstream& stream, LibraryFileResidue* residue)
{
    stream << "!entry." << residue->GetName() << ".unit.residues table  str name  int seq  int childseq  int startatomx  str restype  int imagingx" << endl;
    stream << endl;
}
void LibraryFile::ResolveResiduePdbSequenceNumberSection(ofstream& stream, LibraryFileResidue* residue)
{
    stream << "!entry." << residue->GetName() <<  ".unit.residuesPdbSequenceNumber array int" << endl;
    stream << endl;
}
void LibraryFile::ResolveSolventCapSection(ofstream& stream, LibraryFileResidue* residue)
{
    stream << "!entry." << residue->GetName() << ".unit.solventcap array dbl" << endl;
    stream << endl;
}
void LibraryFile::ResolveVelocitiesSection(ofstream& stream, LibraryFileResidue* residue)
{
    stream << "!entry." << residue->GetName() << ".unit.velocities table  dbl x  dbl y  dbl z" << endl;
    stream << endl;
}

//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
void LibraryFile::Print(std::ostream& out)
{
    for(LibraryFile::ResidueMap::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        //        for(unsigned int i = 0; i < it -> second -> GetAtoms().size(); i++)
        //        {
        //            for(unsigned int j = 0; j < it->second->GetAtomByIndex(i+1)->GetBondedAtomsIndicies().size(); j++)
        //                cout << it -> second -> GetAtomByIndex(i+1)->GetBondedAtomsIndicies()[j] << ", ";
        //            cout << endl;
        //        }
        //        cout << endl;
        it->second->Print(out);
    }
}
