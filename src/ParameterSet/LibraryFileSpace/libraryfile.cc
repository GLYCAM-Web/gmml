#include <fstream>
#include <iostream>
#include <iomanip>

#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfile.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileprocessingexception.hpp"

using LibraryFileSpace::LibraryFile;

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
std::vector<std::string> LibraryFile::GetAllResidueNames()
{
    std::vector<std::string> residue_names;
    for(LibraryFile::ResidueMap::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        std::string residue_name = (*it).first;
        residue_names.push_back(residue_name);
    }
    return residue_names;
}

gmml::ResidueNameMap LibraryFile::GetAllResidueNamesMap()
{
    gmml::ResidueNameMap residue_names = gmml::ResidueNameMap();
    for(LibraryFile::ResidueMap::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        std::string residue_name = (*it).first;
        residue_names[residue_name] = residue_name;
    }
    return residue_names;
}

std::vector<std::string> LibraryFile::GetAllAtomNamesOfResidue(std::string residue_name)
{
    std::vector<std::string> atom_names_of_residue;
    ResidueMap residue_map = GetResidues();
    LibraryFileSpace::LibraryFileResidue* library_file_residue = residue_map[residue_name];
    LibraryFileSpace::LibraryFileResidue::AtomMap atoms = library_file_residue->GetAtoms();
    for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        LibraryFileSpace::LibraryFileAtom* atom = (*it).second;
        atom_names_of_residue.push_back(atom->GetName());
    }
    return atom_names_of_residue;
}

//////////////////////////////////////////////////////////
//                         MUTATORS                     //
//////////////////////////////////////////////////////////
void LibraryFile::SetPath(std::string path)
{
    path_ = path;
}

void LibraryFile::SetResidues(ResidueMap residues)
{
    residues_.clear();
    for(ResidueMap::iterator it = residues.begin(); it != residues.end(); it++)
    {
        LibraryFileSpace::LibraryFileResidue* residue = (*it).second;
        std::string residue_name = (*it).first;
        residues_[residue_name] = residue;
    }
}

//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
/// Read the given file and extract all required fileds into the library file data structure
void LibraryFile::Read(std::ifstream& in_file)
{
    std::string line;

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

    if(line.find("index") != std::string::npos)
    {
        int listing_index = 1;
        getline(in_file, line);
        while(line[0] != '!')
        {
            try
            {
                /// Process index section
                gmml::RemoveQuotes(line);
                gmml::RemoveSpaces(line);
                residues_[line] = new LibraryFileSpace::LibraryFileResidue(line, listing_index);
                listing_index++;
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
        if(line.find("atoms") != std::string::npos)
        {
            int order = 0;
            getline(in_file, line);                 /// Get the first line of the section
            while(line[0] != '!')                   /// Iterate until to the next section that indicates by ! at the begining of the read line
            {
                try
                {
                    /// Process atom section
                    order++;
                    LibraryFileSpace::LibraryFileAtom* newAtom = ProcessAtom(line);
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
        if(line.find("atomspertinfo") != std::string::npos)
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

        if(line.find("boundbox") != std::string::npos)
        {
            getline(in_file, line);                 /// Get the first line of the section
            while(line[0] != '!')                   /// Iterate until to the next section that indicates by ! at the begining of the read line
            {
                try
                {
                    /// Process boundbox section
                    double is_box_set;
                    std::stringstream ss(line);
                    is_box_set = gmml::ConvertString<double>(ss.str());
                    if(is_box_set < 0)              /// If the number written in the first line of the section is negative then boundbox attributes have not been defined in the file
                    {
                        getline(in_file, line);
                        it->second->SetBoxAngle(gmml::dNotSet);
                        getline(in_file, line);
                        it->second->SetBoxLength(gmml::dNotSet);
                        getline(in_file, line);
                        it->second->SetBoxWidth(gmml::dNotSet);
                        getline(in_file, line);
                        it->second->SetBoxHeight(gmml::dNotSet);
                    }
                    else                        /// If the number written in the first line of the section is positive then set the attributes of the boundbox section
                    {
                        getline(in_file, line);
                        double val;

                        line = gmml::Trim(line);
                        std::stringstream angle(line);
                        angle >> val;
                        it->second->SetBoxAngle(val);
                        getline(in_file, line);

                        line = gmml::Trim(line);
                        std::stringstream length(line);
                        length >> val;
                        it->second->SetBoxLength(val);
                        getline(in_file, line);

                        line = gmml::Trim(line);
                        std::stringstream width(line);
                        width >> val;
                        it->second->SetBoxWidth(val);
                        getline(in_file, line);

                        line = gmml::Trim(line);
                        std::stringstream height(line);
                        height >> val;
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
        if(line.find("childsequence") != std::string::npos)
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
        if(line.find("connect") != std::string::npos)
        {
            getline(in_file, line);                 /// Get the first line of the section
            while(line[0] != '!')                   /// Iterate until to the next section that indicates by ! at the begining of the read line
            {
                try
                {
                    /// Process connect section
                    int head_index;
                    std::stringstream ss(line);
                    head_index = gmml::ConvertString<int>(ss.str());
                    getline(in_file,line);

                    int tail_index;
                    std::stringstream sss(line);
                    tail_index = gmml::ConvertString<int>(sss.str());
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
        if(line.find("connectivity") != std::string::npos)
        {
            getline(in_file, line);                 /// Get the first line of the section
            while(line[0] != '!')                   /// Iterate until to the next section that indicates by ! at the begining of the read line
            {
                try
                {
                    /// Process connectivity section
                    std::stringstream ss(line);
                    int from;
                    int to;
                    int t_int;
                    ss >> from >> to >> t_int;
                    it -> second -> GetAtomByOrder(from)->AddBondedAtomIndex(to);
                    it -> second -> GetAtomByOrder(to)->AddBondedAtomIndex(from);
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
        if(line.find("hierarchy") != std::string::npos)
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
        if(line.find("name") != std::string::npos)
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
        if(line.find("positions") != std::string::npos)
        {
            int order = 0;
            getline(in_file, line);                 /// Get the first line of the section
            while(line[0] != '!')                   /// Iterate until to the next section that indicates by ! at the begining of the read line
            {
                try
                {
                    /// Process positions section
                    order ++;
                    std::stringstream ss(line);
                    double x, y, z;
                    ss >> x >> y >> z;
                    GeometryTopology::Coordinate crd(x, y, z);
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
        if(line.find("residueconnect") != std::string::npos)
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
        if(line.find("residues") != std::string::npos)
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
        if(line.find("residuesPdbSequenceNumber") != std::string::npos)
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
        if(line.find("solventcap") != std::string::npos)
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
        if(line.find("velocities") != std::string::npos)
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
LibraryFileSpace::LibraryFileAtom* LibraryFile::ProcessAtom(std::string &line)
{
    std::string name;
    std::string type;
    int residue_index;
    int atom_index;
    int atomic_number;
    double charge;
    int int_t;

    std::stringstream ss(line);             /// Create a stream from the given line
    ss >> name >> type >> int_t >> residue_index >> int_t >> atom_index >> atomic_number >> charge;     /// Split the line by space to extract attributes of the atom

    gmml::RemoveQuotes(name);
    gmml::RemoveQuotes(type);

    LibraryFileSpace::LibraryFileAtom* atom = new LibraryFileSpace::LibraryFileAtom(type, name, residue_index, atom_index, atomic_number, charge);
    return atom;
}

LibraryFileSpace::LibraryFileResidue* LibraryFile::GetLibraryResidueByResidueName(std::string residue_name)
{
    ResidueMap residue_map = this->GetResidues();
    for(ResidueMap::iterator it = residue_map.begin(); it != residue_map.end(); it++)
    {
        std::string name = (*it).first;
        LibraryFileSpace::LibraryFileResidue* residue = (*it).second;
        if(name.compare(residue_name) == 0)
            return residue;
    }
    return NULL;
}
void LibraryFile::Write(const std::string& library_file)
{
    std::ofstream out_file;
    try
    {
        out_file.open(library_file.c_str());
    }
    catch(...)
    {
        throw LibraryFileSpace::LibraryFileProcessingException(__LINE__,"File could not be created");
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
void LibraryFile::BuildLibraryFile(std::ofstream& out_stream)
{
    out_stream << "!!index array str" << std::endl;
    for(ResidueMap::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        LibraryFileSpace::LibraryFileResidue* residue = (*it).second;
        out_stream << "\"" << residue->GetName() << "\"" << std::endl;
    }
    for(ResidueMap::iterator it = residues_.begin(); it != residues_.end(); it++)
    {
        LibraryFileSpace::LibraryFileResidue* residue = (*it).second;
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
void LibraryFile::ResolveAtomSection(std::ofstream &stream, LibraryFileSpace::LibraryFileResidue *residue)
{
    LibraryFileSpace::LibraryFileResidue::AtomMap atoms = residue->GetAtoms();
    const std::string FLAG = "131072";
    stream << "!entry." << residue->GetName() << ".unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg" << std::endl;
    for(unsigned int i = 0; i < atoms.size(); i++)
    {
        LibraryFileSpace::LibraryFileAtom* atom = residue->GetAtomByOrder(i+1);
        stream << "\"" << atom->GetName() << "\" " << "\"" << atom->GetType() << "\" " << "0" << " " << atom->GetResidueIndex() << " " << FLAG << " "
               << atom->GetAtomIndex() << " " << atom->GetAtomicNumber() << " " << std::fixed << atom->GetCharge() << std::endl;
    }
}
void LibraryFile::ResolveAtomPertInfoSection(std::ofstream& stream, LibraryFileSpace::LibraryFileResidue* residue)
{
    stream << "!entry." << residue->GetName() << ".unit.atomspertinfo table  str pname  str ptype  int ptypex  int pelmnt  dbl pchg" << std::endl;
    stream << std::endl;

}
void LibraryFile::ResolveBoundBoxSection(std::ofstream& stream, LibraryFileSpace::LibraryFileResidue* residue)
{
    stream << "!entry." << residue->GetName() << ".unit.boundbox array dbl" << std::endl;
    stream << "-1.000000" << std::endl;
    if(residue->GetBoxAngle() == gmml::dNotSet)
        stream << "0.0" << std::endl;
    else
        stream << residue->GetBoxAngle() << std::endl;
    if(residue->GetBoxLength() == gmml::dNotSet)
        stream << "0.0" << std::endl;
    else
        stream << std::fixed << residue->GetBoxLength() << std::endl;
    if(residue->GetBoxWidth() == gmml::dNotSet)
        stream << "0.0" << std::endl;
    else
        stream << std::fixed << residue->GetBoxWidth() << std::endl;
    if(residue->GetBoxHeight() == gmml::dNotSet)
        stream << "0.0" << std::endl;
    else
        stream << std::fixed << residue->GetBoxHeight() << std::endl;
}
void LibraryFile::ResolveChildSequenceSection(std::ofstream& stream, LibraryFileSpace::LibraryFileResidue* residue)
{
    stream << "!entry." << residue->GetName() << ".unit.childsequence single int" << std::endl;
    stream << std::endl;
}
void LibraryFile::ResolveConnectSection(std::ofstream& stream, LibraryFileSpace::LibraryFileResidue* residue)
{
    stream << "!entry." << residue->GetName() << ".unit.connect array int" << std::endl;
    if(residue->GetHeadAtomIndex() != gmml::iNotSet)
        stream << residue->GetHeadAtomIndex() << std::endl;
    if(residue->GetTailAtomIndex() != gmml::iNotSet)
        stream << residue->GetTailAtomIndex() << std::endl;
}
void LibraryFile::ResolveConnectivitySection(std::ofstream& stream, LibraryFileSpace::LibraryFileResidue* residue)
{
    stream << "!entry." << residue->GetName() << ".unit.connectivity table  int atom1x  int atom2x  int flags" << std::endl;
    LibraryFileSpace::LibraryFileResidue::AtomMap atoms = residue->GetAtoms();
    for(unsigned int i = 0; i < atoms.size(); i++)
    {
        LibraryFileSpace::LibraryFileAtom* atom = residue->GetAtomByOrder(i+1);
        std::vector<int> bonded_atoms_indices = atom->GetBondedAtomsIndices();
//        std::cout << bonded_atoms_indices.size() << std::endl;
        for(std::vector<int>::iterator it = bonded_atoms_indices.begin(); it != bonded_atoms_indices.end(); it++)
        {
            int bonded_atom_index = (*it);
//            std::cout << bonded_atom_index << std::endl;
            if(bonded_atom_index > atom->GetAtomIndex())
            {
                stream << atom->GetAtomIndex() << " " << bonded_atom_index << " " << "1" << std::endl;
            }
        }
    }
}
void LibraryFile::ResolveHierarchySection(std::ofstream& stream, LibraryFileSpace::LibraryFileResidue* residue)
{
    stream << "!entry." << residue->GetName() << ".unit.hierarchy table  str abovetype  int abovex  str belowtype  int belowx" << std::endl;
    stream << std::endl;
}
void LibraryFile::ResolveNameSection(std::ofstream& stream, LibraryFileSpace::LibraryFileResidue* residue)
{
    stream << "!entry." << residue->GetName() << ".unit.name single str" << std::endl;
    stream << "\"" << residue->GetName() << "\"" << std::endl;
}
void LibraryFile::ResolvePositionSection(std::ofstream& stream, LibraryFileSpace::LibraryFileResidue* residue)
{
    stream << "!entry." << residue->GetName() << ".unit.positions table  dbl x  dbl y  dbl z" << std::endl;
    LibraryFileSpace::LibraryFileResidue::AtomMap atoms = residue->GetAtoms();
    for(unsigned int i = 0; i < atoms.size(); i++)
    {
        LibraryFileSpace::LibraryFileAtom* atom = residue->GetAtomByOrder(i+1);
        GeometryTopology::Coordinate coordinate = atom->GetCoordinate();
        stream << std::fixed << coordinate.GetX() << " " << std::fixed << coordinate.GetY() << " " << std::fixed << coordinate.GetZ() << std::endl;
    }
}
void LibraryFile::ResolveResidueConnectSection(std::ofstream& stream, LibraryFileSpace::LibraryFileResidue* residue)
{
    stream << "!entry." << residue->GetName() << ".unit.residueconnect table  int c1x  int c2x  int c3x  int c4x  int c5x  int c6x" << std::endl;
    stream << std::endl;
}
void LibraryFile::ResolveResiduesSection(std::ofstream& stream, LibraryFileSpace::LibraryFileResidue* residue)
{
    stream << "!entry." << residue->GetName() << ".unit.residues table  str name  int seq  int childseq  int startatomx  str restype  int imagingx" << std::endl;
    stream << std::endl;
}
void LibraryFile::ResolveResiduePdbSequenceNumberSection(std::ofstream& stream, LibraryFileSpace::LibraryFileResidue* residue)
{
    stream << "!entry." << residue->GetName() <<  ".unit.residuesPdbSequenceNumber array int" << std::endl;
    stream << std::endl;
}
void LibraryFile::ResolveSolventCapSection(std::ofstream& stream, LibraryFileSpace::LibraryFileResidue* residue)
{
    stream << "!entry." << residue->GetName() << ".unit.solventcap array dbl" << std::endl;
    stream << std::endl;
}
void LibraryFile::ResolveVelocitiesSection(std::ofstream& stream, LibraryFileSpace::LibraryFileResidue* residue)
{
    stream << "!entry." << residue->GetName() << ".unit.velocities table  dbl x  dbl y  dbl z" << std::endl;
    stream << std::endl;
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
        //                std::cout << it -> second -> GetAtomByIndex(i+1)->GetBondedAtomsIndicies()[j] << ", ";
        //            std::cout << std::endl;
        //        }
        //        std::cout << std::endl;
        it->second->Print(out);
    }
}
