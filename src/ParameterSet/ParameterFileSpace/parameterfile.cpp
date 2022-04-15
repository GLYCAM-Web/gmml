#include <fstream>
#include <iostream>
#include <iomanip>

#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfile.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp"
#include "../../../includes/ParameterSet/ParameterFileSpace/parameterfileprocessingexception.hpp"

using ParameterFileSpace::ParameterFile;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
ParameterFile::ParameterFile(std::string param_file, int file_type)
{
    path_ = param_file;
    file_type_ = file_type;
    std::ifstream in_file;
    if(std::ifstream(param_file.c_str()))
        in_file.open(param_file.c_str());
    else
    {
        throw ParameterFileProcessingException(__LINE__, "Parameter file not found");
    }
    switch(file_type_)
    {
        case gmml::MAIN:
            ReadMainParameter(in_file);
            break;
        case gmml::MODIFIED:
            ReadModifiedParameter(in_file);
            break;
        case gmml::IONICMOD:
            ReadIonicModifiedParameter(in_file);
            break;
    }
    in_file.close();            /// Close the parameter files

}

//////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////
const std::string& ParameterFile::GetFilePath()
{
    return path_;
}
const std::string& ParameterFile::GetTitle()
{
    return title_;
}
const ParameterFile::AtomTypeMap& ParameterFile::GetAtomTypes()
{
    return atom_types_;
}
const ParameterFile::BondMap& ParameterFile::GetBonds()
{
    return bonds_;
}
const ParameterFile::AngleMap& ParameterFile::GetAngles()
{
    return angles_;
}
const ParameterFile::DihedralMap& ParameterFile::GetDihedrals()
{
    return dihedrals_;
}
int ParameterFile::GetParameterFileType()
{
    return file_type_;
}
ParameterFile::DihedralMap ParameterFile::GetAllImproperDihedrals()
{
    DihedralMap improper_dihedral_map = DihedralMap();
    for(DihedralMap::iterator it = dihedrals_.begin(); it != dihedrals_.end(); it++)
    {
        ParameterFileDihedral* dihedral = (*it).second;
        if(dihedral->GetIsImproper())
            improper_dihedral_map[(*it).first] = dihedral;///???
    }
    return improper_dihedral_map;
}
ParameterFile::DihedralMap ParameterFile::GetAllproperDihedrals()
{
    DihedralMap proper_dihedral_map = DihedralMap();
    for(DihedralMap::iterator it = dihedrals_.begin(); it != dihedrals_.end(); it++)
    {
        ParameterFileDihedral* dihedral = (*it).second;
        if(!dihedral->GetIsImproper())
            proper_dihedral_map[(*it).first] = dihedral;///???
    }
    return proper_dihedral_map;
}

//////////////////////////////////////////////////////////
//                         MUTATOR                      //
//////////////////////////////////////////////////////////
void ParameterFile::SetParameterFileType(int file_type)
{
    file_type_ = file_type;
}

//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
/// Read from the given stream and extract into the parameter file data structure
/// This function reads the stream line-by-line and uses the whole structure of the file to point to the std::right position for data extraction
/// For more information check a sample parameter file to find out how the function works
void ParameterFile::ReadMainParameter(std::ifstream& in_file)
{
    std::string line;
    int line_number = 0;

    /// Unable to read file
    if (!getline(in_file, line))
    {
        throw ParameterFileProcessingException("Error reading file");
    }

    /// Set the title of the parameter file
    title_ = gmml::Trim(line);
    line_number++;

    /// Atom type section reading
    getline(in_file, line);             /// Read the first line of the atom type section
    line_number++;                      /// Increment line counter
    while (!gmml::Trim(line).empty())         /// Reading until a blank line which is the end of the section
    {
        try
        {
            ProcessAtomType(line);      /// Processing atom type line
            getline(in_file,line);      /// Read the next line
            line_number++;              /// Increment line counter
        } catch(...)
        {
            throw ParameterFileProcessingException(line_number, "Error processing atom type");
        }
    }

    /// Hydrophilic atom type section reading
    getline(in_file, line);             /// Read the first line of the hydrophilic atom type section
    line_number++;                      /// Increment line counter
    while (line.find('-') == std::string::npos)      /// Process the line until the first line of the bond section gets read
    {
        try
        {
            ProcessHydrophilicAtomType(line);   /// Processing hydrophilic atom type line
            getline(in_file, line);             /// Read the next line
            line_number++;                      /// Increment line counter
        } catch(...)
        {
            throw ParameterFileProcessingException(line_number, "Error processing hydrophilic atom list");
        }
    }

    /// Bond section reading
    while (!gmml::Trim(line).empty())         /// Reading until a blank line which is the end of the section
    {
        try
        {
            ProcessBond(line);          /// Processing bond line
            getline(in_file, line);     /// Read the next line
            line_number++;              /// Increment line counter
        } catch(...)
        {
            throw ParameterFileProcessingException(line_number, "Error processing bond");
        }
    }

    /// Angle section reading
    getline(in_file, line);             /// Read the first line of the angle section
    line_number++;                      /// Increment line counter
    while (!gmml::Trim(line).empty())         /// Reading until a blank line which is the end of the section
    {
        try
        {
            ProcessAngle(line);         /// Processing angle line
            getline(in_file, line);     /// Read the next line
            line_number++;              /// Increment line counter
        } catch(...)
        {
            throw ParameterFileProcessingException(line_number, "Error processing angle");
        }
    }

    /// Dihedral section reading
    getline(in_file, line);             /// Read the first line of the dihedral section
    line_number++;                      /// Increment line counter
    while (!gmml::Trim(line).empty())         /// Reading until a blank line which is the end of the section
    {
        try
        {
            ProcessDihedral(line, line_number, in_file);    /// Processing dihedral line
            getline(in_file, line);                         /// Read the next line
            line_number++;                                  /// Increment line counter
        } catch(...)
        {
            throw ParameterFileProcessingException(line_number, "Error processing dihedral");
        }
    }

    /// Improper dihedral section reading
    getline(in_file, line);             /// Read the first line of the improper dihedral section
    line_number++;                      /// Increment line counter
    while (!gmml::Trim(line).empty())         /// Reading until a blank line which is the end of the section
    {
        try
        {
            ProcessImproperDihedral(line, line_number, in_file);        /// Processing improper dihedral line
            getline(in_file, line);                                     /// Read the next line
            line_number++;                                              /// Increment line counter
        } catch(...)
        {
            throw ParameterFileProcessingException(line_number, "Error processing improper dihedral");
        }
    }

    /// Hydrogen bond section reading
    getline(in_file, line);             /// Read the first line of the hydrogen-bond section
    line_number++;                      /// Increment line counter
    while (!gmml::Trim(line).empty())         /// Reading until a blank line which is the end of the section
    {
        try
        {
            ProcessHydrogenBond(line);      /// Processing hydrogen-bond line
            getline(in_file, line);         /// Read the next line
            line_number++;                  /// Increment line counter
        } catch(...)
        {
            throw ParameterFileProcessingException(line_number, "Error processing hydrogen-bond parameters");
        }
    }

    /// Equivalent symbol section reading
    getline(in_file, line);             /// Read the first line of the equivalent symbol section
    line_number++;                      /// Increment line counter
    while (!gmml::Trim(line).empty())         /// Reading until a blank line which is the end of the section
    {
        try
        {
            ProcessEquivalentSymbols(line);     /// Processing equivalent symbol line
            getline(in_file, line);             /// Read the next line
            line_number++;                      /// Increment line counter
        } catch(...)
        {
            throw ParameterFileProcessingException(line_number, "Error processing equivalent symbol list");
        }
    }

    /// Skip MOD4 line
    getline(in_file, line);         /// Read the title of potential parameter section
    line_number++;                  /// Increment line counter

    /// Potential parameter section reading
    getline(in_file, line);         /// Read the first line of the potential parameter section
    line_number++;                  /// Increment line counter
    while (!gmml::Trim(line).empty())         /// Reading until a blank line which is the end of the section
    {
        try
        {
            ProcessPotentialParameter(line);        /// Processing potential parameter line
            getline(in_file, line);                 /// Read the next line
            line_number++;                          /// Increment line counter
        } catch(...)
        {
            throw ParameterFileProcessingException(line_number, "Error processing potential parameters");
        }
    }
}

void ParameterFile::ReadModifiedParameter(std::ifstream& in_file)
{
    std::string line;
    int line_number = 0;

    /// Unable to read file
    if (!getline(in_file, line))
    {
        throw ParameterFileProcessingException("Error reading file");
    }

    /// Set the title of the parameter file
    title_ = gmml::Trim(line);
    line_number++;

    /// Atom type section reading
    getline(in_file, line);             /// Read the first line of the atom type section
    line_number++;                      /// Increment line counter
    if(gmml::Trim(line).find("MASS") != std::string::npos)
    {
        getline(in_file, line);             /// Read the first line of the atom type section
        line_number++;                      /// Increment line counter
        while (!gmml::Trim(line).empty())         /// Reading until a blank line which is the end of the section
        {
            try
            {
                ProcessAtomType(line);      /// Processing atom type line
                getline(in_file,line);      /// Read the next line
                line_number++;              /// Increment line counter
            } catch(...)
            {
                throw ParameterFileProcessingException(line_number, "Error processing atom type");
            }
        }
    }

    /// Bond section reading
    getline(in_file, line);             /// Read the first line of the atom type section
    line_number++;                      /// Increment line counter
    if(gmml::Trim(line).find("BOND") != std::string::npos)
    {
        getline(in_file, line);             /// Read the first line of the atom type section
        line_number++;                      /// Increment line counter
        while (!gmml::Trim(line).empty())         /// Reading until a blank line which is the end of the section
        {
            try
            {
                ProcessBond(line);          /// Processing bond line
                getline(in_file, line);     /// Read the next line
                line_number++;              /// Increment line counter
            } catch(...)
            {
                throw ParameterFileProcessingException(line_number, "Error processing bond");
            }
        }
    }

    /// Angle section reading
    getline(in_file, line);             /// Read the first line of the angle section
    line_number++;                      /// Increment line counter
    if(gmml::Trim(line).find("ANGL") != std::string::npos || gmml::Trim(line).find("ANGLE") != std::string::npos)
    {
        getline(in_file, line);             /// Read the first line of the angle section
        line_number++;                      /// Increment line counter
        while (!gmml::Trim(line).empty())         /// Reading until a blank line which is the end of the section
        {
            try
            {
                ProcessAngle(line);         /// Processing angle line
                getline(in_file, line);     /// Read the next line
                line_number++;              /// Increment line counter
            } catch(...)
            {
                throw ParameterFileProcessingException(line_number, "Error processing angle");
            }
        }
    }

    /// Dihedral section reading
    getline(in_file, line);             /// Read the first line of the dihedral section
    line_number++;                      /// Increment line counter
    if(gmml::Trim(line).find("DIHE") != std::string::npos || gmml::Trim(line).find("DIHEDRAL") != std::string::npos)
    {
        getline(in_file, line);             /// Read the first line of the angle section
        line_number++;                      /// Increment line counter
        while (!gmml::Trim(line).empty())         /// Reading until a blank line which is the end of the section
        {
            try
            {
                ProcessDihedral(line, line_number, in_file);    /// Processing dihedral line
                getline(in_file, line);                         /// Read the next line
                line_number++;                                  /// Increment line counter
            } catch(...)
            {
                throw ParameterFileProcessingException(line_number, "Error processing dihedral");
            }
        }
    }

    /// Improper dihedral section reading
    getline(in_file, line);             /// Read the first line of the improper dihedral section
    line_number++;                      /// Increment line counter
    if(gmml::Trim(line).find("IMPR") != std::string::npos || gmml::Trim(line).find("IMPROPER") != std::string::npos)
    {
        getline(in_file, line);             /// Read the first line of the improper dihedral section
        line_number++;                      /// Increment line counter
        while (!gmml::Trim(line).empty())         /// Reading until a blank line which is the end of the section
        {
            try
            {
                ProcessImproperDihedral(line, line_number, in_file);        /// Processing improper dihedral line
                getline(in_file, line);                                     /// Read the next line
                line_number++;                                              /// Increment line counter
            } catch(...)
            {
                throw ParameterFileProcessingException(line_number, "Error processing improper dihedral");
            }
        }
    }

    /// Potential parameter section reading
    getline(in_file, line);         /// Read the first line of the potential parameter section
    line_number++;                  /// Increment line counter
    if(gmml::Trim(line).find("NONB") != std::string::npos || gmml::Trim(line).find("NONBON") != std::string::npos)
    {
        getline(in_file, line);             /// Read the first line of the improper dihedral section
        line_number++;                      /// Increment line counter
        while (!gmml::Trim(line).empty())         /// Reading until a blank line which is the end of the section
        {
            try
            {
                ProcessPotentialParameter(line);        /// Processing potential parameter line
                getline(in_file, line);                 /// Read the next line
                line_number++;                          /// Increment line counter
            } catch(...)
            {
                throw ParameterFileProcessingException(line_number, "Error processing potential parameters");
            }
        }
    }
}

void ParameterFile::ReadIonicModifiedParameter(std::ifstream& in_file)
{
    std::string line;
    int line_number = 0;

    /// Unable to read file
    if (!getline(in_file, line))
    {
        throw ParameterFileProcessingException("Error reading file");
    }

    /// Set the title of the parameter file
    title_ = gmml::Trim(line);
    line_number++;

    /// Atom type section reading
    getline(in_file, line);             /// Read the first line of the atom type section
    line_number++;                      /// Increment line counter
    if(gmml::Trim(line).find("MASS") != std::string::npos)
    {
        getline(in_file, line);             /// Read the first line of the atom type section
        line_number++;                      /// Increment line counter
        while (!gmml::Trim(line).empty())         /// Reading until a blank line which is the end of the section
        {
            try
            {
                ProcessAtomType(line);      /// Processing atom type line
                getline(in_file,line);      /// Read the next line
                line_number++;              /// Increment line counter
            } catch(...)
            {
                throw ParameterFileProcessingException(line_number, "Error processing atom type");
            }
        }
    }

    /// Potential parameter section reading
    getline(in_file, line);         /// Read the first line of the potential parameter section
    line_number++;                  /// Increment line counter
    if(gmml::Trim(line).find("NONB") != std::string::npos || gmml::Trim(line).find("NONBON") != std::string::npos)
    {
        getline(in_file, line);             /// Read the first line of the improper dihedral section
        line_number++;                      /// Increment line counter
        while (!gmml::Trim(line).empty())         /// Reading until a blank line which is the end of the section
        {
            try
            {
                ProcessPotentialParameter(line);        /// Processing potential parameter line
                getline(in_file, line);                 /// Read the next line
                line_number++;                          /// Increment line counter
            } catch(...)
            {
                throw ParameterFileProcessingException(line_number, "Error processing potential parameters");
            }
        }
    }
}

/// Process the atom type lines of the parameter file
void ParameterFile::ProcessAtomType(const std::string& line)
{
    double mass, polarizability;
    std::string type, dscr;

    std::istringstream in(line);
    switch(this->file_type_)
    {
        case gmml::MAIN:
        case gmml::MODIFIED:
            in >> std::setw(2) >> type                          /// Extract type from the line
               >> std::setw(10) >> mass                         /// Extract mass from the line
               >> std::setw(10) >> polarizability;              /// Extract polarizability from the line
            break;
        case gmml::IONICMOD:
            in >> std::setw(3) >> type                          /// Extract type from the line
               >> std::setw(10) >> mass                         /// Extract mass from the line
               >> std::setw(10) >> polarizability;              /// Extract polarizability from the line
            break;
    }

    if(polarizability == 0)
        polarizability = gmml::dNotSet;
    gmml::Trim(type);
    if (line.size() > 24)               /// Line has description
        dscr = line.substr(24);

    /// atom_types_.insert(std::pair<std::string, ParameterFileAtom*>(type, new ParameterFileAtom(type, mass, polarizability, dscr)));
    atom_types_[type] = new ParameterFileAtom(type, mass, polarizability, dscr);
}

/// Process the hydrophilic atom type line of the parameter file
void ParameterFile::ProcessHydrophilicAtomType(const std::string& line)
{
    std::string type;
    std::istringstream in(line);             /// Create an stream from the read line
    while (in >> std::setw(4) >> type && !gmml::Trim(type).empty())         /// Iterate on the tokens in the read line
    {
        if(distance(atom_types_.begin(), atom_types_.find(type)) >= 0)         /// Check for the existing atom type in the map
        {
            atom_types_[type] -> SetIsHydrophilic(true);            /// Set is_hydrophilic_ attribute to true
        }
    }
}

/// Process the bond lines of the parameter file
void ParameterFile::ProcessBond(const std::string& line)
{
    char c;
    std::vector<std::string> types(2);
    std::string dscr;
    double force_constant, length;

    std::istringstream in(line);                 /// Create an stream from the read bond line
    in >> std::setw(2) >> types[0] >> c
       >> std::setw(2) >> types[1];         /// Tokenize the bond atom types by '-'
    in >> std::setw(10) >> force_constant
       >> length;                           /// Tokenize the rest of the line into the corresponding variables

    if (in.fail())                  /// Invalid template of the read line
        throw std::exception();
    if (line.size() > 26)           /// Line has description
        dscr = line.substr(26);

    for(unsigned int i = 0; i < types.size(); i++)
        gmml::Trim(types[i]);
    bonds_[types] = new ParameterFileBond(types, force_constant, length, dscr);          /// Create a new bond and insert into the bond list
}

/// Process the angle lines of the parameter file
void ParameterFile::ProcessAngle(const std::string& line)
{
    char c;
    std::vector<std::string> types(3);
    std::string dscr;
    double force_constant, angle;

    std::istringstream in(line);                 /// Create an stream from the read bond line
    in >> std::setw(2) >> types[0] >> c
       >> std::setw(2) >> types[1] >> c
       >> std::setw(2) >> types[2];         /// Tokenize the angle atom types by '-'
    in >> std::setw(10) >> force_constant
       >> std::setw(10) >> angle;           /// Tokenize the rest of the line into the corresponding variables

    if (in.fail())                          /// Invalid template of the read line
        throw std::exception();
    if (line.size() > 29)                   /// Line has description
        dscr = line.substr(29);

    for(unsigned int i = 0; i < types.size(); i++)
        gmml::Trim(types[i]);
    angles_[types] = new ParameterFileAngle(types, force_constant, angle, dscr);         /// Create a new angle and insert into the angle list
}

/// Process the dihedral lines of the parameter file
void ParameterFile::ProcessDihedral(std::string &line, int &line_number, std::ifstream &in_file)
{
    char c;
    std::vector<std::string> types(4);
    ParameterFileDihedralTerm t;
    std::vector<ParameterFileDihedralTerm> terms;
    std::string dscr;
    double scee, scnb, temp_force_constant, temp_phase, temp_periodicity;
    int temp_factor;

    std::istringstream in(line);                 /// Create an stream from the read bond line
    in >> std::setw(2) >> types[0] >> c
       >> std::setw(2) >> types[1] >> c
       >> std::setw(2) >> types[2] >> c
       >> std::setw(2) >> types[3];         /// Tokenize the angle atom types by '-'
    in >> std::setw(4) >> temp_factor
       >> std::setw(15) >> temp_force_constant
       >> std::setw(15) >> temp_phase
       >> std::setw(15) >> temp_periodicity;  /// Tokenize the rest of the line into the corresponding variables

    t.SetFactor(temp_factor);
    t.SetForceConstant(temp_force_constant);
    t.SetPhase(temp_phase);
    t.SetPeriodicity(temp_periodicity);

    if (in.fail())
        throw std::exception();

    if (line.size() > 60)                   /// Line has description
        dscr = line.substr(60);

    t.SetDscr(dscr);
    terms.push_back(t);

    scee = ProcessDoubleDihedralDescription(dscr, "SCEE");              /// Extract scee from the description column of the line
    scnb = ProcessDoubleDihedralDescription(dscr, "SCNB");              /// Extract scnb from the description column of the line

    for(unsigned int i = 0; i < types.size(); i++)
        gmml::Trim(types[i]);
    ParameterFileDihedral *dihedral = new ParameterFileDihedral(types, t, scee, scnb);
    while (dihedral->GetTerms().at(dihedral->GetTerms().size() - 1).GetPeriodicity() < 0)       /// Processing the following lines with the same dihedral;
        /// While the periodicity is negative the following lines are the same dihedrals with different attributes
    {
        getline(in_file, line);                     /// Read the following line with the same atom types in the dihedral
        line_number++;                              /// Increment the line counter
        ParameterFileDihedralTerm new_term;

        std::istringstream in2(line.substr(11));         /// Skip the first column of the line which is the same dihedral atom types
        in2 >> std::setw(4) >> temp_factor
            >> std::setw(15) >> temp_force_constant
            >> std::setw(15) >> temp_phase
            >> std::setw(15) >> temp_periodicity;      /// Tokenize the rest of the line into the corresponding variables

        new_term.SetFactor(temp_factor);
        new_term.SetForceConstant(temp_force_constant);
        new_term.SetPhase(temp_phase);
        new_term.SetPeriodicity(temp_periodicity);

        if (in2.fail())                             /// Ivalid entry
        {
            throw ParameterFileProcessingException(line_number, "Error processing dihedral term");
        }

        if (line.size() > 60)                   /// Line has description
            dscr = line.substr(60);

        new_term.SetDscr(dscr);
        terms.push_back(new_term);

        dihedral->AddTerm(new_term);
    }

    std::vector<std::string> inverse_types(4);                /// Create an inverse std::vector of atom types for duplicate checking
    for(unsigned int i = 0; i < types.size(); i++)
    {
        inverse_types[i] = types[types.size() - (i+1)];
    }

    if(dihedrals_.find(types) == dihedrals_.end() && dihedrals_.find(inverse_types) == dihedrals_.end())     /// Duplicate checking
    {
        dihedral->SetIsImproper(false);                                 /// Set is_improper_ attribute to false; improper dihedral doesn't have factor
        if (dihedral->GetTypes().at(0) == "X" || dihedral->GetTypes().at(1) == "X" || dihedral->GetTypes().at(2) == "X" || dihedral->GetTypes().at(3) == "X")   ///  Check for generic dihedral
        {
            dihedral->SetIsGeneric(true);
        }
        else
        {
            dihedral->SetIsGeneric(false);
        }
        dihedrals_[dihedral->GetTypes()] = dihedral;
    }
    else
    {
        throw ParameterFileProcessingException(line_number, "Duplicate dihedral entry");
    }
}

/// Process the improper dihedral lines of the parameter file
void ParameterFile::ProcessImproperDihedral(std::string &line, int &line_number, std::ifstream &in_file)
{
    char c;
    std::vector<std::string> types(4);
    ParameterFileDihedralTerm t;
    std::vector<ParameterFileDihedralTerm> terms;
    std::string dscr;
    double scee, scnb, temp_force_constant, temp_phase, temp_periodicity;

    std::istringstream in(line);                     /// Create an stream from the read bond line
    in >> std::setw(2) >> types[0] >> c
       >> std::setw(2) >> types[1] >> c
       >> std::setw(2) >> types[2] >> c
       >> std::setw(2) >> types[3];             /// Tokenize the angle atom types by '-'
    in >> std::setw(15) >> temp_force_constant
       >> std::setw(15) >> temp_phase
       >> std::setw(15) >> temp_periodicity;      /// Tokenize the rest of the line into the corresponding variables

    t.SetForceConstant(temp_force_constant);
    t.SetPhase(temp_phase);
    t.SetPeriodicity(temp_periodicity);

    if (in.fail())                              /// Invalid entry
        throw std::exception();

    t.SetFactor(gmml::iNotSet);                        /// Improper dihedral doesn't have factor

    if (line.size() > 60)                       /// Line has description
        dscr = line.substr(60);

    t.SetDscr(dscr);
    terms.push_back(t);

    scee = ProcessDoubleDihedralDescription(dscr, "SCEE");              /// Extract scee from the description column of the line
    scnb = ProcessDoubleDihedralDescription(dscr, "SCNB");              /// Extract scnb from the description column of the line

    for(unsigned int i = 0; i < types.size(); i++)
        gmml::Trim(types[i]);
    ParameterFileDihedral *dihedral = new ParameterFileDihedral(types, t, scee, scnb);
    while (dihedral->GetTerms().at(dihedral->GetTerms().size() - 1).GetPeriodicity() < 0)       /// Processing the following lines with the same dihedral;
        /// While the periodicity is negative the following lines are the same dihedrals with different attributes
    {
        getline(in_file, line);                     /// Read the following line with the same atom types in the dihedral
        line_number++;                              /// Increment the line counter
        ParameterFileDihedralTerm new_term;

        std::istringstream in2(line.substr(11));             /// Skip the first column of the line which is the same dihedral atom types
        in2 >> std::setw(15) >> temp_force_constant
            >> std::setw(15) >> temp_phase
            >> std::setw(15) >> temp_periodicity;  /// Tokenize the rest of the line into the corresponding variables

        new_term.SetForceConstant(temp_force_constant);
        new_term.SetPhase(temp_phase);
        new_term.SetPeriodicity(temp_periodicity);

        if (in2.fail())
        {
            throw ParameterFileProcessingException(line_number, "Error processing dihedral term");
        }
        if (line.size() > 60)                   /// Line has description
            dscr = line.substr(60);

        new_term.SetDscr(dscr);
        terms.push_back(new_term);

        dihedral->AddTerm(new_term);
    }

    std::vector<std::string> inverse_types(4);                /// Create an inverse std::vector of atom types for duplicate checking
    for(unsigned int i = 0; i < types.size(); i++)
    {
        inverse_types[i] = types[types.size() - (i+1)];
    }

    if(dihedrals_.find(types) == dihedrals_.end() && dihedrals_.find(inverse_types) == dihedrals_.end())
    {
        dihedral->SetIsImproper(true);                                 /// Set is_improper_ attribute to true; improper dihedral doesn't have factor
        if (dihedral->GetTypes().at(0) == "X" || dihedral->GetTypes().at(1) == "X" || dihedral->GetTypes().at(2) == "X" || dihedral->GetTypes().at(3) == "X")   ///  Check for generic dihedral
        {
            dihedral->SetIsGeneric(true);
        }
        else
        {
            dihedral->SetIsGeneric(false);
        }
        dihedrals_[dihedral->GetTypes()] = dihedral;
    }
}

/// Process the hydrogen-bond lines of the parameter file
void ParameterFile::ProcessHydrogenBond(const std::string& line)
{
    std::vector<std::string> types(2);
    std::vector<double> coefficients(2);
    std::string dscr;

    std::istringstream in(line);                         /// Create an stream from the read bond line
    in >> types[0] >> types[1]
       >> coefficients[0] >> coefficients[1];       /// Tokenize the read line to the types and coefficients

    if (in.fail())                                  /// Invalid entry
        throw std::exception();

    if (line.size() > 58)                           /// Line has description
        dscr = line.substr(58);

    for(unsigned int i = 0; i < types.size(); i++)
        gmml::Trim(types[i]);

    std::vector<std::string> inverse_types(2);                /// Create an inverse std::vector of types
    inverse_types[0] = types[1];
    inverse_types[1] = types[0];

    std::vector<double> inverse_coefficients(2);
    inverse_coefficients[0] = coefficients[1];
    inverse_coefficients[1] = coefficients[0];

    if(bonds_.find(types) != bonds_.end())          /// Check for existing bond
    {
        bonds_[types] ->  SetHbondCoefficients(coefficients);   /// Update hydrogen-bond attribute of the bond
    }
    else if (bonds_.find(inverse_types) != bonds_.end())        /// Check for existing inverse bond
    {
        bonds_[inverse_types] -> SetHbondCoefficients(inverse_coefficients);    /// Update hydrogen-bond attribute of the bond
    }
}

/// Process the lines of equivalent atom types of the parameter file
void ParameterFile::ProcessEquivalentSymbols(const std::string& line)
{
    /// Old style: Create an equivalent list for the first atom in the line
    /**
    std::string type, t;
    std::istringstream in(line);         // Create an stream from the read line
    in >> std::setw(4) >> type;     // Extract the first atom type from the line which the equivalent list is created for
    gmml::Trim(type);
    while (in >> std::setw(4) >> t && !gmml::Trim(t).empty())        // Read until the end of the line
    {
        if(atom_types_.find(type) != atom_types_.end())         // Check for existing atom type in the atom type map
            atom_types_[type] -> equivalent_list_.push_back(t); // Add the equivalent atom type to the corresponding list of the found atom
    }
    **/

    /// New style: Create an equivalent list for all atoms in the line
    std::string type;
    std::istringstream in(line);         /// Create an stream from the read line
    std::vector<std::string> types;

    while (in >> std::setw(4) >> type && !gmml::Trim(type).empty())     /// Read until the end of the line
    {
        types.push_back(type);;     /// Add the read atom type into the temporary atom list
    }

    for(unsigned int i = 0; i < types.size(); i++)       /// For each atom type in the list create a list of equivalent atom types and update the corresponding attribute
    {
        if(distance(atom_types_.begin(), atom_types_.find(types[i])) >= 0)     /// Check for the existing atom type in the atom type map
        {
            std::vector<std::string> equivalent_types = types;
            equivalent_types.erase(equivalent_types.begin() + i);           /// Remove the base atom from the list
            atom_types_[types[i]] -> SetEquivalentList(equivalent_types);   /// Assign the equivalent list to the corresponding atom in the map
        }
    }
}

/// Process the potential parameter lines of the parameter file: Furthur information for atom types
void ParameterFile::ProcessPotentialParameter(const std::string& line)
{
    std::string type, dscr;
    double radius, depth;

    std::istringstream in(line);             /// Create an stream from the read line
    in >> type >> radius >> depth;      /// Tokenize the read line into the type, radius and well depth

    if (in.fail())                      /// Invalid entry
        throw std::exception();

    if (line.size() > 38)               /// Line has description
        dscr = line.substr(38);

    gmml::Trim(type);                         /// Remove spaces from the begining and the end of the std::string
    if (atom_types_.find(type) == atom_types_.end())        /// Check for the existing atom type in the map
    {
        atom_types_[type] = new ParameterFileAtom(type, gmml::dNotSet, gmml::dNotSet, radius, depth, dscr); /// Create a new entry in the map for a non-existing atom type in the map
    }
    else
    {
        atom_types_[type]->SetRadius(radius) ;            /// Update radius attribute of the existing atom type in the map
        atom_types_[type]->SetWellDepth(depth);         /// Update well depth attribute of the existing atom type in the map
        atom_types_[type]->SetMod4Dscr(dscr);           /// Update mod4 description attribute of the existing atom type in the map
    }

    std::vector<std::string>::const_iterator it;
    std::vector<std::string> equivalent_atoms;
    if (atom_types_.find(type) != atom_types_.end())    /// Check if the atom type in the map has equivalent list
        equivalent_atoms = atom_types_[type] -> GetEquivalentList();
    else
        return;

    for (it = equivalent_atoms.begin(); it != equivalent_atoms.end(); ++it)     /// For all atom types in the equivalent list of the atom set the corresponding attributes
    {
        if (atom_types_.find(*it) == atom_types_.end())         /// Check for the non-existing atom type in the map
        {
            atom_types_[*it] = new ParameterFileAtom(*it, gmml::dNotSet, gmml::dNotSet, radius, depth, dscr);       /// Create a new entry in the map for the non-existing atom type in the map
        }
        else
        {
            atom_types_[*it]->SetRadius(radius);         /// Update radius attribute of the existing atom type in the map
            atom_types_[*it]->SetWellDepth(depth);      /// Update well depth attribute of the existing atom type in the map
            atom_types_[*it]->SetMod4Dscr(dscr);        /// Update mod4 description attribute of the existing atom type in the map
        }
    }
}

/// Process the description part of the dihedral line to extract a double value with a key reference
double ParameterFile::ProcessDoubleDihedralDescription(const std::string& dscr, const std::string& key)
{
    double val;
    size_t pos = dscr.find(std::string(key + "="));      /// Find the starting position of the value of the given key in the given description
    if (pos == std::string::npos)                        /// Key is not found
        return gmml::dNotSet;
    std::istringstream ss(dscr.substr(pos + key.size() + 1));
    ss >> val;                                      /// Extract the value of the key from the found position
    if (ss.fail())
        return gmml::dNotSet;
    return val;
}
void ParameterFile::Write(const std::string &parameter_file)
{
    std::ofstream out_file;
    try
    {
        out_file.open(parameter_file.c_str());
    }
    catch(...)
    {
        throw ParameterFileProcessingException(__LINE__,"File could not be created");
    }
    try
    {
        switch(file_type_)
        {
            case gmml::MAIN:
                this->BuildMainParameterFile(out_file);
                break;
            case gmml::MODIFIED:
                this->BuildModifiedParameterFile(out_file);
                break;
            case gmml::IONICMOD:
                this->BuildIonicModifiedParameterFile(out_file);
                break;
        }

    }
    catch(...)
    {
        out_file.close();
    }
}
void ParameterFile::BuildMainParameterFile(std::ofstream& out_stream)
{
    out_stream << GetTitle() << std::endl;
    ResolveAtomTypeSection(out_stream);
    ResolveHydrophilicAtomTypeSection(out_stream);
    ResolveBondSection(out_stream);
    ResolveAngleSection(out_stream);
    ResolveDihedralSection(out_stream);
    ResolveImproperDihedralSection(out_stream);
    ResolveHydrogenBondSection(out_stream);
    ResolveEquivalentSymbolsSection(out_stream);
    ResolvePotentialParameterSection(out_stream);
    out_stream << "END";
}

void ParameterFile::BuildModifiedParameterFile(std::ofstream& out_stream)
{
    out_stream << GetTitle() << std::endl;
    ResolveAtomTypeSection(out_stream);
    ResolveBondSection(out_stream);
    ResolveAngleSection(out_stream);
    ResolveDihedralSection(out_stream);
    ResolveImproperDihedralSection(out_stream);
    ResolvePotentialParameterSection(out_stream);
}

void ParameterFile::BuildIonicModifiedParameterFile(std::ofstream& out_stream)
{
    out_stream << GetTitle() << std::endl;
    ResolveAtomTypeSection(out_stream);
    ResolvePotentialParameterSection(out_stream);
}

void ParameterFile::ResolveAtomTypeSection(std::ofstream& stream)
{
    switch(file_type_)
    {
        case gmml::MODIFIED:
            stream << "MASS" << std::endl;
        case gmml::MAIN:
            for(AtomTypeMap::iterator it = atom_types_.begin(); it != atom_types_.end(); it++)
            {
                ParameterFileAtom* atom = (*it).second;
                stream << std::left << std::setw(2) << atom->GetType() << " " << std::right << std::setw(10) << std::fixed << std::setprecision(2) << atom->GetMass() << " " ;
                if(atom->GetPolarizability() == gmml::dNotSet)
                    stream << std::right << std::setw(10) << " ";
                else
                    stream << std::right << std::setw(10) << std::fixed << std::setprecision(2) << atom->GetPolarizability();
                stream << " " << std::left << atom->GetDscr() << std::endl;
            }
            stream << std::endl;
            break;
        case gmml::IONICMOD:
            for(AtomTypeMap::iterator it = atom_types_.begin(); it != atom_types_.end(); it++)
            {
                ParameterFileAtom* atom = (*it).second;
                stream << std::left << std::setw(3) << atom->GetType() << " " << std::right << std::setw(10) << std::fixed << std::setprecision(2) << atom->GetMass() << " " ;
                if(atom->GetPolarizability() == gmml::dNotSet)
                    stream << std::right << std::setw(10) << " ";
                else
                    stream << std::right << std::setw(10) << std::fixed << std::setprecision(2) << atom->GetPolarizability();
                stream << " " << std::left << atom->GetDscr() << std::endl;
            }
            stream << std::endl;
            break;
    }

}
void ParameterFile::ResolveHydrophilicAtomTypeSection(std::ofstream& stream)
{
    int count = 0;
    const int MAX_IN_LINE = 20;
    for(AtomTypeMap::iterator it1 = atom_types_.begin(); it1 != atom_types_.end(); it1++)
    {
        ParameterFileAtom* atom = (*it1).second;
        if(atom->GetIsHydrophilic())
        {
            stream << std::left << std::setw(2) << atom->GetType() << std::setw(2) << " ";
            count++;
            if(count == MAX_IN_LINE)
            {
                count = 0;
                stream << std::endl;
            }
        }
    }
    if(count < MAX_IN_LINE)
        stream << std::endl;
    stream << std::endl;
}
void ParameterFile::ResolveBondSection(std::ofstream& stream)
{
    switch(file_type_)
    {
        case gmml::MODIFIED:
            stream << "BOND" << std::endl;
        case gmml::MAIN:
            for(BondMap::iterator it2 = bonds_.begin(); it2 != bonds_.end(); it2++)
            {
                std::vector<std::string> atom_types = (*it2).first;
                ParameterFileBond* bond = (*it2).second;
                stream << std::left << std::setw(2) << atom_types.at(0) << "-" << std::left << std::setw(2) << atom_types.at(1) << " " << std::right << std::setw(10) << std::fixed << std::setprecision(2) << bond->GetForceConstant()
                       << " " << std::right << std::setw(10) << std::fixed << std::setprecision(2) << bond->GetLength() << " " << std::left << bond->GetDscr() << std::endl;
            }
            stream << std::endl;
            break;
    }

}
void ParameterFile::ResolveAngleSection(std::ofstream& stream)
{
    switch(file_type_)
    {
        case gmml::MODIFIED:
            stream << "ANGL" << std::endl;
        case gmml::MAIN:
            for(AngleMap::iterator it3 = angles_.begin(); it3 != angles_.end(); it3++)
            {
                std::vector<std::string> angle_types = (*it3).first;
                ParameterFileAngle* angle = (*it3).second;
                stream << std::left << std::setw(2) << angle_types.at(0) << "-" << std::left << std::setw(2) << angle_types.at(1) << "-" << std::left << std::setw(2) << angle_types.at(2)
                       << " " << std::right << std::setw(10) << std::fixed << std::setprecision(2) << angle->GetForceConstant() << " " << std::right << std::setw(10) << std::fixed << std::setprecision(2) << angle->GetAngle()
                       << " " << std::left << angle->GetDscr() << std::endl;
            }
            stream << std::endl;
            break;
    }

}
void ParameterFile::ResolveDihedralSection(std::ofstream& stream)
{
    switch(file_type_)
    {
        case gmml::MODIFIED:
            stream << "DIHE" << std::endl;
        case gmml::MAIN:
            for(DihedralMap::iterator it4 = dihedrals_.begin(); it4 != dihedrals_.end(); it4++)
            {
                std::vector<std::string> dihedral_types = (*it4).first;
                ParameterFileDihedral* dihedral = (*it4).second;
                if(!dihedral->GetIsImproper())
                {
                    std::vector<ParameterFileDihedralTerm> dihedral_terms = dihedral->GetTerms();
                    for(std::vector<ParameterFileDihedralTerm>::iterator it5 = dihedral_terms.begin(); it5 != dihedral_terms.end(); it5++)
                    {
                        ParameterFileDihedralTerm dihedral_term = (*it5);
                        stream << std::left << std::setw(2) << dihedral_types.at(0) << "-" << std::left << std::setw(2) << dihedral_types.at(1) << "-"
                               << std::left << std::setw(2) << dihedral_types.at(2) << "-" << std::left << std::setw(2) << dihedral_types.at(3) << " "
                               << std::right << std::setw(4) << dihedral_term.GetFactor() << " " << std::right << std::setw(10) << std::fixed << std::setprecision(2) << dihedral_term.GetForceConstant()
                               << " " << std::right << std::setw(10) << std::fixed << std::setprecision(2) << dihedral_term.GetPhase() << " " << std::right << std::setw(10) << std::fixed << std::setprecision(2) << dihedral_term.GetPeriodicity()
                               << " " << std::left << dihedral_term.GetDscr() << std::endl;
                    }
                }
            }
            stream << std::endl;
            break;
    }
}
void ParameterFile::ResolveImproperDihedralSection(std::ofstream& stream)
{
    switch(file_type_)
    {
        case gmml::MODIFIED:
            stream << "IMPR" << std::endl;
        case gmml::MAIN:
            for(DihedralMap::iterator it6 = dihedrals_.begin(); it6 != dihedrals_.end(); it6++)
            {
                std::vector<std::string> dihedral_types = (*it6).first;
                ParameterFileDihedral* dihedral = (*it6).second;
                if(dihedral->GetIsImproper())
                {
                    std::vector<ParameterFileDihedralTerm> dihedral_terms = dihedral->GetTerms();
                    for(std::vector<ParameterFileDihedralTerm>::iterator it7 = dihedral_terms.begin(); it7 != dihedral_terms.end(); it7++)
                    {
                        ParameterFileDihedralTerm dihedral_term = (*it7);
                        stream << std::left << std::setw(2) << dihedral_types.at(0) << "-" << std::left << std::setw(2) << dihedral_types.at(1) << "-"
                               << std::left << std::setw(2) << dihedral_types.at(2) << "-" << std::left << std::setw(2) << dihedral_types.at(3)
                               << " " << std::right << std::setw(4) << " " << std::right << std::setw(10) << std::fixed << std::setprecision(2) << dihedral_term.GetForceConstant()
                               << " " << std::right << std::setw(10) << std::fixed << std::setprecision(2) << dihedral_term.GetPhase() << " " << std::right << std::setw(10) << std::fixed << std::setprecision(2) << dihedral_term.GetPeriodicity()
                               << " " << std::left << dihedral_term.GetDscr() << std::endl;
                    }
                }
            }
            stream << std::endl;
            break;
    }
}
void ParameterFile::ResolveHydrogenBondSection(std::ofstream& stream)
{
    for(BondMap::iterator it8 = bonds_.begin(); it8 != bonds_.end(); it8++)
    {
        std::vector<std::string> atom_types = (*it8).first;
        ParameterFileBond* bond = (*it8).second;
        if(bond->GetHbondCoefficients().size() != 0)
        {
            stream << std::left << std::setw(2)  << " " << std::left << std::setw(2) << atom_types.at(0) << std::left << std::setw(2) << " " << std::left << std::setw(2) << atom_types.at(1) << std::left << std::setw(2) << " "
                   << std::right << std::setw(10) << std::fixed << std::setprecision(2) << bond->GetHbondCoefficients().at(0) << " " << std::right << std::setw(10) << std::fixed << std::setprecision(2)
                   << bond->GetHbondCoefficients().at(1) << std::endl;
        }
    }
    stream << std::endl;
}
void ParameterFile::ResolveEquivalentSymbolsSection(std::ofstream& stream)
{
    std::vector<std::string> printed = std::vector<std::string>();
    for(AtomTypeMap::iterator it9 = atom_types_.begin(); it9 != atom_types_.end(); it9++)
    {
        ParameterFileAtom* atom = (*it9).second;
        if(atom->GetEquivalentList().size() != 0)
        {
            if(find(printed.begin(), printed.end(), atom->GetType()) != printed.end())
            {}
            else
            {
                stream << std::left << std::setw(2) << atom->GetType() << std::left << std::setw(2) << " ";
                for(unsigned int i = 0; i < atom->GetEquivalentList().size(); i++)
                {
                    stream << std::left << std::setw(2) << atom->GetEquivalentList().at(i) << std::left << std::setw(2) << " ";
                    printed.push_back(atom->GetEquivalentList().at(i));
                }
                stream << std::endl;
            }
        }
    }
    stream << std::endl;
}
void ParameterFile::ResolvePotentialParameterSection(std::ofstream& stream)
{
    switch(file_type_)
    {
        case gmml::MODIFIED:
            stream << "NONB" << std::endl;
            break;
        case gmml::MAIN:
            stream << "MOD4" << std::endl;
            break;
    }
    for(AtomTypeMap::iterator it10 = atom_types_.begin(); it10 != atom_types_.end(); it10++)
    {
        ParameterFileAtom* atom = (*it10).second;
        if(atom->GetRadius() != gmml::dNotSet || atom->GetWellDepth() != gmml::dNotSet)
        {
            switch(this->file_type_)
            {
                case gmml::MAIN:
                case gmml::MODIFIED:
                    stream << std::left << std::setw(2) << " " << std::left << std::setw(2) << atom->GetType() << std::left << std::setw(6) << " " << std::left << std::setw(2) << " " << " "
                           << std::right << std::setw(10) << std::fixed << std::setprecision(2) << atom->GetRadius() << " "
                           << std::right << std::setw(10) << std::fixed << std::setprecision(2) << atom->GetWellDepth() << std::left << atom->GetMod4Dscr() << std::endl;
                    break;
                case gmml::IONICMOD:
                    stream << std::left << std::setw(2) << " " << std::left << std::setw(3) << atom->GetType() << std::left << std::setw(6) << " " << std::left << std::setw(2) << " " << " "
                           << std::right << std::setw(10) << std::fixed << std::setprecision(2) << atom->GetRadius() << " "
                           << std::right << std::setw(10) << std::fixed << std::setprecision(2) << atom->GetWellDepth() << std::left << atom->GetMod4Dscr() << std::endl;
                    break;
            }

        }
    }
    stream << std::endl;
}

//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
void ParameterFile::Print(std::ostream& out)
{
    out << "Path: " << path_
        << std::endl;
    out << "Title: " << title_
        << std::endl;
    out << std::setw(60) << "******************** Atom Types ********************"
        << std::endl;
    out << std::setw(6) << "Type"
        << std::setw(6) << "Mass"
        << std::setw(16) << "Polarizability"
        << std::setw(8) << "Radius"
        << std::setw(12) << "Well Depth"
        << std::setw(13) << "Hydrophilic"
        << std::setw(60) << "Description"
        << std::setw(60) << "MOD4 Description"
        << std::endl;
    for(ParameterFile::AtomTypeMap::iterator it = atom_types_.begin(); it != atom_types_.end(); it++)
    {
        it->second->Print(out);
    }
    out << std::setw(60) << "******************** Bonds ********************"
        << std::endl;
    out << std::setw(12) << "Bond Types"
        << std::setw(15) << "Force"
        << std::setw(10) << "Length"
        << std::setw(60) << "Description"
        << std::endl;
    for(ParameterFile::BondMap::iterator it = bonds_.begin(); it != bonds_.end(); it ++)
    {
        it->second->Print(out);
    }
    out << std::setw(60) << "******************** Angles ********************"
        << std::endl;
    out << std::setw(16) << "Angle Types"
        << std::setw(15) << "Force"
        << std::setw(10) << "Angle"
        << std::setw(60) << "Description"
        << std::endl;
    for(ParameterFile::AngleMap::iterator it = angles_.begin(); it != angles_.end(); it ++)
    {
        it->second->Print(out);
    }
    out << std::setw(60) << "******************** Dihedrals ********************"
        << std::endl;
    out << std::setw(22) << "Dihedral Types"
        << std::setw(10) << "Generic"
        << std::setw(10) << "Improper"
        << std::setw(6) << "SCEE"
        << std::setw(6) << "SCNB"
        << std::setw(10) << "Factor"
        << std::setw(10) << "Force"
        << std::setw(10) << "Phase"
        << std::setw(13) << "Periodicity"
        << std::setw(60) << "Description"
        << std::endl;
    for(ParameterFile::DihedralMap::iterator it = dihedrals_.begin(); it != dihedrals_.end(); it ++)
    {
        it->second->Print(out);
    }
}
