#include "includes/ParameterSet/parameterManager.hpp"
#include "includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp"
#include "includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp"
#include "includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"

#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/directories.hpp"
#include "includes/CodeUtils/logging.hpp"

using parameters::Manager;

Manager::Manager()
{
    // How exactly this happens can be improved, but the information should only ever be loaded into gmml in one place.
    // Find $GMMLHOME
    std::string gmmlHomeDir = codeUtils::getGmmlHomeDir();
    gmml::log(__LINE__, __FILE__, gmml::INF, "gmmlhome is: " + codeUtils::getGmmlHomeDir());
    // Library files of 3D structures with parameters for simulations.
    // DUCT_TAPE until I get a generic library templates class thingy
    //    std::vector<LibraryFileSpace::LibraryFile> libFiles;
    //    std::vector<PrepFileSpace::PrepFile> prepFiles;

    std::vector<std::string> libFiles, prepFiles;
    prepFiles.push_back(gmmlHomeDir + "/dat/prep/GLYCAM_06j-1_GAGS.prep");
    gmml::log(__LINE__, __FILE__, gmml::INF, "Initializing PrepResidueMap");
    this->InitializePrepResidueMap(prepFiles);
    libFiles.push_back(gmmlHomeDir + "/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_amino_06j_12SB.lib");
    libFiles.push_back(gmmlHomeDir + "/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_aminoct_06j_12SB.lib");
    libFiles.push_back(gmmlHomeDir + "/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_aminont_06j_12SB.lib");
    gmml::log(__LINE__, __FILE__, gmml::INF, "Initializing Glycam Residue names");
    this->SetGlycamResidueNames(libFiles, prepFiles);
    libFiles.push_back(gmmlHomeDir + "/dat/CurrentParams/leaprc.ff12SB_2014-04-24/nucleic12.lib");
    libFiles.push_back(gmmlHomeDir + "/dat/CurrentParams/leaprc.ff12SB_2014-04-24/nucleic12.lib");
    libFiles.push_back(gmmlHomeDir + "/dat/CurrentParams/other/solvents.lib");
    libFiles.push_back(gmmlHomeDir + "/dat/CurrentParams/leaprc.ff12SB_2014-04-24/amino12.lib");
    libFiles.push_back(gmmlHomeDir + "/dat/CurrentParams/leaprc.ff12SB_2014-04-24/aminoct12.lib");
    libFiles.push_back(gmmlHomeDir + "/dat/CurrentParams/leaprc.ff12SB_2014-04-24/aminont12.lib");
    gmml::log(__LINE__, __FILE__, gmml::INF, "Initializing LibResidueMap");
    this->InitializeLibResidueMap(libFiles);
    gmml::log(__LINE__, __FILE__, gmml::INF, "Initializing silly maps");
    this->InitializeMaps();
    gmml::log(__LINE__, __FILE__, gmml::INF, "Finished construction.");
}

// This is real silly now
void Manager::SetGlycamResidueNames(std::vector<std::string>& libFiles, std::vector<std::string>& prepFiles)
{
    for (auto& libFileName : libFiles)
    {
        LibraryFileSpace::LibraryFile libFile(libFileName);
        std::vector<std::string> libResNames = libFile.GetAllResidueNames();
        glycamResidueNames_.insert(glycamResidueNames_.begin(), libResNames.begin(), libResNames.end());
    }
    return;
}

// Kept to maintain compatibility with Preprocessor for now. Need to refactor Preprocessor to not be crazy.
void Manager::InitializeMaps() // This are agnostic of prep and lib I guess. Should not be necessary, but using them for
                               // now until I can change elsewhere.
{
    for (auto& resMap : this->GetLibraryResidueMap())
    {
        residueNamesToSelfMap_[resMap.first] =
            (resMap.first); // This is actually used in the code. If it's still here I gave up.
        LibraryFileSpace::LibraryFileResidue* libResidue = resMap.second;
        for (auto& libAtom : libResidue->GetAtomsVector())
        {
            residueNamesToTheirAtomsNameMap_[resMap.first].push_back(libAtom->GetName());
        }
    }
    for (auto& resMap : this->GetPrepResidueMap())
    {
        residueNamesToSelfMap_[resMap.first] =
            (resMap.first); // This is actually used in the code. If it's still here I gave up.
        PrepFileSpace::PrepFileResidue* prepResidue = resMap.second;
        for (auto& prepAtom : prepResidue->GetAtoms())
        {
            residueNamesToTheirAtomsNameMap_[resMap.first].push_back(prepAtom->GetName());
        }
    }
    return;
}

void Manager::InitializePrepResidueMap(std::vector<std::string>& prepFiles)
{

    for (auto& prepFileName : prepFiles)
    {
        PrepFileSpace::PrepFile prepFile(prepFileName);
        std::map<std::string, PrepFileSpace::PrepFileResidue*> residueMap = prepFile.GetResidues();
        prepResMap_.merge(residueMap); // @suppress("Method cannot be resolved")
    }
    return;
}

void Manager::InitializeLibResidueMap(std::vector<std::string>& libFiles)
{
    for (auto& libFileName : libFiles)
    {
        LibraryFileSpace::LibraryFile libFile(libFileName);
        std::map<std::string, LibraryFileSpace::LibraryFileResidue*> residueMap = libFile.GetResidues();
        libResMap_.merge(residueMap); // @suppress("Method cannot be resolved")
    }
    return;
}

LibraryFileSpace::LibraryFileResidue* Manager::FindLibResidue(std::string residueName)
{
    std::map<std::string, LibraryFileSpace::LibraryFileResidue*>::iterator it = libResMap_.find(residueName);
    if (it != libResMap_.end())
    {
        return it->second;
    }
    //    gmml::log(__LINE__, __FILE__, gmml::WAR, "Did not find a LibraryFileSpace::LibraryFileResidue in my map for "
    //    + residueName );
    return nullptr;
}

PrepFileSpace::PrepFileResidue* Manager::FindPrepResidue(std::string residueName)
{
    std::map<std::string, PrepFileSpace::PrepFileResidue*>::iterator it = prepResMap_.find(residueName);
    if (it != prepResMap_.end())
    {
        return it->second;
    }
    //    gmml::log(__LINE__, __FILE__, gmml::WAR, "Did not find a PrepFileSpace::PrepFileResidue in my map for " +
    //    residueName );
    return nullptr;
}

double Manager::GetChargeForResidue(std::string residueName)
{
    double charge                                = 0.0;
    LibraryFileSpace::LibraryFileResidue* libRes = this->FindLibResidue(residueName);
    if (libRes == nullptr) // not found
    {
        PrepFileSpace::PrepFileResidue* prepRes = this->FindPrepResidue(residueName);
        if (prepRes == nullptr)
        {
            return charge;
        }
        return prepRes->GetCharge();
    }
    return libRes->GetCharge();
}

std::vector<std::string> Manager::GetAtomNamesForResidue(const std::string& residueName)
{
    LibraryFileSpace::LibraryFileResidue* libRes = this->FindLibResidue(residueName);
    if (libRes == nullptr) // not found
    {
        PrepFileSpace::PrepFileResidue* prepRes = this->FindPrepResidue(residueName);
        if (prepRes == nullptr)
        {
            return std::vector<std::string>(); // it will be empty
        }
        return prepRes->GetAtomNames();
    }
    return libRes->GetAtomNames();
}

std::vector<std::string> Manager::GetHeavyAtomNamesForResidue(const std::string& residueName)
{
    LibraryFileSpace::LibraryFileResidue* libRes = this->FindLibResidue(residueName);
    if (libRes == nullptr) // not found
    {
        PrepFileSpace::PrepFileResidue* prepRes = this->FindPrepResidue(residueName);
        if (prepRes == nullptr)
        {
            return std::vector<std::string>(); // it will be empty
        }
        return prepRes->GetHeavyAtomNames();
    }
    return libRes->GetHeavyAtomNames();
}

//
// gmml::ResidueNameMap PdbPreprocessor::GetAllResidueNamesFromMultipleLibFilesMap(std::vector<std::string> lib_files)
//{
//    gmml::ResidueNameMap all_residue_names;
//    std::vector<std::string> residue_names;
//    for(std::vector<std::string>::iterator it = lib_files.begin(); it != lib_files.end(); it++)
//    {
//        LibraryFileSpace::LibraryFile* lib_file = new LibraryFileSpace::LibraryFile(*it);
//        residue_names = lib_file->GetAllResidueNames();
//        for(std::vector<std::string>::iterator it1 = residue_names.begin(); it1 != residue_names.end(); it1++)
//        {
//            std::string residue_name = (*it1);
//            all_residue_names[residue_name] = (residue_name);
//        }
//    }
//    return all_residue_names;
//}
//
//
//
// LibraryFileSpace::LibraryFile::ResidueMap
// PdbPreprocessor::GetAllResiduesFromMultipleLibFilesMap(std::vector<std::string> lib_files)
//{
//    LibraryFileSpace::LibraryFile::ResidueMap all_residues;
//    LibraryFileSpace::LibraryFile::ResidueMap residues;
//    for(std::vector<std::string>::iterator it = lib_files.begin(); it != lib_files.end(); it++)
//    {
//        LibraryFileSpace::LibraryFile* lib_file = new LibraryFileSpace::LibraryFile(*it);
//        residues = lib_file->GetResidues();
//        for(LibraryFileSpace::LibraryFile::ResidueMap::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
//        {
//            std::string lib_residue_name = (*it1).first;
//            all_residues[lib_residue_name] = (*it1).second;
//        }
//    }
//    return all_residues;
//}
// PrepFileSpace::PrepFile::ResidueMap PdbPreprocessor::GetAllResiduesFromMultiplePrepFilesMap(std::vector<std::string>
// prep_files)
//{
//    PrepFileSpace::PrepFile::ResidueMap all_residues;
//    PrepFileSpace::PrepFile::ResidueMap residues;
//    for(std::vector<std::string>::iterator it = prep_files.begin(); it != prep_files.end(); it++)
//    {
//        PrepFileSpace::PrepFile* prep_file = new PrepFileSpace::PrepFile(*it);
//        residues = prep_file->GetResidues();
//        for(PrepFileSpace::PrepFile::ResidueMap::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
//        {
//            std::string prep_residue_name = (*it1).first;
//            all_residues[prep_residue_name] = (*it1).second;
//        }
//    }
//    return all_residues;
//}
//
// std::vector<std::string> PdbPreprocessor::GetAllResidueNamesFromMultiplePrepFiles(std::vector<std::string>
// prep_files)
//{
//    std::vector<std::string> all_residue_names;
//    std::vector<std::string> residue_names;
//    for(std::vector<std::string>::iterator it = prep_files.begin(); it != prep_files.end(); it++)
//    {
//        PrepFileSpace::PrepFile* prep_file = new PrepFileSpace::PrepFile(*it);
//        residue_names = prep_file->GetAllResidueNames();
//        for(std::vector<std::string>::iterator it1 = residue_names.begin(); it1 != residue_names.end(); it1++)
//        {
//            all_residue_names.push_back(*it1);
//        }
//    }
//    return all_residue_names;
//}
//
// gmml::ResidueNameMap PdbPreprocessor::GetAllResidueNamesFromMultiplePrepFilesMap(std::vector<std::string> prep_files)
//{
//    gmml::ResidueNameMap all_residue_names;
//    std::vector<std::string> residue_names;
//    for(std::vector<std::string>::iterator it = prep_files.begin(); it != prep_files.end(); it++)
//    {
//        PrepFileSpace::PrepFile* prep_file = new PrepFileSpace::PrepFile(*it);
//        residue_names = prep_file->GetAllResidueNames();
//        for(std::vector<std::string>::iterator it1 = residue_names.begin(); it1 != residue_names.end(); it1++)
//        {
//            std::string residue_name = (*it1);
//            all_residue_names[residue_name] = residue_name;
//        }
//    }
//    return all_residue_names;
//}
//
// std::vector<std::string> PdbPreprocessor::SetAllResidueNamesFromDatasetFiles(std::vector<std::string> lib_files,
// std::vector<std::string> prep_files)
//{
//    std::vector<std::string> all_lib_residue_names = GetAllResidueNamesFromMultipleLibFiles(lib_files);
//    std::vector<std::string> all_prep_residue_names = GetAllResidueNamesFromMultiplePrepFiles(prep_files);
//    std::vector<std::string> all_residue_names;
//    for(std::vector<std::string>::iterator it1 = all_lib_residue_names.begin(); it1 != all_lib_residue_names.end();
//    it1++)
//    {
//        all_residue_names.push_back(*it1);
//    }
//    for(std::vector<std::string>::iterator it1 = all_prep_residue_names.begin(); it1 != all_prep_residue_names.end();
//    it1++)
//    {
//        all_residue_names.push_back(*it1);
//    }
//    return all_residue_names;
//}
//

//
// gmml::ResidueNameAtomNamesMap
// PdbPreprocessor::GetAllAtomNamesOfResidueNamesFromMultipleLibFiles(std::vector<std::string> lib_files)
//{
//    gmml::ResidueNameAtomNamesMap residue_atom_map = gmml::ResidueNameAtomNamesMap();
//    LibraryFileSpace::LibraryFile::ResidueMap residues;
//    for(std::vector<std::string>::iterator it = lib_files.begin(); it != lib_files.end(); it++)
//    {
//        LibraryFileSpace::LibraryFile* lib_file = new LibraryFileSpace::LibraryFile(*it);
//        residues = lib_file->GetResidues();
//        for(LibraryFileSpace::LibraryFile::ResidueMap::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
//        {
//            std::string residue_name = (*it1).first;
//            LibraryFileSpace::LibraryFileResidue* residue = (*it1).second;
//            LibraryFileSpace::LibraryFileResidue::AtomMap atoms = residue->GetAtoms();
//            for(LibraryFileSpace::LibraryFileResidue::AtomMap::iterator it2 = atoms.begin(); it2 != atoms.end();
//            it2++)
//            {
//                LibraryFileSpace::LibraryFileAtom* atom = (*it2).second;
//                std::string atom_name = atom->GetName();
//                residue_atom_map[residue_name].push_back(atom_name);
//            }
//        }
//    }
//    return residue_atom_map;
//}
//
// gmml::ResidueNameAtomNamesMap
// PdbPreprocessor::GetAllAtomNamesOfResidueNamesFromMultiplePrepFiles(std::vector<std::string> prep_files)
//{
//    gmml::ResidueNameAtomNamesMap residue_atom_map = gmml::ResidueNameAtomNamesMap();
//    for(std::vector<std::string>::iterator it = prep_files.begin(); it != prep_files.end(); it++)
//    {
//        PrepFileSpace::PrepFile* prep_file = new PrepFileSpace::PrepFile(*it);
//        PrepFileSpace::PrepFile::ResidueMap residues = prep_file->GetResidues();
//        for(PrepFileSpace::PrepFile::ResidueMap::iterator it1 = residues.begin(); it1 != residues.end(); it1++)
//        {
//            std::string residue_name = (*it1).first;
//            PrepFileSpace::PrepFileResidue* residue = (*it1).second;
//            PrepFileSpace::PrepFileAtomVector atoms = residue->GetAtoms();
//            for(PrepFileSpace::PrepFileAtomVector::iterator it2 = atoms.begin(); it2 != atoms.end(); it2++)
//            {
//                PrepFileSpace::PrepFileAtom* atom = (*it2);
//                std::string atom_name = atom->GetName();
//                residue_atom_map[residue_name].push_back(atom_name);
//            }
//        }
//    }
//    return residue_atom_map;
//}
//
// gmml::ResidueNameAtomNamesMap PdbPreprocessor::GetAllAtomNamesOfResidueNamesFromDatasetFiles(std::vector<std::string>
// lib_files, std::vector<std::string> prep_files)
//{
//    gmml::ResidueNameAtomNamesMap all_residue_atom_map = gmml::ResidueNameAtomNamesMap();
//    gmml::ResidueNameAtomNamesMap all_residue_atom_map_from_lib =
//    GetAllAtomNamesOfResidueNamesFromMultipleLibFiles(lib_files); gmml::ResidueNameAtomNamesMap
//    all_residue_atom_map_from_prep = GetAllAtomNamesOfResidueNamesFromMultiplePrepFiles(prep_files);
//    for(gmml::ResidueNameAtomNamesMap::iterator it = all_residue_atom_map_from_lib.begin(); it !=
//    all_residue_atom_map_from_lib.end(); it++)
//    {
//        std::string residue_name = (*it).first;
//        std::vector<std::string> atom_names = (*it).second;
//        all_residue_atom_map[residue_name] = atom_names;
//    }
//    for(gmml::ResidueNameAtomNamesMap::iterator it = all_residue_atom_map_from_prep.begin(); it !=
//    all_residue_atom_map_from_prep.end(); it++)
//    {
//        std::string residue_name = (*it).first;
//        std::vector<std::string> atom_names = (*it).second;
//        all_residue_atom_map[residue_name] = atom_names;
//    }
//    return all_residue_atom_map;
//}
//
//
// LibraryFileSpace::LibraryFileResidue* PdbPreprocessor::GetLibraryResidueByNameFromMultipleLibraryFiles(std::string
// residue_name, std::vector<std::string> lib_files)
//{
//    LibraryFileSpace::LibraryFileResidue* library_residue = new LibraryFileSpace::LibraryFileResidue();
//    for(std::vector<std::string>::iterator it = lib_files.begin(); it != lib_files.end(); it++)
//    {
//        LibraryFileSpace::LibraryFile* lib_file = new LibraryFileSpace::LibraryFile(*it);
//        library_residue = lib_file->GetLibraryResidueByResidueName(residue_name);
//        if(library_residue != NULL)
//            break;
//    }
//    return library_residue;
//}
