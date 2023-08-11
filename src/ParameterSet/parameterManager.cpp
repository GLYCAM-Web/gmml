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
