#include "includes/CentralDataStructure/cdsFunctions/parameters.hpp"
#include "includes/CodeUtils/directories.hpp"

using cdsParameters::ParameterManager;

ParameterManager::ParameterManager()
{     // Library files of 3D structures with parameters for simulations.
    // How exactly this happens can be improved, but the information should only ever be loaded into gmml in one place.
    // Find $GMMLHOME
    std::string gmmlHomeDir = codeUtils::getGmmlHomeDir();
    gmml::log(__LINE__, __FILE__, gmml::INF, "gmmlhome is: " + codeUtils::getGmmlHomeDir());
    prepFiles_.emplace_back(gmmlHomeDir + "/dat/prep/GLYCAM_06j-1_GAGS_KDN.prep");
    gmml::log(__LINE__, __FILE__, gmml::INF, "Initializing ResidueMap with lprepFiles");
    for(auto & file : prepFiles_)
    {
        this->InitializeResidueMap(file.getResidues());
    }
    libFiles_.emplace_back(gmmlHomeDir + "/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_amino_06j_12SB.lib");
    libFiles_.emplace_back(gmmlHomeDir + "/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_aminoct_06j_12SB.lib");
    libFiles_.emplace_back(gmmlHomeDir + "/dat/CurrentParams/leaprc_GLYCAM_06j-1_2014-03-14/GLYCAM_aminont_06j_12SB.lib");
    libFiles_.emplace_back(gmmlHomeDir + "/dat/CurrentParams/leaprc.ff12SB_2014-04-24/nucleic12.lib");
    libFiles_.emplace_back(gmmlHomeDir + "/dat/CurrentParams/leaprc.ff12SB_2014-04-24/nucleic12.lib");
    libFiles_.emplace_back(gmmlHomeDir + "/dat/CurrentParams/other/solvents.lib");
    libFiles_.emplace_back(gmmlHomeDir + "/dat/CurrentParams/leaprc.ff12SB_2014-04-24/amino12.lib");
    libFiles_.emplace_back(gmmlHomeDir + "/dat/CurrentParams/leaprc.ff12SB_2014-04-24/aminoct12.lib");
    libFiles_.emplace_back(gmmlHomeDir + "/dat/CurrentParams/leaprc.ff12SB_2014-04-24/aminont12.lib");
    gmml::log(__LINE__, __FILE__, gmml::INF, "Initializing ResidueMap with libFiles");
    for(auto & file : libFiles_)
    {
        this->InitializeResidueMap(file.getResidues());
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "Finished construction of ParameterManager.");
}

void ParameterManager::InitializeResidueMap(std::vector<cds::Residue*> incomingResidues)
{
    for (auto & residue : incomingResidues)
    {
        parameterResidueMap_[residue->getName()] = residue;
    }
}

void ParameterManager::setAtomCharges(std::vector<cds::Residue*> queryResidues)
{
    for (auto & residue : queryResidues)
    {
        this->setAtomCharges(residue);
    }
}


bool ParameterManager::setAtomCharges(cds::Residue* queryResidue)
{
    bool allAtomsPresent = true;
    if (auto search = parameterResidueMap_.find(queryResidue->GetParmName()); search == parameterResidueMap_.end())
    {
        gmml::log(__LINE__,__FILE__,gmml::WAR, "Did not find parameters and so cannot set charges for residue named: " + queryResidue->GetParmName());
        return false;
    }
    std::vector<cds::Atom*> refAtoms = parameterResidueMap_[queryResidue->GetParmName()]->getAtoms();
    for(auto & queryAtom : queryResidue->getAtoms())
    {
        if(!cdsParameters::setChargeForAtom(queryAtom, refAtoms))
        {
            allAtomsPresent = false;
        }
    }
    return allAtomsPresent;
}


bool cdsParameters::setChargeForAtom(cds::Atom* queryAtom, std::vector<cds::Atom*> referenceAtoms)
{
    for(auto & refAtom : referenceAtoms)
    {
        if (queryAtom->getName() == refAtom->getName())
        {
            queryAtom->setCharge(refAtom->getCharge());
            queryAtom->setType(refAtom->getType());
            return true;
        }
    }
    gmml::log(__LINE__,__FILE__,gmml::WAR, "No charges found for atom named " + queryAtom->getName());
    return false;
}
