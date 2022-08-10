#ifndef GMML_INCLUDES_PARAMETERSET_MANAGER_HPP
#define GMML_INCLUDES_PARAMETERSET_MANAGER_HPP

#include <string>
#include "includes/ParameterSet/LibraryFileSpace/libraryfile.hpp"
#include "includes/ParameterSet/PrepFileSpace/prepfile.hpp"

// Created to serve old gmml's setup. I don't want to pass out any classes like prep file or library file, they should only be used internally
// by a 3dTemplate class.

namespace parameters
{
class Manager
{
public:
    Manager();
    double GetChargeForResidue(std::string residueName);
    LibraryFileSpace::LibraryFileResidue* FindLibResidue(std::string residueName);
    PrepFileSpace::PrepFileResidue* FindPrepResidue(std::string residueName);
    inline std::map<std::string, std::string> GetAllResidueNameMap() {return residueNamesToSelfMap_;}
    inline std::map<std::string, std::vector<std::string>> GetResidueNamesToTheirAtomNamesMap() { return residueNamesToTheirAtomsNameMap_;}
    // Ok don't do this, figure out what they want and implement that query here:
    inline std::map<std::string, LibraryFileSpace::LibraryFileResidue*> GetLibraryResidueMap() { return libResMap_;}
    inline std::map<std::string, PrepFileSpace::PrepFileResidue*> GetPrepResidueMap() { return prepResMap_;}
    inline std::vector<std::string> GetGlycamResidueNames() {return glycamResidueNames_;}
    std::vector<std::string> GetAtomNamesForResidue(const std::string& residueName);
    std::vector<std::string> GetHeavyAtomNamesForResidue(const std::string &residueName);
private:
    void InitializeMaps();
    void InitializePrepResidueMap(std::vector<std::string> &prepFiles);
    void InitializeLibResidueMap(std::vector<std::string> &libFiles);
    void SetGlycamResidueNames(std::vector<std::string> &libFiles, std::vector<std::string> &prepFiles);
    // Attributes
    std::map<std::string, std::string> residueNamesToSelfMap_; // Silly
    std::map<std::string, std::vector<std::string>> residueNamesToTheirAtomsNameMap_;
    std::map<std::string, LibraryFileSpace::LibraryFileResidue*> libResMap_;
    std::map<std::string, PrepFileSpace::PrepFileResidue*> prepResMap_;
    std::vector<std::string> glycamResidueNames_;
};
}
#endif //GMML_INCLUDES_PARAMETERSET_MANAGER_HPP
