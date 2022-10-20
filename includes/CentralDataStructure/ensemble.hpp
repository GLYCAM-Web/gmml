#ifndef INCLUDES_CENTRALDATASTRUCTURE_ENSEMBLE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_ENSEMBLE_HPP

#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"
#include <string>
#include <vector>
#include <memory> // unique_ptr

namespace cds
{
class Ensemble : public glygraph::Node<Ensemble>
{
public:
    //////////////////////////////////////////////////////////
    //                    CONSTRUCTOR                       //
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    //                    ACCESSOR                          //
    //////////////////////////////////////////////////////////
    std::vector<Atom*> getAtoms() const;
    std::vector<Residue*> getResidues() const;
    std::vector<Molecule*> getMolecules() const;
    std::vector<Assembly*> getAssemblies() const;
    //////////////////////////////////////////////////////////
    //                    MUTATOR                           //
    //////////////////////////////////////////////////////////
    void addAssembly(std::unique_ptr<Assembly> myAssembly);
    //////////////////////////////////////////////////////////
    //                    FUNCTIONS                         //
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    //                    DISPLAY                           //
    //////////////////////////////////////////////////////////
private:
    //////////////////////////////////////////////////////////
    //                    ATTRIBUTES                        //
    //////////////////////////////////////////////////////////
    std::vector<std::unique_ptr<Assembly>> assemblies_;
};
} // namespace
#endif // INCLUDES_CENTRALDATASTRUCTURE_ENSEMBLE_HPP
