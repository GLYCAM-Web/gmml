#ifndef INCLUDES_CENTRALDATASTRUCTURE_MOLECULE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_MOLECULE_HPP

#include "includes/CentralDataStructure/residue.hpp"
#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"

#include <vector>
#include <memory> // unique_ptr

namespace cds
{
class Molecule : public glygraph::Node<Molecule>
{
public:
    //////////////////////////////////////////////////////////
    //                    CONSTRUCTOR                       //
    //////////////////////////////////////////////////////////
    Molecule() : number_(0) {}//std::cout << "Molecule default ctor\n";}
    virtual ~Molecule() {}//std::cout << "Molecule default dtor for " << this->getName() << "\n";}
    //////////////////////////////////////////////////////////
    //                    ACCESSOR                          //
    //////////////////////////////////////////////////////////
    inline const int& getNumber() {return number_;}
    std::vector<Atom*> getAtoms() const;
    std::vector<Residue*> getResidues() const;
    //////////////////////////////////////////////////////////
    //                    MUTATOR                           //
    //////////////////////////////////////////////////////////
    inline void setNumber(const int i) {number_ = i;}
    //////////////////////////////////////////////////////////
    //                    FUNCTIONS                         //
    //////////////////////////////////////////////////////////
    void addResidue(std::unique_ptr<Residue> myResidue);
    void setResidues(std::vector<std::unique_ptr<Residue>> myResidues);
    Residue* insertNewResidue(std::unique_ptr<Residue> myResidue, const Residue& positionReferenceResidue);
    std::vector<std::unique_ptr<Residue>>::iterator findPositionOfResidue(const Residue* queryResidue);
    std::vector<Residue*> getResidues(std::vector<std::string> queryNames);
    Residue* getResidue(const std::string& queryName);
    void deleteResidue(Residue*);
    void renumberResidues(int newStartNumber = 1);
    //////////////////////////////////////////////////////////
    //                    DISPLAY                           //
    //////////////////////////////////////////////////////////
    void WritePdb(std::ostream& stream) const;
    void WriteOff(std::ostream& stream) const;
private:
    //////////////////////////////////////////////////////////
    //                    FUNCTIONS                         //
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    //                    ATTRIBUTES                        //
    //////////////////////////////////////////////////////////
    std::vector<std::unique_ptr<Residue>> residues_;
    int number_;
};
} // namespace
#endif // INCLUDES_CENTRALDATASTRUCTURE_MOLECULE_HPP
