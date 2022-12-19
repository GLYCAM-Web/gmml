#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSATOM_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CDSATOM_HPP

#include <string>
#include <iostream>
#include <vector>

#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"
#include "includes/Abstract/absAtom.hpp"

namespace cds
{
class Atom : public Abstract::absAtom, public glygraph::Node<Atom>
{
public:
	//////////////////////////////////////////////////////////
	//                       CONSTRUCTORS                   //
	//////////////////////////////////////////////////////////
	Atom() {
	    //std::cout << "Atom default ctor with name_index: " << this->getName() << "_ " << this->getIndex() << "\n";
	}
	//Atom(Atom* refAtom);
	Atom(const std::string& name, const Coordinate& coord);
	virtual ~Atom() {//std::cout << "cds::Atom default dtor for " << this->getName() << "_" << this->getIndex() << "\n";
	}
    //////////////////////////////////////////////////////////
	//                       ACCESSORS                      //
	//////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////
	//                       MUTATOR                        //
	//////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////
	//                       FUNCTIONS                      //
	//////////////////////////////////////////////////////////
	void addBond(Atom* otherAtom);
	void bondIfClose(Atom* otherAtom);
    bool isWithinBondingDistance(const Atom* otherAtom) const;
    std::string getElement() const;
    int getAtomicNumber() const;
    virtual std::string getId() const;
    std::vector<Coordinate*> getCoordinatesOfNeighbors();
    //////////////////////////////////////////////////////////
    //                   OVERLOADED OPERATORS               //
    //////////////////////////////////////////////////////////
    bool operator== (const Atom &otherAtom);
    bool operator!= (const Atom &otherAtom);
	//////////////////////////////////////////////////////////
	//                       DISPLAY FUNCTION               //
	//////////////////////////////////////////////////////////
    virtual void Print(std::ostream& out) const;
private:
	//////////////////////////////////////////////////////////
	//                       ATTRIBUTES                     //
	//////////////////////////////////////////////////////////

};
} // namespace
#endif // CDSATOM_HPP
