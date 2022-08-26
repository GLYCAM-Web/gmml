#ifndef INCLUDES_CENTRALDATASTRUCTURE_CDSATOM_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_CDSATOM_HPP

#include <string>
#include <iostream>
#include <vector>
#include <ctype.h> // isalpha

#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"
#include "includes/MolecularMetadata/atomicBonds.hpp" // bondIfClose
#include "includes/GeometryTopology/coordinate.hpp"
#include "includes/common.hpp"
#include "includes/Abstract/absAtom.hpp"
#include "includes/CodeUtils/logging.hpp"

namespace cds
{
template <class atomT>
class cdsAtom : public Abstract::absAtom, public glygraph::Node<atomT>
{
public:
	//////////////////////////////////////////////////////////
	//                       CONSTRUCTORS                   //
	//////////////////////////////////////////////////////////
	cdsAtom() {};
	cdsAtom(const std::string& name, const Coordinate& coord);
	//////////////////////////////////////////////////////////
	//                       ACCESSORS                      //
	//////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////
	//                       MUTATOR                        //
	//////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////
	//                       FUNCTIONS                      //
	//////////////////////////////////////////////////////////
	void addBond(atomT* otherAtom);
	void bondIfClose(atomT* otherAtom);
    bool isWithinBondingDistance(const atomT* otherAtom) const;
    const std::string getElement() const;
	//////////////////////////////////////////////////////////
	//                       DISPLAY FUNCTION               //
	//////////////////////////////////////////////////////////
private:
	//////////////////////////////////////////////////////////
	//                       ATTRIBUTES                     //
	//////////////////////////////////////////////////////////
};

//////////////////////////////////////////////////////////
//                    FUNCTIONS                         //
//////////////////////////////////////////////////////////

template <class atomT>
cdsAtom<atomT>::cdsAtom(const std::string& name, const Coordinate& coord)
: Abstract::absAtom(coord)
{
	this->setName(name);
}


template <class atomT>
void cdsAtom<atomT>::addBond(atomT *otherAtom)
{
	this->addNeighbor("bondByDistance", otherAtom);
	return;
}

template <class atomT>
void cdsAtom<atomT>::bondIfClose(atomT* otherAtom)
{
	if (this->isWithinBondingDistance(otherAtom))
	{
		this->addBond(otherAtom);
		//std::stringstream ss;
		//std::cout << "Bonded " << this->getName() << "_" << this->getIndex() << " to " << otherAtom->getName() << "_" << otherAtom->getIndex() << "\n";
		//gmml::log(__LINE__,__FILE__,gmml::INF, ss.str());
	}
	return;
}

template <class atomT>
const std::string cdsAtom<atomT>::getElement() const // derived classes should overwrite if more explicit about element.
{
    std::string name = this->getName();
    if (!name.empty())
    {
        if (isalpha(name.at(0))) // if first char is in the alphabet
        {
             return name.substr(0,1); // return first character as string
        }
    }
    gmml::log(__LINE__,__FILE__,gmml::WAR, "Did not find an element for atom named: " + name);
    return "";
}

template <class atomT>
bool cdsAtom<atomT>::isWithinBondingDistance(const atomT* otherAtom) const
{
    double maxLength = atomicBonds::getMaxBondLengthByAtomType(this->getElement(), otherAtom->getElement());
    if (this->getCoordinate()->withinDistance(otherAtom->getCoordinate(), maxLength))
    {
        return true;
    }
    return false;
}

} // namespace
#endif // CDSATOM_HPP
