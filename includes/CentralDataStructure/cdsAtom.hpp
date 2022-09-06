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
#include "includes/MolecularMetadata/elementattributes.hpp"

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
    std::string getElement() const;
    int getAtomicNumber() const;
	//////////////////////////////////////////////////////////
	//                       DISPLAY FUNCTION               //
	//////////////////////////////////////////////////////////
    void WritePdb(std::ostream &stream,
        		std::string residueName = " ",
    			std::string residueNumber = " ",
    			std::string recordName = "ATOM",
    			std::string chainId = " ",
    			std::string insertionCode = " ",
    			std::string alternativeLocation = " ",
    			std::string occupancy = " ",
    			std::string temperatureFactor = " "
    			) const;
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
std::string cdsAtom<atomT>::getElement() const // derived classes should overwrite if more explicit about element.
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
int cdsAtom<atomT>::getAtomicNumber() const
{
	return MolecularMetadata::findElementAtomicNumber(this->getElement());
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
//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
template <class atomT> // This or just construct a pdbAtom class and write from there?
void cdsAtom<atomT>::WritePdb(std::ostream &stream,
		std::string residueName,
		std::string residueNumber,
		std::string recordName,
		std::string chainId,
		std::string insertionCode,
		std::string alternativeLocation,
		std::string occupancy,
		std::string temperatureFactor
		) const
{
	stream << std::left << std::setw(6) << recordName;
	stream << std::right << std::setw(5) << this->getIndex() << std::left << std::setw(1) << " ";
	stream << std::left << std::setw(4) << this->getName();
	stream << std::left << std::setw(1) << alternativeLocation;
    stream << std::right << std::setw(3) << residueName << std::left << std::setw(1) << " ";
    stream << std::left << std::setw(1) << chainId;
    stream << std::right << std::setw(4) << residueNumber;
    stream << std::left << std::setw(1) << insertionCode << std::left << std::setw(3) << " ";
    stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << this->getCoordinate()->GetX();
    stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << this->getCoordinate()->GetY();
    stream << std::right << std::setw(8) << std::fixed << std::setprecision(3) << this->getCoordinate()->GetZ();
    stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << occupancy;
    stream << std::right << std::setw(6) << std::fixed << std::setprecision(2) << temperatureFactor << std::left << std::setw(10) << " ";
    stream << std::right << std::setw(2) << this->getElement();
    if (this->getCharge() != gmml::dNotSet)
    {
    	stream << std::left << std::setw(2) << std::to_string(this->getCharge()) << std::endl;
    }
    return;
}
} // namespace
#endif // CDSATOM_HPP
