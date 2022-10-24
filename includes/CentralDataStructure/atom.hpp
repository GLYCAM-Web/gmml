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
	Atom() {std::cout << "Atom default ctor\n";};
	Atom(const std::string& name, const Coordinate& coord);
	virtual ~Atom() {std::cout << "cds::Atom default dtor for " << this->getName() << "\n";}
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
	//////////////////////////////////////////////////////////
	//                       DISPLAY FUNCTION               //
	//////////////////////////////////////////////////////////
    virtual void Print(std::ostream& out) const;
    // No, paint this like you painted your offWriter
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
} // namespace
#endif // CDSATOM_HPP
