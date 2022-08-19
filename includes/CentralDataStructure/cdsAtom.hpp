#ifndef INCLUDES_CENTRALDATASTRUCTURE_ATOM_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_ATOM_HPP

#include <string>
#include <iostream>
#include <vector>

#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"
//#include "includes/CentralDataStructure/cdsCoordinate.hpp"
#include "includes/GeometryTopology/coordinate.hpp"
#include "includes/common.hpp"

using GeometryTopology::Coordinate;
namespace cds
{
class cdsAtom : public glygraph::Node<cdsAtom>
{
public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTORS                   //
    cdsAtom();
    cdsAtom(const std::string& name, const Coordinate& coord);
    //////////////////////////////////////////////////////////
    //                       ACCESSORS                      //
    //////////////////////////////////////////////////////////
    Coordinate* getCoordinate();
    const Coordinate* getCoordinate() const;
    inline const double& getCharge() const { return charge_;};
    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////
    void setCoordinate(const Coordinate& c);
    void addCoordinate(const Coordinate& c);
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    const std::string getElement() const;
    double calculateDistance(const cdsAtom* otherAtom) const;
    void addBond(cdsAtom* otherAtom);
    bool isWithinBondingDistance(const cdsAtom* otherAtom) const;
    void bondIfClose(cdsAtom* otherAtom);
    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////
private:
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////
    std::vector<std::unique_ptr<Coordinate>> coordinates_;     /*!< Position of the atom >*/
    double charge_ = gmml::dNotSet;
};
}
#endif // ATOM_HPP
