#ifndef INCLUDES_ABSTRACT_ATOM_HPP
#define INCLUDES_ABSTRACT_ATOM_HPP

#include <string>
#include <iostream>
#include <vector>

#include "includes/GeometryTopology/coordinate.hpp"
#include "includes/common.hpp" // dNotSet

using GeometryTopology::Coordinate;
namespace Abstract
{
class absAtom
{
public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTORS                   //
	absAtom() {};
	absAtom(const Coordinate& coord);
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
    double calculateDistance(const absAtom* otherAtom) const;
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
#endif // INCLUDES_ABSTRACT_ATOM_HPP
