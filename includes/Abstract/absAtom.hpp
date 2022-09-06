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
    inline double getCharge() const { return charge_;}
    inline std::string getType() const {return atomType_;}
    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////
    inline void setCharge(const double c) {charge_ = c;}
    inline void setType(const std::string s) {atomType_ = s;}
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
    std::string atomType_ = " ";
};
}
#endif // INCLUDES_ABSTRACT_ATOM_HPP
