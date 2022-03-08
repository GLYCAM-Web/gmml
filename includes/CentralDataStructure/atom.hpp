#ifndef INCLUDES_CENTRALDATASTRUCTURE_ATOM_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_ATOM_HPP

#include <string>
#include <iostream>
#include <vector>

#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"

namespace cds
{
class Coordinate;
class Atom : public glygraph::Node<Atom>
{
public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTORS                   //
    //////////////////////////////////////////////////////////
    Atom();
    Atom(const std::string& name, const Coordinate& coord);
    //////////////////////////////////////////////////////////
    //                       ACCESSORS                      //
    //////////////////////////////////////////////////////////
    Coordinate* getCoordinate();
    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////
    void setCoordinate(const Coordinate& c);
    void addCoordinate(const Coordinate& c);
    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////
private:
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////
    std::vector<std::unique_ptr<Coordinate>> coordinates_;     /*!< Position of the atom >*/
};
}
#endif // ATOM_HPP
