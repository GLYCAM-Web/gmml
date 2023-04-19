#ifndef GMML_INCLUDES_ABSTRACT_RESIDUE_HPP
#define GMML_INCLUDES_ABSTRACT_RESIDUE_HPP

#include "includes/CodeUtils/constants.hpp" // iNotSet

namespace Abstract
{
enum ResidueType {Protein, Sugar, Aglycone, Derivative, Solvent, Deoxy, Undefined};
class absResidue
{
public:
    //enum Type {Aglycone, Sugar, Derivative, Solvent, Protein, Deoxy, Undefined};
    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////
    inline ResidueType GetType() const {return type_;}
    inline unsigned int getNumber() const {return number_;}
    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////
    inline void SetType(ResidueType type) {type_ = type;}
    inline void setNumber(unsigned int i) {number_ = i;}
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    ResidueType determineType(const std::string &residueName);
private:
    //////////////////////////////////////////////////////////
    //                       ATTRRIBUTES                    //
    //////////////////////////////////////////////////////////
    ResidueType type_ = Undefined;  // enum Type. See enum above.
    unsigned int number_ = constants::iNotSet;
};
}
#endif
