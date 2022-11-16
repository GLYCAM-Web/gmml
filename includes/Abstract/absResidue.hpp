#ifndef GMML_INCLUDES_ABSTRACT_RESIDUE_HPP
#define GMML_INCLUDES_ABSTRACT_RESIDUE_HPP

#include <string>

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
    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////
    inline void SetType(ResidueType type) {type_ = type;}
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    ResidueType determineType(const std::string &residueName) const;
private:
    //////////////////////////////////////////////////////////
    //                       ATTRRIBUTES                    //
    //////////////////////////////////////////////////////////
    ResidueType type_ = Undefined;  // enum Type. See enum above.
};
}
#endif
