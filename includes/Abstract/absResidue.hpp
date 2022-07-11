#ifndef GMML_INCLUDES_ABSTRACT_RESIDUE_HPP
#define GMML_INCLUDES_ABSTRACT_RESIDUE_HPP

#include <string>

namespace Abstract
{
	class absResidue 
	{
	public:
        enum Type {Aglycone, Sugar, Derivative, Solvent, Protein, Deoxy, Undefined};  
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        inline absResidue::Type GetType() const {return type_;}
        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        inline void SetType(absResidue::Type type) {type_ = type;}
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        Type determineType(const std::string &residueName) const;
    private:      
        //////////////////////////////////////////////////////////
        //                       ATTRRIBUTES                    //
        //////////////////////////////////////////////////////////
        absResidue::Type type_;  // enum Type. See enum above.   
	};
}
#endif
