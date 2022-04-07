#ifndef GMML_INCLUDES_ABSTRACT_RESIDUE_HPP
#define GMML_INCLUDES_ABSTRACT_RESIDUE_HPP

namespace Abstract
{
	class absResidue 
	{
	public:
        enum Type {Aglycone, Sugar, Derivative, Solvent, Protein, Deoxy, Undefined};  
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        inline absResidue::Type GetType() {return type_;}
        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        inline void SetType(absResidue::Type type) {type_ = type;}
    private:      
        //////////////////////////////////////////////////////////
        //                       ATTRRIBUTES                    //
        //////////////////////////////////////////////////////////
        absResidue::Type type_;  // enum Type. See enum above.   
	};
}
#endif
