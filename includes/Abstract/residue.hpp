#ifndef GMML_INCLUDES_ABSTRACT_RESIDUE_HPP
#define GMML_INCLUDES_ABSTRACT_RESIDUE_HPP

namespace Abstract
{
	class Residue 
	{
	public:
        enum Type {Aglycone, Sugar, Derivative, Solvent, Protein, Deoxy, Undefined};  
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        inline Residue::Type GetType() {return type_;}
        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        inline void SetType(Residue::Type type) {type_ = type;}
    private:      
        //////////////////////////////////////////////////////////
        //                       ATTRRIBUTES                    //
        //////////////////////////////////////////////////////////
        Residue::Type type_;  // enum Type. See enum above.   
	};
}
#endif // GMML_INCLUDES_ABSTRACT_RESIDUE_HPP