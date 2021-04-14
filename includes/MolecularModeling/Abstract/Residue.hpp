#ifndef ABSTRACT_RESIDUE_HPP
#define ABSTRACT_RESIDUE_HPP

// Can I #include nothing and still compile? Nah...

namespace Abstract
{
	class Residue 
	{
	public:
        enum Type {Aglycone, Sugar, Derivative, Solvent, Protein, Undefined};    
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
#endif