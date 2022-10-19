#ifndef GMML_INPUTSET_CONDENSED_SEQUENCE_PARSED_RESIDUE_HPP
#define GMML_INPUTSET_CONDENSED_SEQUENCE_PARSED_RESIDUE_HPP

#include <string>
//#include "includes/Abstract/absResidue.hpp"
#include "includes/CentralDataStructure/cdsResidue.hpp"
#include "includes/CentralDataStructure/atom.hpp"
#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp" // TemplateGraph


namespace CondensedSequence
{
//	class ParsedResidue : public Abstract::absResidue , public glygraph::Node<ParsedResidue>
class ParsedResidue : public cds::cdsResidue<cds::Atom>
	{
	public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        ParsedResidue(std::string residueString, absResidue::Type specifiedType = absResidue::Type::Undefined);
        ParsedResidue(std::string residueString, ParsedResidue* neighbor, absResidue::Type specifiedType = absResidue::Type::Undefined);
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        std::string GetName(const bool withLabels = false) const;
        std::string GetLinkageName(const bool withLabels = false) const;
        inline std::string GetInputString() const {return fullResidueString_;}
        inline std::string GetIsomer() const {return isomer_;}
        inline std::string GetResidueName() const {return residueName_;}
        inline std::string GetRingType() const {return ringType_;}
        inline std::string GetRingShape() const {return ringShape_;}
        inline std::string GetResidueModifier() const {return residueModifier_;}
        inline std::string GetConfiguration() const {return configuration_;}
        inline std::string GetLinkage() const {return linkage_;}
        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
    	void AddLinkage(ParsedResidue* otherRes);
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        std::string GetLink() const;
        std::vector<ParsedResidue*> GetChildren() const;
        std::vector<ParsedResidue*> GetParents() const;
        std::string GetChildLinkagesForGlycamResidueNaming() const;
        std::string Print() const;
        std::string GetGraphVizLine(std::string SnfgFilePath = "") const;
        std::string GetMonosaccharideName() const;
        //////////////////////////////////////////////////////////
        //                  OPERATOR OVERLOADING                //
        //////////////////////////////////////////////////////////
        bool operator== ( ParsedResidue& rhs)  { return (this->GetLink() == rhs.GetLink());}
        bool operator!= ( ParsedResidue& rhs)  { return (this->GetLink() != rhs.GetLink());}
        bool operator> ( ParsedResidue& rhs)  { return (this->GetLink() > rhs.GetLink());}
        bool operator< ( ParsedResidue& rhs)  { return (this->GetLink() < rhs.GetLink());}
	private:
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        void ParseResidueStringIntoComponents(std::string residueString, ParsedResidue::Type specifiedType = ParsedResidue::Type::Undefined);
        void ExciseRingShapeFromModifier();
        //////////////////////////////////////////////////////////
        //                       MUTATORS                       //
        //////////////////////////////////////////////////////////
        inline void SetIsomer(std::string isomer) {isomer_ = isomer;}
        inline void SetResidueName(std::string name) {residueName_ = name;}
        inline void SetRingType(std::string type) {ringType_ = type;}
        inline void SetRingShape(std::string shape) {ringShape_ = shape;}
        inline void SetResidueModifier(std::string modifier) {residueModifier_ = modifier;}
        inline void SetConfiguration(std::string config) {configuration_ = config;}
        inline void SetLinkage(std::string label) {linkage_ = label;}
        //////////////////////////////////////////////////////////
        //                       ATTRRIBUTES                    //
        //////////////////////////////////////////////////////////
		std::string fullResidueString_;           // DManpNAca1-4, etc
		std::string isomer_;                      // D or L
        std::string residueName_;                 // Man, Neu, Ido etc
        std::string ringType_;                    // f or p
        std::string ringShape_;                   // 2SO, 4C1, 1C4 etc
        std::string residueModifier_;             // NAc, Gc, A (IdoA) etc
    	std::string configuration_;               // a or b
        std::string linkage_;                     // 1-4, 2-6, 1- (when connected to OH) etc
	};
}
#endif
