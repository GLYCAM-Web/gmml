#ifndef GMML_INPUTSET_CONDENSED_SEQUENCE_PARSED_RESIDUE_HPP
#define GMML_INPUTSET_CONDENSED_SEQUENCE_PARSED_RESIDUE_HPP

#include <string>
#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp" // TemplateGraph
#include "includes/MolecularModeling/Abstract/Residue.hpp"

namespace CondensedSequence
{
	class ParsedResidue : public Abstract::Residue , public glygraph::Node<ParsedResidue>
	{  
	public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        ParsedResidue(std::string residueString, ParsedResidue::Type specifiedType = ParsedResidue::Type::Undefined);
        ParsedResidue(std::string residueString, ParsedResidue* neighbor, ParsedResidue::Type specifiedType = ParsedResidue::Type::Undefined);
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        std::string GetName(const bool withLabels = false);
        std::string GetLinkageName(const bool withLabels = false);
        inline std::string GetInputString() {return fullResidueString_;}
        inline std::string GetIsomer() {return isomer_;}
        inline std::string GetResidueName() {return residueName_;}
        inline std::string GetRingType() {return ringType_;}
        inline std::string GetRingShape() {return ringShape_;}
        inline std::string GetResidueModifier() {return residueModifier_;}
        inline std::string GetConfiguration() {return configuration_;}
        inline std::string GetLinkage() {return linkage_;}
        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
    	void AddLinkage(ParsedResidue* otherRes);
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        std::string GetLink();
        std::vector<ParsedResidue*> GetChildren();
        std::vector<ParsedResidue*> GetParents();
        std::string GetChildLinkagesForGlycamResidueNaming();
        std::string Print();
        std::string GetGlycamResidueName();
        std::string GetGraphVizLine(std::string SnfgFilePath = "");
        std::string GetMonosaccharideName();
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
