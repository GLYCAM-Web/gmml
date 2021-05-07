#ifndef PARSED_RESIDUE_HPP
#define PARSED_RESIDUE_HPP

#include <string>
#include "includes/MolecularModeling/Graph/Node.hpp" // TemplateGraph
#include "includes/MolecularModeling/Abstract/Residue.hpp"

namespace CondensedSequence
{
	class ParsedResidue : public Abstract::Residue , public TemplateGraph::Node<ParsedResidue>
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
        inline char GetIsomer() {return isomer_;}
        inline std::string GetResidueName() {return residueName_;}
        inline char GetRingType() {return ringType_;}
        inline std::string GetResidueModifier() {return residueModifier_;}
        inline char GetConfiguration() {return configuration_;}
        inline std::string GetLinkage() {return linkage_;}
        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
    	void AddLinkage(ParsedResidue* otherRes);
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        char GetLink();
        std::vector<ParsedResidue*> GetChildren();
        std::vector<ParsedResidue*> GetParents();
        std::string GetChildLinkages();
        std::string Print();
        std::string GetGlycamResidueName();
        std::string GetGraphVizLine();
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
        std::string GetSimpleName();
        std::string GetImageFileName();
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        void ParseResidueStringIntoComponents(std::string residueString, ParsedResidue::Type specifiedType = ParsedResidue::Type::Undefined);
        //////////////////////////////////////////////////////////
        //                       MUTATORS                       //
        //////////////////////////////////////////////////////////
        inline void SetIsomer(char isomer) {isomer_ = isomer;}
        inline void SetResidueName(std::string name) {residueName_ = name;}
        inline void SetRingType(char type) {ringType_ = type;}
        inline void SetResidueModifier(std::string modifier) {residueModifier_ = modifier;}
        inline void SetConfiguration(char config) {configuration_ = config;}
        inline void SetLinkage(std::string label) {linkage_ = label;}
        //////////////////////////////////////////////////////////
        //                       ATTRRIBUTES                    //
        //////////////////////////////////////////////////////////
		std::string fullResidueString_;           // DManpNAca1-4.
		char isomer_;                             // D or L
        std::string residueName_;                 // Man, Neu, Ido etc
        char ringType_;                           // f or p
        std::string residueModifier_;             // NAc, Gc, A (IdoA)
    	char configuration_;                      // a or b
        std::string linkage_;                // 1-4, 2-6, 1- (when connected to OH)
	};
}
#endif