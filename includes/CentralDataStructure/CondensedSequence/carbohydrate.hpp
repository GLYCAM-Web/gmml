#ifndef GMML_INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_CARBOHYDRATE_HPP
#define GMML_INCLUDES_CENTRALDATASTRUCTURE_CONDENSEDSEQUENCE_CARBOHYDRATE_HPP

#include "includes/CentralDataStructure/CondensedSequence/sequenceManipulator.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/Abstract/absBuilder.hpp"
#include "includes/CentralDataStructure/Readers/Prep/prepFile.hpp"
#include <vector>

namespace cdsCondensedSequence
{
    class Carbohydrate : public SequenceManipulator, public Abstract::absBuilder
    {
    public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
    	Carbohydrate(std::string inputSequence = "DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeup5Aca2-6DGalpb1-4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH", std::string prepFilePath = "../dat/prep/GLYCAM_06j-1_GAGS_KDN.prep");
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
    	inline std::string GetInputSequenceString() const {return inputSequenceString_;}
        inline std::vector<cds::ResidueLinkage>& GetGlycosidicLinkages() {return glycosidicLinkages_;}
        inline unsigned long int GetResidueCount() const {return this->getResidues().size();}
    	//////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        void deleteResidue(cds::Residue* byeBye);
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
    	void Generate3DStructureFiles(std::string fileOutputDirectory = "unspecified", std::string outputFileNaming = "structure");
    	void ResolveOverlaps();
        void SetDefaultShapeUsingMetadata();
        unsigned long int CountShapes(bool likelyShapesOnly = false) const;
        std::string GetNumberOfShapes(bool likelyShapesOnly = false) const; // This one is for gems. ToDo try to deprecate and use CountShapes.
        cds::Residue* GetReducingResidue();
        cds::Residue* GetAglycone();
        cds::Atom* GetAnomericAtom();
    private:
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        void deleteLinkage(cds::ResidueLinkage* linkage);
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
    	void ApplyDeoxy(ParsedResidue* deoxyResidue);
    	void MoveAtomsFromPrepResidueToParsedResidue(prep::PrepFile& prepResidues, ParsedResidue* parsedResidue);
    	void DerivativeChargeAdjustment(ParsedResidue* parsedResidue);
        void EnsureIntegralCharge(double charge);
//        void RecurveGenerateResidues(ParsedResidue *currentChild, MolecularModeling::Residue &parent, MolecularModeling::Assembly* assembly);
        void ConnectAndSetGeometry(cds::Residue* parentResidue, cds::Residue* childResidue);
        std::vector<std::string> GetGlycamNamesOfResidues() const;
        std::string GetGlycamResidueName(ParsedResidue *residue) const;
        void DepthFirstSetConnectivityAndGeometry(cds::Residue* currentParent);
        //////////////////////////////////////////////////////////
        //                 PRIVATE MEMBERS                      //
        //////////////////////////////////////////////////////////
        //std::map<std::string, PrepFileSpace::PrepFileResidue*> prepResidueMap_;
        std::string inputSequenceString_;
        std::vector<cds::ResidueLinkage> glycosidicLinkages_;
    };
}
#endif
