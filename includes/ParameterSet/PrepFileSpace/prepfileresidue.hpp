#ifndef PREPFILERESIDUE_HPP
#define PREPFILERESIDUE_HPP

#include <string>
#include <vector>
#include <map>

namespace PrepFileSpace
{
    enum CoordinateType { kINT, kXYZ };
    enum OutputFormat { kFormatted, kBinary };
    enum GeometryType { kGeometryCorrect, kGeometryChange };
    enum DummyAtomPosition { kPositionAll, kPositionBeg };
    enum DummyAtomOmission { kOmit, kNomit };
    enum SectionType { kSectionLoop, kSectionImproper, kSectionDone, kSectionOther };
    class PrepFileAtom;

    class PrepFileResidue
    {
        public:
            ///////////////////////////// TYPE DEFINITION /////////////////////////////
            typedef std::map<std::string, std::string> Loop;
            typedef std::vector<std::string> Dihedral;

            /////////////////////////////// CONSTRUCTOR ///////////////////////////////
            PrepFileResidue();

            ///////////////////////////// FUNCTION ////////////////////////////////////
            PrepFileResidue* LoadFromStream(std::ifstream& in_file);
            std::string ExtractResidueName(std::istream& ss);
            CoordinateType ExtractResidueCoordinateType(std::istream& ss);
            OutputFormat ExtractResidueOutputFormat(std::istream& ss);
            GeometryType ExtractResidueGeometryType(std::istream& ss);
            DummyAtomOmission ExtractResidueDummyAtomOmission(std::istream& ss);
            DummyAtomPosition ExtractResidueDummyAtomPosition(std::istream& ss);
            SectionType ExtractSectionType(std::string& line);

            ////////////////////////// DISPLAY FUNCTION ///////////////////////////////
            void Print(std::ostream& out);

            ///////////////////////////// ATTRIBUTES //////////////////////////////////
            std::string title_;
            std::string name_;
            CoordinateType coordinate_type_;
            OutputFormat output_format_;
            GeometryType geometry_type_;
            DummyAtomOmission dummy_atom_omission_;
            std::string dummy_atom_type_;
            DummyAtomPosition dummy_atom_position_;
            double charge_;
            std::vector<PrepFileAtom*> atoms_;
            std::vector<Dihedral> improper_dihedrals_;
            std::vector<Loop> loops_;
    };
}

#endif // PREPFILERESIDUE_HPP
