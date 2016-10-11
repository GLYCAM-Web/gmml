#ifndef COMMON_HPP
#define COMMON_HPP

#include <string>
#include <vector>
#include <math.h>
#include <map>
#include <set>

#include "GeometryTopology/coordinate.hpp"
#include "Glycan/sugarname.hpp"

namespace gmml
{

    //*******************************************
    typedef std::map<std::string, std::string> ResidueNameMap;
    typedef ResidueNameMap GlycamAtomNameMap;
    typedef std::map<std::string, std::vector<std::string> > AtomMatchingMap;
    typedef std::map<std::string, std::vector<std::string> > GlycamResidueNamingMap;
    typedef std::map<std::string, std::vector<std::string> > ResidueNameAtomNamesMap;
    typedef GeometryTopology::Coordinate Vector;    

    //*******************************************    

    struct AtomTypesInfo{
            std::string atom_type_;
            std::string element_symbol_;
            std::string hybridization_;
    };

    struct ResidueCodeName{
            std::string name_;
            std::string code_;
            int index_;
    };

    struct AmberGlycamMap{
            std::string amber_name_;
            std::string glycam_name_;
    };

    const double dNotSet = 123456789.0;
    const int iNotSet = -123456;
    const int iPdbLineLength = 80;
    const double dSulfurCutoff = 2.5;
    const double PI_RADIAN = 4.0*atan(1.0);
    const double PI_DEGREE = 180.0;
    const double EPSILON = 0.001;
    const double dCutOff = 1.65;
    const double BOND_LENGTH = 1.4;
    const double ROTATION_ANGLE = 109.5;
    const int PdbResidueThreshold = 500;
    const int DEFAULT_DUMMY_ATOMS = 3;
    const char BLANK_SPACE = '?';
    const double MINIMUM_RADIUS = 0.6;
    const double DEFAULT_RADIUS = 1.35;
    const double DIST_EPSILON = 0.000001;
    const int MAX_PDB_ATOM = 99999;
    const double CHARGE_DIVIDER = 18.2223;
    const double CARBON_SURFACE_AREA = 36.31681103;

    // Ionizing
    const double DEFAULT_GRID_LENGTH = 2;//0.5;
    const double DEFAULT_GRID_WIDTH = 2;//0.5;
    const double DEFAULT_GRID_HEIGHT = 2;//0.5;
    const double DEFAULT_BOX_WIDTH = INFINITY;//20;
    const double DEFAULT_BOX_LENGTH = INFINITY;//20;
    const double DEFAULT_BOX_HEIGHT = INFINITY;//20;
    const double THRESHOLD = 0.001;//0.001
    const double THRESHOLD_PARTITIONING = 0.001;//0.0001;
    const double CHARGE_TOLERANCE = 0.001;//0.0001;
    const double GRID_OFFSET = 1.0;
    const double MARGIN = 0.0;//10.0;
    const double CRITICAL_RADIOUS = 1.0;        

    const double EXTERNAL28LINKAGEROTAMERS[][6] = {
        {-66.0, 0.0, -66.0, 80.0, -167.0},
        {-68.0, 15.0, -67.0, -56.0, -59.0},
        {-53.0, 44.0, -68.0, 172.0, 62.0},
        {-73.0, 7.0, -66.0, -58.0, -162.0},
        {-165.0, -15.0, -60.0, 79.0, -169.0},
        {-66.0, 11.0, -68.0, 75.0, 61.0}
    };

    const double INTERNAL28LINKAGEROTAMERS[][5] = {
        {-80.0, -30.0, -68.0, 70.0, -174.0},
        {-74.0, -9.0, -68.0, 71.0, 59.0},
        {-172.0, -9.0, -59.0, 76.0, -172.0},
        {-72.0, 9.0, -72.0, -62.0, -61.0},
        {-43.0, 44.0, -62.0, -178.0, 62.0}
    };

    const AtomTypesInfo ATOMTYPESINFOLOOKUP[] = {
        {"","",""},
        {"2C", "C", "sp3"},
        {"3C", "C", "sp3"},
        {"AC", "C", "sp3"},
        {"BR", "Br", "sp3"},
        {"Br", "Br", "sp3"},
        {"br", "Br", "sp3"},
        {"Br-", "Br", "sp3"},
        {"C", "C", "sp2"},
        {"c", "C", "sp2"},
        {"C*", "C", "sp2"},
        {"C0", "C", "sp2"},
        {"C1", "C", "sp2"},
        {"c1", "C", "sp2"},
        {"C2", "C", "sp2"},
        {"c2", "C", "sp2"},
        {"C3", "C", "sp2"},
        {"c3", "C", "sp3"},
        {"C4", "C", "sp2"},
        {"C5", "C", "sp2"},
        {"C6", "C", "sp2"},
        {"C8", "C", "sp3"},
        {"CA", "C", "sp2"},
        {"cA", "C", "sp3"},
        {"ca", "C", "sp2"},
        {"CB", "C", "sp2"},
        {"cB", "C", "sp2"},
        {"CC", "C", "sp2"},
        {"cC", "C", "sp2"},
        {"cc", "C", "sp2"},
        {"CD", "C", "sp2"},
        {"cd", "C", "sp2"},
        {"CE", "C", "sp2"},
        {"ce", "C", "sp2"},
        {"CF", "C", "sp2"},
        {"cf", "C", "sp2"},
        {"CG", "C", "sp3"},
        {"Cg", "C", "sp3"},
        {"cg", "C", "sp2"},
        {"CH", "C", "sp3"},
        {"ch", "C", "sp2"},
        {"CI", "C", "sp2"},
        {"CJ", "C", "sp2"},
        {"Cj", "C", "sp2"},
        {"CK", "C", "sp2"},
        {"Ck", "C", "sp2"},
        {"CL", "Cl", "sp3"},
        {"Cl", "Cl", "sp3"},
        {"cl", "Cl", "sp3"},
        {"Cl-", "Cl", "sp3"},
        {"CM", "C", "sp2"},
        {"CN", "C", "sp2"},
        {"CO", "C", "sp2"},
        {"CP", "C", "sp2"},
        {"Cp", "C", "sp3"},
        {"cP", "C", "sp3"},
        {"cp", "C", "sp2"},
        {"CQ", "C", "sp2"},
        {"cq", "C", "sp2"},
        {"CR", "C", "sp2"},
        {"cR", "C", "sp3"},
        {"CS", "C", "sp2"},
        {"Cs", "Cs", "sp3"},
        {"Cs+", "Cs", "sp3"},
        {"CT", "C", "sp3"},
        {"cu", "C", "sp2"},
        {"CV", "C", "sp2"},
        {"cv", "C", "sp2"},
        {"CW", "C", "sp2"},
        {"CX", "C", "sp3"},
        {"cx", "C", "sp2"},
        {"CY", "C", "sp2"},
        {"Cy", "C", "sp3"},
        {"cy", "C", "sp2"},
        {"EC", "C", "sp3"},
        {"EP", "", "sp3"},
        {"F", "F", "sp3"},
        {"f", "F", "sp3"},
        {"F-", "F", "sp3"},
        {"FE", "Fe", "sp3"},
        {"H", "H", "sp3"},
        {"H0", "H", "sp3"},
        {"H1", "H", "sp3"},
        {"h1", "H", "sp3"},
        {"H2", "H", "sp3"},
        {"h2", "H", "sp3"},
        {"H3", "H", "sp3"},
        {"h3", "H", "sp3"},
        {"H4", "H", "sp3"},
        {"h4", "H", "sp3"},
        {"H5", "H", "sp3"},
        {"h5", "H", "sp3"},
        {"HA", "H", "sp3"},
        {"Ha", "H", "sp3"},
        {"hA", "H", "sp3"},
        {"ha", "H", "sp3"},
        {"hB", "H", "sp3"},
        {"HC", "H", "sp3"},
        {"Hc", "H", "sp3"},
        {"hc", "H", "sp3"},
        {"hE", "H", "sp3"},
        {"hN", "H", "sp3"},
        {"hn", "H", "sp3"},
        {"HO", "H", "sp3"},
        {"Ho", "H", "sp3"},
        {"hO", "H", "sp3"},
        {"ho", "H", "sp3"},
        {"HP", "H", "sp3"},
        {"Hp", "H", "sp3"},
        {"hp", "H", "sp3"},
        {"hR", "H", "sp3"},
        {"HS", "H", "sp3"},
        {"hS", "H", "sp3"},
        {"hs", "H", "sp3"},
        {"HW", "H", "sp3"},
        {"hw", "H", "sp3"},
        {"hX", "H", "sp3"},
        {"hx", "H", "sp3"},
        {"HZ", "H", "sp3"},
        {"I", "I", "sp3"},
        {"i", "I", "sp3"},
        {"I-", "I", "sp3"},
        {"IB", "Na", "sp3"},
        {"IM", "Cl", "sp3"},
        {"IP", "Na", "sp3"},
        {"K", "K", "sp3"},
        {"K+", "K", "sp3"},
        {"Li", "Li", "sp3"},
        {"Li+", "Li", "sp3"},
        {"LP", "", "sp3"},
        {"MG", "Mg", "sp3"},
        {"Mg+", "Mg", "sp3"},
        {"N", "N", "sp2"},
        {"n", "N", "sp2"},
        {"N*", "N", "sp2"},
        {"n1", "N", "sp2"},
        {"N2", "N", "sp2"},
        {"n2", "N", "sp2"},
        {"N3", "N", "sp3"},
        {"n3", "N", "sp3"},
        {"n4", "N", "sp3"},
        {"NA", "N", "sp2"},
        {"nA", "N", "sp3"},
        {"na", "N", "sp2"},
        {"Na+", "Na", "sp3"},
        {"NB", "N", "sp2"},
        {"nb", "N", "sp2"},
        {"NC", "N", "sp2"},
        {"nc", "N", "sp2"},
        {"nd", "N", "sp2"},
        {"ne", "N", "sp2"},
        {"nf", "N", "sp2"},
        {"Ng", "N", "sp2"},
        {"nh", "N", "sp2"},
        {"no", "N", "sp2"},
        {"NP", "N", "sp2"},
        {"NQ", "N", "sp2"},
        {"NT", "N", "sp3"},
        {"NY", "N", "sp2"},
        {"O", "O", "sp2"},
        {"o", "O", "sp2"},
        {"O2", "O", "sp2"},
        {"o2", "O", "sp2"},
        {"Oa", "N", "sp2"},
        {"oC", "O", "sp2"},
        {"Oe", "N", "sp2"},
        {"OG", "O", "sp3"},
        {"OH", "O", "sp3"},
        {"Oh", "O", "sp3"},
        {"oH", "O", "sp3"},
        {"oh", "O", "sp3"},
        {"OL", "O", "sp3"},
        {"OM", "O", "sp2"},
        {"oO", "O", "sp2"},
        {"OP", "O", "sp2"},
        {"oP", "O", "sp2"},
        {"oR", "O", "sp3"},
        {"OS", "O", "sp3"},
        {"Os", "O", "sp3"},
        {"oS", "O", "sp3"},
        {"os", "O", "sp3"},
        {"oT", "O", "sp3"},
        {"OW", "O", "sp3"},
        {"ow", "O", "sp3"},
        {"OY", "O", "sp3"},
        {"Oy", "O", "sp3"},
        {"P", "P", "sp3"},
        {"p2", "P", "sp2"},
        {"p3", "P", "sp3"},
        {"p4", "P", "sp3"},
        {"p5", "P", "sp3"},
        {"pA", "P", "sp3"},
        {"pb", "P", "sp3"},
        {"pd", "P", "sp3"},
        {"px", "P", "sp3"},
        {"py", "P", "sp3"},
        {"Rb", "Rb", "sp3"},
        {"Rb+", "Rb", "sp3"},
        {"S", "S", "sp3"},
        {"s", "S", "sp2"},
        {"s2", "S", "sp2"},
        {"s3", "S", "sp3"},
        {"s4", "S", "sp3"},
        {"s6", "S", "sp3"},
        {"SH", "S", "sp3"},
        {"sh", "S", "sp3"},
        {"SM", "S", "sp3"},
        {"Sm", "S", "sp3"},
        {"ss", "S", "sp3"},
        {"sx", "S", "sp3"},
        {"sy", "S", "sp3"},
        {"Zn", "Zn", "sp3"}

    };

    const int ATOMTYPESINFOLOOKUPSIZE = (sizeof(ATOMTYPESINFOLOOKUP)/sizeof(ATOMTYPESINFOLOOKUP[0]));

    const ResidueCodeName RESIDUENAMECODELOOKUP[] = {
        {"", "", 0},
        {"All", "N", 1},
        {"Alt", "E", 2},
        {"Ara", "A", 3},
        {"Bac", "BC", 4},
        {"Fru", "C", 5},
        {"Fuc", "F", 6},
        {"Gal", "L", 7},
        {"GalA", "O", 8},
        {"GalNAc", "V", 9},
        {"Glc", "G", 10},
        {"GlcA", "Z", 11},
        {"GlcNAc", "Y", 12},
        {"Gul", "K", 13},
        {"Ido", "I", 14},
        {"IdoA", "U", 15},
        {"IdoA(1C4)", "U1", 16},
        {"IdoA(2SO)", "U2", 17},
        {"IdoA(4C1)", "U3", 18},
        {"Lyx", "D", 19},
        {"Man", "M", 20},
        {"ManNAc", "W", 21},
        {"Neu", "", 22},
        {"Neu5Ac", "S", 23},
        {"Neu5Gc", "GL", 24},
        {"NeuNAc", "S", 25},
        {"NeuNGc", "GL", 26},
        {"Psi", "P", 27},
        {"Qui", "Q", 28},
        {"Rha", "H", 29},
        {"Rib", "R", 30},
        {"Sor", "B", 31},
        {"Tag", "J", 32},
        {"Tal", "T", 33},
        {"Tyv", "TV", 34},
        {"Xyl", "X", 35},
        {"GlcNS", "Y", 36},
        {"U", "045", 37},
        {"S", "245", 38},
    };

    const int RESIDUENAMECODELOOKUPSIZE = (sizeof(RESIDUENAMECODELOOKUP)/sizeof(RESIDUENAMECODELOOKUP[0]));

    const AmberGlycamMap AMBERGLYCAMLOOKUP[] = {
        {"", ""},
        {"ASN", "NLN"},
        {"THR", "OLT"},
        {"SER", "OLS"}
    };

    const int AMBERGLYCAMLOOKUPSIZE = (sizeof(AMBERGLYCAMLOOKUP)/sizeof(AMBERGLYCAMLOOKUP[0]));

    const Glycan::SugarName SUGARNAMELOOKUP[] = {

        {"", "", "", "", "", "", "", "", "", ""},
        ///Alpha D Aldohexapyranoses
        {"_2_3_4P_a", "a-D-ribopyranose", "DRibpa", "D", "", "P", "a", "", "", "RIB"},
        {"_3_4^2P_a", "a-D-arabinopyranose" , "DArapa", "D", "", "P", "a", "", "", ""},
        {"_2_4^3P_a", "a-D-xylopyranose",	"DXylpa", "D", "", "P", "a", "", "", "XYZ"},
        {"_4^2^3P_a",	"a-D-lyxopyranose", "DLyxpa", "D", "", "P", "a", "", "", "LDY"},
        {"_2_3_4P_a^+1", "a-D-allopyranose", "DAllpa", "D", "", "P", "a", "", "", "AFD"},
        {"_3_4^2P_a^+1", "a-D-altropyranose", "DAltpa", "D", "", "P", "a", "", "", ""},
        {"_2_4^3P_a^+1", "a-D-glucopyranose",	"DGlcpa", "D", "", "P", "a", "", "", "GLC"},
        {"_4^2^3P_a^+1", "a-D-mannopyranose", "DManpa", "D", "", "P", "a", "", "", "MAN"},
        {"_2_3^4P_a^+1", "a-D-gulopyranose", "DGulpa", "D", "", "P", "a", "", "", ""},
        {"_3^2^4P_a^+1", "a-D-idopyranose", "DIdopa", "D", "", "P", "a", "", "", ""},
        {"_2^3^4P_a^+1", "a-D-galactopyranose", "DGalpa", "D", "", "P", "a", "", "", "GLA"},
        {"^2^3^4P_a^+1", "a-D-talopyranose", "DTalpa", "D", "", "P", "a", "", "", ""},
        ///Beta D Aldohexapyranoses
        {"_2_3_4P^a", "b-D-ribopyranose", "DRibpb", "D", "", "P", "b", "", "", "RIP"},
        {"_3_4^2P^a", "b-D-arabinopyranose","DArapb", "D", "", "P", "b", "", "", ""},
        {"_2_4^3P^a",	"b-D-xylopyranose",	"DXylpb", "D", "", "P", "b", "", "", "XYP"},
        {"_4^2^3P_a", "b-D-lyxopyranose", "DLyxpb", "D", "", "P", "b", "", "", ""},
        {"_2_3_4P^a^+1", "b-D-allopyranose",	"DAllpb", "D", "", "P", "b", "", "", "ALL"},
        {"_3_4^2P^a^+1", "b-D-altropyranose", "DAltpb", "D", "", "P", "b", "", "", ""},
        {"_2_4^3P^a^+1",	"b-D-glucopyranose", "DGlcpb", "D", "", "P", "b", "", "", "BGC"},
        {"_4^2^3P^a^+1", "b-D-mannopyranose",	"DManpb", "D", "", "P", "b", "", "", "BMA"},
        {"_2_3^4P^a^+1", "b-D-gulopyranose", "DGulpb", "D", "", "P", "b", "", "", "GL0"},
        {"_3^2^4P^a^+1",	"b-D-idopyranose",	"DIdopb", "D", "", "P", "b", "", "", ""},
        {"_2^3^4P^a^+1", "b-D-galactopyranose", "DGalpb", "D", "", "P", "b", "", "", "GAL,GLB"},
        {"^2^3^4P^a^+1",	"b-D-talopyranose",	"DTalpb", "D", "", "P", "b", "", "", ""},
        ///Alpha D Ketohexapyranoses
        {"_2_3_4P_a^-1", "a-D-psicopyranose", "DPsipa", "D", "", "P", "a", "", "", "PSJ"},
        {"_3_4^2P_a^-1", "a-D-fructopyranose", "DFrupa", "D", "", "P", "a", "", "", "FRU,FUD"},
        {"_2_4^3P_a^-1", "a-D-sorbopyranose", "DSorpa", "D", "", "P", "a", "", "", "SDD"},
        {"_4^2^3P_a^-1", "a-D-tagatopyranose", "DTagpa", "D", "", "P", "a", "", "", "T6T,TAG"},
        ///Beta D Ketohexapyranoses
        {"_2_3_4P_-1^a", "b-D-psicopyranose", "DPsipb", "D", "", "P", "b", "", "", "PSJ"},
        {"_3_4^2P_-1^a", "b-D-fructopyranose", "DFrupb", "D", "", "P", "b", "", "", "BDF"},
        {"_2_4^3P_-1^a", "b-D-sorbopyranose", "DSorpb", "D", "", "P", "b", "", "", "SDD"},
        {"_4^2^3P_-1^a", "b-D-tagatopyranose", "DTagpb", "D", "", "P", "b", "", "", "TAG"},
        ///D aldotetrafuranoses
        {"_2_3F_a", "a-D-erythrofuranose", "DEryfa", "D", "", "F", "a", "", "", ""},
        {"_3^2F_a", "a-D-threofuranose", "DThrfa", "D", "", "F", "a", "", "", ""},
        {"_2_3F^a", "b-D-erythrofuranose", "DEryfb", "D", "", "F", "b", "", "", ""},
        {"_3^2F^a", "b-D-threofuranose", "DThrfb", "D", "", "F", "b", "", "", ""},
        ///D aldopentafuranoses
        {"_2_3F_a^+1", "a-D-ribofuranose", "DRibfa", "D", "", "F", "a", "", "", ""},
        {"_3^2F_a^+1", "a-D-arabinofuranose", "DArafa", "D", "", "F", "a", "", "", "BXY"},
        {"_2^3F_a^+1", "a-D-xylofuranose", "DXylfa", "D", "", "F", "a", "", "", ""},
        {"^2^3F_a^+1", "a-D-lyxofuranose", "DLyxfa", "D", "", "F", "a", "", "", ""},
        {"_2_3F^a^+1", "b-D-ribofuranose", "DRibfb", "D", "", "F", "b", "", "", ""},
        {"_3^2F^a^+1", "b-D-arabinofuranose", "DArafb", "D", "", "F", "b", "", "", ""},
        {"_2^3F^a^+1", "b-D-xylofuranose", "DXylfb", "D", "", "F", "b", "", "", ""},
        {"^2^3F^a^+1", "b-D-lyxofuranose", "DLyxfb", "D", "", "F", "b", "", "", ""},
        ///D ketopentafuranoses
        {"_2_3F_a^-1", "a-D-ribulofuranose", "DRulfa", "D", "", "F", "a", "", "", ""},
        {"_3^2F_a^-1", "a-D-Xylulofuranose", "DXulfa", "D", "", "F", "a", "", "", ""},
        {"_2_3F_-1^a", "b-D-ribulofuranose", "DRulfb", "D", "", "F", "b", "", "", ""},
        {"_3^2F_-1^a", "b-D-Xylulofuranose", "DXulfb", "D", "", "F", "b", "", "", ""},
        ///D Ketohexafuranoses
        {"_2_3F_a^-1^+1", "a-D-psicofuranose", "DPsifa", "D", "", "F", "a", "", "", "PSV"},
        {"_3^2F_a^-1^+1", "a-D-fructofuranose", "DFrufa", "D", "", "F", "a", "", "", ""},
        {"_2^3F_a^-1^+1", "a-D-sorbofuranose", "DSorfa", "D", "", "F", "a", "", "", "SDD"},
        {"^2^3F_a^-1^+1", "a-D-tagatofuranose", "DTagfa", "D", "", "F", "a", "", "", "TAG"},
        {"_2_3F_-1^a^+1", "b-D-psicofuranose", "DPsifb", "D", "", "F", "b", "", "", "PSJ"},
        {"_3^2F_-1^a^+1", "b-D-fructofuranose", "DFrufb", "D", "", "F", "b", "", "", ""},
        {"_2^3F_-1^a^+1", "b-D-sorbofuranose", "DSorfb", "D", "", "F", "b", "", "", "SDD"},
        {"^2^3F_-1^a^+1", "b-D-tagatofuranose", "DTagfb", "D", "", "F", "b", "", "", "TAG"},
        ///Other
        ///Deoxy
        {"_2^3^4P^a^+1d", "b-D-fucopyranose", "DFucpb", "D", "", "P", "b", "", "", "FCB"},
        {"_2^3^4P_a^+1d", "a-D-fucopyranose", "DFucpa", "D", "", "P", "a", "", "", "FCA"},
        {"_2^43dP^a^+1d", "b-D-abequopyranose", "DAbepb", "D", "", "P", "b", "", "", ""},
        {"_2^43dP_a^+1d", "a-D-abequopyranose",	"DAbepa", "D", "", "P", "a", "", "", ""},
        {"_4^2^3P^a^+1d", "b-D-rhamnoopyranose", "DRhapb", "D", "", "P", "b", "", "", ""},
        {"_4^2^3P_a^+1d", "a-D-rhamnoopyranose", "DRhapa", "D", "", "P", "a", "", "", ""},
        {"_4^23dP^a^+1d", "b-D-ascarylopyranose", "", "D", "", "P", "b", "", "", ""},
        {"_4^23dP_a^+1d", "a-D-ascarylopyranose", "","D", "", "P", "a", "", "", ""},
        {"_32dF^a^+1", "b-D-deoxyribofuranose", "", "D", "", "P", "b", "", "", ""},
        {"_32dF_a^+1", "a-D-deoxyribofuranose", "", "D", "", "P", "a", "", "", ""},
        {"_23dF^a^+1" ,"b-D-cordycepofuranose", "", "D", "", "P", "b", "", "", ""},
        {"_23dF_a^+1", "a-D-cordycepofuranose", "", "D", "", "P", "a", "", "", ""},

        ///Alpha L Aldohexapyranoses
        {"^2^3^4P^a", "a-L-ribopyranose", "LRibpa", "L", "", "P", "a", "", "", ""},
        {"_2^3^4P^a", "a-L-arabinopyranose" , "LArapa", "L", "", "P", "a", "", "", "ARA"},
        {"_3^2^4P^a", "a-L-xylopyranose",	"LXylpa", "L", "", "P", "a", "", "", "HSY"},
        {"_2_3^4P^a", "a-L-lyxopyranose", "LLyxpa", "L", "", "P", "a", "", "", ""},
        {"^2^3^4P_+1^a", "a-L-allopyranose", "LAllpa", "L", "", "P", "a", "", "", ""},
        {"_2^3^4P_+1^a", "a-D-altropyranose", "LAltpa", "L", "", "P", "a", "", "", ""},
        {"_3^2^4P_+1^a", "a-L-glucopyranose",	"LGlcpa", "L", "", "P", "a", "", "", ""},
        {"_2_3^4P_+1^a", "a-L-mannopyranose", "LManpa", "L", "", "P", "a", "", "", ""},
        {"_4^2^3P_+1^a", "a-L-gulopyranose", "LGulpa", "L", "", "P", "a", "", "", "GUP"},
        {"_2_4^3P_+1^a", "a-L-idopyranose", "LIdopa", "L", "", "P", "a", "", "", ""},
        {"_3_4^2P_+1^a", "a-L-galactopyranose", "LGalpa", "L", "", "P", "a", "", "", ""},
        {"_2_3_4P_+1^a", "a-L-talopyranose", "LTalpa", "L", "", "P", "a", "", "", ""},
        ///Beta L Aldohexapyranoses
        {"^2^3^4P_a", "b-L-ribopyranose", "LRibpb", "L", "", "P", "b", "", "", ""},
        {"_2^3^4P_a", "b-L-arabinopyranose","LArapb", "L", "", "P", "b", "", "", "ARB"},
        {"_3^2^4P_a", "b-L-xylopyranose",	"LXylpb", "L", "", "P", "b", "", "", "LXC"},
        {"_2_3^4P^a", "b-L-lyxopyranose", "LLyxpb", "L", "", "P", "b", "", "", ""},
        {"^2^3^4P_a_+1", "b-L-allopyranose",	"LAllpb", "L", "", "P", "b", "", "", ""},
        {"_2^3^4P_a_+1", "b-L-altropyranose", "LAltpb", "L", "", "P", "b", "", "", ""},
        {"_3^2^4P_a_+1", "b-L-glucopyranose", "LGlcpb", "L", "", "P", "b", "", "", ""},
        {"_2_3^4P_a_+1", "b-L-mannopyranose",	"LManpb", "L", "", "P", "b", "", "", ""},
        {"_4^2^3P_a_+1", "b-L-gulopyranose", "LGulpb", "L", "", "P", "b", "", "", ""},
        {"_2_4^3P_a_+1", "b-L-idopyranose",	"LIdopb", "L", "", "P", "b", "", "", ""},
        {"_3_4^2P_a_+1", "b-L-galactopyranose", "LGalpb", "L", "", "P", "b", "", "", ""},
        {"_2_3_4P_a_+1", "b-L-talopyranose",	"LTalpb", "L", "", "P", "b", "", "", ""},
        ///Alpha L Ketohexapyranoses
        {"^2^3^4P_-1^a", "a-L-psicopyranose", "LPsipa", "L", "", "P", "a", "", "", ""},
        {"_2^3^4P_-1^a", "a-L-fructopyranose", "LFrupa", "L", "", "P", "a", "", "", "FRU,FUD"},
        {"_3^2^4P_-1^a", "a-L-sorbopyranose", "LSorpa", "L", "", "P", "a", "", "", "SDD,SOE,SOL"},
        {"_2_3^4P_-1^a", "a-L-tagatopyranose", "LTagpa", "L", "", "P", "a", "", "", "TAG"},
        ///Beta L Ketohexapyranoses
        {"^2^3^4P_a^-1", "b-L-psicopyranose", "LPsipb", "L", "", "P", "b", "", "", ""},
        {"_2^3^4P_a^-1", "b-L-fructopyranose", "LFrupb", "L", "", "P", "b", "", "", ""},
        {"_3^2^4P_a^-1", "b-L-sorbopyranose", "LSorpb", "L", "", "P", "b", "", "", "SDD"},
        {"_2_3^4P_a^-1", "b-L-tagatopyranose", "LTagpb", "L", "", "P", "b", "", "", "TAG"},
        ///L aldotetrafuranoses
        {"^2^3F^a", "a-L-erythrofuranose", "LEryfa", "L", "", "F", "a", "", "", ""},
        {"_2^3F^a", "a-L-threofuranose", "LThrfa", "L", "", "F", "a", "", "", ""},
        {"^2^3F_a", "b-L-erythrofuranose", "LEryfb", "L", "", "F", "b", "", "", ""},
        {"_2^3F_a", "b-L-threofuranose", "LThrfb", "L", "", "F", "b", "", "", ""},
        ///L aldopentafuranoses
        {"^2^3F_+1^a", "a-L-ribofuranose", "LRibfa", "L", "", "F", "a", "", "", "Z6J"},
        {"_2^3F_+1^a", "a-L-arabinofuranose", "LArafa", "L", "", "F", "a", "", "", "AHR"},
        {"_3^2F_+1^a", "a-L-xylofuranose", "LXylfa", "L", "", "F", "a", "", "", ""},
        {"_2_3F_+1^a", "a-L-lyxofuranose", "LLyxfa", "L", "", "F", "a", "", "", ""},
        {"^2^3F_a_+1", "b-L-ribofuranose", "LRibfb", "L", "", "F", "b", "", "", "32O"},
        {"_2^3F_a_+1", "b-L-arabinofuranose", "LArafb", "L", "", "F", "b", "", "", ""},
        {"_3^2F_a_+1", "b-L-xylofuranose", "LXylfb", "L", "", "F", "b", "", "", ""},
        {"_2_3F_a_+1", "b-L-lyxofuranose", "LLyxfb", "L", "", "F", "b", "", "", ""},
        ///L ketopentafuranoses
        {"^2^3F_-1^a", "a-L-ribulofuranose", "LRulfa", "L", "", "F", "a", "", "", ""},
        {"_2^3F_-1^a", "a-L-Xylulofuranose", "LXulfa", "L", "", "F", "a", "", "", ""},
        {"^2^3F_a^-1", "b-L-ribulofuranose", "LRulfb", "L", "", "F", "b", "", "", ""},
        {"_2^3F_a^-1", "b-L-Xylulofuranose", "LXulfb", "L", "", "F", "b", "", "", ""},
        ///L Ketohexafuranoses
        {"^2^3F_-1_+1^a", "a-L-psicofuranose", "LPsifa", "L", "", "F", "a", "", "", "SF6"},
        {"_2^3F_-1_+1^a", "a-L-fructofuranose", "LFrufa", "L", "", "F", "a", "", "", ""},
        {"_3^2F_-1_+1^a", "a-L-sorbofuranose", "LSorfa", "L", "", "F", "a", "", "", "SDD"},
        {"_2_3F_-1_+1^a", "a-L-tagatofuranose", "LTagfa", "L", "", "F", "a", "", "", "TAG"},
        {"^2^3F_a_+1^-1", "b-L-psicofuranose", "LPsifb", "L", "", "F", "b", "", "", "SF9"},
        {"_2^3F_a_+1^-1", "b-L-fructofuranose", "LFrufb", "L", "", "F", "b", "", "", "LFR"},
        {"_3^2F_a_+1^-1", "b-L-sorbofuranose", "LSorfb", "L", "", "F", "b", "", "", "SDD"},
        {"_2_3F_a_+1^-1", "b-L-tagatofuranose", "LTagfb", "L", "", "F", "b", "", "", "TAG"},
        ///Other
        ///Deoxy
        {"_3_4^2P_a_+1d", "b-L-fucopyranose", "LFucpb", "L", "", "P", "b", "", "", ""},
        {"_3_4^2P_+1d^a", "a-L-fucopyranose", "LFucpa", "L", "", "P", "a", "", "", "AFL,FUC"},
        {"_4^23dP_a_+1d", "b-L-abequopyranose", "LAbepb", "L", "", "P", "b", "", "", ""},
        {"_4^23dP_a^+1d", "a-L-abequopyranose",	"LAbepa", "L", "", "P", "a", "", "", ""},
        {"_2_3^4P_a_+1d", "b-L-rhamnoopyranose", "LRhapb", "L", "", "P", "b", "", "", ""},
        {"_2_3^4P_+1d^a", "a-L-rhamnoopyranose", "LRhapa", "L", "", "P", "a", "", "", ""},
        {"_2^43dP_+1d_a", "b-L-ascarylopyranose", "", "L", "", "P", "b", "", "", ""},
        {"_2^43dP_+1d^a", "a-L-ascarylopyranose", "", "L", "", "P", "a", "", "", ""},
        {"^32dF_a_+1", "b-L-deoxyribofuranose", "", "L", "", "F", "b", "", "", ""},
        {"^32dF_+1^a", "a-L-deoxyribofuranose", "", "L", "", "F", "a", "", "", ""},
        {"^23dF_a_+1" ,"b-L-cordycepofuranose", "", "L", "", "F", "b", "", "", ""},
        {"^23dF_+1^a", "a-L-cordycepofuranose", "", "L", "", "F", "a", "", "", ""},

        ///Indeterminate Anomeric Oxygen
        ///X D Aldohexapyranoses
        {"_2_3_4P", "x-D-ribopyranose", "DRibpx", "D", "", "P", "", "", "", ""},
        {"_3_4^2P", "x-D-arabinopyranose" , "DArapx", "D", "", "P", "", "", "", ""},
        {"_2_4^3P", "x-D-xylopyranose",	"DXylpx", "D", "", "P", "", "", "", "XYS"},
        {"_4^2^3P",	"x-D-lyxopyranose", "DLyxpx", "D", "", "P", "", "", "", ""},
        {"_2_3_4P^+1", "x-D-allopyranose", "DAllpx", "D", "", "P", "", "", "", ""},
        {"_3_4^2P^+1", "x-D-altropyranose", "DAltpx", "D", "", "P", "", "", "", ""},
        {"_2_4^3P^+1", "x-D-glucopyranose",	"DGlcpx", "D", "", "P", "", "", "", ""},
        {"_4^2^3P^+1", "x-D-mannopyranose", "DManpx", "D", "", "P", "", "", "", ""},
        {"_2_3^4P^+1", "x-D-gulopyranose", "DGulpx", "D", "", "P", "", "", "", ""},
        {"_3^2^4P^+1", "x-D-idopyranose", "DIdopx", "D", "", "P", "", "", "", ""},
        {"_2^3^4P^+1", "x-D-galactopyranose", "DGalpx", "D", "", "P", "", "", "", ""},
        {"^2^3^4P^+1", "x-D-talopyranose", "DTalpx", "D", "", "P", "", "", "", ""},
        ///Alpha D Ketohexapyranoses
        {"_2_3_4P^-1", "a-D-psicopyranose", "DPsipa", "D", "", "P", "a", "", "", "PSJ"},
        {"_3_4^2P^-1", "a-D-fructopyranose", "DFrupa", "D", "", "P", "a", "", "", "FRU,FUD"},
        {"_2_4^3P^-1", "a-D-sorbopyranose", "DSorpa", "D", "", "P", "a", "", "", "SDD"},
        {"_4^2^3P^-1", "a-D-tagatopyranose", "DTagpa", "D", "", "P", "a", "", "", "TAG"},
        ///Beta D Ketohexapyranoses
        {"_2_3_4P_-1", "b-D-psicopyranose", "DPsipb", "D", "", "P", "b", "", "", "PSJ"},
        {"_3_4^2P_-1", "b-D-fructopyranose", "DFrupb", "D", "", "P", "b", "", "", "FRU,FUD"},
        {"_2_4^3P_-1", "b-D-sorbopyranose", "DSorpb", "D", "", "P", "b", "", "", "SDD"},
        {"_4^2^3P_-1", "b-D-tagatopyranose", "DTagpb", "D", "", "P", "b", "", "", "TAG"},
        ///X D aldotetrafuranoses
        {"_2_3F", "x-D-erythrofuranose", "DEryfx", "D", "", "F", "", "", "", ""},
        {"_3^2F", "x-threofuranose", "DThrfx", "D", "", "F", "", "", "", ""},
        ///X D aldopentafuranoses
        {"_2_3F^+1", "x-D-ribofuranose", "DRibfx", "D", "", "F", "", "", "", ""},
        {"_3^2F^+1", "x-D-arabinofuranose", "DArafx", "D", "", "F", "", "", "", ""},
        {"_2^3F^+1", "x-D-xylofuranose", "DXylfx", "D", "", "F", "", "", "", ""},
        {"^2^3F^+1", "x-D-lyxofuranose", "DLyxfx", "D", "", "F", "", "", "", ""},
        ///D ketopentafuranoses
        {"_2_3F^-1", "a-D-ribulofuranose", "DRulfa", "D", "", "F", "a", "", "", ""},
        {"_3^2F^-1", "a-D-Xylulofuranose", "DXulfa", "D", "", "F", "a", "", "", ""},
        {"_2_3F_-1", "b-D-ribulofuranose", "DRulfb", "D", "", "F", "b", "", "", ""},
        {"_3^2F_-1", "b-D-Xylulofuranose", "DXulfb", "D", "", "F", "b", "", "", ""},
        ///D Ketohexafuranoses
        {"_2_3F^-1^+1", "a-D-psicofuranose", "DPsifa", "D", "", "F", "a", "", "", ""},
        {"_3^2F^-1^+1", "a-D-fructofuranose", "DFrufa", "D", "", "F", "a", "", "", ""},
        {"_2^3F^-1^+1", "a-D-sorbofuranose", "DSorfa", "D", "", "F", "a", "", "", "SDD"},
        {"^2^3F^-1^+1", "a-D-tagatofuranose", "DTagfa", "D", "", "F", "a", "", "", "TAG"},
        {"_2_3F_-1^+1", "b-D-psicofuranose", "DPsifb", "D", "", "F", "b", "", "", ""},
        {"_3^2F_-1^+1", "b-D-fructofuranose", "DFrufb", "D", "", "F", "b", "", "", ""},
        {"_2^3F_-1^+1", "b-D-sorbofuranose", "DSorfb", "D", "", "F", "b", "", "", "SDD"},
        {"^2^3F_-1^+1", "b-D-tagatofuranose", "DTagfb", "D", "", "F", "b", "", "", "TAG"},
        ///Other
        ///Deoxy
        {"_2^3^4P^+1d", "x-D-fucopyranose", "DFucpx", "D", "", "P", "", "", "", ""},
        {"_2^43dP^+1d", "x-D-abequopyranose", "DAbepx", "D", "", "P", "", "", "", ""},
        {"_4^2^3P^+1d", "x-D-rhamnoopyranose", "DRhapx", "D", "", "P", "", "", "", ""},
        {"_4^23dP^+1d", "x-D-ascarylopyranose", "", "D", "", "P", "", "", "", ""},
        {"_32dF^+1", "x-D-deoxyribofuranose", "", "D", "", "P", "", "", "", ""},
        {"_23dF^+1", "x-D-cordycepofuranose", "", "D", "", "P", "", "", "", ""},

        ///X L Aldohexapyranoses
        {"^2^3^4P", "x-L-ribopyranose", "LRibpx", "L", "", "P", "", "", "", ""},
        {"_2^3^4P", "x-L-arabinopyranose" , "LArapx", "L", "", "P", "", "", "", ""},
        {"_3^2^4P", "x-L-xylopyranose",	"LXylpx", "L", "", "P", "", "", "", ""},
        {"_2_3^4P",	"x-L-lyxopyranose", "LLyxpx", "L", "", "P", "", "", "", ""},
        {"^2^3^4P_+1", "x-L-allopyranose", "LAllpx", "L", "", "P", "", "", "", ""},
        {"_2^3^4P_+1", "x-L-altropyranose", "LAltpx", "L", "", "P", "", "", "", ""},
        {"_3^2^4P_+1", "x-L-glucopyranose",	"LGlcpx", "L", "", "P", "", "", "", ""},
        {"_2_3^4P_+1", "x-L-mannopyranose", "LManpx", "L", "", "P", "", "", "", ""},
        {"_4^2^3P_+1", "x-L-gulopyranose", "LGulpx", "L", "", "P", "", "", "", ""},
        {"_2_4^3P_+1", "x-L-idopyranose", "LIdopx", "L", "", "P", "", "", "", ""},
        {"_3_4^2P_+1", "x-L-galactopyranose", "LGalpx", "L", "", "P", "", "", "", ""},
        {"_2_3_4P_+1", "x-L-talopyranose", "LTalpx", "L", "", "P", "", "", "", ""},
        ///Alpha L Ketohexapyranoses
        {"^2^3^4P_-1", "a-L-psicopyranose", "LPsipa", "L", "", "P", "a", "", "", ""},
        {"_2^3^4P_-1", "a-L-fructopyranose", "LFrupa", "L", "", "P", "a", "", "", "FRU,FUD"},
        {"_3^2^4P_-1", "a-L-sorbopyranose", "LSorpa", "L", "", "P", "a", "", "", "SDD,SOL,SOE"},
        {"_2_3^4P_-1", "a-L-tagatopyranose", "LTagpa", "L", "", "P", "a", "", "", "TAG"},
        ///Beta L Ketohexapyranoses
        {"^2^3^4P^-1", "b-L-psicopyranose", "LPsipb", "L", "", "P", "b", "", "", ""},
        {"_2^3^4P^-1", "b-L-fructopyranose", "LFrupb", "L", "", "P", "b", "", "", ""},
        {"_3^2^4P^-1", "b-L-sorbopyranose", "LSorpb", "L", "", "P", "b", "", "", "SDD"},
        {"_2_3^4P^-1", "b-L-tagatopyranose", "LTagpb", "L", "", "P", "b", "", "", "TAG"},
        ///X L aldotetrafuranoses
        {"^2^3F", "x-L-erythrofuranose", "LEryfx", "L", "", "F", "", "", "", ""},
        {"_2^3F", "x-threofuranose", "LThrfx", "L", "", "F", "", "", "", ""},
        ///X L aldopentafuranoses
        {"^2^3F_+1", "x-L-ribofuranose", "LRibfx", "L", "", "F", "", "", "", ""},
        {"_2^3F_+1", "x-L-arabinofuranose", "LArafx", "L", "", "F", "", "", "", ""},
        {"_3^2F_+1", "x-L-xylofuranose", "LXylfx", "L", "", "F", "", "", "", ""},
        {"_2_3F_+1", "x-L-lyxofuranose", "LLyxfx", "L", "", "F", "", "", "", ""},
        ///L ketopentafuranoses
        {"^2^3F_-1", "a-L-ribulofuranose", "LRulfa", "L", "", "F", "a", "", "", ""},
        {"_2^3F_-1", "a-L-Xylulofuranose", "LXulfa", "L", "", "F", "a", "", "", ""},
        {"^2^3F^-1", "b-L-ribulofuranose", "LRulfb", "L", "", "F", "b", "", "", ""},
        {"_2^3F^-1", "b-L-Xylulofuranose", "LXulfb", "L", "", "F", "b", "", "", ""},
        ///L Ketohexafuranoses
        {"^2^3F_-1_+1", "a-L-psicofuranose", "LPsifa", "L", "", "F", "a", "", "", ""},
        {"_2^3F_-1_+1", "a-L-fructofuranose", "LFrufa", "L", "", "F", "a", "", "", ""},
        {"_3^2F_-1_+1", "a-L-sorbofuranose", "LSorfa", "L", "", "F", "a", "", "", "SDD"},
        {"_2_3F_-1_+1", "a-L-tagatofuranose", "LTagfa", "L", "", "F", "a", "", "", "TAG"},
        {"^2^3F_+1^-1", "b-L-psicofuranose", "LPsifb", "L", "", "F", "b", "", "", ""},
        {"_2^3F_+1^-1", "b-L-fructofuranose", "LFrufb", "L", "", "F", "b", "", "", ""},
        {"_3^2F_+1^-1", "b-L-sorbofuranose", "LSorfb", "L", "", "F", "b", "", "", "SDD"},
        {"_2_3F_+1^-1", "b-L-tagatofuranose", "LTagfb", "L", "", "F", "b", "", "", "TAG"},
        ///Other
        ///Deoxy
        {"_3_4^2P_+1d", "x-L-fucopyranose", "LFucpx", "L", "", "P", "", "", "", ""},
        {"_4^23dP_+1d", "x-L-abequopyranose", "LAbepx", "L", "", "P", "", "", "", ""},
        {"_2_3^4P_+1d", "x-L-rhamnoopyranose", "LRhapx", "L", "", "P", "", "", "", ""},
        {"_2^43dP_+1d", "x-L-ascarylopyranose", "", "L", "", "P", "", "", "", ""},
        {"^32dF_+1", "x-L-deoxyribofuranose", "", "L", "", "P", "", "", "", ""},
        {"^23dF_+1", "x-L-cordycepofuranose", "", "L", "", "P", "", "", "", ""}
    };

    const int SUGARNAMELOOKUPSIZE = (sizeof(SUGARNAMELOOKUP)/sizeof(SUGARNAMELOOKUP[0]));

    const Glycan::SugarName COMPLEXSUGARNAMELOOKUP[] = {

        {"", "", "", "", "", "", "", "", "", ""},
        ///Alpha D aldohexafuranoses
        {"_2_3F_a^+1R^+2", "a-D-allofuranose", "DAllfa", "D", "", "F", "a", "", "", ""},
        {"_3^2F_a^+1R^+2", "a-D-altrofuranose", "DAltfa", "D", "", "F", "a", "", "", ""},
        {"_2^3F_a^+1R^+2", "a-D-glucofuranose", "DGlcfa", "D", "", "F", "a", "", "", ""},
        {"^2^3F_a^+1R^+2", "a-D-mannofuranose", "DManfa", "D", "", "F", "a", "", "", ""},
        {"_2_3F_a_+1R_+2", "a-D-gulofuranose", "DGulfa", "D", "", "F", "a", "", "", ""},
        {"_3^2F_a_+1R_+2", "a-D-idofuranose", "DIdofa", "D", "", "F", "a", "", "", ""},
        {"_2^3F_a_+1R_+2", "a-D-galactofuranose", "DGalfa", "D", "", "F", "a", "", "", ""},
        {"^2^3F_a_+1R_+2", "a-D-talofuranose", "DTalfa", "D", "", "F", "a", "", "", ""},
        ///Beta D aldohexafuranoses
        {"_2_3F^a^+1R^+2", "b-D-allofuranose", "DAllfb", "D", "", "F", "b", "", "", ""},
        {"_3^2F^a^+1R^+2", "b-D-altrofuranose", "DAltfb", "D", "", "F", "b", "", "", ""},
        {"_2^3F^a^+1R^+2", "b-D-glucofuranose", "DGlcfb", "D", "", "F", "b", "", "", ""},
        {"^2^3F^a^+1R^+2", "b-D-mannofuranose", "DManfb", "D", "", "F", "b", "", "", ""},
        {"_2_3F_+1R_+2^a", "b-D-gulofuranose", "DGulfb", "D", "", "F", "b", "", "", ""},
        {"_3^2F_+1R_+2^a", "b-D-idofuranose", "DIdofb", "D", "", "F", "b", "", "", ""},
        {"_2^3F_+1R_+2^a", "b-D-galactofuranose", "DGalfb", "D", "", "P", "b", "", "", ""},
        {"^2^3F_+1R_+2^a", "b-D-talofuranose", "DTalfb", "D", "", "F", "b", "", "", ""},
        ///Eight carbons
        {"^3^42dP_a^-1A^+1R^+2", "", "", "D", "", "P", "a", "2-keto-3-deoxy-a-D-octulopyranosate", "DKDOpa", ""},
        {"^3^42dP_-1A^a^+1R^+2", "", "", "D", "", "P", "b", "2-keto-3-deoxy-b-D-octulopyranosate", "DKDOpb", ""},
        {"^32dF_a_+1R_+2R_+3^-1A", "", "", "D", "", "F", "a", "2-keto-3-deoxy-a-D-octulofuranosate", "DKDOfa", ""},
        {"^32dF_-1A_+1R_+2R_+3^a", "", "", "D", "", "F", "b", "2-keto-3-deoxy-b-D-octulofuranosate", "DKDOfb", ""},
        ///Nine carbons
        {"_3^42dP_a_+1R_+2R_+3^-1A", "", "", "D", "", "P", "a", "2-keto-3-deoxy-a-D-nonulopyranosate", "DKDNpa", ""},
        {"_3^42dP_-1A_+1R_+2R_+3^a", "", "", "D", "", "P", "b", "2-keto-3-deoxy-b-D-nonulopyranosate", "DKDNpb", ""},
        {"_3^4NAc2dP_a_+1R_+2R_+3^-1A", "", "", "D", "", "P", "a", "N-acetyl-a-D-neuraminate", "DNeupNAca", "SIA"},
        {"_3^4NAc2dP_-1A_+1R_+2R_+3^a", "", "", "D", "", "P", "b", "N-acetyl-b-D-neuraminate", "DNeupNAcb", ""},
        {"_3Ac^4N2dP_+1R_+2R_+3^-1A_a", "", "", "D", "", "P", "a", "4-O-acetyl-a-D-neuraminate", "DNeup4Aca", ""},
        {"_3Ac^4N2dP_-1A_+1R_+2R_+3^a", "", "", "D", "", "P", "b", "4-O-acetyl-b-D-neuraminate", "DNeup4Acb", ""},
        {"_3^4NGc2dP_a_+1R_+2R_+3^-1A", "", "", "D", "", "P", "a", "N-glycolyl-a-D-neuraminate", "DNeupNGca", ""},
        {"_3^4NGc2dP_-1A_+1R_+2R_+3^a", "", "", "D", "", "P", "b", "N-glycolyl-b-D-neuraminate", "DNeupNGcb", ""},
        {"_3^4N2dP_a_+1R_+2R_+3^-1A", "", "", "D", "", "P", "a", "a-D-neuraminate", "DNeupa", ""},
        {"_3^4N2dP_-1A_+1R_+2R_+3^a", "", "", "D", "", "P", "b", "b-D-neuraminate", "DNeupb", ""},

        ///Alpha L aldohexafuranoses
        {"^2^3F_+1S_+2^a", "a-L-allofuranose", "LAllfa", "L", "", "F", "a", "", "", ""},
        {"_2^3F_+1S_+2^a", "a-L-altrofuranose", "LAltfa", "L", "", "F", "a", "", "", ""},
        {"_3^2F_+1S_+2^a", "a-L-glucofuranose", "LGlcfa", "L", "", "F", "a", "", "", ""},
        {"_2_3F_+1S_+2^a", "a-L-mannofuranose", "LManfa", "L", "", "F", "a", "", "", ""},
        {"^2^3F^a^+1S^+2", "a-L-gulofuranose", "LGulfa", "L", "", "F", "a", "", "", ""},
        {"_2^3F^a^+1S^+2", "a-L-idofuranose", "LIdofa", "L", "", "F", "a", "", "", ""},
        {"_3^2F^a^+1S^+2", "a-L-galactofuranose", "LGalfa", "L", "", "F", "a", "", "", ""},
        {"_2_3F^a^+1S^+2", "a-L-talofuranose", "LTalfa", "L", "", "F", "a", "", "", ""},
        ///Beta L aldohexafuranoses
        {"^2^3F_a_+1S_+2", "b-L-allofuranose", "LAllfb", "L", "", "F", "b", "", "", ""},
        {"_2^3F_a_+1S_+2", "b-L-altrofuranose", "LAltfb", "L", "", "F", "b", "", "", ""},
        {"_3^2F_a_+1S_+2", "b-L-glucofuranose", "LGlcfb", "L", "", "F", "b", "", "", ""},
        {"_2_3F_a_+1S_+2", "b-L-mannofuranose", "LManfb", "L", "", "F", "b", "", "", ""},
        {"^2^3F_a^+1S^+2", "b-L-gulofuranose", "LGulfb", "L", "", "F", "b", "", "", ""},
        {"_2^3F_a^+1S^+2", "b-L-idofuranose", "LIdofb", "L", "", "F", "b", "", "", ""},
        {"_3^2F_a^+1S^+2", "b-L-galactofuranose", "LGalfb", "L", "", "F", "b", "", "", ""},
        {"_2_3F_a^+1S^+2", "b-L-talofuranose", "LTalfb", "L", "", "F", "b", "", "", ""},
        ///Eight carbons
        {"_3_42dP_-1A_+1S_+2^a", "", "", "L", "", "P", "a", "2-keto-3-deoxy-a-L-octulopyranosate", "LKLOpa", ""},
        {"_3_42dP_a_+1S_+2^-1A", "", "", "L", "", "P", "b", "2-keto-3-deoxy-b-L-octulopyranosate", "LKLOpb", ""},
        {"_32dF_-1A^a^+1S^+2S^+3", "", "", "L", "", "F", "a", "2-keto-3-deoxy-a-L-octulofuranosate", "LKLOfa", ""},
        {"_32dF_a^-1A^+1S^+2S^+3", "", "", "L", "", "F", "b", "2-keto-3-deoxy-b-L-octulofuranosate", "LKLOfb", ""},
        ///Nine carbons
        {"_4^32dP_-1A^a^+1S^+2S^+3", "", "", "L", "", "P", "a", "2-keto-3-deoxy-a-L-nonulopyranosate", "LKLNpa", ""},
        {"_4^32dP_a^-1A^+1S^+2S^+3", "", "", "L", "", "P", "b", "2-keto-3-deoxy-b-L-nonulopyranosate", "LKLNpb", ""},
        {"_4NAc^32dP_-1A^a^+1S^+2S^+3", "", "", "L", "", "P", "a", "N-acetyl-a-L-neuraminate",	"LNeupNAca", ""},
        {"_4NAc^32dP_a^-1A^+1S^+2S^+3", "", "", "L", "", "P", "b", "N-acetyl-b-L-neuraminate",	"LNeupNAcb", ""},
        {"_4N^3Ac2dP_-1A^a^+1S^+2S^+3", "", "", "L", "", "P", "a", "4-O-acetyl-a-L-neuraminate", "LNeup4Aca", ""},
        {"_4N^3Ac2dP_a^-1A^+1S^+2S^+3", "", "", "L", "", "P", "b", "4-O-acetyl-b-L-neuraminate", "LNeup4Acb", ""},
        {"_4NGc^32dP_-1A^a^+1S^+2S^+3", "", "", "L", "", "P", "a", "N-glycolyl-a-L-neuraminate", "LNeupNGca", ""},
        {"_4NGc^32dP_a^-1A^+1S^+2S^+3", "", "", "L", "", "P", "b", "N-glycolyl-b-L-neuraminate", "LNeupNGcb", ""},
        {"_4N^32dP_-1A^a^+1S^+2S^+3", "", "", "L", "", "P", "a", "a-L-neuraminate", "LNeupa", ""},
        {"_4N^32dP_a^-1A^+1S^+2S^+3", "", "", "L", "", "P", "b", "b-L-neuraminate", "LNeupb", ""},

        ///Indeterminate Anomeric Oxygen
        ///Alpha D aldohexafuranoses
        {"_2_3F^+1R^+2", "x-D-allofuranose", "DAllfx", "D", "", "F", "", "", "", ""},
        {"_3^2F^+1R^+2", "x-D-altrofuranose", "DAltfx", "D", "", "F", "", "", "", ""},
        {"_2^3F^+1R^+2", "x-D-glucofuranose", "DGlcfx", "D", "", "F", "", "", "", ""},
        {"^2^3F^+1R^+2", "x-D-mannofuranose", "DManfx", "D", "", "F", "", "", "", ""},
        {"_2_3F_+1R_+2", "x-D-gulofuranose", "DGulfx", "D", "", "F", "", "", "", ""},
        {"_3^2F_+1R_+2", "x-D-idofuranose", "DIdofx", "D", "", "F", "", "", "", ""},
        {"_2^3F_+1R_+2", "x-D-galactofuranose", "DGalfx", "D", "", "F", "", "", "", ""},
        {"^2^3F_+1R_+2", "x-D-talofuranose", "DTalfx", "D", "", "F", "", "", "", ""},
        ///Eight carbons
        {"^3^42dP^-1A^+1R^+2", "", "", "D", "", "P", "a", "2-keto-3-deoxy-a-D-octulopyranosate", "DKDOpa", ""},
        {"^3^42dP_-1A^+1R^+2", "", "", "D", "", "P", "b", "2-keto-3-deoxy-b-D-octulopyranosate", "DKDOpb", ""},
        {"^32dF_+1R_+2R_+3^-1A", "", "", "D", "", "F", "a", "2-keto-3-deoxy-a-D-octulofuranosate", "DKDOfa", ""},
        {"^32dF_-1A_+1R_+2R_+3", "", "", "D", "", "F", "b", "2-keto-3-deoxy-b-D-octulofuranosate", "DKDOfb", ""},
        ///Nine carbons
        {"_3^42dP_+1R_+2R_+3^-1A", "", "", "D", "", "P", "a", "2-keto-3-deoxy-a-D-nonulopyranosate", "DKDNpa", ""},
        {"_3^42dP_-1A_+1R_+2R_+3", "", "", "D", "", "P", "b", "2-keto-3-deoxy-b-D-nonulopyranosate", "DKDNpb", ""},
        {"_3^4NAc2dP_+1R_+2R_+3^-1A", "", "", "D", "", "P", "a", "N-acetyl-a-D-neuraminate", "DNeupNAca", ""},
        {"_3^4NAc2dP_-1A_+1R_+2R_+3", "", "", "D", "", "P", "b", "N-acetyl-b-D-neuraminate", "DNeupNAcb", ""},
        {"_3Ac^4N2dP_+1R_+2R_+3^-1A", "", "", "D", "", "P", "a", "4-O-acetyl-a-D-neuraminate", "DNeup4Aca", ""},
        {"_3Ac^4N2dP_-1A_+1R_+2R_+3", "", "", "D", "", "P", "b", "4-O-acetyl-b-D-neuraminate", "DNeup4Acb", ""},
        {"_3^4NGc2dP_+1R_+2R_+3^-1A", "", "", "D", "", "P", "a", "N-glycolyl-a-D-neuraminate", "DNeupNGca", ""},
        {"_3^4NGc2dP_-1A_+1R_+2R_+3", "", "", "D", "", "P", "b", "N-glycolyl-b-D-neuraminate", "DNeupNGcb", ""},
        {"_3^4N2dP_+1R_+2R_+3^-1A", "", "", "D", "", "P", "a", "a-D-neuraminate", "Neupa", ""},
        {"_3^4N2dP_-1A_+1R_+2R_+3", "", "", "D", "", "P", "b", "b-D-neuraminate", "Neupb", ""},

        ///Alpha D aldohexafuranoses
        {"^2^3F_+1S_+2", "x-L-allofuranose", "LAllfx", "L", "", "F", "", "", "", ""},
        {"_2^3F_+1S_+2", "x-L-altrofuranose", "LAltfx", "L", "", "F", "", "", "", ""},
        {"_3^2F_+1S_+2", "x-L-glucofuranose", "LGlcfx", "L", "", "F", "", "", "", ""},
        {"_2_3F_+1S_+2", "x-L-mannofuranose", "LManfx", "L", "", "F", "", "", "", ""},
        {"^2^3F^+1S^+2", "x-L-gulofuranose", "LGulfx", "L", "", "F", "", "", "", ""},
        {"_2^3F^+1S^+2", "x-L-idofuranose", "LIdofx", "L", "", "F", "", "", "", ""},
        {"_3^2F^+1S^+2", "x-L-galactofuranose", "LGalfx", "L", "", "F", "", "", "", ""},
        {"_2_3F^+1S^+2", "x-L-talofuranose", "LTalfx", "L", "", "F", "", "", "", ""},
        ///Eight carbons
        {"_3_42dP_-1A_+1S_+2", "", "", "L", "", "P", "a", "2-keto-3-deoxy-a-L-octulopyranosate", "LKDOpa", ""},
        {"_3_42dP_+1S_+2^-1A", "", "", "L", "", "P", "b", "2-keto-3-deoxy-b-L-octulopyranosate",	"LKDOpb", ""},
        {"_32dF_-1A^+1S^+2S^+3", "", "", "L", "", "F", "a", "2-keto-3-deoxy-a-L-octulofuranosate", "LKDOfa", ""},
        {"_32dF^-1A^+1S^+2S^+3", "", "", "L", "", "F", "b", "2-keto-3-deoxy-b-L-octulofuranosate", "LKDOfb", ""},
        ///Nine carbons
        {"_4^32dP_-1A^+1S^+2S^+3", "", "", "L", "", "P", "a", "2-keto-3-deoxy-a-L-nonulopyranosate", "LKDNpa", ""},
        {"_4^32dP^-1A^+1S^+2S^+3", "", "", "L", "", "P", "b", "2-keto-3-deoxy-b-L-nonulopyranosate", "LKDNpb", ""},
        {"_4NAc^32dP_-1A^+1S^+2S^+3", "", "", "L", "", "P", "a", "N-acetyl-a-L-neuraminate", "LNeupNAca", ""},
        {"_4NAc^32dP^-1A^+1S^+2S^+3", "", "", "L", "", "P", "b", "N-acetyl-b-L-neuraminate", "LNeupNAcb", ""},
        {"_4N^3Ac2dP_-1A^+1S^+2S^+3", "", "", "L", "", "P", "a", "4-O-acetyl-a-L-neuraminate", "LNeup4Aca", ""},
        {"_4N^3Ac2dP^-1A^+1S^+2S^+3", "", "", "L", "", "P", "b", "4-O-acetyl-b-L-neuraminate", "LNeup4Acb", ""},
        {"_4NGc^32dP_-1A^+1S^+2S^+3", "", "", "L", "", "P", "a", "N-glycolyl-a-L-neuraminate", "LNeupNGca", ""},
        {"_4NGc^32dP^-1A^+1S^+2S^+3", "", "", "L", "", "P", "b", "N-glycolyl-b-L-neuraminate", "LNeupNGcb", ""},
        {"_4N^32dP_-1A^+1S^+2S^+3", "", "", "L", "", "P", "a", "a-L-neuraminate", "LNeupa", ""},
        {"_4N^32dP^-1A^+1S^+2S^+3", "", "", "L", "", "P", "b", "b-L-neuraminate", "LNeupb", ""}
    };

        const int COMPLEXSUGARNAMELOOKUPSIZE = (sizeof(COMPLEXSUGARNAMELOOKUP)/sizeof(COMPLEXSUGARNAMELOOKUP[0]));

    /*! \enum
      * Enumerator to possible n chain termination
      */
    enum PossibleNChainTermination
    {
        COCH3 = 1,
        NH3 = 2
    };

    /*! \enum
      * Enumerator to possible c chain termination
      */
    enum PossibleCChainTermination
    {
        NH2 = 1,
        NHCH3 = 2,
        CO2 = 3
    };

    /*! \enum
      * Enumerator to possible HIS mapping in pdb preprocessor
      */
    enum PdbPreprocessorHISMapping
    {
        HIE = 1,
        HIP = 2,
        HID = 3
    };

    /*! \enum
      * Enumerator to possible input file type to central data structure
      */
    enum InputFileType
    {
        PDB,
        PDBQT,
        LIB,
        PREP,
        TOP,
        TOP_CRD,
        MULTIPLE,
        UNKNOWN
    };

    /*! \enum
      * Enumerator to possible options for building a graph for a central data structure
      */
    enum BuildingStructureOption
    {
        DISTANCE,
        ORIGINAL,
        DATABASE
    };

    /*! \enum
      * Enumerator to parameter file type
      */
    enum ParameterFileType
    {
        MAIN = 0,           /*!< Main parameter file >*/
        MODIFIED = 1,        /*!< Force modified parameter file >*/
        IONICMOD = 2
    };

    /*! \enum
      * Topological type enumerator
      */
    enum TopologicalType
    {
        kTopTypeE,
        kTopTypeS,
        kTopTypeB,
        kTopType3,
        kTopType4,
        kTopTypeM
    };

    /*! \enum
      * Enumerator to topological type stack
      */
    enum TopologicalTypeStackElement
    {
        EMPTY = 0,
        M1 = 1,
        M2 = 2,
        M3 = 3,
        M4 = 4,
        S1 = 5,
        B1 = 6,
        B2 = 7,
        T1 = 8,
        T2 = 9,
        T3 = 10
    };

    /*! \enum
      * Enumerator to topological type stack
      */
    enum GraphSearchNodeStatus
    {
        UNVISITED = 0,
        VISITED = 1,
        DONE = 2
    };
    /*! \enum
      * LogLevel enumerator
      */
    enum LogLevel
    {
        INF,
        ERR,
        WAR
    };
    /*! \enum
      * Condensed sequence token type
      */
    enum CondensedSequenceTokenType
    {
        CONDENSED_SEQUENCE_LEFT_BRACKET,
        CONDENSED_SEQUENCE_RIGHT_BRACKET,
        CONDENSED_SEQUENCE_RESIDUE
    };

    /*! \enum
      * Enumerator to possible resource type URIs for ontology
      */
    enum URIType
    {
        OntPDB,
        OntResidue,
        OntOligosaccharide,
        OntNote,
        OntLinkage,
        OntMonosaccharide,
        OntSugarName,
        OntAtom
    };    
}

#endif // COMMON_HPP
