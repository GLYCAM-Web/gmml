#ifndef COMMON_HPP
#define COMMON_HPP

#include <string>
#include <vector>
#include <math.h>
#include <map>

#include "Geometry/coordinate.hpp"
#include "Glycam/sugarname.hpp"

namespace gmml
{

    //*******************************************
    typedef std::map<std::string, std::string> ResidueNameMap;
    typedef std::map<std::string, std::vector<std::string> > ResidueNameAtomNamesMap;
    typedef Geometry::Coordinate Vector;

    //*******************************************

    const double dNotSet = 123456789.0;
    const int iNotSet = -123456;
    const int iPdbLineLength = 80;
    const double dSulfurCutoff = 2.5;
    const double PI_RADIAN = 4.0*atan(1.0);
    const double PI_DEGREE = 180.0;
    const double EPSILON = 0.001;
    const double dCutOff = 1.65;
    const int PdbResidueThreshold = 500;
    const int DEFAULT_DUMMY_ATOMS = 3;
    const char BLANK_SPACE = '?';

    // Ionizing
    const double DEFAULT_GRID_LENGTH = 2;//0.5;
    const double DEFAULT_GRID_WIDTH = 2;//0.5;
    const double DEFAULT_GRID_HEIGHT = 2;//0.5;
    const double DEFAULT_BOX_WIDTH = INFINITY;//20;
    const double DEFAULT_BOX_LENGTH = INFINITY;//20;
    const double DEFAULT_BOX_HEIGHT = INFINITY;//20;
    const double THRESHOLD = 0.001;//0.001
    const double THRESHOLD_PARTITIONING = 0.001;//0.0001;


    const Glycam::SugarName SUGARNAMELOOKUP[] = {

        {"", "", "", "", "", "", "", "", ""},
        ///Alpha D Aldohexapyranoses
        {"_2_3_4P_a", "a-D-ribopyranose", "DRibpa", "D", "", "P", "a", "", ""},
        {"_3_4^2P_a", "a-D-arabinopyranose" , "DArapa", "D", "", "P", "a", "", ""},
        {"_2_4^3P_a", "a-D-xylopyranose",	"DXylpa", "D", "", "P", "a", "", ""},
        {"_4^2^3P_a",	"a-D-lyxopyranose", "DLyxpa", "D", "", "P", "a", "", ""},
        {"_2_3_4P_a^+1", "a-D-allopyranose", "DAllpa", "D", "", "P", "a", "", ""},
        {"_3_4^2P_a^+1", "a-D-altropyranose", "DAltpa", "D", "", "P", "a", "", ""},
        {"_2_4^3P_a^+1", "a-D-glucopyranose",	"DGlcpa", "D", "", "P", "a", "", ""},
        {"_4^2^3P_a^+1", "a-D-mannopyranose", "DManpa", "D", "", "P", "a", "", ""},
        {"_2_3^4P_a^+1", "a-D-gulopyranose", "DGulpa", "D", "", "P", "a", "", ""},
        {"_3^2^4P_a^+1", "a-D-idopyranose", "DIdopa", "D", "", "P", "a", "", ""},
        {"_2^3^4P_a^+1", "a-D-galactopyranose", "DGalpa", "D", "", "P", "a", "", ""},
        {"^2^3^4P_a^+1", "a-D-talopyranose", "DTalpa", "D", "", "P", "a", "", ""},
        ///Beta D Aldohexapyranoses
        {"_2_3_4P^a", "b-D-ribopyranose", "DRibpb", "D", "", "P", "b", "", ""},
        {"_3_4^2P^a", "b-D-arabinopyranose","DArapb", "D", "", "P", "b", "", ""},
        {"_2_4^3P^a",	"b-D-xylopyranose",	"DXylpb", "D", "", "P", "b", "", ""},
        {"_4^2^3P_a", "b-D-lyxopyranose", "DLyxpb", "D", "", "P", "b", "", ""},
        {"_2_3_4P^a^+1", "b-D-allopyranose",	"DAllpb", "D", "", "P", "b", "", ""},
        {"_3_4^2P^a^+1", "b-D-altropyranose", "DAltpb", "D", "", "P", "b", "", ""},
        {"_2_4^3P^a^+1",	"b-D-glucopyranose", "DGlcpb", "D", "", "P", "b", "", ""},
        {"_4^2^3P^a^+1", "b-D-mannopyranose",	"DManpb", "D", "", "P", "b", "", ""},
        {"_2_3^4P^a^+1", "b-D-gulopyranose", "DGulpb", "D", "", "P", "b", "", ""},
        {"_3^2^4P^a^+1",	"b-D-idopyranose",	"DIdopb", "D", "", "P", "b", "", ""},
        {"_2^3^4P^a^+1", "b-D-galactopyranose", "DGalpb", "D", "", "P", "b", "", ""},
        {"^2^3^4P^a^+1",	"b-D-talopyranose",	"DTalpb", "D", "", "P", "b", "", ""},
        ///Alpha D Ketohexapyranoses
        {"_2_3_4P_a^-1", "a-D-psicopyranose", "DPsipa", "D", "", "P", "a", "", ""},
        {"_3_4^2P_a^-1", "a-D-fructopyranose", "DFrupa", "D", "", "P", "a", "", ""},
        {"_2_4^3P_a^-1", "a-D-sorbopyranose", "DSorpa", "D", "", "P", "a", "", ""},
        {"_4^2^3P_a^-1", "a-D-tagatopyranose", "DTagpa", "D", "", "P", "a", "", ""},
        ///Beta D Ketohexapyranoses
        {"_2_3_4P_-1^a", "b-D-psicopyranose", "DPsipb", "D", "", "P", "b", "", ""},
        {"_3_4^2P_-1^a", "b-D-fructopyranose", "DFrupb", "D", "", "P", "b", "", ""},
        {"_2_4^3P_-1^a", "b-D-sorbopyranose", "DSorpb", "D", "", "P", "b", "", ""},
        {"_4^2^3P_-1^a", "b-D-tagatopyranose", "DTagpb", "D", "", "P", "b", "", ""},
        ///D aldotetrafuranoses
        {"_2_3F_a", "a-D-erythrofuranose", "DEryfa", "D", "", "F", "a", "", ""},
        {"_3^2F_a", "a-D-threofuranose", "DThrfa", "D", "", "F", "a", "", ""},
        {"_2_3F^a", "b-D-erythrofuranose", "DEryfb", "D", "", "F", "b", "", ""},
        {"_3^2F^a", "b-D-threofuranose", "DThrfb", "D", "", "F", "b", "", ""},
        ///D aldopentafuranoses
        {"_2_3F_a^+1", "a-D-ribofuranose", "DRibfa", "D", "", "F", "a", "", ""},
        {"_3^2F_a^+1", "a-D-arabinofuranose", "DArafa", "D", "", "F", "a", "", ""},
        {"_2^3F_a^+1", "a-D-xylofuranose", "DXylfa", "D", "", "F", "a", "", ""},
        {"^2^3F_a^+1", "a-D-lyxofuranose", "DLyxfa", "D", "", "F", "a", "", ""},
        {"_2_3F^a^+1", "b-D-ribofuranose", "DRibfb", "D", "", "F", "b", "", ""},
        {"_3^2F^a^+1", "b-D-arabinofuranose", "DArafb", "D", "", "F", "b", "", ""},
        {"_2^3F^a^+1", "b-D-xylofuranose", "DXylfb", "D", "", "F", "b", "", ""},
        {"^2^3F^a^+1", "b-D-lyxofuranose", "DLyxfb", "D", "", "F", "b", "", ""},
//        ///Alpha D aldohexafuranoses
//        {"_3F^+1R^+2_a", "a-D-allofuranose", "DAllfa", "D", "", "F", "a", "", ""},
//        {"_3F^+1R^+2_a", "a-D-altrofuranose", "DAltfa", "D", "", "F", "a", "", ""},
//        {"^3F^+1R^+2_a", "a-D-glucofuranose", "DGlcfa", "D", "", "F", "a", "", ""},
//        {"^3F^+1R^+2_a", "a-D-mannofuranose", "DManfa", "D", "", "F", "a", "", ""},
//        {"_3F_+1R_+2_a", "a-D-gulofuranose", "DGulfa", "D", "", "F", "a", "", ""},
//        {"_3F_+1R_+2_a", "a-D-idofuranose", "DIdofa", "D", "", "F", "a", "", ""},
//        {"^3F_+1R_+2_a", "a-D-galactofuranose",	"DGalfa", "D", "", "F", "a", "", ""},
//        {"^3F_+1R_+2_a", "a-D-talofuranose", "DTalfa", "D", "", "F", "a", "", ""},
//        ///Beta D aldohexafuranoses
//        {"_3F^+1R^+2^a", "b-D-allofuranose", "DAllfb", "D", "", "F", "b", "", ""},
//        {"_3F^+1R^+2^a", "b-D-altrofuranose", "DAltfb", "D", "", "F", "b", "", ""},
//        {"^3F^+1R^+2^a", "b-D-glucofuranose", "DGlcfb", "D", "", "F", "b", "", ""},
//        {"^3F^+1R^+2^a", "b-D-mannofuranose", "DManfb", "D", "", "F", "b", "", ""},
//        {"_3F_+1R_+2^a", "b-D-gulofuranose", "DGulfb", "D", "", "F", "b", "", ""},
//        {"_3F_+1R_+2^a", "b-D-idofuranose", "DIdofb", "D", "", "F", "b", "", ""},
//        {"^3F_+1R_+2^a", "b-D-galactofuranose", "DGalfb", "D", "", "P", "b", "", ""},
//        {"^3F_+1R_+2^a", "b-D-talofuranose", "DTalfb", "D", "", "F", "b", "", ""},
        ///D ketopentafuranoses
        {"_2_3F_a^-1", "a-D-ribulofuranose", "DRulfa", "D", "", "F", "a", "", ""},
        {"_3^2F_a^-1", "a-D-Xylulofuranose", "DXulfa", "D", "", "F", "a", "", ""},
        {"_2_3F_-1^a", "b-D-ribulofuranose", "DRulfb", "D", "", "F", "b", "", ""},
        {"_3^2F_-1^a", "b-D-Xylulofuranose", "DXulfb", "D", "", "F", "b", "", ""},
        ///D Ketohexafuranoses
        {"_2_3F_a^-1^+1", "a-D-psicofuranose", "DPsifa", "D", "", "F", "a", "", ""},
        {"_3^2F_a^-1^+1", "a-D-fructofuranose", "DFrufa", "D", "", "F", "a", "", ""},
        {"_2^3F_a^-1^+1", "a-D-sorbofuranose", "DSorfa", "D", "", "F", "a", "", ""},
        {"^2^3F_a^-1^+1", "a-D-tagatofuranose", "DTagfa", "D", "", "F", "a", "", ""},
        {"_2_3F_-1^a^+1", "b-D-psicofuranose", "DPsifb", "D", "", "F", "b", "", ""},
        {"_3^2F_-1^a^+1", "b-D-fructofuranose", "DFrufb", "D", "", "F", "b", "", ""},
        {"_2^3F_-1^a^+1", "b-D-sorbofuranose", "DSorfb", "D", "", "F", "b", "", ""},
        {"^2^3F_-1^a^+1", "b-D-tagatofuranose", "DTagfb", "D", "", "F", "b", "", ""},
        ///Other
        ///Deoxy
        {"_2^3^4P^a^+1d", "b-D-fucopyranose", "DFucpb", "D", "", "P", "b", "", ""},
        {"^2^3^4P_a^+1d", "a-D-fucopyranose", "DFucpa", "D", "", "P", "a", "", ""},
        {"_2^43dP^a^+1d", "b-D-abequopyranose", "DAbepb", "D", "", "P", "b", "", ""},
        {"^2^43dP_a^+1d", "a-D-abequopyranose",	"DAbepa", "D", "", "P", "a", "", ""},
        {"_2_4^3P^a^+1d", "b-D-rhamnoopyranose", "DRhapb", "D", "", "P", "b", "", ""},
        {"_4^2^3P_a^+1d", "a-D-rhamnoopyranose", "DRhapa", "D", "", "P", "a", "", ""},
        {"_2_43dP^a^+1d", "b-D-ascarylopyranose", "", "D", "", "P", "b", "", ""},
        {"_4^23dP_a^+1d", "a-D-ascarylopyranose", "","D", "", "P", "a", "", ""},
        {"_32dF^a^+1", "b-D-deoxyribofuranose", "", "D", "", "P", "b", "", ""},
        {"_32dF_a^+1", "a-D-deoxyribofuranose", "", "D", "", "P", "a", "", ""},
        {"_23dF^a^+1" ,"b-D-cordycepofuranose", "", "D", "", "P", "b", "", ""},
        {"_23dF_a^+1", "a-D-cordycepofuranose", "", "D", "", "P", "a", "", ""},

        ///Alpha L Aldohexapyranoses
        {"^2^3^4P^a", "a-L-ribopyranose", "LRibpa", "L", "", "P", "a", "", ""},
        {"_2^3^4P^a", "a-L-arabinopyranose" , "LArapa", "L", "", "P", "a", "", ""},
        {"_3^2^4P^a", "a-L-xylopyranose",	"LXylpa", "L", "", "P", "a", "", ""},
        {"_2_3^4P^a", "a-L-lyxopyranose", "LLyxpa", "L", "", "P", "a", "", ""},
        {"^2^3^4P_+1^a", "a-L-allopyranose", "LAllpa", "L", "", "P", "a", "", ""},
        {"_2^3^4P_+1^a", "a-D-altropyranose", "LAltpa", "L", "", "P", "a", "", ""},
        {"_3^2^4P_+1^a", "a-L-glucopyranose",	"LGlcpa", "L", "", "P", "a", "", ""},
        {"_2_3^4P_+1^a", "a-L-mannopyranose", "LManpa", "L", "", "P", "a", "", ""},
        {"_4^2^3P_+1^a", "a-L-gulopyranose", "LGulpa", "L", "", "P", "a", "", ""},
        {"_2_4^3P_+1^a", "a-L-idopyranose", "LIdopa", "L", "", "P", "a", "", ""},
        {"_3_4^2P_+1^a", "a-L-galactopyranose", "LGalpa", "L", "", "P", "a", "", ""},
        {"_2_3_4P_+1^a", "a-L-talopyranose", "LTalpa", "L", "", "P", "a", "", ""},
        ///Beta L Aldohexapyranoses
        {"^2^3^4P_a", "b-L-ribopyranose", "LRibpb", "L", "", "P", "b", "", ""},
        {"_2^3^4P_a", "b-L-arabinopyranose","LArapb", "L", "", "P", "b", "", ""},
        {"_3^2^4P_a", "b-L-xylopyranose",	"LXylpb", "L", "", "P", "b", "", ""},
        {"_2_3^4P^a", "b-L-lyxopyranose", "LLyxpb", "L", "", "P", "b", "", ""},
        {"^2^3^4P_a_+1", "b-L-allopyranose",	"LAllpb", "L", "", "P", "b", "", ""},
        {"_2^3^4P_a_+1", "b-L-altropyranose", "LAltpb", "L", "", "P", "b", "", ""},
        {"_3^2^4P_a_+1", "b-L-glucopyranose", "LGlcpb", "L", "", "P", "b", "", ""},
        {"_2_3^4P_a_+1", "b-L-mannopyranose",	"LManpb", "L", "", "P", "b", "", ""},
        {"_4^2^3P_a_+1", "b-L-gulopyranose", "LGulpb", "L", "", "P", "b", "", ""},
        {"_2_4^3P_a_+1", "b-L-idopyranose",	"LIdopb", "L", "", "P", "b", "", ""},
        {"_3_4^2P_a_+1", "b-L-galactopyranose", "LGalpb", "L", "", "P", "b", "", ""},
        {"_2_3_4P_a_+1", "b-L-talopyranose",	"LTalpb", "L", "", "P", "b", "", ""},
        ///Alpha L Ketohexapyranoses
        {"^2^3^4P_-1^a", "a-L-psicopyranose", "LPsipa", "L", "", "P", "a", "", ""},
        {"_2^3^4P_-1^a", "a-L-fructopyranose", "LFrupa", "L", "", "P", "a", "", ""},
        {"_3^2^4P_-1^a", "a-L-sorbopyranose", "LSorpa", "L", "", "P", "a", "", ""},
        {"_2_3^4P_-1^a", "a-L-tagatopyranose", "LTagpa", "L", "", "P", "a", "", ""},
        ///Beta L Ketohexapyranoses
        {"^2^3^4P_a^-1", "b-L-psicopyranose", "LPsipb", "L", "", "P", "b", "", ""},
        {"_2^3^4P_a^-1", "b-L-fructopyranose", "LFrupb", "L", "", "P", "b", "", ""},
        {"_3^2^4P_a^-1", "b-L-sorbopyranose", "LSorpb", "L", "", "P", "b", "", ""},
        {"_2_3^4P_a^-1", "b-L-tagatopyranose", "LTagpb", "L", "", "P", "b", "", ""},
        ///L aldotetrafuranoses
        {"^2^3F^a", "a-L-erythrofuranose", "LEryfa", "L", "", "F", "a", "", ""},
        {"_2^3F^a", "a-L-threofuranose", "LThrfa", "L", "", "F", "a", "", ""},
        {"^2^3F_a", "b-L-erythrofuranose", "LEryfb", "L", "", "F", "b", "", ""},
        {"_2^3F_a", "b-L-threofuranose", "LThrfb", "L", "", "F", "b", "", ""},
        ///L aldopentafuranoses
        {"^2^3F_+1^a", "a-L-ribofuranose", "LRibfa", "L", "", "F", "a", "", ""},
        {"_2^3F_+1^a", "a-L-arabinofuranose", "LArafa", "L", "", "F", "a", "", ""},
        {"_3^2F_+1^a", "a-L-xylofuranose", "LXylfa", "L", "", "F", "a", "", ""},
        {"_2_3F_+1^a", "a-L-lyxofuranose", "LLyxfa", "L", "", "F", "a", "", ""},
        {"^2^3F_a_+1", "b-L-ribofuranose", "LRibfb", "L", "", "F", "b", "", ""},
        {"_2^3F_a_+1", "b-L-arabinofuranose", "LArafb", "L", "", "F", "b", "", ""},
        {"_3^2F_a_+1", "b-L-xylofuranose", "LXylfb", "L", "", "F", "b", "", ""},
        {"_2_3F_a_+1", "b-L-lyxofuranose", "LLyxfb", "L", "", "F", "b", "", ""},
//        ///Alpha L aldohexafuranoses
//        {"^3F_+1R_+2^a", "a-L-allofuranose", "LAllfa", "L", "", "F", "a", "", ""},
//        {"^3F_+1R_+2^a", "a-L-altrofuranose", "LAltfa", "L", "", "F", "a", "", ""},
//        {"_3F_+1R_+2^a", "a-L-glucofuranose", "LGlcfa", "L", "", "F", "a", "", ""},
//        {"_3F_+1R_+2^a", "a-L-mannofuranose", "LManfa", "L", "", "F", "a", "", ""},
//        {"^3F^+1R^+2^a", "a-L-gulofuranose", "LGulfa", "L", "", "F", "a", "", ""},
//        {"^3F^+1R^+2^a", "a-L-idofuranose", "LIdofa", "L", "", "F", "a", "", ""},
//        {"_3F^+1R^+2^a", "a-L-galactofuranose",	"LGalfa", "L", "", "F", "a", "", ""},
//        {"_3F^+1R^+2^a", "a-L-talofuranose", "LTalfa", "L", "", "F", "a", "", ""},
//        ///Beta L aldohexafuranoses
//        {"^3F_+1R_+2_a", "b-L-allofuranose", "LAllfb", "L", "", "F", "b", "", ""},
//        {"^3F_+1R_+2_a", "b-L-altrofuranose", "LAltfb", "L", "", "F", "b", "", ""},
//        {"_3F_+1R_+2_a", "b-L-glucofuranose", "LGlcfb", "L", "", "F", "b", "", ""},
//        {"_3F_+1R_+2_a", "b-L-mannofuranose", "LManfb", "L", "", "F", "b", "", ""},
//        {"^3F^+1R^+2_a", "b-L-gulofuranose", "LGulfb", "L", "", "F", "b", "", ""},
//        {"^3F^+1R^+2_a", "b-L-idofuranose", "LIdofb", "L", "", "F", "b", "", ""},
//        {"_3F^+1R^+2_a", "b-L-galactofuranose", "LGalfb", "L", "", "F", "b", "", ""},
//        {"_3F^+1R^+2_a", "b-L-talofuranose", "LTalfb", "L", "", "F", "b", "", ""},
        ///L ketopentafuranoses
        {"^2^3F_-1^a", "a-L-ribulofuranose", "LRulfa", "L", "", "F", "a", "", ""},
        {"_2^3F_-1^a", "a-L-Xylulofuranose", "LXulfa", "L", "", "F", "a", "", ""},
        {"^2^3F_a^-1", "b-L-ribulofuranose", "LRulfb", "L", "", "F", "b", "", ""},
        {"_2^3F_a^-1", "b-L-Xylulofuranose", "LXulfb", "L", "", "F", "b", "", ""},
        ///L Ketohexafuranoses
        {"^2^3F_-1_+1^a", "a-L-psicofuranose", "LPsifa", "L", "", "F", "a", "", ""},
        {"_2^3F_-1_+1^a", "a-L-fructofuranose", "LFrufa", "L", "", "F", "a", "", ""},
        {"_3^2F_-1_+1^a", "a-L-sorbofuranose", "LSorfa", "L", "", "F", "a", "", ""},
        {"_2_3F_-1_+1^a", "a-L-tagatofuranose", "LTagfa", "L", "", "F", "a", "", ""},
        {"^2^3F_a_+1^-1", "b-L-psicofuranose", "LPsifb", "L", "", "F", "b", "", ""},
        {"_2^3F_a_+1^-1", "b-L-fructofuranose", "LFrufb", "L", "", "F", "b", "", ""},
        {"_3^2F_a_+1^-1", "b-L-sorbofuranose", "LSorfb", "L", "", "F", "b", "", ""},
        {"_2_3F_a_+1^-1", "b-L-tagatofuranose", "LTagfb", "L", "", "F", "b", "", ""},
        ///Other
        ///Deoxy
        {"_3_4^2P_a_+1d", "b-L-fucopyranose", "LFucpb", "L", "", "P", "b", "", ""},
        {"_2_3_4P_+1d^a", "a-L-fucopyranose", "LFucpa", "L", "", "P", "a", "", ""},
        {"_4^23dP_a_+1d", "b-L-abequopyranose", "LAbepb", "L", "", "P", "b", "", ""},
        {"_2_43dP_+1d^a", "a-L-abequopyranose",	"LAbepa", "L", "", "P", "a", "", ""},
        {"_3^2^4P_a_+1d", "b-L-rhamnoopyranose", "LRhapb", "L", "", "P", "b", "", ""},
        {"_2_3^4P_+1d^a", "a-L-rhamnoopyranose", "LRhapa", "L", "", "P", "a", "", ""},
        {"^2^43dP_a_+1d", "b-L-ascarylopyranose", "", "L", "", "P", "b", "", ""},
        {"_2^43dP_+1d^a", "a-L-ascarylopyranose", "", "L", "", "P", "a", "", ""},
        {"^32dF_a_+1", "b-L-deoxyribofuranose", "", "L", "", "F", "b", "", ""},
        {"^32dF_+1^a", "a-L-deoxyribofuranose", "", "L", "", "F", "a", "", ""},
        {"^23dF_a_+1" ,"b-L-cordycepofuranose", "", "L", "", "F", "b", "", ""},
        {"^23dF_+1^a", "a-L-cordycepofuranose", "", "L", "", "F", "a", "", ""}

    };

    const int SUGARNAMELOOKUPSIZE = (sizeof(SUGARNAMELOOKUP)/sizeof(SUGARNAMELOOKUP[0]));

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
        MODIFIED = 1        /*!< Force modified parameter file >*/
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

}

#endif // COMMON_HPP
