#ifndef GMML_INCLUDES_CODEUTILS_BIOLOGY_HPP
#define GMML_INCLUDES_CODEUTILS_BIOLOGY_HPP

#include <string>
#include <vector>

namespace biology
{
    static const std::vector<std::string> proteinResidueNames = {
        "ALA",  "ASP",  "ASN",  "ARG",  "GLY",  "GLU",   "GLN",   "PRO",  "HIS",  "HIP",  "CYS",  "VAL",
        "LEU",  "THR",  "SER",  "LYS",  "MET",  "MSE",   "TYR",   "TRP",  "PHE",  "SEC",  "ILE",  "CYX",
        "CYM",  "HID",  "HIE",  "NLN",  "OLY",  "OLS",   "OLT",   "ARN",  "ASH",  "GLH",  "HYP",  "LYN",
        "NALA", "NASP", "NASN", "NARG", "NGLY", "NGLU",  "NGLN",  "NPRO", "NHIS", "NCYS", "NVAL", "NLEU",
        "NTHR", "NSER", "NLYS", "NMET", "NTYR", "NTRP",  "NPHE",  "NSEC", "NILE", "NCYX", "NCYM", "NHID",
        "NHIE", "NASH", "NGLH", "NHYP", "NLYN", "CNALA", "CNASP", "CASN", "CARG", "CGLY", "CGLU", "CGLN",
        "CPRO", "CHIS", "CCYS", "CVAL", "CLEU", "CTHR",  "CSER",  "CLYS", "CMET", "CTYR", "CTRP", "CPHE",
        "CSEC", "CILE", "CCYX", "CCYM", "CHID", "CHIE",  "CASH",  "CGLH", "CHYP", "CLYN"};
} // namespace biology
#endif
