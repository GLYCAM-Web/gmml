#ifndef PARAMETERFILEDIHEDRAL_HPP
#define PARAMETERFILEDIHEDRAL_HPP

#include <string>
#include <vector>
#include "../../common.hpp"

namespace ParameterFileSpace
{
    class ParameterFileDihedralTerm;
    class ParameterFileDihedral
    {
        public:
            /////////////////////////////// CONSTRUCTOR ////////////////////////////////
            ParameterFileDihedral();
            ParameterFileDihedral(const std::vector<std::string>& types, const ParameterFileDihedralTerm& term,
                                  double scee = gmml::kNotSet, double scnb = gmml::kNotSet,
                                  bool is_generic = false, bool is_improper = false);

            /////////////////////////// DISPLAY FUNCTION ///////////////////////////////
            void Print(std::ostream& out);

            /////////////////////////////// ATTRIBUTES /////////////////////////////////
            std::vector<std::string> types_;                // Atom types involved in the dihedral; Fill by the first column of the 5th or 6th section of the parameter file split by '-'
            std::vector<ParameterFileDihedralTerm> terms_;
            double scee_;                                   // scee coefficient; Extract from the 6th column of the 5th section or 5th column of the 6th section of the parameter file
            double scnb_;                                   // scnb coefficient; Extract from the 6th column of the 5th section or 5th column of the 6th section of the parameter file
            bool is_generic_;                               // Set based on the first column of the 5th or 6th section of the parameter file (if first or last atom type is 'X')
            bool is_improper_;                              // Set for 6th section of the parameter file
            // A sample of the 5th (dihedral) section of the parameter file : X -C -C -X    4   14.50        180.0             2.         Junmei et al, 1999
            // A sample of the 6th (improper dihedral) section of the parameter file : X -X -C -O          10.5         180.          2.           JCC,7,(1986),230
    };
}

#endif // PARAMETERFILEDIHEDRAL_HPP
