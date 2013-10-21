#ifndef PARAMETERFILEDIHEDRAL_HPP
#define PARAMETERFILEDIHEDRAL_HPP

#include <string>
#include <vector>
#include <iostream>
#include "../../common.hpp"

namespace ParameterFileSpace
{
    class ParameterFileDihedralTerm;
    class ParameterFileDihedral
    {
        public:
            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            ParameterFileDihedral();
            /*! \fn
              * Constructor with essential parameters
              * @param types Vector of four atom types involving in a dihedral
              * @param term Set of parameters according to a dihedral
              * @param scee
              * @param scnb
              * @param is_generic Boolean value to distinguish between generic and non-generic dihedrals
              * @param is_improper Boolean value to distinguish between proper and improper dihedrals
              */
            ParameterFileDihedral(const std::vector<std::string>& types, const ParameterFileDihedralTerm& term,
                                  double scee = gmml::kNotSet, double scnb = gmml::kNotSet,
                                  bool is_generic = false, bool is_improper = false);

            //////////////////////////////////////////////////////////
            //                     DISPLAY FUNCTIONS                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the file contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

            //////////////////////////////////////////////////////////
            //                         ATTRIBUTES                   //
            //////////////////////////////////////////////////////////
            std::vector<std::string> types_;                /*!< Atom types involved in the dihedral; Fill by the first column of the 5th or 6th section of the parameter file split by '-'*/
            std::vector<ParameterFileDihedralTerm> terms_;  /*!< Vector of terms that are assigned to a single dihedral*/
            double scee_;                                   /*!< scee coefficient; Extract from the 6th column of the 5th section or 5th column of the 6th section of the parameter file*/
            double scnb_;                                   /*!< scnb coefficient; Extract from the 6th column of the 5th section or 5th column of the 6th section of the parameter file*/
            bool is_generic_;                               /*!< Set based on the first column of the 5th or 6th section of the parameter file (if first or last atom type is 'X')*/
            bool is_improper_;                              /*!< Set for 6th section of the parameter file*/
            /*!< A sample of the 5th (dihedral) section of the parameter file : X -C -C -X    4   14.50        180.0             2.         Junmei et al, 1999*/
            /*!< A sample of the 6th (improper dihedral) section of the parameter file : X -X -C -O          10.5         180.          2.           JCC,7,(1986),230*/
    };
}

#endif // PARAMETERFILEDIHEDRAL_HPP
