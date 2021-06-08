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
              * Constructor to initialize the attributes of a dihedral that exists in a prep file
              * @param types Vector of four atom types involving in a dihedral
              * @param term Set of parameters according to a dihedral
              * @param scee A float number that indicates scee coefficient of a dihedral
              * @param scnb A float number that indicates scnb coefficient of a dihedral
              * @param is_generic Boolean value to distinguish between generic and non-generic dihedrals
              * @param is_improper Boolean value to distinguish between proper and improper dihedrals
              */
            ParameterFileDihedral(const std::vector<std::string>& types, const ParameterFileDihedralTerm& term,
                                  double scee = gmml::dNotSet, double scnb = gmml::dNotSet,
                                  bool is_generic = false, bool is_improper = false);

            //////////////////////////////////////////////////////////
            //                           ACCESSOR                   //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to the atom types involved in an dihedral in the current object
              * @return types_ attribute of the current object of this class
              */
            std::vector<std::string> GetTypes();
            /*! \fn
              * An accessor function in order to access to the terms attribute of the current object
              * The attribute is set by the contents of the given file
              * @return terms_ of the current object of this class
              */
            std::vector<ParameterFileDihedralTerm> GetTerms();
            /*! \fn
              * An accessor function in order to access to scee coefficient of the current object
              * The attribute is set by the contents of the given file
              * @return scee_ of the current object of this class
              */
            double GetScee();
            /*! \fn
              * An accessor function in order to access to scnb coefficient of the current object
              * The attribute is set by the contents of the given file
              * @return scnb_ of the current object of this class
              */
            double GetScnb();
            /*! \fn
              * An accessor function in order to access to is_generic_ attribute of the current object
              * The attribute is set by the contents of the given file
              * @return is_generic_ of the current object of this class
              */
            bool GetIsGeneric();
            /*! \fn
              * An accessor function in order to access to is_improper_ attribute of the current object
              * The attribute is set by the contents of the given file
              * @return is_improper_ of the current object of this class
              */
            bool GetIsImproper();
/**@}*/
            //////////////////////////////////////////////////////////
            //                           MUTATOR                    //
            //////////////////////////////////////////////////////////
/** @addtogroup Manipulators
*  @{
*/
            /*! \fn
              * A mutator function in order to set the list of atom types of the current object
              * Set the types_ attribute of the current dihedral
              * @param types A list of types of atoms of the current object
              */
            void SetTypes(const std::vector<std::string> types);
            /*! \fn
              * A mutator function in order to set the list of terms of the current object
              * Set the terms_ attribute of the current dihedral
              * @param terms A list of terms of the current object
              */
            void SetTerms(const std::vector<ParameterFileDihedralTerm> terms);
            /*! \fn
              * A mutator function in order to add a new term belonging to the current object
              * Add a new entry to the terms_ attribute of the current dihedral
              * @param term A new term belonging to the current dihedral
              */
            void AddTerm(ParameterFileDihedralTerm term);
            /*! \fn
              * A mutator function in order to set the scee coefficient of the current object
              * Set the scee_ attribute of the current dihedral
              * @param scee The scee of the current object
              */
            void SetScee(double scee);
            /*! \fn
              * A mutator function in order to set the scnb coefficient of the current object
              * Set the scnb_ attribute of the current dihedral
              * @param scnb The scnb of the current object
              */
            void SetScnb(double scnb);
            /*! \fn
              * A mutator function in order to set the is_generic attribute of the current object
              * Set the is_generic_ attribute of the current dihedral
              * @param is_generic The is_generic attribute of the current object
              */
            void SetIsGeneric(bool is_generic);
            /*! \fn
              * A mutator function in order to set the is_improper attribute of the current object
              * Set the is_improper_ attribute of the current dihedral
              * @param is_improper The is_improper attribute of the current object
              */
            void SetIsImproper(bool is_improper);
/**@}*/
            //////////////////////////////////////////////////////////
            //                     DISPLAY FUNCTIONS                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the file contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
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
