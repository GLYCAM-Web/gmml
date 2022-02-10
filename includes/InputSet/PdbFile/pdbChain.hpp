//#ifndef GMML_INCLUDES_INPUTSET_PDBFILE_PDBCHAIN_HPP
//#define GMML_INCLUDES_INPUTSET_PDBFILE_PDBCHAIN_HPP
//
//#include <string>
//#include <iostream>
//#include <functional>
//#include "includes/InputSet/PdbFile/pdbResidue.hpp"
//
//namespace pdb
//{
//    class PdbChain
//    {
//        public:
//            //////////////////////////////////////////////////////////
//            //                       CONSTRUCTOR                    //
//            //////////////////////////////////////////////////////////
//            PdbChain(PdbResidue pdbResidue);
//            PdbChain(std::vector<PdbResidue> pdbResidues);
//            //////////////////////////////////////////////////////////
//            //                       ACCESSOR                       //
//            //////////////////////////////////////////////////////////
//            const std::string& GetChainId() const;
//            //////////////////////////////////////////////////////////
//            //                       MUTATOR                        //
//            //////////////////////////////////////////////////////////
//            void AddResidue(PdbResidue &pdbResidue);
//            //////////////////////////////////////////////////////////
//            //                       FUNCTIONS                      //
//            //////////////////////////////////////////////////////////
//
//            //////////////////////////////////////////////////////////
//            //                       DISPLAY FUNCTION               //
//            //////////////////////////////////////////////////////////
//            void Print(std::ostream& out = std::cerr) const;
//        private:
//            PdbResidue& GetFirstResidue() const;
//            //////////////////////////////////////////////////////////
//            //                       ATTRIBUTES                     //
//            //////////////////////////////////////////////////////////
//            std::vector<PdbResidue> pdbResidues_;
//    };
//}
//#endif // GMML_INCLUDES_INPUTSET_PDBFILE_PDBCHAIN_HPP
