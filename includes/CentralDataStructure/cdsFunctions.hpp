//#ifndef INCLUDES_CENTRALDATASTRUCTURE_FUNCTIONS_HPP
//#define INCLUDES_CENTRALDATASTRUCTURE_FUNCTIONS_HPP
//
//#include <iostream>
//#include <string>
//#include <utility>
//
//namespace cds
//{
//
//// From https://www.modernescpp.com/index.php/perfect-forwarding.
//// More background here https://eli.thegreenplace.net/2014/perfect-forwarding-and-universal-references-in-c
//// Or find info on perfect forwarding.
//// I'm using it to forward constructor arguments into cds classes so e.g. when I want to construct a pdbAtom in the pdbResidue<pdbAtom> class
//// the vector of pdbAtom in pdbResidue is inaccessible to pdbResidue as it's inherited from cdsResidue. But I can use this method to have
//// the cdsResidue (base class) create the pdbAtom, even though it doesn't "know" what it is until compile time when the templates get resolved.
//// This issue is that I want to use it from within the cdsResidue class and making it a member function seems fraught with issues. Feck a deck.
//
//template <typename T, typename ... Args>
//T create(Args&& ... args)
//{
//    return T(std::forward<Args>(args)...);
//}
//
//} // namespace
//
//#endif // CDS Functions
//
