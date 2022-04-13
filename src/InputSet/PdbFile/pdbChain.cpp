#include "includes/InputSet/PdbFile/pdbChain.hpp"
#include "includes/InputSet/PdbFile/pdbResidue.hpp"
//#include "includes/InputSet/PdbFile/atomRecord.hpp"
#include "includes/GeometryTopology/coordinate.hpp"
#include "includes/GeometryTopology/geometrytopology.hpp" // get_cartesian_point_from_internal_coords
#include "includes/CodeUtils/logging.hpp"

using pdb::PdbChain;
////////////////////////////////////////////////////////////
////                       CONSTRUCTOR                    //
////////////////////////////////////////////////////////////
//PdbChain::PdbChain(PdbResidue* pdbResidue)
//{
//    this->AddResidue(pdbResidue);
//}
//PdbChain::PdbChain(std::vector<PdbResidue*> pdbResidues)
//{
//    pdbResidues_ = pdbResidues;
//}
////////////////////////////////////////////////////////////
////                         ACCESSOR                     //
////////////////////////////////////////////////////////////
//pdb::PdbResidue* PdbChain::GetFirstResidue() const
//{
//    return pdbResidues_.front();
//}
//const std::string& PdbChain::GetChainId() const
//{
//    return this->GetFirstResidue()->GetChainId();
//}
////////////////////////////////////////////////////////////
////                    FUNCTIONS                         //
////////////////////////////////////////////////////////////
void PdbChain::InsertCap(const PdbResidue& refResidue, const std::string& type)
{
    // This approach is bad, should be using templates. When parameter manager is good we can use that to remove the get_carestian stuff
    using GeometryTopology::Coordinate;
    if (type == "NHCH3") // NME
    {
//        int sequenceNumber = refResidue.GetSequenceNumber() + 1; // Single gaps will end up with the same ACE NME resid numbers. Otherwise good.
        const Coordinate& cCoordProtein = refResidue.FindAtom("C")->getCoordinate();
        const Coordinate& caCoordProtein = refResidue.FindAtom("CA")->GetCoordinate();
        const Coordinate& oCoordProtein = refResidue.FindAtom("O")->GetCoordinate();
        Coordinate nCoordNME = GeometryTopology::get_cartesian_point_from_internal_coords(oCoordProtein, caCoordProtein, cCoordProtein, 120.0, 180.0, 1.4);
        Coordinate hCoordNME = GeometryTopology::get_cartesian_point_from_internal_coords(oCoordProtein, caCoordProtein, nCoordNME, 109.0, 180.0, 1.0);
        Coordinate ch3CoordNME = GeometryTopology::get_cartesian_point_from_internal_coords(caCoordProtein, cCoordProtein, nCoordNME, 125.0, 180.0, 1.48);
        Coordinate hh31CoordNME = GeometryTopology::get_cartesian_point_from_internal_coords(hCoordNME, nCoordNME, ch3CoordNME, 109.0, 180.0, 1.09);
        Coordinate hh32CoordNME = GeometryTopology::get_cartesian_point_from_internal_coords(hCoordNME, nCoordNME, ch3CoordNME, 109.0, 60.0, 1.09);
        Coordinate hh33CoordNME = GeometryTopology::get_cartesian_point_from_internal_coords(hCoordNME, nCoordNME, ch3CoordNME, 109.0, -60.0, 1.09);
        //AtomRecordIterator atomPosition = this->GetCoordinateSection().FindPositionOfAtom(refResidue.GetLastAtom());
        PdbResidue *newNMEResidue = this->createNewResidue("NME", refResidue);
        newNMEResidue->createAtom("N", nCoordNME);
        newNMEResidue->createAtom("H", hCoordNME);
        newNMEResidue->createAtom("CH3", ch3CoordNME);
        newNMEResidue->createAtom("HH31", hh31CoordNME);
        newNMEResidue->createAtom("HH32", hh32CoordNME);
        newNMEResidue->createAtom("HH33", hh33CoordNME);
        newNMEResidue->AddTerCard();
//        atomPosition = this->GetCoordinateSection().CreateNewAtomRecord("N", "NME", sequenceNumber, nCoordNME, refResidue.GetChainId(), refResidue.GetModelNumber(), atomPosition);
//        atomPosition = this->GetCoordinateSection().CreateNewAtomRecord("H", "NME", sequenceNumber, hCoordNME, refResidue.GetChainId(), refResidue.GetModelNumber(), atomPosition);
//        atomPosition = this->GetCoordinateSection().CreateNewAtomRecord("CH3", "NME", sequenceNumber, ch3CoordNME, refResidue.GetChainId(), refResidue.GetModelNumber(), atomPosition);
//        atomPosition = this->GetCoordinateSection().CreateNewAtomRecord("HH31", "NME", sequenceNumber, hh31CoordNME, refResidue.GetChainId(), refResidue.GetModelNumber(), atomPosition);
//        atomPosition = this->GetCoordinateSection().CreateNewAtomRecord("HH32", "NME", sequenceNumber, hh32CoordNME, refResidue.GetChainId(), refResidue.GetModelNumber(), atomPosition);
//        atomPosition = this->GetCoordinateSection().CreateNewAtomRecord("HH33", "NME", sequenceNumber, hh33CoordNME, refResidue.GetChainId(), refResidue.GetModelNumber(), atomPosition);
    }
    else if (type == "COCH3") // ACE
    {
//        int sequenceNumber = refResidue.GetSequenceNumber() - 1; // Single gaps will end up with the same ACE NME resid numbers. Otherwise good.
        // These are the atoms in residue that I use to build the ACE out from.
        const Coordinate& cCoordProtein = refResidue.FindAtom("C")->GetCoordinate();
        const Coordinate& caCoordProtein = refResidue.FindAtom("CA")->GetCoordinate();
        const Coordinate& nCoordProtein = refResidue.FindAtom("N")->GetCoordinate();
        // This is bad, should use templates loaded from lib/prep file instead.
        Coordinate cCoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(cCoordProtein, caCoordProtein, nCoordProtein, 120.0, -130.0, 1.4);
        Coordinate oCoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(caCoordProtein, nCoordProtein, cCoordACE, 120.0, 0.0, 1.23);
        Coordinate ch3CoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(caCoordProtein, nCoordProtein, cCoordACE, 125.0, 180.0, 1.48);
        Coordinate hh31CoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(oCoordACE, cCoordACE, ch3CoordACE, 109.0, 180.0, 1.09);
        Coordinate hh32CoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(oCoordACE, cCoordACE, ch3CoordACE, 109.0, 60.0, 1.09);
        Coordinate hh33CoordACE = GeometryTopology::get_cartesian_point_from_internal_coords(oCoordACE, cCoordACE, ch3CoordACE, 109.0, -60.0, 1.09);
        // Ok this next bit is convoluted, but I look up the position of the first atom in the protein residue and insert the new Atom before it, and get passed back the position of the newly created atom, so I can use that when creating the next one and so on.
      //  AtomRecordIterator atomPosition = this->GetCoordinateSection().FindPositionOfAtom(refResidue.GetFirstAtom());

        // With ACE we want to insert before the residue, so I'm finding the residue before here:
        auto refPosition = this->GetCoordinateSection().FindPositionOfResidue(&refResidue);
        --refPosition;
        PdbResidue* previousResidue = (*refPosition).get(); // Its an iterator to a unique ptr, so deref and get the raw. Ugh.
        PdbResidue *newACEResidue = this->GetCoordinateSection().CreateNewResidue("ACE", "C", cCoordACE, *previousResidue);
        newACEResidue->CreateAtom("O", oCoordACE);
        newACEResidue->CreateAtom("CH3", ch3CoordACE);
        newACEResidue->CreateAtom("HH31", hh31CoordACE);
        newACEResidue->CreateAtom("HH32", hh32CoordACE);
        newACEResidue->CreateAtom("HH33", hh33CoordACE);

//        --atomPosition; // Want to insert before the first atom, inserting at begin() position is fine.
//        atomPosition = this->GetCoordinateSection().CreateNewAtomRecord("C", "ACE", sequenceNumber, cCoordACE, refResidue.GetChainId(), refResidue.GetModelNumber(), atomPosition);
//        atomPosition = this->GetCoordinateSection().CreateNewAtomRecord("O", "ACE", sequenceNumber, oCoordACE, refResidue.GetChainId(), refResidue.GetModelNumber(), atomPosition);
//        atomPosition = this->GetCoordinateSection().CreateNewAtomRecord("CH3", "ACE", sequenceNumber, ch3CoordACE, refResidue.GetChainId(), refResidue.GetModelNumber(), atomPosition);
//        atomPosition = this->GetCoordinateSection().CreateNewAtomRecord("HH31", "ACE", sequenceNumber, hh31CoordACE, refResidue.GetChainId(), refResidue.GetModelNumber(), atomPosition);
//        atomPosition = this->GetCoordinateSection().CreateNewAtomRecord("HH32", "ACE", sequenceNumber, hh32CoordACE, refResidue.GetChainId(), refResidue.GetModelNumber(), atomPosition);
//        atomPosition = this->GetCoordinateSection().CreateNewAtomRecord("HH33", "ACE", sequenceNumber, hh33CoordACE, refResidue.GetChainId(), refResidue.GetModelNumber(), atomPosition);
        gmml::log(__LINE__, __FILE__, gmml::INF, "Created ACE residue: " + newACEResidue->GetId());
    }
}

////////////////////////////////////////////////////////////
////                          MUTATOR                     //
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
////                      DISPLAY FUNCTION                //
////////////////////////////////////////////////////////////
//void PdbChain::Print(std::ostream &out) const
//{
//    out << "Printing chain " << this->GetChainId() << "\n";
//    for(auto &residue : pdbResidues_)
//    {
//        residue->Print();
//    }
//}
