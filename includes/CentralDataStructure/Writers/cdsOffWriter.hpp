#ifndef INCLUDES_CENTRALDATASTRUCTURE_WRITERS_CDSOFFWRITER_HPP_
#define INCLUDES_CENTRALDATASTRUCTURE_WRITERS_CDSOFFWRITER_HPP_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "includes/Abstract/absResidue.hpp"

namespace cds
{

inline std::string getOffType(const Abstract::ResidueType queryType)
{
	if ( queryType == Abstract::ResidueType::Protein )
	{
		return "p";
	}
	if ( queryType == Abstract::ResidueType::Solvent )
	{
		return "w";
	}
	return "?";
}

template <typename residueT>
void WriteOffFileUnit(std::vector<residueT*> residues, std::ostream& stream, const std::string unitName)
{
	//WriteAtomSection
	const std::string FLAG = "131072";
	stream << "!entry." << unitName << ".unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg" << std::endl;
	for(auto &residue : residues)
	{
		unsigned int atomNumberInResidue = 1;
		for(auto &atom : residue->getAtoms())
		{
			stream << " \"" << atom->getName() << "\" " << "\"" << atom->getType() << "\" " << "0" << " " << residue->getIndex() << " " << FLAG << " "
					<< atomNumberInResidue << " " << atom->getAtomicNumber() << " " << std::fixed << atom->getCharge() << std::endl;

			atomNumberInResidue++;
		}
	}
	//WriteAtomPertInfoSection
	stream << "!entry." << unitName << ".unit.atomspertinfo table  str pname  str ptype  int ptypex  int pelmnt  dbl pchg" << std::endl;
	for(auto &residue : residues)
	{
		for(auto &atom : residue->getAtoms())
		{
			stream << " \"" << atom->getName() << "\" " << "\"" << atom->getType() << "\" " << 0 << " " << -1 << " " << 0.0 << std::endl;
		}
	}
	//WriteBoundBoxSection
	stream << "!entry." << unitName << ".unit.boundbox array dbl" << std::endl;
	stream << " " << "-1.000000" << std::endl;
	stream << " " << 0.0 << std::endl;
	stream << " " << 0.0 << std::endl;
	stream << " " << 0.0 << std::endl;
	stream << " " << 0.0 << std::endl;
	//WriteChildSequenceSection
	stream << "!entry." << unitName << ".unit.childsequence single int" << std::endl;
	stream << " " << residues.size() + 1 << std::endl;

	//WriteConnectSection
	// ToDo: this should be the atom number of the atoms that are bonded to the unit before and after. They were set to zero in the original, as below. Is this ok?
	stream << "!entry." << unitName << ".unit.connect array int" << std::endl;
	stream << " " << 0 << std::endl;
	stream << " " << 0 << std::endl;
	//WriteConnectivitySection
	stream << "!entry." << unitName << ".unit.connectivity table  int atom1x  int atom2x  int flags" << std::endl;
	for(auto &residue : residues)
	{
		for(auto &atom : residue->getAtoms())
		{
			for(auto &neighbor : atom->getChildren())
			{
				stream << " " << atom->getIndex() << " " << neighbor->getIndex() << " " << 1 << std::endl;
			}
		}
	}
	//WriteHierarchySection
	stream << "!entry." << unitName << ".unit.hierarchy table  str abovetype  int abovex  str belowtype  int belowx" << std::endl;
	for(auto &residue : residues)
	{
		stream << " \"" << "U" << "\"" << " " << 0 << " " << "\"" << "R" << "\"" << " " << residue->getIndex() << std::endl;
		for (int atomCounter = 0; atomCounter < residue->getAtoms().size(); ++atomCounter)
		{
			stream << " \"" << "R" << "\"" << " " << residue->getIndex() << " " << "\"" << "A" << "\"" << " " << atomCounter << std::endl;
		}
	}
	//WriteNameSection
	stream << "!entry." << unitName << ".unit.name single str" << std::endl;
	stream << " \"" << unitName << "\"" << std::endl;
	//WritePositionSection
	stream << "!entry." << unitName << ".unit.positions table  dbl x  dbl y  dbl z" << std::endl;
	for(auto &residue : residues)
	{
		for(auto &atom : residue->getAtoms())
		{
			stream << " " << std::fixed << atom->getCoordinate()->GetX() << " "<< std::fixed << atom->getCoordinate()->GetY() << " " << std::fixed << atom->getCoordinate()->GetZ() << std::endl;
		}
	}
	//WriteResidueConnectSection
	stream << "!entry." << unitName << ".unit.residueconnect table  int c1x  int c2x  int c3x  int c4x  int c5x  int c6x" << std::endl;
	for(auto &residue : residues)
	{
		for(auto &atom : residue->getAtomsConnectedToOtherResidues())
		{
			stream << " " << atom->getIndex();
		}
		int columnsWithZero = 6 - residue->getAtomsConnectedToOtherResidues().size();
		for (int i = 0; i < columnsWithZero; ++i)
		{
			stream << " " << "0";
		}
		stream << std::endl;
	}
	//WriteResiduesSection
	stream << "!entry." << unitName << ".unit.residues table  str name  int seq  int childseq  int startatomx  str restype  int imagingx" << std::endl;
	for(auto &residue : residues)
	{
		unsigned int childseq = residue->getAtoms().size() + 1;
		unsigned int startatomx = residue->getAtoms().front()->getIndex();
		std::string restype = cds::getOffType(residue->GetType());
		unsigned int imagingx = 0;
		stream << " \"" << residue->getName() << "\"" << " " << residue->getIndex() << " " << childseq << " " << startatomx << " " << "\""<<restype <<"\""<< " " << imagingx << std::endl;
	}
	//WriteSolventCapSection
	stream << "!entry." << unitName << ".unit.solventcap array dbl" << std::endl;
	stream << " " << "-1.000000" << std::endl;
	stream << " " <<"0.0" << std::endl;
	stream << " " <<"0.0" << std::endl;
	stream << " " <<"0.0" << std::endl;
	stream << " " <<"0.0" << std::endl;
	//WriteVelocitiesSection
	stream << "!entry." << unitName << ".unit.velocities table  dbl x  dbl y  dbl z" << std::endl;
	int totalNumberOfAtoms = 0;
	for(auto &residue : residues)
	{
		for (auto &atom : residue->getAtoms())
		{ // Maybe later we'll want to deal with atom velocities...
			stream << " " << "0.0" << " " << "0.0" << " " << "0.0" << std::endl;
		}
	}
	return;
}

template <typename residueT>
void WriteResiduesToOffFile(std::vector<residueT*> residues, std::ostream& stream)
{ // For writing each residue separately
	stream << "!!index array str" << std::endl;
	for(auto &residue : residues)
	{
		stream << " \"" << residue->getName() << "\"" << std::endl;
	}
	for(auto &residue : residues)
	{
		cds::WriteOffFileUnit(std::vector<residueT*>{residue}, stream, residue->getName());
	}
	return;
}

template <typename residueT>
void WriteMoleculeToOffFile(std::vector<residueT*> residues, std::ostream& stream, const std::string unitName)
{ // For writing residues together as a molecule
	stream << "!!index array str" << std::endl;
	stream << " \"" << unitName << "\"" << std::endl;
	cds::WriteOffFileUnit(residues, stream, unitName);
	return;
}

} // namespace

#endif
