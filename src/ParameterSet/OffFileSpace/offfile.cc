#include <fstream>
#include <iostream>
#include <iomanip>

#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/ParameterSet/OffFileSpace/offfile.hpp"

using OffFileSpace::OffFile;


    //////////////////////////////////////////////////////////
    //                       Constructor                    //
    //////////////////////////////////////////////////////////
    OffFileSpace::OffFile::OffFile(){}

    OffFileSpace::OffFile::~OffFile(){}

    //////////////////////////////////////////////////////////
    //                           ACCESSOR                   //
    //////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////
    //                         MUTATORS                     //
    //////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////
    //                         FUNCTIONS                    //
    //////////////////////////////////////////////////////////

    void OffFileSpace::OffFile::Write(std::string file_name, int CoordinateIndex, MolecularModeling::Assembly* assembly)
    {

            std::ofstream out_file;
            out_file.open(file_name.c_str());
            ResidueVector residues = assembly->GetResidues();
            unit_name_=assembly->GetName();
            out_file << "!!index array str" << std::endl;
            out_file  << unit_name_ << std::endl;
            WriteAtomSection(out_file,residues);
            WriteAtomPertInfoSection(out_file,residues);
            WriteBoundBoxSection(out_file,assembly);
            WriteChildSequenceSection(out_file,residues);
            WriteConnectSection(out_file,residues);
            WriteConnectivitySection(out_file,residues);
            WriteHierarchySection(out_file,residues);
            WriteNameSection(out_file);
            WritePositionSection(out_file,residues, CoordinateIndex);
            WriteResidueConnectSection(out_file,residues);
            WriteResiduesSection(out_file,residues);
            WriteSolventCapSection(out_file);
            WriteVelocitiesSection(out_file,residues);

            out_file.close();
    }

    void OffFileSpace::OffFile::WriteAtomSection(std::ofstream &stream, ResidueVector assembly_residues)
    {

        const std::string FLAG = "131072";
        int residue_count=0;
        int atom_count=0;
        stream << "!entry." << unit_name_<< ".unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg" << std::endl;
           
           for(ResidueVector::iterator it = assembly_residues.begin(); it != assembly_residues.end(); it++)
            {

                MolecularModeling::Residue* residue = (*it);
                residue_count++;
                AtomVector all_atoms_of_residue = residue->GetAtoms();
                for(AtomVector::iterator it = all_atoms_of_residue.begin(); it != all_atoms_of_residue.end(); it++)
                {
                    MolecularModeling::Atom* atom = *it;
                    atom_count++;
                stream << "\"" << atom->GetName() << "\" " << "\"" << atom->GetAtomType() << "\" " << "0" << " " << residue_count << " " << FLAG << " "
                << atom->GetIndex() << " " << " 0 " << " " << std::fixed << atom->GetCharge() << std::endl;
                }
            }
    }
    void OffFileSpace::OffFile::WriteAtomPertInfoSection(std::ofstream& stream, ResidueVector assembly_residues)
    {

        stream << "!entry." << unit_name_ << ".unit.atomspertinfo table  str pname  str ptype  int ptypex  int pelmnt  dbl pchg" << std::endl;

        for(ResidueVector::iterator it = assembly_residues.begin(); it != assembly_residues.end(); it++)
            {
                MolecularModeling::Residue* residue = (*it);
                AtomVector all_atoms_of_residue = residue->GetAtoms();
                for(AtomVector::iterator it = all_atoms_of_residue.begin(); it != all_atoms_of_residue.end(); it++)
                {
                     MolecularModeling::Atom* current_atom = *it;
                     stream<<current_atom->GetName()<<" "<<current_atom->GetAtomType()<<" "<<current_atom->GetIndex()<<" "<<current_atom->GetElementSymbol()<<" "<<current_atom->GetCharge()<<std::endl;
                }
            }
    }

    void OffFileSpace::OffFile::WriteBoundBoxSection(std::ofstream& stream, MolecularModeling::Assembly* assembly)
    {
            stream << "!entry." << unit_name_ << ".unit.boundbox array dbl" << std::endl;
            stream << "-1.000000" << std::endl;
            stream << "0.0" << std::endl;
            stream << "0.0" << std::endl;
            stream << "0.0" << std::endl;
            stream << "0.0" << std::endl;
   
    }
    void OffFileSpace::OffFile::WriteChildSequenceSection(std::ofstream& stream, ResidueVector assembly_residues)
    {

            stream << "!entry." << unit_name_ << ".unit.childsequence single int" << std::endl;
            unsigned int residue_count = assembly_residues.size();
            stream<<residue_count+1<<std::endl;

    }
    void OffFileSpace::OffFile::WriteConnectSection(std::ofstream& stream, ResidueVector assembly_residues)
    {
            stream << "!entry." << unit_name_ << ".unit.connect array int" << std::endl;
            stream << "0" << std::endl;
            stream << "0" << std::endl;
    }
    void OffFileSpace::OffFile::WriteConnectivitySection(std::ofstream& stream, ResidueVector assembly_residues)
    {
         stream << "!entry." << unit_name_ << ".unit.connectivity table  int atom1x  int atom2x  int flags" << std::endl;

         for(ResidueVector::iterator it = assembly_residues.begin(); it != assembly_residues.end(); it++)
            {
                MolecularModeling::Residue* residue = (*it);
                AtomVector all_atoms_of_residue = residue->GetAtoms();
                for(AtomVector::iterator it = all_atoms_of_residue.begin(); it != all_atoms_of_residue.end(); it++)
                {
                    MolecularModeling::Atom* atom = *it;
                    AtomVector bonded_atoms = atom->GetNode()->GetNodeNeighbors();
                   for(AtomVector::iterator it = bonded_atoms.begin(); it != bonded_atoms.end(); it++)
                    {
                        MolecularModeling::Atom* bonded_atom = *it;
                        stream << atom->GetIndex() << " " << bonded_atom->GetIndex() << " " << "1" << std::endl;
                    }

                }
            }
    }
    void OffFileSpace::OffFile::WriteHierarchySection(std::ofstream& stream, ResidueVector assembly_residues)
    {
        stream << "!entry." << unit_name_ << ".unit.hierarchy table  str abovetype  int abovex  str belowtype  int belowx" << std::endl;
        int residue_count=0;
        for(ResidueVector::iterator it = assembly_residues.begin(); it != assembly_residues.end(); it++)
            {   
                residue_count++;
                MolecularModeling::Residue* residue = (*it);
                stream<<"U"<<" "<<0<<" R"<<" "<<residue_count<<std::endl;
                AtomVector all_atoms_of_residue = residue->GetAtoms();
                for(AtomVector::iterator it = all_atoms_of_residue.begin(); it != all_atoms_of_residue.end(); it++)
                {

                    MolecularModeling::Atom* atom = *it;
                    stream<<"R"<<" "<<residue_count<<" A"<<" "<<atom->GetIndex() <<std::endl;
                }
            }
        stream << std::endl;
    }
    void OffFileSpace::OffFile::WriteNameSection(std::ofstream& stream)
    {

            stream << "!entry." << unit_name_ << ".unit.name single str" << std::endl;
            stream << "\"" << unit_name_ << "\"" << std::endl;
    }
    void OffFileSpace::OffFile::WritePositionSection(std::ofstream& stream, ResidueVector assembly_residues, int CoordinateIndex)
    {
        stream << "!entry." << unit_name_<< ".unit.positions table  dbl x  dbl y  dbl z" << std::endl;

            for(ResidueVector::iterator it = assembly_residues.begin(); it != assembly_residues.end(); it++)
            {

                MolecularModeling::Residue* residue = (*it);
                AtomVector all_atoms_of_residue = residue->GetAtoms();
                for(AtomVector::iterator it = all_atoms_of_residue.begin(); it != all_atoms_of_residue.end(); it++)
                {
                    MolecularModeling::Atom* atom = *it;
                    GeometryTopology::Coordinate coordinate = atom->GetCoordinates().at(CoordinateIndex);
                    stream << std::fixed << coordinate.GetX() << " " << std::fixed << coordinate.GetY() << " " << std::fixed << coordinate.GetZ() << std::endl;
                }
            }

    }
    void OffFileSpace::OffFile::WriteResidueConnectSection(std::ofstream& stream, ResidueVector assembly_residues)
    {
            
            stream << "!entry." << unit_name_ << ".unit.residueconnect table  int c1x  int c2x  int c3x  int c4x  int c5x  int c6x" << std::endl;
        for(ResidueVector::iterator it = assembly_residues.begin(); it != assembly_residues.end(); it++)
            {

                MolecularModeling::Residue* residue = (*it);

                AtomVector head_atoms_of_residue= residue->GetHeadAtoms();
                AtomVector tail_atoms_of_residue= residue->GetTailAtoms();

                int c1x[head_atoms_of_residue.size()];
                int c2x[tail_atoms_of_residue.size()];

                int c1x_count=0;
                int c2x_count=0;
                for(AtomVector::iterator it = head_atoms_of_residue.begin(); it != head_atoms_of_residue.end(); it++)
                {
                    MolecularModeling::Atom* atom = *it;
                    c1x[c1x_count]=atom->GetIndex();
                    c1x_count++;
                   
                }

                 for(AtomVector::iterator it = tail_atoms_of_residue.begin(); it != tail_atoms_of_residue.end(); it++)
                {
                   MolecularModeling::Atom* atom = *it;
                    c2x[c2x_count]=atom->GetIndex();
                    c2x_count++;
                }

                for(int i=0;i<c1x_count || i<c2x_count;i++)
                {
                    stream<<c1x[i]<<" "<<c2x[i]<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<std::endl;
                }
                        
            }
    }
    

    void OffFileSpace::OffFile::WriteResiduesSection(std::ofstream& stream, ResidueVector assembly_residues)
    {
            std::string name;
            unsigned int seq=0;
            unsigned int childseq;
            unsigned int startatomx;
            std::string restype;
            unsigned int imagingx;

            stream << "!entry." <<unit_name_<<".unit.residues table  str name  int seq  int childseq  int startatomx  str restype  int imagingx" << std::endl;
            for(ResidueVector::iterator it = assembly_residues.begin(); it != assembly_residues.end(); it++)
                {
                    MolecularModeling::Residue* residue = (*it);
                    name= residue->GetName();
                    seq++;
                    childseq= residue->GetAtoms().size()+1;
                    AtomVector head_atoms_of_residue = residue->GetHeadAtoms();
                    MolecularModeling::Atom* current_atom=head_atoms_of_residue[0];
                    startatomx= current_atom->GetIndex();
                    restype="?";
                    imagingx=0;

                    stream<<name<<" "<<seq<<" "<<childseq<<" "<<startatomx<<" "<<restype<<" "<<imagingx<<std::endl;
                }
    }

    void OffFileSpace::OffFile::WriteSolventCapSection(std::ofstream& stream)
    {
        stream <<"!entry."<<unit_name_<<".unit.solventcap array dbl" <<std::endl;
        stream<<"-1.000000"<<std::endl;
        stream<<"0.0"<<std::endl;
        stream<<"0.0"<<std::endl;
        stream<<"0.0"<<std::endl;
        stream<<"0.0"<<std::endl;
    }

    void OffFileSpace::OffFile::WriteVelocitiesSection(std::ofstream& stream, ResidueVector assembly_residues)
    {
        stream << "!entry." <<unit_name_ <<".unit.velocities table  dbl x  dbl y  dbl z" << std::endl;
        for(ResidueVector::iterator it = assembly_residues.begin(); it != assembly_residues.end(); it++)
            {
                MolecularModeling::Residue* residue = (*it);
                AtomVector all_atoms_of_residue = residue->GetAtoms();
                for(AtomVector::iterator it = all_atoms_of_residue.begin(); it != all_atoms_of_residue.end(); it++)
                {
                     stream<<"0.0"<<"0.0"<<"0.0"<<std::endl;
                }
            }
    }

    //////////////////////////////////////////////////////////
    //                     DISPLAY FUNCTIONS                //
    //////////////////////////////////////////////////////////
