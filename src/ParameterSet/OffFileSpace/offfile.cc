#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/ParameterSet/OffFileSpace/offfile.hpp"
#include "../../../includes/ParameterSet/OffFileSpace/offfileresidue.hpp"
#include "../../../includes/ParameterSet/OffFileSpace/offfileatom.hpp"

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
             std::cout<<"aaa"<<std::endl;
            unit_name_=assembly->GetName();
              std::cout<<"bb"<<std::endl;
            out_file << "!!index array str" << std::endl;
            out_file << " \"" << unit_name_ << "\"" << std::endl;

             std::cout<<"Before PopulateOffFileResiduesFromAssembly"<<std::endl;
            OffFileResidueVector off_file_residues=PopulateOffFileResiduesFromAssembly(residues,CoordinateIndex);
            std::cout<<"After PopulateOffFileResiduesFromAssembly"<<std::endl;
            WriteAtomSection(out_file,this->off_file_residues_);
//            WriteAtomPertInfoSection(out_file,residues);
//            WriteBoundBoxSection(out_file,assembly);
//            WriteChildSequenceSection(out_file,residues);
//            WriteConnectSection(out_file,residues);
//            WriteConnectivitySection(out_file,residues);
//            WriteHierarchySection(out_file,residues);
//            WriteNameSection(out_file);
//            WritePositionSection(out_file,residues, CoordinateIndex);
//            WriteResidueConnectSection(out_file,residues);
//            WriteResiduesSection(out_file,residues);
//            WriteSolventCapSection(out_file);
//            WriteVelocitiesSection(out_file,residues);
            out_file.close();
    }

    OffFileSpace::OffFile::OffFileResidueVector OffFileSpace::OffFile::PopulateOffFileResiduesFromAssembly(ResidueVector assembly_residues,int CoordinateIndex)
    {
            int ResidueIndex=0;
            for(ResidueVector::iterator it = assembly_residues.begin(); it != assembly_residues.end(); it++)
            {     std::cout<<"000000"<<std::endl;
                ResidueIndex++;
                OffFileResidue* off_file_residue = new OffFileResidue();

                MolecularModeling::Residue* assembly_residue = (*it);
                std::string name=assembly_residue->GetName();
                off_file_residue->SetName(name);
                off_file_residue->SetListingIndex(ResidueIndex);
                
                  AtomVector all_atoms_of_residue = assembly_residue->GetAtoms();
                std::cout<<"1111111"<<std::endl;
                int AtomIndex=0;
                for(AtomVector::iterator it = all_atoms_of_residue.begin(); it != all_atoms_of_residue.end(); it++)
                {
                    AtomIndex++;
                    MolecularModeling::Atom* residue_atom = *it;
                    OffFileAtom* off_file_atom= new OffFileAtom();
                    std::cout<<"222222"<<std::endl;
                    //setting values from residue atom to off file atom
                    off_file_atom->SetName(residue_atom->GetName());
                    off_file_atom->SetType(residue_atom->MolecularDynamicAtom::GetAtomType());
                    off_file_atom->SetResidueIndex(ResidueIndex);
                    off_file_atom->SetAtomIndex(AtomIndex);
                    std::cout<<"333333"<<std::endl;
                    int AtomicNumber;
                    std::string atom_element_symbol= residue_atom->GetElementSymbol();
                   int size_of_lookup_map = sizeof(gmml::ElementAtributes::Elements) / sizeof(gmml::ElementAtributes::Elements[0]);
                       for (int i = 0; i < size_of_lookup_map; i++){
                                  gmml::ElementAtributes::ElementAttributeInfo entry = gmml::ElementAtributes::Elements[i];
                                  if (atom_element_symbol.compare(entry.elment_type_) == 0){
                                      AtomicNumber = entry.atomic_number_;
                                  }
                              }

                     std::cout<<"44444"<<std::endl;
                    off_file_atom->SetAtomicNumber(AtomicNumber);
                    off_file_atom->SetCoordinate(residue_atom->GetCoordinates().at(CoordinateIndex));
                    off_file_atom->SetAtomCharge(residue_atom->GetCharge());
                    this->atom_index_map_[residue_atom->GetIndex()]=AtomIndex;

                    off_file_residue->AddAtom(off_file_atom);
                }

                //for adding head and tail atoms of residues
                AtomVector head_atoms_of_residue= assembly_residue->GetHeadAtoms();
                AtomVector tail_atoms_of_residue= assembly_residue->GetTailAtoms();

                std::map<int,int>::iterator it1,it2;
                int head_temp_index= head_atoms_of_residue.at(0)->GetIndex();
                it1 = atom_index_map_.find(head_temp_index);
                    if(it1 != atom_index_map_.end()){
                          off_file_residue->SetHeadAtomIndex(atom_index_map_[head_temp_index]);
                    }


                int tail_temp_index=tail_atoms_of_residue.at(0)->GetIndex();
                  it2 = atom_index_map_.find(tail_temp_index);
                    if(it2 != atom_index_map_.end()){
                          off_file_residue->SetTailAtomIndex(atom_index_map_[tail_temp_index]);
                    }

                    //adding an off file residue to the offfileresidue vector
                    off_file_residues_.push_back(off_file_residue);
            }


            //For populating the bounding information 

       // gmml::AtomVector center_atoms_visited = gmml::AtomVector();
        for(ResidueVector::iterator it = assembly_residues.begin(); it != assembly_residues.end(); it++)
        {
            MolecularModeling::Residue* residue = (*it);
            AtomVector all_atoms_of_residue = residue->GetAtoms();
            for(AtomVector::iterator it2 = all_atoms_of_residue.begin(); it2 != all_atoms_of_residue.end(); it2++)
            {
                MolecularModeling::Atom* atom = (*it2);

                  int main_atom_index;
                    //adding the bonding atom index to 
                    std::map<int,int>::iterator it = atom_index_map_.find(atom->GetIndex());
                    if(it != atom_index_map_.end()){
                            main_atom_index = atom_index_map_[atom->GetIndex()];
                    }


                std::string bonded_index;
                AtomVector bonded_atoms = atom->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it3 = bonded_atoms.begin(); it3 != bonded_atoms.end(); it3++)
                {
                    MolecularModeling::Atom* bonded_atom = (*it3);

                    int bonded_atom_index= bonded_atom->GetIndex();
                    std::map<int,int>::iterator it = atom_index_map_.find(bonded_atom_index);
                    if(it != atom_index_map_.end()){
                      std::size_t found = bonded_index.find_first_of(atom_index_map_[bonded_atom_index]);
                      if(found==std::string::npos)
                      {
                        bonded_index.append(" "+atom_index_map_[bonded_atom_index]);
                      }
                    }
                }
                atom_bonding_map_[main_atom_index]=bonded_index;
            }
        }
    }
    void OffFileSpace::OffFile::WriteAtomSection(std::ofstream &stream, OffFileResidueVector off_file_residues)
    {

        const std::string FLAG = "131072";
        stream << "!entry." << unit_name_ << ".unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg" << std::endl;
           
           for(OffFileResidueVector::iterator it = off_file_residues.begin(); it != off_file_residues.end(); it++)
            {
                    OffFileSpace::OffFileResidue* residue =(*it);
                    OffFileAtomVector residue_atoms=residue->GetAtoms();
                    for(OffFileAtomVector::iterator it1 = residue_atoms.begin(); it1 != residue_atoms.end(); it1++)
                    {
                        OffFileAtom* atom=(*it1);
                        stream << " \"" << atom->GetName() << "\" " << "\"" << atom->GetType() << "\" " << "0" << " " << atom->GetResidueIndex() << " " << FLAG << " "
                         << atom->GetAtomIndex() << " " << atom->GetAtomicNumber() << " " << std::fixed << atom->GetCharge() << std::endl;
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
                     stream << " \"" << current_atom->GetName() << "\" " << "\"" << current_atom->MolecularDynamicAtom::GetAtomType() << "\" " << 0
                    << " " << -1 << " " << "0.0" << std::endl;
                }
            }
    }

    void OffFileSpace::OffFile::WriteBoundBoxSection(std::ofstream& stream, MolecularModeling::Assembly* assembly)
    {
            stream << "!entry." << unit_name_ << ".unit.boundbox array dbl" << std::endl;
            stream << " " << "-1.000000" << std::endl;
            stream << " " << "0.0" << std::endl;
            stream << " " << "0.0" << std::endl;
            stream << " " << "0.0" << std::endl;
            stream << " " << "0.0" << std::endl;
    }
    void OffFileSpace::OffFile::WriteChildSequenceSection(std::ofstream& stream, ResidueVector assembly_residues)
    {

            stream << "!entry." << unit_name_ << ".unit.childsequence single int" << std::endl;
            unsigned int residue_count = assembly_residues.size();
            stream << " " << residue_count+1 << std::endl;

    }
    void OffFileSpace::OffFile::WriteConnectSection(std::ofstream& stream, ResidueVector assembly_residues)
    {
            stream << "!entry." << unit_name_ << ".unit.connect array int" << std::endl;
            stream << " " << "0" << std::endl;
            stream << " " << "0" << std::endl;
    }
    void OffFileSpace::OffFile::WriteConnectivitySection(std::ofstream& stream, ResidueVector assembly_residues)
    {
        stream << "!entry." << unit_name_ << ".unit.connectivity table  int atom1x  int atom2x  int flags" << std::endl;
	//Hi Ayush: In this section, each bond appears only once. If you already find that 1 bonds to 2, you can't later say 2 bonds to 1.
	//Below I've add some logic to your code to avoid duplication.
	gmml::AtomVector center_atoms_visited = gmml::AtomVector();
        for(ResidueVector::iterator it = assembly_residues.begin(); it != assembly_residues.end(); it++)
        {
            MolecularModeling::Residue* residue = (*it);
            AtomVector all_atoms_of_residue = residue->GetAtoms();
            for(AtomVector::iterator it2 = all_atoms_of_residue.begin(); it2 != all_atoms_of_residue.end(); it2++)
            {
                MolecularModeling::Atom* atom = (*it2);
		center_atoms_visited.push_back(atom);
                AtomVector bonded_atoms = atom->GetNode()->GetNodeNeighbors();
                for(AtomVector::iterator it3 = bonded_atoms.begin(); it3 != bonded_atoms.end(); it3++)
                {
                    MolecularModeling::Atom* bonded_atom = (*it3);
		    if (std::find(center_atoms_visited.begin(), center_atoms_visited.end(), bonded_atom) == center_atoms_visited.end()){
                        stream << " " << atom->GetIndex() << " " << bonded_atom->GetIndex() << " " << "1" << std::endl;
		    }
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
                stream << " \"" << "U" << "\"" << " " << 0 << " " << "\"" << "R" << "\"" << " " << residue_count << std::endl;
                AtomVector all_atoms_of_residue = residue->GetAtoms();
                for(AtomVector::iterator it = all_atoms_of_residue.begin(); it != all_atoms_of_residue.end(); it++)
                {

                    MolecularModeling::Atom* atom = *it;
                    stream << " \"" << "R" << "\"" << " " << residue_count << " " << "\"" << "A" << "\"" << " " << atom->GetIndex() << std::endl;
                }
            }
    }
    void OffFileSpace::OffFile::WriteNameSection(std::ofstream& stream)
    {
            stream << "!entry." << unit_name_ << ".unit.name single str" << std::endl;
            stream << " \"" << unit_name_ << "\"" << std::endl;
    }
    void OffFileSpace::OffFile::WritePositionSection(std::ofstream& stream, ResidueVector assembly_residues, int CoordinateIndex)
    {
        stream << "!entry." << unit_name_ << ".unit.positions table  dbl x  dbl y  dbl z" << std::endl;

            for(ResidueVector::iterator it = assembly_residues.begin(); it != assembly_residues.end(); it++)
            {

                MolecularModeling::Residue* residue = (*it);
                AtomVector all_atoms_of_residue = residue->GetAtoms();
                for(AtomVector::iterator it = all_atoms_of_residue.begin(); it != all_atoms_of_residue.end(); it++)
                {
                    MolecularModeling::Atom* atom = *it;
                    GeometryTopology::Coordinate coordinate = atom->GetCoordinates().at(CoordinateIndex);
                    stream << " " << std::fixed << coordinate.GetX() << " "<< std::fixed << coordinate.GetY() << " " << std::fixed << coordinate.GetZ() << std::endl;
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

                int c1x[head_atoms_of_residue.size()]={0};
                int c2x[tail_atoms_of_residue.size()]={0};

                int c1x_count=0;
                int c2x_count=0;
                for(AtomVector::iterator it1 = head_atoms_of_residue.begin(); it1 != head_atoms_of_residue.end(); it1++)
                {
                    MolecularModeling::Atom* atom = *it1;
                    c1x[c1x_count]=atom->GetIndex();
                    c1x_count++;
                }

                 for(AtomVector::iterator it2 = tail_atoms_of_residue.begin(); it2!= tail_atoms_of_residue.end(); it2++)
                {
                   MolecularModeling::Atom* atom = *it2;
                    c2x[c2x_count]=atom->GetIndex();
                    c2x_count++;
                }

                for(int i=0;i<c1x_count || i<c2x_count;i++)
                {
                    if(i<c1x_count && i<c2x_count)
                    {
                        stream << " " << c1x[i] << " " << c2x[i] << " " << 0 << " " << 0 << " " << 0 << " " << 0 << std::endl;
                    }
                    else if(i<c1x_count && i>c2x_count )
                    {
                        stream << " " << c1x[i] << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << std::endl;
                    }
                    else if (i>c1x_count && i<c2x_count)
                    {
                        stream << " " << 0 << " " << c2x[i]  << " " << 0 << " " << 0 << " " << 0 << " " << 0 << std::endl;
                    }
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

            stream << "!entry." << unit_name_ << ".unit.residues table  str name  int seq  int childseq  int startatomx  str restype  int imagingx" << std::endl;
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

                    stream << " \"" << name << "\"" << " " << seq << " " << childseq << " " << startatomx << " " << restype << " " << imagingx << std::endl;
                }
    }

    void OffFileSpace::OffFile::WriteSolventCapSection(std::ofstream& stream)
    {
        stream << "!entry." << unit_name_ << ".unit.solventcap array dbl" << std::endl;
        stream << " " << "-1.000000" << std::endl;
        stream << " " <<"0.0" << std::endl;
        stream << " " <<"0.0" << std::endl;
        stream << " " <<"0.0" << std::endl;
        stream << " " <<"0.0" << std::endl;
    }

    void OffFileSpace::OffFile::WriteVelocitiesSection(std::ofstream& stream, ResidueVector assembly_residues)
    {
        stream << "!entry." << unit_name_ << ".unit.velocities table  dbl x  dbl y  dbl z" << std::endl;
        for(ResidueVector::iterator it = assembly_residues.begin(); it != assembly_residues.end(); it++)
            {
                MolecularModeling::Residue* residue = (*it);
                AtomVector all_atoms_of_residue = residue->GetAtoms();
                for(AtomVector::iterator it = all_atoms_of_residue.begin(); it != all_atoms_of_residue.end(); it++)
                {
                     stream << " " << "0.0" << " " << "0.0" << " " << "0.0" << std::endl;
                }
            }
    }

    //////////////////////////////////////////////////////////
    //                     DISPLAY FUNCTIONS                //
    //////////////////////////////////////////////////////////
