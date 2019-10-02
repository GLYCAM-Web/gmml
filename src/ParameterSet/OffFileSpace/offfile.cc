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
            MolecularModeling::ResidueVector residues = assembly->GetResidues();
            unit_name_=assembly->GetName();
            out_file << "!!index array str" << std::endl;
            out_file << " \"" << unit_name_ << "\"" << std::endl;
            PopulateOffFileResiduesFromAssembly(residues,CoordinateIndex);
            WriteAtomSection(out_file,this->off_file_residues_);
            WriteAtomPertInfoSection(out_file,this->off_file_residues_);
            WriteBoundBoxSection(out_file);
            WriteChildSequenceSection(out_file,this->off_file_residues_);
            WriteConnectSection(out_file);
            WriteConnectivitySection(out_file,residues);
            WriteHierarchySection(out_file,this->off_file_residues_);
            WriteNameSection(out_file);
            WritePositionSection(out_file,residues, CoordinateIndex);
            WriteResidueConnectSection(out_file,residues);
            WriteResiduesSection(out_file,residues);
            WriteSolventCapSection(out_file);
            WriteVelocitiesSection(out_file,this->off_file_residues_);
              out_file.close();

    }

    void OffFileSpace::OffFile::PopulateOffFileResiduesFromAssembly(MolecularModeling::ResidueVector assembly_residues,int CoordinateIndex)
    {
            int ResidueIndex=0;
             int BoundingAtomIndex=0;
            for(MolecularModeling::ResidueVector::iterator it = assembly_residues.begin(); it != assembly_residues.end(); it++)
            {
                ResidueIndex++;
                OffFileResidue* off_file_residue = new OffFileResidue();

                MolecularModeling::Residue* assembly_residue = (*it);
                std::string name=assembly_residue->GetName();
                off_file_residue->SetName(name);
                off_file_residue->SetListingIndex(ResidueIndex);
                
                  MolecularModeling::AtomVector all_atoms_of_residue = assembly_residue->GetAtoms();

                int AtomIndex=0;

                for(MolecularModeling::AtomVector::iterator it = all_atoms_of_residue.begin(); it != all_atoms_of_residue.end(); it++)
                {
                    AtomIndex++;
                    BoundingAtomIndex++;
                    MolecularModeling::Atom* residue_atom = *it;
                    OffFileAtom* off_file_atom= new OffFileAtom();

                    //setting values from residue atom to off file atom
                    off_file_atom->SetName(residue_atom->GetName());
                    off_file_atom->SetType(residue_atom->MolecularDynamicAtom::GetAtomType());
                    off_file_atom->SetResidueIndex(ResidueIndex);
                    off_file_atom->SetAtomIndex(AtomIndex);

                    int AtomicNumber=0;
                    std::string atom_element_symbol= residue_atom->GetElementSymbol();
                   int size_of_lookup_map = sizeof(gmml::ElementAtributes::Elements) / sizeof(gmml::ElementAtributes::Elements[0]);
                       for (int i = 0; i < size_of_lookup_map; i++){
                                  gmml::ElementAtributes::ElementAttributeInfo entry = gmml::ElementAtributes::Elements[i];
                                  if (atom_element_symbol.compare(entry.elment_type_) == 0){
                                      AtomicNumber = entry.atomic_number_;
                                  }
                              }


                    off_file_atom->SetAtomicNumber(AtomicNumber);
                    off_file_atom->SetCoordinate(residue_atom->GetCoordinates().at(CoordinateIndex));
                    off_file_atom->SetAtomCharge(residue_atom->GetCharge());
                    this->atom_index_map_[residue_atom->GetIndex()]=AtomIndex;
                     this->atom_bonding_map_[residue_atom->GetIndex()]=BoundingAtomIndex;
                    off_file_residue->AddAtom(off_file_atom);

                }

                /*
                //for adding head and tail atoms of residues
                MolecularModeling::AtomVector head_atoms_of_residue= assembly_residue->GetHeadAtoms();
                MolecularModeling::AtomVector tail_atoms_of_residue= assembly_residue->GetTailAtoms();

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
                */
                    //adding an off file residue to the offfileresidue vector
                    off_file_residues_.push_back(off_file_residue);
            }


            //For populating the bounding information 
        //for(MolecularModeling::ResidueVector::iterator it = assembly_residues.begin(); it != assembly_residues.end(); it++)
        //{
            //MolecularModeling::Residue* residue = (*it);
            //MolecularModeling::AtomVector all_atoms_of_residue = residue->GetAtoms();
/*
This whole loop seems to be unused
            for(MolecularModeling::AtomVector::iterator it2 = all_atoms_of_residue.begin(); it2 != all_atoms_of_residue.end(); ++it2)
            {
                MolecularModeling::Atom* atom = (*it2);
                  // int main_atom_index=0; // commented because not used
                    //adding the bonding atom index to 
                    std::map<int,int>::iterator it = atom_index_map_.find(atom->GetIndex());
                    if(it != atom_index_map_.end()){
                        main_atom_index = atom_index_map_[atom->GetIndex()];
                        } 
            }
*/
         // }
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

    void OffFileSpace::OffFile::WriteAtomPertInfoSection(std::ofstream& stream, OffFileResidueVector off_file_residues)
    {

        stream << "!entry." << unit_name_ << ".unit.atomspertinfo table  str pname  str ptype  int ptypex  int pelmnt  dbl pchg" << std::endl;
        for(OffFileResidueVector::iterator it = off_file_residues.begin(); it != off_file_residues.end(); it++)
         {
            OffFileSpace::OffFileResidue* residue =(*it);
            OffFileAtomVector residue_atoms=residue->GetAtoms();
            for(OffFileAtomVector::iterator it1 = residue_atoms.begin(); it1 != residue_atoms.end(); it1++)
            {
                OffFileAtom* atom=(*it1);
                stream << " \"" << atom->GetName() << "\" " << "\"" << atom->GetType() << "\" " << 0
                << " " << -1 << " " << 0.0 << std::endl;
            }
        }
    }

    void OffFileSpace::OffFile::WriteBoundBoxSection(std::ofstream& stream)
    {
            stream << "!entry." << unit_name_ << ".unit.boundbox array dbl" << std::endl;
            stream << " " << "-1.000000" << std::endl;
            stream << " " << 0.0 << std::endl;
            stream << " " << 0.0 << std::endl;
            stream << " " << 0.0 << std::endl;
            stream << " " << 0.0 << std::endl;
    }

    void OffFileSpace::OffFile::WriteChildSequenceSection(std::ofstream& stream, OffFileResidueVector off_file_residues)
    {

            stream << "!entry." << unit_name_ << ".unit.childsequence single int" << std::endl;
            unsigned int residue_count = off_file_residues.size();
            stream << " " << residue_count+1 << std::endl;
    }

    void OffFileSpace::OffFile::WriteConnectSection(std::ofstream& stream)
    {
            stream << "!entry." << unit_name_ << ".unit.connect array int" << std::endl;
            stream << " " << 0 << std::endl;
            stream << " " << 0 << std::endl;
    }
    void OffFileSpace::OffFile::WriteConnectivitySection(std::ofstream& stream, MolecularModeling::ResidueVector assembly_residues)
    {
        stream << "!entry." << unit_name_ << ".unit.connectivity table  int atom1x  int atom2x  int flags" << std::endl;

        MolecularModeling::AtomVector center_atoms_visited = MolecularModeling::AtomVector();

                for(MolecularModeling::ResidueVector::iterator it = assembly_residues.begin(); it != assembly_residues.end(); it++)
                {
                    MolecularModeling::Residue* residue = (*it);
                    MolecularModeling::AtomVector all_atoms_of_residue = residue->GetAtoms();
                    for(MolecularModeling::AtomVector::iterator it2 = all_atoms_of_residue.begin(); it2 != all_atoms_of_residue.end(); it2++)
                    {
                        MolecularModeling::Atom* atom = (*it2);

                        center_atoms_visited.push_back(atom);
                        MolecularModeling::AtomVector bonded_atoms = atom->GetNode()->GetNodeNeighbors();
                        for(MolecularModeling::AtomVector::iterator it3 = bonded_atoms.begin(); it3 != bonded_atoms.end(); it3++)
                        {
                            MolecularModeling::Atom* bonded_atom = (*it3);
                                if (std::find(center_atoms_visited.begin(), center_atoms_visited.end(), bonded_atom) == center_atoms_visited.end())
                                {
                                    stream << " " << atom_bonding_map_[atom->GetIndex()] << " " << atom_bonding_map_[bonded_atom->GetIndex()] << " " << 1 << std::endl;
                                }
                        }
                    }
                }
    }

    void OffFileSpace::OffFile::WriteHierarchySection(std::ofstream& stream, OffFileResidueVector off_file_residues)
    {
        stream << "!entry." << unit_name_ << ".unit.hierarchy table  str abovetype  int abovex  str belowtype  int belowx" << std::endl;

        int atom_counter=0;
        for(OffFileResidueVector::iterator it = off_file_residues.begin(); it != off_file_residues.end(); it++)
            {   

                OffFileSpace::OffFileResidue* residue = (*it);
                stream << " \"" << "U" << "\"" << " " << 0 << " " << "\"" << "R" << "\"" << " " << residue->GetListingIndex() << std::endl;
                OffFileAtomVector all_atoms_of_residue = residue->GetAtoms();
                for(OffFileAtomVector::iterator it = all_atoms_of_residue.begin(); it != all_atoms_of_residue.end(); it++)
                {
                    atom_counter++;
//  unused                    OffFileSpace::OffFileAtom* atom = *it;
                    stream << " \"" << "R" << "\"" << " " << residue->GetListingIndex() << " " << "\"" << "A" << "\"" << " " << atom_counter << std::endl;
                }
            }
    }
    void OffFileSpace::OffFile::WriteNameSection(std::ofstream& stream)
    {
            stream << "!entry." << unit_name_ << ".unit.name single str" << std::endl;
            stream << " \"" << unit_name_ << "\"" << std::endl;
    }


    void OffFileSpace::OffFile::WritePositionSection(std::ofstream& stream, MolecularModeling::ResidueVector assembly_residues, int CoordinateIndex)
    {
        stream << "!entry." << unit_name_ << ".unit.positions table  dbl x  dbl y  dbl z" << std::endl;

            for(MolecularModeling::ResidueVector::iterator it = assembly_residues.begin(); it != assembly_residues.end(); it++)
            {

                MolecularModeling::Residue* residue = (*it);
                MolecularModeling::AtomVector all_atoms_of_residue = residue->GetAtoms();
                for(MolecularModeling::AtomVector::iterator it = all_atoms_of_residue.begin(); it != all_atoms_of_residue.end(); it++)
                {
                    MolecularModeling::Atom* atom = *it;
                    GeometryTopology::Coordinate coordinate = atom->GetCoordinates().at(CoordinateIndex);
                    stream << " " << std::fixed << coordinate.GetX() << " "<< std::fixed << coordinate.GetY() << " " << std::fixed << coordinate.GetZ() << std::endl;
                }
            }

    }


    void OffFileSpace::OffFile::WriteResidueConnectSection(std::ofstream& stream, MolecularModeling::ResidueVector assembly_residues)
    {    
        stream << "!entry." << unit_name_ << ".unit.residueconnect table  int c1x  int c2x  int c3x  int c4x  int c5x  int c6x" << std::endl;
        for(MolecularModeling::ResidueVector::iterator it = assembly_residues.begin(); it != assembly_residues.end(); it++)
            {

                MolecularModeling::Residue* residue = (*it);

                MolecularModeling::AtomVector head_atoms_of_residue= residue->GetHeadAtoms();
                MolecularModeling::AtomVector tail_atoms_of_residue= residue->GetTailAtoms();

                //int c1x[head_atoms_of_residue.size()]={0};
                //int c2x[tail_atoms_of_residue.size()]={0};
		std::vector<int> c1x(head_atoms_of_residue.size(), 0);
		std::vector<int> c2x(tail_atoms_of_residue.size(), 0);

                int c1x_count=0;
                int c2x_count=0;
                for(MolecularModeling::AtomVector::iterator it1 = head_atoms_of_residue.begin(); it1 != head_atoms_of_residue.end(); it1++)
                {
                    MolecularModeling::Atom* atom = *it1;
                    c1x[c1x_count]=atom_bonding_map_[atom->GetIndex()];
		    //std::cout << "c1x[" << c1x_count << "] = " << atom_bonding_map_[atom->GetIndex()] << std::endl;
                    c1x_count++;
                }

                 for(MolecularModeling::AtomVector::iterator it2 = tail_atoms_of_residue.begin(); it2!= tail_atoms_of_residue.end(); it2++)
                {
                   MolecularModeling::Atom* atom = *it2;
                    c2x[c2x_count]=atom_bonding_map_[atom->GetIndex()];
		    //std::cout << "c2x[" << c2x_count << "] = " << atom_bonding_map_[atom->GetIndex()] << std::endl;
                    c2x_count++;
                }

		stream << " " << c1x[0] << " " << c2x[0];
		for (int i=1; i<c1x_count; i++){
		    stream << " " << c1x[i];
		}
		for (int i=1; i<c2x_count; i++){
		    stream << " " << c2x[i];
		}

		int zero_column_count = 6 - head_atoms_of_residue.size() - tail_atoms_of_residue.size();
		for (int i=0; i<zero_column_count; i++){
		    stream << " " << "0";
		}
		stream << std::endl;
                /*for(int i=0;i<c1x_count || i<c2x_count;i++)
                {
                    if(i<c1x_count && i<c2x_count)
                    {
                        stream << " " << c1x[i] << " " << c2x[i] << " " << 0 << " " << 0 << " " << 0 << " " << 0 << std::endl;
                    }
                    else if(i<c1x_count && i>c2x_count )//if more head than tail
                    {
                        stream << " " << c1x[i] << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << std::endl;
                    }
                    else if (i>c1x_count && i<c2x_count)
                    {
                        stream << " " << 0 << " " << c2x[i]  << " " << 0 << " " << 0 << " " << 0 << " " << 0 << std::endl;
                    }
                }*/
		
            }
    }
    

    void OffFileSpace::OffFile::WriteResiduesSection(std::ofstream& stream, MolecularModeling::ResidueVector assembly_residues)
    {
            std::string name;
            unsigned int seq=0;
            unsigned int childseq;
            unsigned int startatomx;
            std::string restype;
            unsigned int imagingx;

            stream << "!entry." << unit_name_ << ".unit.residues table  str name  int seq  int childseq  int startatomx  str restype  int imagingx" << std::endl;
            for(MolecularModeling::ResidueVector::iterator it = assembly_residues.begin(); it != assembly_residues.end(); it++)
                {
                    MolecularModeling::Residue* residue = (*it);
                    name= residue->GetName();
                    seq++;
                    childseq= residue->GetAtoms().size()+1;
                    //MolecularModeling::AtomVector head_atoms_of_residue = residue->GetHeadAtoms();
                    //MolecularModeling::Atom* current_atom =head_atoms_of_residue[0];
                    MolecularModeling::Atom* current_atom = residue->GetAtoms().at(0);
		    //Yao: Don't determine starting atom index from head atoms. Instead, search for AtomIndex = 1 in an offfile atom object, or use the index of an element in this ResidueVector.
		    //Since there are offfile residue/atom objects populated, it might be better to rely on these objects, rather than assembly residues.
		    //In one case,head_atoms_of_residue[0] is indeed the 2nd atom in residue.I'm talking about the residue ROH as extracted from prep file.
                    startatomx= atom_bonding_map_[current_atom->GetIndex()];
                    restype="?";
                    imagingx=0;

                    stream << " \"" << name << "\"" << " " << seq << " " << childseq << " " << startatomx << " " << "\""<<restype <<"\""<< " " << imagingx << std::endl;
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

    void OffFileSpace::OffFile::WriteVelocitiesSection(std::ofstream& stream, OffFileResidueVector off_file_residues)
    {
        stream << "!entry." << unit_name_ << ".unit.velocities table  dbl x  dbl y  dbl z" << std::endl;
        for(OffFileResidueVector::iterator it = off_file_residues.begin(); it != off_file_residues.end(); it++)
            {
                OffFileSpace::OffFileResidue * residue = (*it);
                OffFileAtomVector all_atoms_of_residue = residue->GetAtoms();
                for(OffFileAtomVector::iterator it = all_atoms_of_residue.begin(); it != all_atoms_of_residue.end(); it++)
                {
                     stream << " " << "0.0" << " " << "0.0" << " " << "0.0" << std::endl;
                }
            }
    }

    //////////////////////////////////////////////////////////
    //                     DISPLAY FUNCTIONS                //
    //////////////////////////////////////////////////////////


    void OffFile::Print(std::ostream& out)
    {
        out << "Printing OffFile Details"<<std::endl;
    }
