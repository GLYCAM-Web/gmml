#include <math.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <queue>
#include <stack>

#include "../../../../includes/MolecularModeling/assembly.hpp"
#include "../../../../includes/MolecularModeling/residue.hpp"
#include "../../../../includes/MolecularModeling/atom.hpp"
#include "../../../../includes/MolecularModeling/atomnode.hpp"
#include "../../../../includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp"
#include "../../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp"
#include "../../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyfile.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyassembly.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyresidue.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyatom.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologybond.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologybondtype.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyangle.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyangletype.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologydihedral.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologydihedraltype.hpp"
#include "../../../../includes/InputSet/TopologyFileSpace/topologyatompair.hpp"
#include "../../../../includes/InputSet/CoordinateFileSpace/coordinatefile.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbfile.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbremarksection.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbtitlesection.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbmodelcard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbmodelsection.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbatomsection.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbconnectsection.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdblinkcard.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdblinksection.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdblinkcardresidue.hpp"
#include "../../../../includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtfile.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtmodel.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelcard.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../../../includes/InputSet/PdbqtFileSpace/pdbqtmodelresidueset.hpp"
#include "../../../../includes/ParameterSet/LibraryFileSpace/libraryfile.hpp"
#include "../../../../includes/ParameterSet/LibraryFileSpace/libraryfileatom.hpp"
#include "../../../../includes/ParameterSet/LibraryFileSpace/libraryfileresidue.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfile.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfilebond.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfileangle.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedral.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfiledihedralterm.hpp"
#include "../../../../includes/ParameterSet/ParameterFileSpace/parameterfileatom.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfile.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfileresidue.hpp"
#include "../../../../includes/ParameterSet/PrepFileSpace/prepfileatom.hpp"
#include "../../../../includes/utils.hpp"
#include "../../../../includes/common.hpp"
#include "../../../../includes/GeometryTopology/grid.hpp"
#include "../../../../includes/GeometryTopology/cell.hpp"

#include "../../../../includes/InputSet/PdbFileSpace/inputfile.hpp"

#include <unistd.h>
#include <errno.h>
#include <string.h>

using MolecularModeling::Assembly;

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
void Assembly::PopulateOntology(std::ofstream& main_stream, OligosaccharideVector oligos)
{
    std::stringstream pdb_stream;
    //Move to input file
    std::string pdb_resource = CreateURIResource(gmml::OntPDB, 0, "", "");
    //    CreateTitle(pdb_resource, pdb_stream);
    std::transform(pdb_resource.begin(), pdb_resource.end(),pdb_resource.begin(), ::tolower);
    std::stringstream ss;
    ss << pdb_resource << "_";
    std::string id_prefix = ss.str();
    std::string pdb_uri = CreateURI(pdb_resource);
    //
    //
    //
    // //    pdb_stream << Ontology::ENTITY_COMMENT << pdb_resource << std::endl;
    // gmml::AddTriple(pdb_uri, Ontology::TYPE, Ontology::PDB, pdb_stream);
    // gmml::AddLiteral(pdb_uri, Ontology::id, pdb_resource, pdb_stream);
    // //    gmml::AddLiteral(pdb_uri, Ontology::LABEL, pdb_resource, pdb_stream);
    // //    gmml::AddLiteral(pdb_uri, Ontology::input_file_path, source_file_, pdb_stream);
    //
    // //Add Inputfile Ontology function
    // // std::stringstream out_stream;
    // // this->input_file_->PrintOntology(out_stream);
    // // std::cout << out_stream.str();
    // gmml::log(__LINE__, __FILE__,  gmml::INF, "About to get pdb info");
    this->input_file_->PrintOntology(pdb_stream);
    // gmml::log(__LINE__, __FILE__,  gmml::INF, "Done getting pdb info");

    int link_id = 1;
    std::stringstream oligo_stream;
    std::stringstream oligo_sequence_stream;
    std::stringstream mono_stream;
    std::stringstream linkage_stream;
    std::vector<std::string> side_or_ring_atoms = std::vector<std::string>();
    std::vector<int> visited_oligos = std::vector<int>();

    // gmml::log(__LINE__, __FILE__,  gmml::INF, "Getting notes");
    NoteVector notes = this->GetNotes();
    std::stringstream note_stream;
    if(notes.size() != 0)
    {
        int note_id = 1;
        // gmml::log(__LINE__, __FILE__,  gmml::INF, "Poulating notes");
        PopulateNotes(pdb_stream, note_stream, pdb_uri, notes, id_prefix, note_id);
    }

    std::map<std::string, std::string> mono_to_short_name_map;
    std::map<std::string, std::string> oligo_to_res_uri_map;
    int root_oligo_id = 0;
    // gmml::log(__LINE__, __FILE__,  gmml::INF, "Populating oligosaccharides");
    //TODO segfaults populating oligos
    PopulateOligosaccharide(pdb_stream, oligo_stream, oligo_sequence_stream, mono_stream, linkage_stream, pdb_uri, id_prefix, link_id, oligos, side_or_ring_atoms, visited_oligos, mono_to_short_name_map, oligo_to_res_uri_map, root_oligo_id);

    // gmml::log(__LINE__, __FILE__,  gmml::INF, "Getting residues");
    ResidueVector residues = this->GetResidues();

		// Need to remove all Residues that are Waters or Proteins(except the Protein at the end of a Glycan::Oligosaccharide)
		for( ResidueVector::iterator it = residues.begin(); it != residues.end(); ) {
			Residue* residue = ( *it );
			if( residue->CheckIfProtein() || residue->CheckIfWater() ) {
				delete residue;
				it = residues.erase( it );
			} else {
				it++;
			}
		}

    std::stringstream residue_stream;
    // if(!residues.empty())
    // {
    //   gmml::log(__LINE__, __FILE__,  gmml::INF, "Populating residues");
      // PopulateResidue(pdb_stream, residue_stream, pdb_uri, id_prefix, residues, side_or_ring_atoms);
    // }


    // gmml::log(__LINE__, __FILE__,  gmml::INF, "PDB Stream");
    // gmml::log(__LINE__, __FILE__,  gmml::INF, pdb_stream.str());
    main_stream << pdb_stream.str() << note_stream.str() << oligo_stream.str() << oligo_sequence_stream.str() << mono_stream.str() << linkage_stream.str() << residue_stream.str() << std::endl;
}

void Assembly::PopulateNotes(std::stringstream& pdb_stream, std::stringstream& note_stream, std::string pdb_uri, NoteVector notes, std::string id_prefix, int note_id)
{
    std::string note_resource = "";
    std::string note_uri = "";
    for(NoteVector::iterator it = notes.begin(); it != notes.end(); it++)
    {
        Glycan::Note* note = (*it);
        note_resource = CreateURIResource(gmml::OntNote, note_id, id_prefix, "");
        note_uri = CreateURI(note_resource);
        gmml::AddTriple(pdb_uri, Ontology::hasNote, note_uri, pdb_stream);

        //        note_stream << Ontology::ENTITY_COMMENT << note_resource << std::endl;
        gmml::AddTriple(note_uri, Ontology::TYPE, Ontology::Note, note_stream);
        //        gmml::AddLiteral(note_uri, Ontology::LABEL, note_resource, note_stream);
        gmml::AddLiteral(note_uri, Ontology::note_type, note->ConvertGlycanNoteType2String(note->type_), note_stream);
        gmml::AddLiteral(note_uri, Ontology::note_category, note->ConvertGlycanNoteCat2String(note->category_), note_stream);
        gmml::AddLiteral(note_uri, Ontology::note_description, note->description_, note_stream);
        note_id++;
    }
}

void Assembly::PopulateOligosaccharide(std::stringstream& pdb_stream, std::stringstream& oligo_stream, std::stringstream& oligo_sequence_stream, std::stringstream& mono_stream, std::stringstream& linkage_stream, std::string pdb_uri, std::string id_prefix,
                                       int& link_id, OligosaccharideVector oligos, std::vector<std::string>& side_or_ring_atoms, std::vector<int>& visited_oligos,
                                       std::map<std::string, std::string>& mono_to_short_name_map, std::map<std::string, std::string>& oligo_to_res_uri_map, int& root_oligo_id)
{
    std::string oligo_resource = "";
    std::string oligo_uri = "";
    std::string child_oligo_resource = "";
    std::string child_oligo_uri = "";
    std::string child_res_resource = "";
    std::string child_res_uri = "";
    std::string child_mono_resource = "";
    std::string child_mono_uri = "";
    // std::string child_res_resource = "";
    // std::string child_res_uri = "";
    std::string parent_mono_resource = "";
    std::string parent_mono_uri = "";
    std::string parent_res_resource = "";
    std::string parent_res_uri = "";
    std::string root_oligo_resource = "";
    std::string root_oligo_uri = "";
    std::string term_resource = "";
    std::string term_uri = "";
    if(oligos.size() != 0) //Earlier it was  if(oligos.size() != NULL) but throws  warning: NULL used in arithmetic [-Wpointer-arith], hence changed NULL to 0 by Ayush on 06/22/2017
    {
      int oligoNum = 1;
        for(OligosaccharideVector::iterator it = oligos.begin(); it != oligos.end(); it++, oligoNum++)
        {
          std::vector<MolecularModeling::Residue*> residueVector;
          residueVector.clear();
          int linkNum = 0;
            Glycan::Oligosaccharide* oligo = (*it);

            oligo_resource = CreateURIResource(gmml::OntOligosaccharide, oligoNum, id_prefix, "");
            oligo_uri = CreateURI(oligo_resource);

            gmml::AddTriple(pdb_uri, Ontology::hasOligo, oligo_uri, pdb_stream);
            //            oligo_stream << Ontology::ENTITY_COMMENT << oligo_resource << std::endl;
            gmml::AddTriple(oligo_uri, Ontology::TYPE, Ontology::Oligosaccharide, oligo_stream);
            //            gmml::AddLiteral(oligo_uri, Ontology::LABEL, oligo_resource, oligo_stream);

            std::string o_name = oligo->oligosaccharide_name_;
            if(o_name.compare("") != 0)
            {
                root_oligo_id = oligoNum;
                gmml::AddLiteral(oligo_uri, Ontology::oligo_name, o_name, oligo_stream);
                gmml::AddLiteral(oligo_uri, Ontology::oligo_sequence_name, o_name, oligo_sequence_stream);

                mono_to_short_name_map.clear();
                oligo_to_res_uri_map.clear();
            }
            std::string oligo_iupac = oligo->IUPAC_name_;
            if(oligo_iupac.compare("") != 0)
            {
              gmml::AddLiteral(oligo_uri, Ontology::oligo_IUPAC_name, oligo_iupac, oligo_sequence_stream);
            }
            std::string author_oligo_iupac = oligo->author_IUPAC_name_;
            if(oligo_iupac.compare("") != 0)
            {
              gmml::AddLiteral(oligo_uri, Ontology::author_oligo_name, author_oligo_iupac, oligo_sequence_stream);
            }
            if(oligo_iupac.find("<R"))
            {
              int numR = 0;
              for(std::vector<Glycan::Monosaccharide*>::reverse_iterator rit = oligo->mono_nodes_.rbegin(); rit != oligo->mono_nodes_.rend(); rit++)
              {
                Glycan::Monosaccharide* thisMono = *rit;
                if(thisMono->on_R_ > 0)
                {
                  for(std::vector<std::pair<std::string, std::string> >::iterator derivative = thisMono->unknown_derivatives_.begin(); derivative != thisMono->unknown_derivatives_.end(); derivative++)
                  {
                    numR++;
                    if((*derivative).second != "" && thisMono->sugar_name_.monosaccharide_short_name_ != thisMono->cycle_atoms_[0]->GetResidue()->GetName())
                    {
                      std::stringstream RgroupStream, RnumStream;
                      RgroupStream << oligo_uri << "_R" << numR;
                      RnumStream << "gmmo:has" << "R" << numR;
                      gmml::AddTriple(oligo_uri, RnumStream.str(), RgroupStream.str(), oligo_sequence_stream); 
                      gmml::AddLiteral(RgroupStream.str(), Ontology::hasFormula,(*derivative).second, oligo_sequence_stream); 
                    }
                  }
                }
              }
            }
            // gmml::log(__LINE__, __FILE__,  gmml::INF, " ");
            float o_b_factor = oligo->oligosaccharide_b_factor_;
            // std::stringstream bfss;
            // bfss << o_b_factor;
            // std::string o_b_factor_str = bfss.str();
            // if(o_b_factor !=NULL)
            gmml::AddDecimal(oligo_uri, Ontology::oligo_b_factor, o_b_factor, oligo_stream);
            // gmml::log(__LINE__, __FILE__,  gmml::INF, " ");

            NoteVector notes = oligo->oligo_notes_;
            if(notes.size() != 0)
            {
                int note_id = 1;
                // gmml::log(__LINE__, __FILE__,  gmml::INF, "Populating notes");
                std::string id_prefix = oligo_uri + "_";
                PopulateNotes(oligo_stream, oligo_stream, oligo_uri, notes, id_prefix, note_id);
            }

            std::string o_residue_links = oligo->oligosaccharide_residue_linkages_;
            // gmml::log(__LINE__, __FILE__,  gmml::INF, " ");
            if(o_residue_links.compare("") != 0)
                gmml::AddLiteral(oligo_uri, Ontology::oligo_residue_linkages, o_residue_links, oligo_stream);
              // gmml::log(__LINE__, __FILE__,  gmml::INF, " ");
            int MonoNum, MonoNeighborNum;
            root_oligo_id = oligoNum;
            for(std::vector<Glycan::Monosaccharide*>::reverse_iterator rit = oligo->mono_nodes_.rbegin(); rit != oligo->mono_nodes_.rend(); rit++)
            {
              Glycan::Monosaccharide* thisMono = *rit;
              residueVector.push_back(thisMono->cycle_atoms_[0]->GetResidue());
              if(thisMono->is_root_)
              {
                MonoNum = thisMono->IUPAC_index_;//change to IUPAC_index_
                // root_oligo_id = thisMono->mono_id_;
                // gmml::log(__LINE__, __FILE__,  gmml::INF, " ");
                // PopulateLinkage(linkage_stream, oligo, oligo_uri, id_prefix, link_id, visited_oligos);
                // gmml::log(__LINE__, __FILE__,  gmml::INF, "About to populate sequence linkages");
                parent_mono_resource = CreateURIResource(gmml::OntOligosaccharide, MonoNum, id_prefix, "");
                parent_mono_uri = CreateURI(parent_mono_resource);
                std::string resID = std::to_string(MonoNum);
                std::string monoSNFG = thisMono->SNFG_name_;
                std::string monoShortName = thisMono->sugar_name_.monosaccharide_short_name_;
                parent_res_resource = CreateURIResource(gmml::OntSequenceResidue, root_oligo_id, id_prefix, resID);
                parent_res_uri = CreateURI(parent_res_resource);
                CheckDerivativesAndPopulate(oligo_sequence_stream, monoShortName, oligo_uri, parent_res_uri, monoSNFG, thisMono);
                gmml::AddLiteral(parent_res_uri, Ontology::hasNameIndex, std::to_string(thisMono->oligosaccharide_index_), oligo_sequence_stream);
                for(std::vector<std::pair<Glycan::GlycosidicLinkage*, Glycan::Monosaccharide*> >::iterator it = thisMono->mono_neighbors_.begin(); it!=thisMono->mono_neighbors_.end(); it++)
                {
                  Glycan::Monosaccharide* thisMonoNeighbor = (*it).second;
                  MonoNeighborNum = thisMonoNeighbor->IUPAC_index_;//change to IUPAC_index_
                  std::string neighborResID = std::to_string(MonoNeighborNum);
                  std::string monoSNFG = thisMonoNeighbor->SNFG_name_;
                  std::string monoShortName = thisMonoNeighbor->sugar_name_.monosaccharide_short_name_;
                  child_res_resource = CreateURIResource(gmml::OntSequenceResidue, root_oligo_id, id_prefix, neighborResID);
                  child_res_uri = CreateURI(child_res_resource);
                  //TODO replace isConnectedTo w/ isx-nLinkedTo (IE is1-3LinkedTo)
                  std::stringstream connectionInfo;
                  connectionInfo << "gmmo:is" << (*it).first->linkage_type_ << "ConnectedTo";
                  gmml::AddTriple(parent_res_uri, connectionInfo.str(), child_res_uri, oligo_sequence_stream);
                  PopulateLinkage(linkage_stream, oligo_uri, parent_res_uri, child_res_uri, linkNum, (*it).first, thisMono, thisMonoNeighbor);
                  
                }
                // PopulateSequenceLinkage(oligo_sequence_stream, oligo, oligo_uri, id_prefix, visited_oligos, mono_to_short_name_map, oligo_to_res_uri_map, root_oligo_id);
                // gmml::log(__LINE__, __FILE__,  gmml::INF, "Done populating sequence linkages");
                PopulateMonosaccharide(mono_stream, oligo_stream, oligo_uri, id_prefix, thisMono, side_or_ring_atoms, pdb_uri);
                std::string term_resource = "";
                std::string term_uri = "";
                term_resource = CreateURIResource(gmml::OntTerminal, root_oligo_id, id_prefix, "");
                term_uri = CreateURI(term_resource);
                //this is the terminal so it's okay to not have the linkage type, as it is in the terminal name (IE 1-OH)
                gmml::AddTriple(parent_res_uri, Ontology::isConnectedTo, term_uri, oligo_sequence_stream);
                gmml::AddTriple(oligo_uri, Ontology::hasTerminal, term_uri, oligo_sequence_stream);
                gmml::AddTriple(term_uri, Ontology::TYPE, Ontology::Terminal, oligo_sequence_stream);
                gmml::AddLiteral(term_uri, Ontology::id, oligo->oligosaccharide_terminal_, oligo_sequence_stream);
              }
              else
              {
                // root_oligo_id = thisMono->mono_id_;
              // gmml::log(__LINE__, __FILE__,  gmml::INF, " ");
                // PopulateLinkage(linkage_stream, oligo, oligo_uri, id_prefix, link_id, visited_oligos);
                // gmml::log(__LINE__, __FILE__,  gmml::INF, "About to populate sequence linkages");
                MonoNum = thisMono->IUPAC_index_;//change to IUPAC_index_
                parent_mono_resource = CreateURIResource(gmml::OntOligosaccharide, MonoNum, id_prefix, "");
                parent_mono_uri = CreateURI(parent_mono_resource);
                std::string resID = std::to_string(MonoNum);
                std::string monoSNFG = thisMono->SNFG_name_;
                std::string monoShortName = thisMono->sugar_name_.monosaccharide_short_name_;
                parent_res_resource = CreateURIResource(gmml::OntSequenceResidue, root_oligo_id, id_prefix, resID);
                parent_res_uri = CreateURI(parent_res_resource);
                CheckDerivativesAndPopulate(oligo_sequence_stream, monoShortName, oligo_uri, parent_res_uri, monoSNFG, thisMono);
                gmml::AddLiteral(parent_res_uri, Ontology::hasNameIndex, std::to_string(thisMono->oligosaccharide_index_), oligo_sequence_stream);
                for(std::vector<std::pair<Glycan::GlycosidicLinkage*, Glycan::Monosaccharide*> >::iterator it = thisMono->mono_neighbors_.begin(); it!=thisMono->mono_neighbors_.end(); it++)
                {
                  Glycan::Monosaccharide* thisMonoNeighbor = (*it).second;
                  MonoNeighborNum = thisMonoNeighbor->IUPAC_index_;//change to IUPAC_index_
                  std::string neighborResID = std::to_string(MonoNeighborNum);
                  std::string monoSNFG = thisMonoNeighbor->SNFG_name_;
                  std::string monoShortName = thisMonoNeighbor->sugar_name_.monosaccharide_short_name_;
                  child_res_resource = CreateURIResource(gmml::OntSequenceResidue, root_oligo_id, id_prefix, neighborResID);
                  child_res_uri = CreateURI(child_res_resource);
                  std::stringstream connectionInfo;
                  connectionInfo << "gmmo:is" << (*it).first->linkage_type_ << "ConnectedTo";
                  gmml::AddTriple(parent_res_uri, connectionInfo.str(), child_res_uri, oligo_sequence_stream);
                  PopulateLinkage(linkage_stream, oligo_uri, parent_res_uri, child_res_uri, linkNum, (*it).first, thisMono, thisMonoNeighbor);
                  
                }
                PopulateMonosaccharide(mono_stream, oligo_stream, oligo_uri, id_prefix, thisMono, side_or_ring_atoms, pdb_uri);
              }
              
            }
            
            MolecularModeling::Assembly subAssembly(residueVector);
            subAssembly.SetModelIndex(0);
            PdbFileSpace::PdbFile* thisPDB = subAssembly.BuildPdbFileStructureFromAssembly();
            std::ostringstream PDBstringstream;
            thisPDB->WriteToStringstream(PDBstringstream);
            gmml::AddLiteral(oligo_uri, "gmmo:PDBfile", PDBstringstream.str(), oligo_stream);
            
            // if(oligo->child_oligos_.size() != 0 && (find(visited_oligos.begin(), visited_oligos.end(), oligo->root_->mono_id_) == visited_oligos.end()))
            // {
            //   // gmml::log(__LINE__, __FILE__,  gmml::INF, " ");
            //     PopulateLinkage(linkage_stream, oligo, oligo_uri, id_prefix, link_id, visited_oligos);
            //     // gmml::log(__LINE__, __FILE__,  gmml::INF, "About to populate sequence linkages");
            //     PopulateSequenceLinkage(oligo_sequence_stream, oligo, oligo_uri, id_prefix, visited_oligos, mono_to_short_name_map, oligo_to_res_uri_map, root_oligo_id);
            //     // gmml::log(__LINE__, __FILE__,  gmml::INF, "Done populating sequence linkages");
            // }
            // else if(oligo->child_oligos_.size() == 0 && o_name.compare("") != 0 && oligo->oligosaccharide_terminal_.compare("") != 0)
            // {
            //   // gmml::log(__LINE__, __FILE__,  gmml::INF, " ");
            //     std::string term_resource = "";
            //     std::string term_uri = "";
            //     term_resource = CreateURIResource(gmml::OntTerminal, root_oligo_id, id_prefix, "");
            //     term_uri = CreateURI(term_resource);
            //     gmml::AddTriple(oligo_uri, Ontology::hasTerminal, term_uri, oligo_sequence_stream);
            //     gmml::AddTriple(term_uri, Ontology::TYPE, Ontology::Terminal, oligo_sequence_stream);
            //     gmml::AddLiteral(term_uri, Ontology::id, oligo->oligosaccharide_terminal_, oligo_sequence_stream);
            //     //std::cout << "Terminalll " << o_term_name << std::endl;
            //     std::stringstream res_id;
            //     std::string res_resource;
            //     std::string res_uri;
            //     res_id << "1";
            //     res_resource = CreateURIResource(gmml::OntSequenceResidue, root_oligo_id, id_prefix, res_id.str());
            //     res_uri = CreateURI(res_resource);
            //     std::string mono_short_name = oligo->root_->sugar_name_.monosaccharide_short_name_;
            //     std::string monoSNFG = oligo->root_->SNFG_name_;
            //     // gmml::log(__LINE__, __FILE__,  gmml::INF, "About to check derivative map");
            //     CheckDerivativesAndPopulate(oligo_sequence_stream, mono_short_name, oligo_uri, res_uri, monoSNFG);
            //     gmml::AddTriple(res_uri, Ontology::isConnectedTo, term_uri, oligo_sequence_stream);
            // }
            // else if(oligo->child_oligos_.size() == 0 && o_name.compare("") != 0)
            // {
            //   // gmml::log(__LINE__, __FILE__,  gmml::INF, " ");
            //     std::stringstream res_id;
            //     std::string res_resource;
            //     std::string res_uri;
            //     res_id << "1";
            //     res_resource = CreateURIResource(gmml::OntSequenceResidue, root_oligo_id, id_prefix, res_id.str());
            //     res_uri = CreateURI(res_resource);
            //     std::string mono_short_name = oligo->root_->sugar_name_.monosaccharide_short_name_;
            //     std::string monoSNFG = oligo->root_->SNFG_name_;
            //     // gmml::log(__LINE__, __FILE__,  gmml::INF, "About to check derivative map");
            //     CheckDerivativesAndPopulate(oligo_sequence_stream, mono_short_name, oligo_uri, res_uri, monoSNFG);
            // }
            // // gmml::log(__LINE__, __FILE__,  gmml::INF, " ");
            // Glycan::Monosaccharide* mono = oligo->root_;
            // PopulateMonosaccharide(mono_stream, oligo_stream, oligo_uri, id_prefix, mono, side_or_ring_atoms, pdb_uri);

            // std::vector<Glycan::Oligosaccharide*> child_oligos = oligo->child_oligos_;
            // PopulateOligosaccharide(pdb_stream, oligo_stream, oligo_sequence_stream, mono_stream, linkage_stream, pdb_uri, id_prefix, link_id, child_oligos, side_or_ring_atoms, visited_oligos, mono_to_short_name_map, oligo_to_res_uri_map, root_oligo_id);
        }
        
    }
}

void Assembly::CheckDerivativesAndPopulate(std::stringstream& oligo_sequence_stream, std::string mono_short_name, std::string oligo_uri, std::string res_uri, std::string monoSNFG, Glycan::Monosaccharide* mono)
{
    if(hasDerivative(mono_short_name)){
        std::vector<std::string> derivatives;
        getDerivatives(mono_short_name, derivatives);
        gmml::AddTriple(oligo_uri, Ontology::hasSequenceResidue, res_uri, oligo_sequence_stream);
        gmml::AddTriple(res_uri, Ontology::TYPE, Ontology::SequenceResidue, oligo_sequence_stream);
        gmml::AddLiteral(res_uri, Ontology::id, mono->cycle_atoms_[0]->GetResidue()->GetId(), oligo_sequence_stream);
        gmml::AddLiteral(res_uri, Ontology::mono_short_name, mono_short_name, oligo_sequence_stream);
        gmml::AddLiteral(res_uri, Ontology::hasSNFGName, monoSNFG, oligo_sequence_stream);
        for (std::vector<std::string>::iterator t=derivatives.begin(); t!=derivatives.end(); ++t)
        {
            gmml::AddLiteral(res_uri, Ontology::seq_derivative, *t, oligo_sequence_stream);
        }

    }else{
        gmml::AddTriple(oligo_uri, Ontology::hasSequenceResidue, res_uri, oligo_sequence_stream);
        gmml::AddTriple(res_uri, Ontology::TYPE, Ontology::SequenceResidue, oligo_sequence_stream);
        gmml::AddLiteral(res_uri, Ontology::id, mono->cycle_atoms_[0]->GetResidue()->GetId(), oligo_sequence_stream);
        gmml::AddLiteral(res_uri, Ontology::mono_short_name, mono_short_name, oligo_sequence_stream);
        gmml::AddLiteral(res_uri, Ontology::hasSNFGName, monoSNFG, oligo_sequence_stream);
    }
}

bool Assembly::hasDerivative(std::string mono_short_name)
{
    std::string::size_type loc1 = mono_short_name.find( "[", 0 );
    std::string::size_type loc2 = mono_short_name.find( "]", 0 );
    if (loc1 != std::string::npos && loc2 != std::string::npos && loc1 < loc2)
        return true;
    else
        return false;
}

void Assembly::getDerivatives(std::string& mono_short_name, std::vector<std::string>& derivatives){
    std::string::size_type loc1 = mono_short_name.find( "[", 0 );
    std::string::size_type loc2 = mono_short_name.find( "]", 0 );
    std::string deriv = "";
    if (loc1 != std::string::npos && loc2 != std::string::npos && loc1 < loc2) {
        deriv = mono_short_name.substr(loc1, loc2-loc1+1);
        gmml::FindReplaceString(mono_short_name, deriv, "");
        deriv = deriv.substr(1, deriv.length()-2);
        derivatives = gmml::Split(deriv, ",");
    }
}

void Assembly::PopulateLinkage(std::stringstream& linkage_stream, std::string oligo_uri, std::string parent_res_uri, std::string child_res_uri, int& linkNum, Glycan::GlycosidicLinkage* thisLinkage, Glycan::Monosaccharide* thisMono, Glycan::Monosaccharide* thisMonoNeighbor)
{
  int local_debug = -1;
    std::string linkage_resource = "";
    std::string linkage_uri = "";
    std::string child_oligo_resource = "";
    std::string child_oligo_uri = "";
    std::string child_atom_resource = "";
    std::string child_atom_uri = "";
    std::string glycosidic_atom_resource = "";
    std::string glycosidic_atom_uri = "";
    std::string parent_atom_resource = "";
    std::string parent_atom_uri = "";
    std::stringstream linkage_str;
    std::stringstream glycosidic_linkage_str;
    
    if(thisMonoNeighbor == thisLinkage->non_reducing_mono_)
    {
      linkage_str << oligo_uri << "_link_" << linkNum;
      linkage_uri = linkage_str.str();
      gmml::AddTriple(oligo_uri, Ontology::hasGlycosidicLinkage, linkage_uri, linkage_stream);
      gmml::AddTriple(linkage_uri,Ontology::TYPE, Ontology::Linkage, linkage_stream);
      gmml::AddTriple(linkage_uri, Ontology::hasParent, parent_res_uri, linkage_stream);
      gmml::AddTriple(linkage_uri, Ontology::hasChild, child_res_uri, linkage_stream);
      gmml::AddLiteral(linkage_uri, Ontology::hasParentMono, thisMono->sugar_name_.monosaccharide_short_name_, linkage_stream);
      gmml::AddLiteral(linkage_uri, Ontology::hasChildMono, thisMonoNeighbor->sugar_name_.monosaccharide_short_name_, linkage_stream);
      gmml::AddLiteral(linkage_uri, Ontology::linkageIndeces, thisLinkage->inverse_linkage_type_, linkage_stream);
      gmml::AddDecimal(linkage_uri, Ontology::hasGlycosidicPhiAngle, thisLinkage->phi_angle_, linkage_stream);
      gmml::AddDecimal(linkage_uri, Ontology::hasGlycosidicPsiAngle, thisLinkage->psi_angle_, linkage_stream);
      if(thisLinkage->omega_angle_ != -9999)
        gmml::AddDecimal(linkage_uri, Ontology::hasGlycosidicOmegaAngle, thisLinkage->omega_angle_, linkage_stream);
      linkNum++;
    }
    
    // for(OligosaccharideVector::iterator it = oligo->child_oligos_.begin(); it != oligo->child_oligos_.end(); it++)
    // {
    //     int index = distance(oligo->child_oligos_.begin(), it);
    // 
    //     Glycan::Oligosaccharide* child_oligo = (*it);
    //     // OligosaccharideVector::iterator it2 = (std::next(it,1));
    //     // Glycan::Oligosaccharide* parent_oligo = (*it2);
    //     //        visited_oligos.push_back(child_oligo->root_->mono_id_);
    // 
    //     linkage_resource = CreateURIResource(gmml::OntLinkage, link_id, id_prefix, "");
    //     linkage_uri = CreateURI(linkage_resource);
    //     //        linkage_stream << Ontology::ENTITY_COMMENT << linkage_resource << std::endl;
    //     gmml::AddTriple(linkage_uri, Ontology::TYPE, Ontology::Linkage, linkage_stream);
    //     //        gmml::AddLiteral(linkage_uri, Ontology::LABEL, linkage_resource, linkage_stream);
    //     link_id++;
    // 
    //     gmml::AddTriple(linkage_uri, Ontology::hasParent, oligo_uri, linkage_stream);
    //     child_oligo_resource = CreateURIResource(gmml::OntOligosaccharide, child_oligo->root_->mono_id_, id_prefix, "");
    //     child_oligo_uri = CreateURI(child_oligo_resource);
    //     gmml::AddTriple(linkage_uri, Ontology::hasChild, child_oligo_uri, linkage_stream);
    // 
    //     std::vector<std::string> linkage_tokens = gmml::Split(oligo->child_oligos_linkages_.at(index), "-");
    //     std::string parent_atom_id = linkage_tokens.at(0);
    //     std::string glycosidic_atom_id = linkage_tokens.at(1);
    //     std::string child_atom_id = linkage_tokens.at(2);
    // 
    //     int parent_c_index = ExtractLinkageCarbonIndex(oligo, parent_atom_id);
    //     int child_c_index = ExtractLinkageCarbonIndex(child_oligo, child_atom_id);
    // 
    //     if(child_c_index != 0 && parent_c_index != 0)
    //     {
    //         std::stringstream link_indeces_str;
    //         link_indeces_str << child_c_index << "-" << parent_c_index;
    //         gmml::AddLiteral(linkage_uri, Ontology::linkageIndeces, link_indeces_str.str(), linkage_stream);
    //     }
    // 
    //     child_atom_resource = CreateURIResource(gmml::OntAtom, 0, id_prefix, child_atom_id);
    //     child_atom_uri = CreateURI(child_atom_resource);
    //     gmml::AddTriple(linkage_uri, Ontology::hasChildAtomLinkage, child_atom_uri, linkage_stream);
    // 
    //     glycosidic_atom_resource = CreateURIResource(gmml::OntAtom, 0, id_prefix, glycosidic_atom_id);
    //     glycosidic_atom_uri = CreateURI(glycosidic_atom_resource);
    //     gmml::AddTriple(linkage_uri, Ontology::hasGlycosidicLinkage, glycosidic_atom_uri, linkage_stream);
    // 
    // 
    //     gmml::AddLiteral(linkage_uri, Ontology::hasChildMono, child_oligo->root_->sugar_name_.monosaccharide_short_name_, linkage_stream);
    //     gmml::AddLiteral(linkage_uri, Ontology::hasParentMono, oligo->root_->sugar_name_.monosaccharide_short_name_, linkage_stream);
    //     double glycosidic_phi_angle = CalculatePhiAngle(child_oligo, parent_atom_id, child_atom_id, glycosidic_atom_id);
    //     glycosidic_phi_angle = gmml::ConvertRadian2Degree(glycosidic_phi_angle);
    //     if(local_debug > 0)
    //     {
    //       gmml::log(__LINE__, __FILE__, gmml::INF, "Phi");
    //       gmml::log(__LINE__, __FILE__, gmml::INF, std::to_string(glycosidic_phi_angle));
    //     }
    //     gmml::AddTriple(linkage_uri, Ontology::hasGlycosidicPhiAngle, std::to_string(glycosidic_phi_angle), linkage_stream);
    // 
    // 
    //     double glycosidic_psi_angle = CalculatePsiAngle(child_oligo, parent_atom_id, child_atom_id, glycosidic_atom_id);
    //     glycosidic_psi_angle = gmml::ConvertRadian2Degree(glycosidic_psi_angle);
    //     if(local_debug > 0)
    //     {
    //       gmml::log(__LINE__, __FILE__, gmml::INF, "Psi");
    //       gmml::log(__LINE__, __FILE__, gmml::INF, std::to_string(glycosidic_psi_angle));
    //     }
    //     gmml::AddTriple(linkage_uri, Ontology::hasGlycosidicPsiAngle, std::to_string(glycosidic_psi_angle), linkage_stream);
    // 
    //     if (parent_c_index == 6)
    //     {
    //       gmml::log(__LINE__, __FILE__,  gmml::INF, "About to calculate Omega Angle");
    //       double glycosidic_omega_angle = CalculateOmegaAngle(oligo, parent_atom_id, glycosidic_atom_id);
    //       glycosidic_omega_angle = gmml::ConvertRadian2Degree(glycosidic_omega_angle);
    //       if(local_debug > 0)
    //       {
    //         gmml::log(__LINE__, __FILE__, gmml::INF, "Omega");
    //         gmml::log(__LINE__, __FILE__, gmml::INF, std::to_string(glycosidic_omega_angle));
    //       }
    //       gmml::AddTriple(linkage_uri, Ontology::hasGlycosidicOmegaAngle, std::to_string(glycosidic_omega_angle), linkage_stream);
    // 
    //     }
    //     // gmml::log(__LINE__, __FILE__,  gmml::INF, "Done with angles");
    //     parent_atom_resource = CreateURIResource(gmml::OntAtom, 0, id_prefix, parent_atom_id);
    //     parent_atom_uri = CreateURI(parent_atom_resource);
    //     // gmml::log(__LINE__, __FILE__,  gmml::INF, "About to add parent linkage");
    //     gmml::AddTriple(linkage_uri, Ontology::hasParentAtomLinkage, parent_atom_uri, linkage_stream);
    // 
    //     std::vector<std::string> child_atom_id_tokens = gmml::Split(child_atom_id, "_");
    //     if(child_atom_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
    //         linkage_str << child_atom_id_tokens.at(2) << "(" << child_atom_id_tokens.at(4) << ")" << child_atom_id_tokens.at(0);
    //     else
    //         linkage_str << child_atom_id_tokens.at(2) << "(" << child_atom_id_tokens.at(4) << "_" << child_atom_id_tokens.at(3) << ")" << child_atom_id_tokens.at(0);
    // 
    //     std::vector<std::string> parent_atom_id_tokens = gmml::Split(parent_atom_id, "_");
    //     if(parent_atom_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
    //         linkage_str << "-" << parent_atom_id_tokens.at(2) << "(" << parent_atom_id_tokens.at(4) << ")"  << parent_atom_id_tokens.at(0);
    //     else
    //         linkage_str << "-" << parent_atom_id_tokens.at(2) << "(" << parent_atom_id_tokens.at(4) <<  "_" << parent_atom_id_tokens.at(3) << ")"  << parent_atom_id_tokens.at(0);
    // 
    //     //        gmml::AddLiteral(linkage_uri, Ontology::linkage_str, linkage_str.str(), linkage_stream);
    // 
    //     std::vector<std::string> glycosidic_atom_id_tokens = gmml::Split(glycosidic_atom_id, "_");
    //     if(glycosidic_atom_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
    //         glycosidic_linkage_str << glycosidic_atom_id_tokens.at(2) << "(" << glycosidic_atom_id_tokens.at(4) << ")" << glycosidic_atom_id_tokens.at(0);
    //     else
    //         glycosidic_linkage_str << glycosidic_atom_id_tokens.at(2) << "(" << glycosidic_atom_id_tokens.at(4) << "_" << glycosidic_atom_id_tokens.at(3)
    //                                << ")"  << glycosidic_atom_id_tokens.at(0);
    //     // gmml::log(__LINE__, __FILE__,  gmml::INF, "About to add linkage info");
    //     gmml::AddLiteral(linkage_uri, Ontology::glycosidic_linkage, glycosidic_linkage_str.str(), linkage_stream);
    // }
    // gmml::log(__LINE__, __FILE__,  gmml::INF, "Done populating linkages");
}

void Assembly::PopulateSequenceLinkage(std::stringstream& oligo_sequence_stream, Glycan::Oligosaccharide* oligo, std::string oligo_uri, std::string id_prefix, std::vector<int>& visited_oligos,
                               std::map<std::string, std::string>& mono_to_short_name_map, std::map<std::string, std::string>& oligo_to_res_uri_map, int& root_oligo_id)
{
    // std::string child_oligo_resource = "";
    // std::string child_oligo_uri = "";
    // std::string child_res_resource = "";
    // std::string child_res_uri = "";
    // std::string root_oligo_resource = "";
    // std::string root_oligo_uri = "";
    // std::string term_resource = "";
    // std::string term_uri = "";
    // std::string o_name = oligo->oligosaccharide_name_;
    // std::string o_term_name = oligo->oligosaccharide_terminal_;
    // bool terminalAdded = false;
    // 
    // 
    // root_oligo_resource = CreateURIResource(gmml::OntOligosaccharide, root_oligo_id, id_prefix, "");
    // root_oligo_uri = CreateURI(root_oligo_resource);
    // 
    // 
    // 
    // visited_oligos.push_back(oligo->root_->mono_id_);
    // 
    // for(std::vector<Glycan::Monosaccharide*>::reverse_iterator rit = oligo->mono_nodes_.rbegin(); rit != oligo->mono_nodes.rend(); rit++)
    // {
    //   thisMono = *rit;
    //   int MonoNum = oligo->mono_nodes_.size() - thisMono->oligosaccharide_index_;
    //   parent_mono_resource = CreateURIResource(gmml::OntOligosaccharide, thisMono->mono_id_, id_prefix, "");
    //   parent_mono_uri = CreateURI(parent_mono_resource);
    //   std::string resID = std::to_string(MonoNum);
    //   std::string monoSNFG = thisMono->SNFG_name_;
    //   std::string monoShortName = thisMono->sugar_name_.monosaccharide_short_name_;
    //   parent_res_resource = CreateURIResource(gmml::OntSequenceResidue, root_oligo_id, id_prefix, resID);
    //   parent_res_uri = CreateURI(parent_res_resource);
    //   CheckDerivativesAndPopulate(oligo_sequence_stream, monoShortName, root_oligo_uri, parent_res_uri, monoSNFG);
    //   for(std::vector<std::pair<Glycan::GlycosidicLinkage*, Glycan::Monosaccharide*> >::iterator it = thisMono->mono_neighbors_.begin(); it!=thisMono->mono_neighbors_.end(); it++)
    //   {
    //     thisMonoNeighbor = *it.second;
    //     int MonoNeighborNum = oligo->mono_nodes_.size() - thisMonoNeighbor->oligosaccharide_index_;
    //     std::string resID = std::to_string(MonoNeighborNum);
    //     std::string monoSNFG = thisMonoNeighbor->SNFG_name_;
    //     std::string monoShortName = thisMonoNeighbor->sugar_name_.monosaccharide_short_name_;
    //     parent_res_resource = CreateURIResource(gmml::OntSequenceResidue, root_oligo_id, id_prefix, resID);
    //     parent_res_uri = CreateURI(parent_res_resource);
    //     gmml::AddTriple(child_res_uri, Ontology::isConnectedTo, parent_res_uri, oligo_sequence_stream);
    //   }
    // 
    // }
    // for(OligosaccharideVector::iterator it = oligo->child_oligos_.begin(); it != oligo->child_oligos_.end(); it++)
    // {
    //     int index = distance(oligo->child_oligos_.begin(), it);
    // 
    //     Glycan::Oligosaccharide* child_oligo = (*it);
    // 
    //     child_oligo_resource = CreateURIResource(gmml::OntOligosaccharide, child_oligo->root_->mono_id_, id_prefix, "");
    //     child_oligo_uri = CreateURI(child_oligo_resource);
    // 
    //     if (mono_to_short_name_map.find(child_oligo_uri) == mono_to_short_name_map.end())
    //     {
    //         mono_to_short_name_map[child_oligo_uri] = child_oligo->root_->sugar_name_.monosaccharide_short_name_;
    //         std::string monoSNFG = child_oligo->root_->SNFG_name_;
    //         std::stringstream res_id;
    //         res_id << mono_to_short_name_map.size();
    //         child_res_resource = CreateURIResource(gmml::OntSequenceResidue, root_oligo_id, id_prefix, res_id.str());
    //         child_res_uri = CreateURI(child_res_resource);
    // 
    //         CheckDerivativesAndPopulate(oligo_sequence_stream, mono_to_short_name_map[child_oligo_uri], root_oligo_uri, child_res_uri, monoSNFG);
    // 
    //         //gmml::AddTriple(root_oligo_uri, Ontology::hasSequenceResidue, child_res_uri, oligo_sequence_stream);
    //         //gmml::AddTriple(child_res_uri, Ontology::TYPE, Ontology::SequenceResidue, oligo_sequence_stream);
    //         //gmml::AddLiteral(child_res_uri, Ontology::id, mono_to_short_name_map[child_oligo_uri], oligo_sequence_stream);
    //         oligo_to_res_uri_map[child_oligo_uri] = child_res_uri;
    //     }
    //     if (mono_to_short_name_map.find(oligo_uri) == mono_to_short_name_map.end())
    //     {
    //         mono_to_short_name_map[oligo_uri] = oligo->root_->sugar_name_.monosaccharide_short_name_;
    //         std::string monoSNFG = oligo->root_->SNFG_name_;
    //         std::stringstream res_id;
    //         res_id << mono_to_short_name_map.size();
    //         child_res_resource = CreateURIResource(gmml::OntSequenceResidue, root_oligo_id, id_prefix, res_id.str());
    //         child_res_uri = CreateURI(child_res_resource);
    // 
    //         CheckDerivativesAndPopulate(oligo_sequence_stream, mono_to_short_name_map[oligo_uri], root_oligo_uri, child_res_uri, monoSNFG);
    // 
    //         //gmml::AddTriple(root_oligo_uri, Ontology::hasSequenceResidue, child_res_uri, oligo_sequence_stream);
    //         //gmml::AddTriple(child_res_uri, Ontology::TYPE, Ontology::SequenceResidue, oligo_sequence_stream);
    //         //gmml::AddLiteral(child_res_uri, Ontology::id, mono_to_short_name_map[oligo_uri], oligo_sequence_stream);
    //         oligo_to_res_uri_map[oligo_uri] = child_res_uri;
    //     }
    // 
    //     std::string child_uri = oligo_to_res_uri_map[child_oligo_uri];
    //     std::string parent_uri = oligo_to_res_uri_map[oligo_uri];
    //     gmml::AddTriple(child_uri, Ontology::isConnectedTo, parent_uri, oligo_sequence_stream);
    // 
    //     if(o_name.compare("") != 0 && o_term_name.compare("") != 0 && !terminalAdded)
    //     {
    //         gmml::AddTriple(parent_uri, Ontology::isConnectedTo, term_uri, oligo_sequence_stream);
    //         terminalAdded = true;
    //     }
    // 
    //     std::vector<std::string> linkage_tokens = gmml::Split(oligo->child_oligos_linkages_.at(index), "-");
    //     std::string parent_atom_id = linkage_tokens.at(0);
    //     std::string child_atom_id = linkage_tokens.at(2);
    // 
    //     int parent_c_index = ExtractLinkageCarbonIndex(oligo, parent_atom_id);
    //     int child_c_index = ExtractLinkageCarbonIndex(child_oligo, child_atom_id);
    // 
    //     if(child_c_index != 0 && parent_c_index != 0)
    //     {
    //         std::stringstream link_indeces_str;
    //         link_indeces_str << child_c_index << "-" << parent_c_index;
    //         gmml::AddLiteral(oligo_to_res_uri_map[child_oligo_uri], Ontology::sequence_linkage, link_indeces_str.str(), oligo_sequence_stream);
    //     }
    // }
    // if(o_name.compare("") != 0 && o_term_name.compare("") != 0)
    // {
    //     term_resource = CreateURIResource(gmml::OntTerminal, root_oligo_id, id_prefix, "");
    //     term_uri = CreateURI(term_resource);
    //     gmml::AddTriple(root_oligo_uri, Ontology::hasTerminal, term_uri, oligo_sequence_stream);
    //     gmml::AddTriple(term_uri, Ontology::TYPE, Ontology::Terminal, oligo_sequence_stream);
    //     gmml::AddLiteral(term_uri, Ontology::id, o_term_name, oligo_sequence_stream);
    //     //std::cout << "Terminalll " << o_term_name << std::endl;
    // }
}

int Assembly::ExtractLinkageCarbonIndex(Glycan::Oligosaccharide* oligo, std::string linkage_carbon_id)
{
    int c_index = 0;
    std::vector<std::string> cycle_atom_tokens = gmml::Split(oligo->root_->cycle_atoms_str_, "-");

    if(oligo->root_->side_atoms_.at(0).at(0) != NULL)
    {
        c_index++;
        Atom* anomeric_side_carbon = oligo->root_->side_atoms_.at(0).at(0);
        if(anomeric_side_carbon->GetId().compare(linkage_carbon_id) == 0)
            return c_index;
    }
    if(!cycle_atom_tokens.empty())
    {
      for(unsigned int i = 0; i < cycle_atom_tokens.size() - 1; i++) /// cycle_atom_tokens.size() - 1 > because the ring oxygen is not considered
      {
          c_index++;
          if(cycle_atom_tokens.at(i).compare(linkage_carbon_id) == 0)
              return c_index;
      }
    }
    AtomVector side_atoms_of_last_ring_carbon = oligo->root_->side_atoms_.at(oligo->root_->side_atoms_.size() - 1);
    if (!side_atoms_of_last_ring_carbon.empty()){
    for(AtomVector::iterator it1 = side_atoms_of_last_ring_carbon.begin(); it1 != side_atoms_of_last_ring_carbon.end(); it1++)
    {
        Atom* side_atom = (*it1);
        c_index++;
        if(side_atom != NULL)
        {
          if(side_atom->GetId().compare(linkage_carbon_id) == 0)
            return c_index;
        }
    }
  }
    return c_index;
}

void Assembly::PopulateMonosaccharide(std::stringstream& mono_stream, std::stringstream& oligo_stream, std::string oligo_uri, std::string id_prefix, Glycan::Monosaccharide* mono,
                                      std::vector<std::string>& side_or_ring_atoms, std::string pdb_uri)
{
    std::stringstream object;
    std::string mono_resource = "";
    std::string mono_uri = "";
    std::string ring_resource = "";
    std::string ring_uri = "";

    mono_resource = CreateURIResource(gmml::OntMonosaccharide, mono->mono_id_, id_prefix, "");
    mono_uri = CreateURI(mono_resource);
    gmml::AddTriple(pdb_uri, Ontology::hasMono, mono_uri, mono_stream);
    gmml::AddTriple(oligo_uri, Ontology::hasMono, mono_uri, oligo_stream);
    //    mono_stream << Ontology::ENTITY_COMMENT << mono_resource << std::endl;
    gmml::AddTriple(mono_uri, Ontology::TYPE, Ontology::monosaccharide, mono_stream);
    gmml::AddLiteral(mono_uri, Ontology::id, mono->cycle_atoms_[0]->GetResidue()->GetId(), mono_stream);
    gmml::AddLiteral(mono_uri, Ontology::hasOligoParent, oligo_uri, mono_stream);
    int Index = mono->IUPAC_index_; //change to IUPAC_index_
    gmml::AddLiteral(mono_uri, Ontology::hasIndex, std::to_string(Index), mono_stream);
    gmml::AddLiteral(mono_uri, Ontology::hasNameIndex, std::to_string(mono->oligosaccharide_index_), mono_stream);
    if(checkIfNucleotide(mono))
    {
      gmml::AddTriple(mono_uri, Ontology::isNucleotide, "true", mono_stream);
    }
    else
    {
      gmml::AddTriple(mono_uri, Ontology::isNucleotide, "false", mono_stream);
    }
    // gmml::AddLiteral(mono_uri, Ontology::hasSNFGName, mono->SNFG_name_, mono_stream);
    //    gmml::AddLiteral(mono_uri, Ontology::LABEL, mono_resource, mono_stream);

    AtomVector ring_atoms = mono->cycle_atoms_;
    object.str(std::string());
    int ring_index = 1;
    std::stringstream ring_atom_stream;
    if(ring_atoms.size() > 0)
    {
      for(AtomVector::iterator it = ring_atoms.begin(); it != ring_atoms.end(); it++)
      {
          Atom* ring_atom = (*it);
          if (ring_atom != NULL)
          {
            ring_resource = CreateURIResource(gmml::OntAtom, 0, id_prefix, ring_atom->GetId());
            ring_uri = CreateURI(ring_resource);
            gmml::AddTriple(mono_uri, Ontology::hasRingAtom, ring_uri, mono_stream);

            PopulateRingAtom(ring_atom_stream, id_prefix, ring_uri, ring_resource, ring_index, ring_atom, mono, side_or_ring_atoms);
            ring_index++;

            if(it == ring_atoms.end() - 1)
                object << ring_resource;
            else
                object << ring_resource << "-";
          }
      }
    }
    gmml::AddLiteral(mono_uri, Ontology::ring_atoms, object.str(), mono_stream);

    object.str(std::string());
    object << mono->anomeric_status_ << " " << CreateURIResource(gmml::OntAtom, 0, id_prefix, mono->cycle_atoms_.at(0)->GetId());
    gmml::AddLiteral(mono_uri, Ontology::anomeric_status, object.str(), mono_stream);

    gmml::AddLiteral(mono_uri, Ontology::stereochemistry_chemical_code, mono->sugar_name_.chemical_code_string_, mono_stream);
    if(mono->bfmp_ring_conformation_.compare("") != 0) {
      std::string bfmp = mono->bfmp_ring_conformation_;
      std::size_t index;
      // This is to only output two decimal places for BFMP values which have a decimal value.
      if( ( index = mono->bfmp_ring_conformation_.find( "." ) ) != std::string::npos ) {
        bfmp = bfmp.substr( 0, index + 3 ) + ")";
      }
      gmml::AddLiteral(mono_uri, Ontology::bfmp_ring_conformation, bfmp, mono_stream);
    }
    gmml::AddLiteral(mono_uri, Ontology::hasSNFGName, mono->SNFG_name_, mono_stream);
    gmml::AddLiteral(mono_uri, Ontology::hasAuthorSNFGName, mono->author_SNFG_name_, mono_stream);
    std::size_t offset = 0;
    int numR;
    if(mono->on_R_ > 0)
    {
      for(std::vector<std::pair<std::string, std::string> >::iterator derivative = mono->unknown_derivatives_.begin(); derivative != mono->unknown_derivatives_.end(); derivative++)
      {
        std::size_t found = mono->sugar_name_.monosaccharide_short_name_.find("<", offset);
        if(found != std::string::npos)
        {
          found+=2;
          if(mono->on_R_ < 10)
            numR = std::stoi(mono->sugar_name_.monosaccharide_short_name_.substr(found, 1));
          else
            numR = std::stoi(mono->sugar_name_.monosaccharide_short_name_.substr(found,2));
            
          if((*derivative).second != "")
          {
            std::stringstream RgroupStream, RnumStream;
            RgroupStream << mono_uri << "_R" << numR;
            RnumStream << "gmmo:hasR" << numR;
            gmml::AddTriple(mono_uri, RnumStream.str(), RgroupStream.str(), mono_stream); 
            gmml::AddLiteral(RgroupStream.str(), Ontology::hasFormula,(*derivative).second, mono_stream); 
          }
          offset = found;
        }
        else
          break;
      }
    }
    gmml::AddLiteral(mono_uri, Ontology::author_mono_name, mono->author_sugar_name_.monosaccharide_name_, mono_stream);
    NoteVector notes = mono->mono_notes_;
    if(notes.size() != 0)
    {
        int note_id = 1;
        // gmml::log(__LINE__, __FILE__,  gmml::INF, "Poulating notes");
        std::string id_prefix = mono_uri.substr(5,mono_uri.length()) + "_";
        PopulateNotes(mono_stream, mono_stream, mono_uri, notes, id_prefix, note_id);
    }
    Glycan::SugarName sugar_name = mono->sugar_name_;
    PopulateSugarName(mono_stream, id_prefix, mono_uri, mono->mono_id_, sugar_name);
    mono_stream << ring_atom_stream.str();

}

void Assembly::PopulateRingAtom(std::stringstream& ring_atom_stream, std::string id_prefix, std::string ring_uri, std::string ring_resource, unsigned int ring_index, Atom* ring_atom, Glycan::Monosaccharide* mono,
                                std::vector<std::string>& side_or_ring_atoms)
{
  if(ring_index <= mono->side_atoms_.size())
  {
    std::stringstream object;
    //    ring_atom_stream << Ontology::ENTITY_COMMENT << ring_resource << std::endl;
    gmml::AddTriple(ring_uri, Ontology::TYPE, Ontology::RingAtom, ring_atom_stream);
    object << ring_index;
    gmml::AddLiteral(ring_uri, Ontology::ring_index, object.str(), ring_atom_stream);
    gmml::AddLiteral(ring_uri, Ontology::id, ring_resource, ring_atom_stream);
    //    gmml::AddLiteral(ring_uri, Ontology::LABEL, ring_resource, ring_atom_stream);
    GeometryTopology::Coordinate* coords = ring_atom->GetCoordinates().at(model_index_);
    /*
    gmml::AddLiteral(ring_uri, Ontology::x, gmml::ConvertT<double>(coords->GetX()), ring_atom_stream);
    gmml::AddLiteral(ring_uri, Ontology::y, gmml::ConvertT<double>(coords->GetY()), ring_atom_stream);
    gmml::AddLiteral(ring_uri, Ontology::z, gmml::ConvertT<double>(coords->GetZ()), ring_atom_stream);
    */
    std::stringstream coord_stream;
    coord_stream << gmml::ConvertT<double>(coords->GetX()) << ", " << gmml::ConvertT<double>(coords->GetY()) << ", " << gmml::ConvertT<double>(coords->GetZ());
    gmml::AddLiteral(ring_uri, Ontology::coordinate, coord_stream.str(), ring_atom_stream);

    side_or_ring_atoms.push_back(ring_atom->GetId());

    // std::string neighbor_resource = "";
    // std::string neighbor_uri = "";
    // AtomVector neighbors = ring_atom->GetNode()->GetNodeNeighbors();
    // for(AtomVector::iterator it = neighbors.begin(); it != neighbors.end(); it++)
    // {
    //     Atom* neighbor = (*it);
    //     neighbor_resource = CreateURIResource(gmml::OntAtom, 0, id_prefix, neighbor->GetId());
    //     neighbor_uri = CreateURI(neighbor_resource);
    //     gmml::AddTriple(ring_uri, Ontology::hasNeighbor, neighbor_uri, ring_atom_stream);
    // }

    std::stringstream side_atom_stream;
    if((ring_atom->GetName().substr(0,1).compare("O") != 0 )) ///side atoms for the oxygen of the ring are not saved
    {
        std::vector<AtomVector> all_sides = mono->side_atoms_;
        AtomVector sides = all_sides.at(ring_index - 1);
        std::string side_resource = "";
        std::string side_uri = "";
        int side_index = 1;
        for(AtomVector::iterator it = sides.begin(); it != sides.end(); it++)
        {
            Atom* side_atom = (*it);
            if(side_atom != NULL)
            {
                side_resource = CreateURIResource(gmml::OntAtom, 0, id_prefix, side_atom->GetId());
                side_uri = CreateURI(side_resource);
                gmml::AddTriple(ring_uri, Ontology::hasSideAtom, side_uri, ring_atom_stream);

                PopulateSideAtom(side_atom_stream, id_prefix, side_uri, side_resource, ring_index, side_index, side_atom, mono, side_or_ring_atoms);
                side_index++;
            }
        }
    }
    ring_atom_stream << side_atom_stream.str();
  }
}

void Assembly::PopulateSideAtom(std::stringstream& side_atom_stream, std::string id_prefix, std::string side_uri, std::string side_resource, int ring_index, int side_index, Atom* side_atom, Glycan::Monosaccharide* mono,
                                std::vector<std::string>& side_or_ring_atoms)
{
    std::stringstream object;

    std::string chemical_code_str = mono->sugar_name_.chemical_code_string_;
    if((mono->sugar_name_.ring_type_.compare("P") == 0 && ring_index == 5) || (mono->sugar_name_.ring_type_.compare("F") == 0 && ring_index == 4))
        object << "+" << side_index;
    else if(ring_index == 1 && (side_atom->GetName().substr(0,1).compare("C") != 0))
        object << "1";
    else if(ring_index == 1 && (side_atom->GetName().substr(0,1).compare("C") == 0 ))
        object << "-1";
    else
        object << ring_index;

    if(find(side_or_ring_atoms.begin(), side_or_ring_atoms.end(), side_atom->GetId()) == side_or_ring_atoms.end())///if this side atom has not been added to the ontology as side atom of another mono
    {
        //    side_atom_stream << Ontology::ENTITY_COMMENT << side_resource << std::endl;
        gmml::AddTriple(side_uri, Ontology::TYPE, Ontology::SideAtom, side_atom_stream);
        gmml::AddLiteral(side_uri, Ontology::id, side_resource, side_atom_stream);
        //        gmml::AddLiteral(side_uri, Ontology::LABEL, side_resource, side_atom_stream);
        GeometryTopology::Coordinate* coords = side_atom->GetCoordinates().at(model_index_);
        /*gmml::AddLiteral(side_uri, Ontology::x, gmml::ConvertT<double>(coords->GetX()), side_atom_stream);
        gmml::AddLiteral(side_uri, Ontology::y, gmml::ConvertT<double>(coords->GetY()), side_atom_stream);
        gmml::AddLiteral(side_uri, Ontology::z, gmml::ConvertT<double>(coords->GetZ()), side_atom_stream);*/

        std::stringstream coord_stream;
        coord_stream << gmml::ConvertT<double>(coords->GetX()) << ", " << gmml::ConvertT<double>(coords->GetY()) << ", " << gmml::ConvertT<double>(coords->GetZ());
        gmml::AddLiteral(side_uri, Ontology::coordinate, coord_stream.str(), side_atom_stream);

        if(object.str().compare("1") == 0)
        {
            // if(mono->derivatives_map_.find("a") != mono->derivatives_map_.end())
            std::vector<std::pair<std::string, std::string> >::iterator thisPosition = std::find_if( mono->derivatives_map_.begin(), mono->derivatives_map_.end(),
              [&](const std::pair<std::string, std::string>& element){ return element.first == "a";} );
            if(thisPosition != mono->derivatives_map_.end())
                gmml::AddLiteral(side_uri, Ontology::derivative, (*thisPosition).second, side_atom_stream);
        }
        else
        {
          std::vector<std::pair<std::string, std::string> >::iterator thisPosition = std::find_if( mono->derivatives_map_.begin(), mono->derivatives_map_.end(),
              [&](const std::pair<std::string, std::string>& element){ return element.first == object.str() ;} );
            // if(mono->derivatives_map_.find(object.str()) != mono->derivatives_map_.end())
            if(thisPosition != mono->derivatives_map_.end())
                gmml::AddLiteral(side_uri, Ontology::derivative, (*thisPosition).second, side_atom_stream);
        }

        // std::string neighbor_resource = "";
        // std::string neighbor_uri = "";
        // AtomVector neighbors = side_atom->GetNode()->GetNodeNeighbors();
        // for(AtomVector::iterator it = neighbors.begin(); it != neighbors.end(); it++)
        // {
        //     Atom* neighbor = (*it);
        //     neighbor_resource = CreateURIResource(gmml::OntAtom, 0, id_prefix, neighbor->GetId());
        //     neighbor_uri = CreateURI(neighbor_resource);
        //     gmml::AddTriple(side_uri, Ontology::hasNeighbor, neighbor_uri, side_atom_stream);
        // }

        side_or_ring_atoms.push_back(side_atom->GetId());
    }

    size_t index = std::string::npos;
    if(object.str().compare("1") == 0)
        index = chemical_code_str.find("a");
    else
        index = chemical_code_str.find(object.str().c_str());
    if (index != std::string::npos)
    {
        index--;
        if(chemical_code_str.at(index) == '^')
        {
            gmml::AddLiteral(side_uri, Ontology::orientation, "Up", side_atom_stream);
        }
        else if(chemical_code_str.at(index) == '_')
        {
            gmml::AddLiteral(side_uri, Ontology::orientation, "Down", side_atom_stream);
        }
        else
        {
            gmml::AddLiteral(side_uri, Ontology::orientation, "Deoxy", side_atom_stream);
        }
        gmml::AddLiteral(side_uri, Ontology::side_index, object.str(), side_atom_stream);
    }

}

void Assembly::PopulateSugarName(std::stringstream& mono_stream, std::string id_prefix, std::string mono_uri, int mono_id_, Glycan::SugarName sugar_name)
{
    std::stringstream sugar_name_stream;

    std::string sugar_name_resource = "";
    std::string sugar_name_uri = "";

    sugar_name_resource = CreateURIResource(gmml::OntSugarName, mono_id_, id_prefix, "");
    sugar_name_uri = CreateURI(sugar_name_resource);

    gmml::AddTriple(mono_uri, Ontology::hasSugarName, sugar_name_uri, mono_stream);

    //    sugar_name_stream << Ontology::ENTITY_COMMENT << sugar_name_resource << std::endl;
    gmml::AddTriple(sugar_name_uri, Ontology::TYPE, Ontology::SugarName, sugar_name_stream);
    gmml::AddLiteral(sugar_name_uri, Ontology::id, sugar_name_resource, sugar_name_stream);
    //    gmml::AddLiteral(sugar_name_uri, Ontology::LABEL, sugar_name_resource, sugar_name_stream);
    gmml::AddLiteral(sugar_name_uri, Ontology::mono_stereo_name, sugar_name.monosaccharide_stereochemistry_name_, sugar_name_stream);
    gmml::AddLiteral(sugar_name_uri, Ontology::mono_stereo_short_name, sugar_name.monosaccharide_stereochemistry_short_name_, sugar_name_stream);
    gmml::AddLiteral(sugar_name_uri, Ontology::mono_name, sugar_name.monosaccharide_name_, sugar_name_stream);
    
    gmml::AddLiteral(sugar_name_uri, Ontology::mono_short_name, sugar_name.monosaccharide_short_name_, sugar_name_stream);
    gmml::AddLiteral(sugar_name_uri, Ontology::isomer, sugar_name.isomer_, sugar_name_stream);
    if(sugar_name.configuration_.compare("a") == 0)
        gmml::AddLiteral(sugar_name_uri, Ontology::configuration, "alpha", sugar_name_stream);
    else if(sugar_name.configuration_.compare("b") == 0)
        gmml::AddLiteral(sugar_name_uri, Ontology::configuration, "beta", sugar_name_stream);
    if(sugar_name.ring_type_.compare("") != 0)
        gmml::AddLiteral(sugar_name_uri, Ontology::ring_type, sugar_name.ring_type_, sugar_name_stream);

    mono_stream << sugar_name_stream.str();
}

void Assembly::PopulateResidue(std::stringstream& pdb_stream, std::stringstream& residue_stream, std::string pdb_uri, std::string id_prefix, ResidueVector residues, std::vector<std::string> side_or_ring_atoms)
{
    std::string res_resource = "";
    std::string res_uri = "";
    std::string atom_resource = "";
    std::string atom_uri = "";

    std::stringstream atom_stream;
    for(ResidueVector::iterator it = residues.begin(); it != residues.end(); it++)
    {
        Residue* residue = (*it);
        res_resource = CreateURIResource(gmml::OntResidue, 0, id_prefix, residue->GetId());
        res_uri = CreateURI(res_resource);

        //        residue_stream << Ontology::ENTITY_COMMENT << res_resource << std::endl;
        gmml::AddTriple(res_uri, Ontology::TYPE, Ontology::Residue, residue_stream);
        gmml::AddLiteral(res_uri, Ontology::id, res_resource, residue_stream);
        //        gmml::AddLiteral(res_uri, Ontology::LABEL, res_resource, residue_stream);

        AtomVector res_atoms = residue->GetAtoms();
        for(AtomVector::iterator it1 = res_atoms.begin(); it1 != res_atoms.end(); it1++)
        {
            Atom* atom = (*it1);
            atom_resource = CreateURIResource(gmml::OntAtom, 0, id_prefix, atom->GetId());
            atom_uri = CreateURI(atom_resource);

            gmml::AddTriple(res_uri, Ontology::hasAtom, atom_uri, residue_stream);

            if(find(side_or_ring_atoms.begin(), side_or_ring_atoms.end(), atom->GetId()) == side_or_ring_atoms.end())
            {
                PopulateAtom(atom_stream, atom_uri, atom_resource, id_prefix, atom);
            }
        }
        gmml::AddTriple(pdb_uri, Ontology::hasResidue, res_uri, pdb_stream);
    }
    residue_stream << atom_stream.str();
}

void Assembly::PopulateAtom(std::stringstream& atom_stream, std::string atom_uri, std::string atom_resource, std::string id_prefix, Atom* atom)
{
    //    atom_stream << Ontology::ENTITY_COMMENT << atom_resource << std::endl;
    gmml::AddTriple(atom_uri, Ontology::TYPE, Ontology::Atom, atom_stream);
    gmml::AddLiteral(atom_uri, Ontology::id, atom_resource, atom_stream);
    //    gmml::AddLiteral(atom_uri, Ontology::LABEL, atom_resource, atom_stream);
    GeometryTopology::Coordinate* coords = atom->GetCoordinates().at(model_index_);
    /*gmml::AddLiteral(atom_uri, Ontology::x, gmml::ConvertT<double>(coords->GetX()), atom_stream);
    gmml::AddLiteral(atom_uri, Ontology::y, gmml::ConvertT<double>(coords->GetY()), atom_stream);
    gmml::AddLiteral(atom_uri, Ontology::z, gmml::ConvertT<double>(coords->GetZ()), atom_stream);*/

    std::stringstream coord_stream;
    coord_stream << gmml::ConvertT<double>(coords->GetX()) << ", " << gmml::ConvertT<double>(coords->GetY()) << ", " << gmml::ConvertT<double>(coords->GetZ());
    gmml::AddLiteral(atom_uri, Ontology::coordinate, coord_stream.str(), atom_stream);

    // std::string neighbor_resource = "";
    // std::string neighbor_uri = "";
    // AtomVector neighbors = atom->GetNode()->GetNodeNeighbors();
    // for(AtomVector::iterator it = neighbors.begin(); it != neighbors.end(); it++)
    // {
    //     Atom* neighbor = (*it);
    //     neighbor_resource = CreateURIResource(gmml::OntAtom, 0, id_prefix, neighbor->GetId());
    //     neighbor_uri = CreateURI(neighbor_resource);
    //     gmml::AddTriple(atom_uri, Ontology::hasNeighbor, neighbor_uri, atom_stream);
    // }
}

void Assembly::CreateTitle(std::string pdb_resource, std::stringstream& pdb_stream)
{
    pdb_stream << std::endl << "####################################" << std::endl;
    pdb_stream << "#" << std::setw(9) << " " << pdb_resource << " Individuals" << std::setw(9) << " " << "#" << std::endl;
    pdb_stream << "####################################" << std::endl;
}

// void Assembly::gmml::AddTriple(std::string s, std::string p, std::string o, std::stringstream& stream)
// {
//     stream << s << " " << p << " " << o << "." << std::endl;
// }
//
// void Assembly::gmml::AddLiteral(std::string s, std::string p, std::string o, std::stringstream& stream)
// {
//     stream << s << " " << p << " \"" << o << "\"." << std::endl;
// }
//
// void Assembly::gmml::AddDecimal(std::string s, std::string p, float o, std::stringstream& stream)
// {
//   stream << s << " " << p << " \"" << o << "\"^^xsd:decimal." << std::endl;
// }

std::string Assembly::CreateURI(std::string uri_resource)
{
    std::stringstream uri;
    uri << Ontology::ONT_PREFIX << uri_resource;
    return uri.str();
}

std::string Assembly::CreateURIResource(gmml::URIType resource , int number, std::string id_prefix, std::string id )
{
    std::stringstream uri_resource;
    switch(resource)
    {
        case gmml::OntPDB:
            uri_resource << (gmml::Split(this->GetSourceFile().substr(this->GetSourceFile().find_last_of('/') + 1), ".").at(0));
            break;
        case gmml::OntOligosaccharide:
            uri_resource << id_prefix << "oligo" << number;
            break;
        case gmml::OntTerminal:
            uri_resource << id_prefix << "oligo" << number << "_terminal";
            break;
        case gmml::OntSequenceResidue:
            uri_resource << id_prefix << "oligo" << number << "_res" << id;
            break;
        case gmml::OntNote:
            uri_resource << id_prefix << "note" << number;
            break;
        case gmml::OntLinkage:
            uri_resource << id_prefix << "link" << number;
            break;
        case gmml::OntMonosaccharide:
            uri_resource << id_prefix << "mono" << number;
            break;
        case gmml::OntSugarName:
            uri_resource << id_prefix << "mono" << number << "_sn";
            break;
        case gmml::OntAtom:
            replace( id.begin(), id.end(), '?', 'n'); // replace all '?' with 'n'
            gmml::FindReplaceString(id, "\'", "q");
            gmml::FindReplaceString(id, ",", "c");
            //gmml::FindReplaceString(id, "*", "s");
            replace( id.begin(), id.end(), '*', 's'); // replace all '*' with 's'
            uri_resource << id_prefix << id;
            break;
        case gmml::OntResidue:
            replace( id.begin(), id.end(), '?', 'n'); // replace all '?' with 'n'
            replace( id.begin(), id.end(), '*', 's'); // replace all '*' with 's'
            uri_resource << id_prefix << id;
            break;
        case gmml::OntSequenceLinkage: // Added to surpress warning.
            break;
    }
    return uri_resource.str();
}

void Assembly::FormulateCURL(std::string output_file_type, std::string query)
{
//    std::cout << "GENERATED QUERY:" << std::endl;
//    std::cout << query << std::endl;
    std::stringstream curl;
    curl << Ontology::CURL_PREFIX;

    if(output_file_type.compare("csv") == 0)
        curl << Ontology::CSV_OUTPUT_FORMAT;
    else if(output_file_type.compare("xml") == 0)
        curl << Ontology::XML_OUTPUT_FORMAT;
    else if(output_file_type.compare("json") == 0)
        curl << Ontology::JSON_OUTPUT_FORMAT;

    curl << Ontology::DATA_STORE_ADDRESS << Ontology::QUERY_PREFIX << query << Ontology::QUERY_POSTFIX;
    std::string tmp = curl.str();
//    std::cout << std::endl << "RESULTS: " << std::endl;
    const char* cstr = tmp.c_str();
    system(cstr);
//    std::cout << std::endl;
}

std::string Assembly::FormulateCURLGF(std::string output_file_type, std::string query, std::string url)
{
//    std::cout << "GENERATED QUERY:" << std::endl;
//    std::cout << query << std::endl;
    std::stringstream curl;
    curl << Ontology::CURL_PREFIX;

    if(output_file_type.compare("csv") == 0)
        curl << Ontology::CSV_OUTPUT_FORMAT;
    else if(output_file_type.compare("xml") == 0)
        curl << Ontology::XML_OUTPUT_FORMAT;
    else if(output_file_type.compare("json") == 0)
        curl << Ontology::JSON_OUTPUT_FORMAT;

	// Had to change this from DATA_STORE_ADDRESS_GF to the url parameter to allow for Docker to work.
	// TODO Functions like this should not be in GMML C++. They should be at the GEMS scripting level.
    curl << url << Ontology::QUERY_PREFIX << query << Ontology::QUERY_POSTFIX;
    std::string tmp = curl.str();
    //std::cout << std::endl << "RESULTS: " << std::endl;
    // const char* cstr = tmp.c_str();
    //system(cstr);
    //std::cout << std::endl;
    return tmp;
}
