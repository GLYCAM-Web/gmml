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
#include "../../../../includes/CodeUtils/logging.hpp"
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
    // This function will populate the ontology file for the PDB

    int local_debug = -1;

    // All of the stringstreams we will need to write to the ontology file
    std::stringstream pdb_stream, note_stream, oligo_stream, mono_stream, linkage_stream, residue_stream,
        oligo_structure_stream;

    // Add title for every section
    std::string pdb_resource = CreateURIResource(gmml::OntPDB, 0, "", "");
    std::transform(pdb_resource.begin(), pdb_resource.end(), pdb_resource.begin(), ::tolower);
    CreateTitle(pdb_resource, pdb_stream);
    CreateTitle("Notes", note_stream);
    CreateTitle("Oligosaccharides", oligo_stream);
    CreateTitle("Monosaccharides", mono_stream);
    CreateTitle("Linkages", linkage_stream);
    CreateTitle("Residues", residue_stream);
    CreateTitle("Oligosaccharide Structures", oligo_structure_stream);

    ///////////////////////////////
    //    Get the PDB info       //
    ///////////////////////////////
    std::string id_prefix = pdb_resource + "_";
    std::string pdb_uri   = CreateURI(pdb_resource);

    // The file needs to either be formatted
    // s rdf:type Subject.
    // s p o.
    // s p o.
    // Or
    // s
    //   p o;
    //   p o1, o2, o3;
    //   rdf:type Subject.
    // The second is so much easier to read so as I'm refactoring I
    // will be changing to the second format
    // rdf:type will always be the last property as it is always present
    // and will have a period instead of a semicolon

    // Print the PDB info
    pdb_stream << pdb_uri << "\n";

    // When we add more file types, we need to make sure they follow the format
    // of the PDB file printing
    this->input_file_->PrintOntology(pdb_stream);

    // Get and print the oligo, mono, and linkage info
    int link_id = 1;
    std::map<std::string, std::string> mono_to_short_name_map;
    std::map<std::string, std::string> oligo_to_res_uri_map;
    std::vector<std::string> side_or_ring_atoms = std::vector<std::string>();
    std::vector<int> visited_oligos             = std::vector<int>();
    int root_oligo_id                           = 0;
    if (oligos.size() != 0)
    {
        PopulateOligosaccharide(pdb_stream, oligo_stream, oligo_structure_stream, residue_stream, mono_stream,
                                linkage_stream, pdb_uri, id_prefix, link_id, oligos, side_or_ring_atoms, visited_oligos,
                                mono_to_short_name_map, oligo_to_res_uri_map, root_oligo_id);
    }

    // todo - get notes working properly all around; refactor

    if (local_debug > 0)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Getting notes");
    }
    NoteVector notes = this->GetNotes();

    if (notes.size() != 0)
    {
        int note_id = 1;
        // gmml::log(__LINE__, __FILE__,  gmml::INF, "Poulating notes");
        PopulateNotes(pdb_stream, note_stream, pdb_uri, notes, id_prefix, note_id);
        pdb_stream << "\t"
                   << "rdf:type"
                   << "\t" << Ontology::PDB << "."
                   << "\n\n";
        pdb_stream << note_stream.str();
    }
    else
    {
        pdb_stream << "\t"
                   << "rdf:type"
                   << "\t" << Ontology::PDB << "."
                   << "\n\n";
    }

    main_stream << pdb_stream.str() << oligo_stream.str() << mono_stream.str() << linkage_stream.str()
                << residue_stream.str() << oligo_structure_stream.str() << std::endl;
}

void Assembly::PopulateOligosaccharide(std::stringstream& pdb_stream, std::stringstream& oligo_stream,
                                       std::stringstream& oligo_structure_stream, std::stringstream& residue_stream,
                                       std::stringstream& mono_stream, std::stringstream& linkage_stream,
                                       std::string pdb_uri, std::string id_prefix, int& link_id,
                                       OligosaccharideVector oligos, std::vector<std::string>& side_or_ring_atoms,
                                       std::vector<int>& visited_oligos,
                                       std::map<std::string, std::string>& mono_to_short_name_map,
                                       std::map<std::string, std::string>& oligo_to_res_uri_map, int& root_oligo_id)
{
    int local_debug = -1;

    std::string oligo_resource = "", oligo_uri = "";
    // std::string child_oligo_resource = "";
    // std::string child_oligo_uri = "";
    // std::string child_res_resource = "";
    // std::string child_res_uri = "";
    // std::string child_mono_resource = "";
    // std::string child_mono_uri = "";
    // std::string child_res_resource = "";
    // std::string child_res_uri = "";
    // std::string parent_mono_resource = "";
    // std::string parent_mono_uri = "";
    // std::string parent_res_resource = "";
    // std::string parent_res_uri = "";
    // std::string root_oligo_resource = "";
    // std::string root_oligo_uri = "";
    // std::string term_resource = "";
    // std::string term_uri = "";

    int oligoNum = 1;
    for (OligosaccharideVector::iterator it = oligos.begin(); it != oligos.end(); it++, oligoNum++)
    {
        Glycan::Oligosaccharide* oligo = (*it);

        // Add oligo info to pdb stream
        oligo_resource = CreateURIResource(gmml::OntOligosaccharide, oligoNum, id_prefix, "");
        oligo_uri      = CreateURI(oligo_resource);
        pdb_stream << "\t"
                   << ":hasOligo "
                   << "\t" << oligo_uri << ";"
                   << "\n";

        // Add oligo info to oligo stream
        CreateSubtitle(oligo_resource, oligo_stream);
        oligo_stream << oligo_uri << "\n";

        // Different oligo names
        std::string o_name      = oligo->oligosaccharide_name_;
        std::string iupac_name  = oligo->IUPAC_name_;
        std::string author_name = oligo->author_IUPAC_name_;
        if (o_name.compare("") != 0)
        {
            gmml::CreateLiteral(o_name);
            oligo_stream << "\t" << Ontology::oligo_sequence_name << "\t" << o_name << ";"
                         << "\n";
        }
        if (iupac_name.compare("") != 0)
        {
            gmml::CreateLiteral(iupac_name);
            oligo_stream << "\t" << Ontology::oligo_IUPAC_name << "\t\t" << iupac_name << ";"
                         << "\n";
        }
        if (author_name.compare("") != 0)
        {
            gmml::CreateLiteral(author_name);
            oligo_stream << "\t" << Ontology::author_IUPAC_name << "\t" << author_name << ";"
                         << "\n";
        }

        // Residue info
        std::string o_residue_links = oligo->oligosaccharide_residue_linkages_;
        if (o_residue_links.compare("") != 0)
        {
            gmml::CreateLiteral(o_residue_links);
            oligo_stream << "\t" << Ontology::oligo_residue_linkages << "\t" << o_residue_links << ";"
                         << "\n";
        }

        // Chemical modification info
        if (oligo->is_chemically_modified_)
        {
            oligo_stream << "\t" << Ontology::isChemicallyModified << "\t"
                         << "true"
                         << ";"
                         << "\n";
        }
        else
        {
            oligo_stream << "\t" << Ontology::isChemicallyModified << "\t"
                         << "false"
                         << ";"
                         << "\n";
        }

        if (oligo->chemically_modified_terminal_)
        {
            oligo_stream << "\t" << Ontology::hasModifiedTerminal << "\t"
                         << "true"
                         << ";"
                         << "\n";
            // todo - add chemical modification info
        }
        else
        {
            oligo_stream << "\t" << Ontology::hasModifiedTerminal << "\t"
                         << "false"
                         << ";"
                         << "\n";
        }

        // Glycosylation info
        if (oligo->is_attached_to_protein_)
        {
            oligo_stream << "\t" << Ontology::isAttachedToProtein << "\t"
                         << "true"
                         << ";"
                         << "\n";
            std::string glycosylation_type    = oligo->glycosylation_type_;
            std::string glycosylation_residue = oligo->glycosylation_residue_;
            std::string glycosylation_pair    = oligo->glycosylation_pair_;
            gmml::CreateLiteral(glycosylation_type);
            gmml::CreateLiteral(glycosylation_residue);
            gmml::CreateLiteral(glycosylation_pair);

            oligo_stream << "\t" << Ontology::glycosylationType << "\t" << glycosylation_type << ";"
                         << "\n";
            oligo_stream << "\t" << Ontology::glycosylationResidue << "\t" << glycosylation_residue << ";"
                         << "\n";
            oligo_stream << "\t" << Ontology::glycosylationPair << "\t" << glycosylation_pair << ";"
                         << "\n";

            if (oligo->is_N_Glycan_)
            {
                oligo_stream << "\t" << Ontology::isNGlycan << "\t\t"
                             << "true"
                             << ";"
                             << "\n";
            }
            else
            {
                oligo_stream << "\t" << Ontology::isNGlycan << "\t\t"
                             << "false"
                             << ";"
                             << "\n";
            }
            if (oligo->is_O_Glycan_)
            {
                oligo_stream << "\t" << Ontology::isOGlycan << "\t\t"
                             << "true"
                             << ";"
                             << "\n";
            }
            else
            {
                oligo_stream << "\t" << Ontology::isOGlycan << "\t\t"
                             << "false"
                             << ";"
                             << "\n";
            }
            if (oligo->is_C_Glycan_)
            {
                oligo_stream << "\t" << Ontology::isCGlycan << "\t\t"
                             << "true"
                             << ";"
                             << "\n";
            }
            else
            {
                oligo_stream << "\t" << Ontology::isCGlycan << "\t\t"
                             << "false"
                             << ";"
                             << "\n";
            }
            if (oligo->is_S_Glycan_)
            {
                oligo_stream << "\t" << Ontology::isSGlycan << "\t\t"
                             << "true"
                             << ";"
                             << "\n";
            }
            else
            {
                oligo_stream << "\t" << Ontology::isSGlycan << "\t\t"
                             << "false"
                             << ";"
                             << "\n";
            }
        }
        else
        {
            oligo_stream << "\t" << Ontology::isAttachedToProtein << "\t"
                         << "false"
                         << ";"
                         << "\n";
        }

        // B factor and eventually other info/stats
        std::string oligo_b_factor = gmml::ConvertT(oligo->oligosaccharide_b_factor_);
        gmml::CreateDecimal(oligo_b_factor);
        oligo_stream << "\t" << Ontology::hasBFactor << "\t\t" << oligo_b_factor << ";"
                     << "\n";

        // These are reset for each oligosaccharide
        int linkNum = 0;
        int numR    = 0;

        int MonoNum, MonoNeighborNum;
        std::string child_res_resource, child_res_uri, parent_res_resource, parent_res_uri, parent_mono_resource,
            parent_mono_uri, child_mono_resource, child_mono_uri;
        std::vector<MolecularModeling::Residue*> residueVector;
        // Monosaaccharide info
        for (std::vector<Glycan::Monosaccharide*>::iterator it = oligo->mono_nodes_.begin();
             it != oligo->mono_nodes_.end(); it++)
        {
            // This loop should really go in the order of the name, but that indexing
            // is full of bugs.

            Glycan::Monosaccharide* thisMono = *it;
            root_oligo_id                    = oligoNum;
            MonoNum                          = thisMono->mono_id_;
            residueVector.push_back(thisMono->cycle_atoms_[0]->GetResidue());

            parent_mono_resource = CreateURIResource(gmml::OntMonosaccharide, MonoNum, id_prefix, "");
            parent_mono_uri      = CreateURI(parent_mono_resource);
            std::string resID    = std::to_string(MonoNum);
            parent_res_resource  = CreateURIResource(gmml::OntSequenceResidue, root_oligo_id, id_prefix, resID);
            parent_res_uri       = CreateURI(parent_res_resource);
            residue_stream << parent_res_uri << "\n";

            PopulateMonosaccharide(mono_stream, oligo_stream, pdb_stream, oligo_uri, id_prefix, thisMono,
                                   side_or_ring_atoms, pdb_uri, numR);
            CheckDerivativesAndPopulate(oligo_stream, residue_stream, oligo_uri, parent_res_uri, thisMono);

            // Add residue and linkage info
            for (std::vector<std::pair<Glycan::GlycosidicLinkage*, Glycan::Monosaccharide*>>::iterator it =
                     thisMono->mono_neighbors_.begin();
                 it != thisMono->mono_neighbors_.end(); it++)
            {
                Glycan::GlycosidicLinkage* thisLink = (*it).first;
                if (thisLink->non_reducing_mono_ == thisMono)
                {
                    Glycan::Monosaccharide* thisMonoNeighbor = (*it).second;
                    MonoNeighborNum                          = thisMonoNeighbor->mono_id_;
                    // change to IUPAC_index_ when bug is fixed
                    child_mono_resource = CreateURIResource(gmml::OntMonosaccharide, MonoNeighborNum, id_prefix, "");
                    child_mono_uri      = CreateURI(child_mono_resource);
                    std::string neighborResID = std::to_string(MonoNeighborNum);
                    child_res_resource =
                        CreateURIResource(gmml::OntSequenceResidue, root_oligo_id, id_prefix, neighborResID);
                    child_res_uri = CreateURI(child_res_resource);

                    std::stringstream connectionInfo;
                    connectionInfo << ":is" << thisLink->linkage_type_ << "ConnectedTo";

                    residue_stream << "\t" << connectionInfo.str() << "\t" << child_res_uri << ";"
                                   << "\n";
                    residue_stream << "\t" << Ontology::isConnectedTo << "\t\t" << child_res_uri << ";"
                                   << "\n";

                    // Linkage info
                    PopulateLinkage(linkage_stream, oligo_stream, oligo_uri, parent_mono_uri, child_mono_uri, linkNum,
                                    thisLink, thisMono, thisMonoNeighbor);
                    linkNum++;
                }
            }
            if (thisMono->is_root_)
            {
                std::string term_resource = "";
                std::string term_uri      = "";
                term_resource             = CreateURIResource(gmml::OntTerminal, root_oligo_id, id_prefix, "");
                term_uri                  = CreateURI(term_resource);
                // this is the terminal so it's okay to not have the linkage type, as it is in the terminal name (IE
                // 1-OH)
                //  gmml::AddTriple(parent_res_uri, Ontology::isConnectedTo, term_uri, oligo_stream);
                //  gmml::AddTriple(oligo_uri, Ontology::hasTerminal, term_uri, oligo_stream);
                //  gmml::AddTriple(term_uri, Ontology::TYPE, Ontology::Terminal, oligo_stream);
                //  gmml::AddLiteral(term_uri, Ontology::id, oligo->oligosaccharide_terminal_, oligo_stream);
                std::string terminal_name = oligo->oligosaccharide_terminal_;
                gmml::CreateLiteral(terminal_name);
                oligo_stream << "\t" << Ontology::hasTerminal << "\t\t" << terminal_name << ";"
                             << "\n";
                residue_stream << "\t" << Ontology::isConnectedTo << "\t\t" << term_uri << ";"
                               << "\n";
            }
            residue_stream << "\t" << Ontology::TYPE << "\t\t" << Ontology::SequenceResidue << "."
                           << "\n\n";
        }

        MolecularModeling::Assembly subAssembly(residueVector);
        subAssembly.SetModelIndex(0);
        PdbFileSpace::PdbFile* thisPDB = subAssembly.BuildPdbFileStructureFromAssembly();
        std::ostringstream PDBstringstream;
        thisPDB->WriteToStringstream(PDBstringstream);
        gmml::AddLiteral(oligo_uri, "gmmo:PDBfile", PDBstringstream.str(), oligo_structure_stream);

        // Always end with rdf:type
        oligo_stream << "\t"
                     << "rdf:type"
                     << "\t\t" << Ontology::Oligosaccharide << "."
                     << "\n\n";
    }

    // The below code is currently being refactored and replaced by the above code

    //   int MonoNum, MonoNeighborNum;
    //   root_oligo_id = oligoNum;
    //   int tempIndex = 0;
    //   for(std::vector<Glycan::Monosaccharide*>::iterator it = oligo->mono_nodes_.begin(); it !=
    //   oligo->mono_nodes_.end(); it++)
    //   {
    //     Glycan::Monosaccharide* thisMono = *it;
    //     thisMono->oligosaccharide_index_ = tempIndex;
    //     tempIndex++;
    //   }
    //   for(std::vector<Glycan::Monosaccharide*>::iterator it = oligo->mono_nodes_.begin(); it !=
    //   oligo->mono_nodes_.end(); it++)
    //   {
    //     Glycan::Monosaccharide* thisMono = *it;
    //     residueVector.push_back(thisMono->cycle_atoms_[0]->GetResidue());
    //     if(thisMono->is_root_)
    //     {
    //       // MonoNum = thisMono->IUPAC_index_;
    //       // change to IUPAC_index_ when bug is fixed
    //
    //
    //       MonoNum = thisMono->oligosaccharide_index_;
    //
    //
    //       if(local_debug > 0)
    //       {
    //         gmml::log(__LINE__, __FILE__,  gmml::INF, "Residue ID: " +
    //         thisMono->cycle_atoms_[0]->GetResidue()->GetId()); gmml::log(__LINE__, __FILE__,  gmml::INF, "Oligo
    //         Index: " + std::to_string(MonoNum)); gmml::log(__LINE__, __FILE__,  gmml::INF, "IUPAC Index: " +
    //         std::to_string(thisMono->IUPAC_index_));
    //       }
    //       // root_oligo_id = thisMono->mono_id_;
    //       // gmml::log(__LINE__, __FILE__,  gmml::INF, " ");
    //       // PopulateLinkage(linkage_stream, oligo, oligo_uri, id_prefix, link_id, visited_oligos);
    //       // gmml::log(__LINE__, __FILE__,  gmml::INF, "About to populate sequence linkages");
    //       parent_mono_resource = CreateURIResource(gmml::OntMonosaccharide, MonoNum, id_prefix, "");
    //       parent_mono_uri = CreateURI(parent_mono_resource);
    //       std::string resID = std::to_string(MonoNum);
    //       std::string monoSNFG = thisMono->SNFG_name_;
    //       std::string monoShortName = thisMono->sugar_name_.monosaccharide_short_name_;
    //       parent_res_resource = CreateURIResource(gmml::OntSequenceResidue, root_oligo_id, id_prefix, resID);
    //       parent_res_uri = CreateURI(parent_res_resource);
    //       if(local_debug > 0)
    //       {
    //         gmml::log(__LINE__, __FILE__,  gmml::INF, "Mono URI: " + parent_mono_uri);
    //         gmml::log(__LINE__, __FILE__,  gmml::INF, "Res URI: " + parent_res_uri);
    //       }
    //       CheckDerivativesAndPopulate(oligo_stream, residue_stream, monoShortName, oligo_uri, parent_res_uri,
    //       monoSNFG, thisMono);
    //
    //       gmml::AddLiteral(parent_res_uri, Ontology::hasNameIndex, std::to_string(tempIndex), oligo_stream);
    //       gmml::AddLiteral(parent_res_uri, "gmmo:hasIUPACIndex", std::to_string(thisMono->IUPAC_index_),
    //       oligo_stream);
    //   for(std::vector<std::pair<Glycan::GlycosidicLinkage*, Glycan::Monosaccharide*> >::iterator it =
    //   thisMono->mono_neighbors_.begin(); it!=thisMono->mono_neighbors_.end(); it++)
    //   {
    //     Glycan::GlycosidicLinkage* thisLink = (*it).first;
    //     if(thisLink->non_reducing_mono_ == thisMono)
    //     {
    //       Glycan::Monosaccharide* thisMonoNeighbor = (*it).second;
    //
    //           // MonoNeighborNum = thisMonoNeighbor->IUPAC_index_;
    //           //change to IUPAC_index_ when bug is fixed
    //
    //           MonoNeighborNum = thisMonoNeighbor->oligosaccharide_index_;
    //           std::string neighborResID = std::to_string(MonoNeighborNum);
    //           std::string monoSNFG = thisMonoNeighbor->SNFG_name_;
    //           std::string monoShortName = thisMonoNeighbor->sugar_name_.monosaccharide_short_name_;
    //           child_res_resource = CreateURIResource(gmml::OntSequenceResidue, root_oligo_id, id_prefix,
    //           neighborResID); child_res_uri = CreateURI(child_res_resource);
    //           //TODO replace isConnectedTo w/ isx-nLinkedTo (IE is1-3LinkedTo)
    //           std::stringstream connectionInfo;
    //           connectionInfo << "gmmo:is" << (*it).first->linkage_type_ << "ConnectedTo";
    //           gmml::AddTriple(parent_res_uri, connectionInfo.str(), child_res_uri, oligo_stream);
    //           PopulateLinkage(linkage_stream, oligo_uri, parent_res_uri, child_res_uri, linkNum, (*it).first,
    //           thisMono, thisMonoNeighbor); linkNum++;
    //         }
    //       }
    //       // PopulateSequenceLinkage(oligo_stream, oligo, oligo_uri, id_prefix, visited_oligos,
    //       mono_to_short_name_map, oligo_to_res_uri_map, root_oligo_id);
    //       // gmml::log(__LINE__, __FILE__,  gmml::INF, "Done populating sequence linkages");
    //
    //       PopulateMonosaccharide(mono_stream, oligo_stream, oligo_uri, id_prefix, thisMono, side_or_ring_atoms,
    //       pdb_uri); std::string term_resource = ""; std::string term_uri = ""; term_resource =
    //       CreateURIResource(gmml::OntTerminal, root_oligo_id, id_prefix, ""); term_uri = CreateURI(term_resource);
    //       //this is the terminal so it's okay to not have the linkage type, as it is in the terminal name (IE 1-OH)
    //       gmml::AddTriple(parent_res_uri, Ontology::isConnectedTo, term_uri, oligo_stream);
    //       gmml::AddTriple(oligo_uri, Ontology::hasTerminal, term_uri, oligo_stream);
    //       gmml::AddTriple(term_uri, Ontology::TYPE, Ontology::Terminal, oligo_stream);
    //       gmml::AddLiteral(term_uri, Ontology::id, oligo->oligosaccharide_terminal_, oligo_stream);
    //
    //
    //     }
    //     else
    //     {
    //       // root_oligo_id = thisMono->mono_id_;
    //       // gmml::log(__LINE__, __FILE__,  gmml::INF, " ");
    //       // PopulateLinkage(linkage_stream, oligo, oligo_uri, id_prefix, link_id, visited_oligos);
    //       // gmml::log(__LINE__, __FILE__,  gmml::INF, "About to populate sequence linkages");
    //       // MonoNum = thisMono->IUPAC_index_;
    //       // change to IUPAC_index_ when bug is fixed
    //       MonoNum = tempIndex;
    //       //this is so ugly but it needs to work asap
    //       thisMono->oligosaccharide_index_ = MonoNum;
    //       if(local_debug > 0)
    //       {
    //         gmml::log(__LINE__, __FILE__,  gmml::INF, "Residue ID: " +
    //         thisMono->cycle_atoms_[0]->GetResidue()->GetId()); gmml::log(__LINE__, __FILE__,  gmml::INF, "Oligo
    //         Index: " + std::to_string(MonoNum)); gmml::log(__LINE__, __FILE__,  gmml::INF, "IUPAC Index: " +
    //         std::to_string(thisMono->IUPAC_index_));
    //       }
    //       parent_mono_resource = CreateURIResource(gmml::OntMonosaccharide, MonoNum, id_prefix, "");
    //       parent_mono_uri = CreateURI(parent_mono_resource);
    //       std::string resID = std::to_string(MonoNum);
    //       std::string monoSNFG = thisMono->SNFG_name_;
    //       std::string monoShortName = thisMono->sugar_name_.monosaccharide_short_name_;
    //       parent_res_resource = CreateURIResource(gmml::OntSequenceResidue, root_oligo_id, id_prefix, resID);
    //       parent_res_uri = CreateURI(parent_res_resource);
    //
    //       CheckDerivativesAndPopulate(oligo_stream, residue_stream, monoShortName, oligo_uri, parent_res_uri,
    //       monoSNFG, thisMono);
    //
    //       gmml::AddLiteral(parent_res_uri, Ontology::hasNameIndex, std::to_string(tempIndex), oligo_stream);
    //       gmml::AddLiteral(parent_res_uri, "gmmo:hasIUPACIndex", std::to_string(thisMono->IUPAC_index_),
    //       oligo_stream); gmml::AddTriple(parent_res_uri, Ontology::hasMono, parent_mono_uri, oligo_stream);
    //       gmml::AddTriple(parent_res_uri, Ontology::TYPE, Ontology::SequenceResidue, oligo_stream);
    //       for(std::vector<std::pair<Glycan::GlycosidicLinkage*, Glycan::Monosaccharide*> >::iterator it =
    //       thisMono->mono_neighbors_.begin(); it!=thisMono->mono_neighbors_.end(); it++)
    //       {
    //         Glycan::GlycosidicLinkage* thisLink = (*it).first;
    //         if(thisLink->non_reducing_mono_ == thisMono)
    //         {
    //           Glycan::Monosaccharide* thisMonoNeighbor = (*it).second;
    //           MonoNeighborNum = thisMonoNeighbor->IUPAC_index_;//change to IUPAC_index_
    //           std::string neighborResID = std::to_string(MonoNeighborNum);
    //           std::string monoSNFG = thisMonoNeighbor->SNFG_name_;
    //           std::string monoShortName = thisMonoNeighbor->sugar_name_.monosaccharide_short_name_;
    //           child_res_resource = CreateURIResource(gmml::OntSequenceResidue, root_oligo_id, id_prefix,
    //           neighborResID); child_res_uri = CreateURI(child_res_resource); std::stringstream connectionInfo;
    //           connectionInfo << "gmmo:is" << (*it).first->linkage_type_ << "ConnectedTo";
    //           gmml::AddTriple(parent_res_uri, connectionInfo.str(), child_res_uri, oligo_stream);
    //           PopulateLinkage(linkage_stream, oligo_uri, parent_res_uri, child_res_uri, linkNum, (*it).first,
    //           thisMono, thisMonoNeighbor);
    //         }
    //       }
    //
    //       PopulateMonosaccharide(mono_stream, oligo_stream, oligo_uri, id_prefix, thisMono, side_or_ring_atoms,
    //       pdb_uri);
    //     }
    //     tempIndex++;
    //   }
    //
    //
    //
    //   MolecularModeling::Assembly subAssembly(residueVector);
    //   subAssembly.SetModelIndex(0);
    //   PdbFileSpace::PdbFile* thisPDB = subAssembly.BuildPdbFileStructureFromAssembly();
    //   std::ostringstream PDBstringstream;
    //   thisPDB->WriteToStringstream(PDBstringstream);
    //   gmml::AddLiteral(oligo_uri, "gmmo:PDBfile", PDBstringstream.str(), oligo_structure_stream);
    //
    //
    //
    //
    //   /*Keeping for reference for now
    //   THis all needs a good cleaning and documenting
    //   if(oligo->child_oligos_.size() != 0 && (find(visited_oligos.begin(), visited_oligos.end(),
    //   oligo->root_->mono_id_) == visited_oligos.end()))
    //   {
    //     // gmml::log(__LINE__, __FILE__,  gmml::INF, " ");
    //       PopulateLinkage(linkage_stream, oligo, oligo_uri, id_prefix, link_id, visited_oligos);
    //       // gmml::log(__LINE__, __FILE__,  gmml::INF, "About to populate sequence linkages");
    //       PopulateSequenceLinkage(oligo_stream, oligo, oligo_uri, id_prefix, visited_oligos, mono_to_short_name_map,
    //       oligo_to_res_uri_map, root_oligo_id);
    //       // gmml::log(__LINE__, __FILE__,  gmml::INF, "Done populating sequence linkages");
    //   }
    //   else if(oligo->child_oligos_.size() == 0 && o_name.compare("") != 0 &&
    //   oligo->oligosaccharide_terminal_.compare("") != 0)
    //   {
    //     // gmml::log(__LINE__, __FILE__,  gmml::INF, " ");
    //       std::string term_resource = "";
    //       std::string term_uri = "";
    //       term_resource = CreateURIResource(gmml::OntTerminal, root_oligo_id, id_prefix, "");
    //       term_uri = CreateURI(term_resource);
    //       gmml::AddTriple(oligo_uri, Ontology::hasTerminal, term_uri, oligo_stream);
    //       gmml::AddTriple(term_uri, Ontology::TYPE, Ontology::Terminal, oligo_stream);
    //       gmml::AddLiteral(term_uri, Ontology::id, oligo->oligosaccharide_terminal_, oligo_stream);
    //       //std::cout << "Terminalll " << o_term_name << std::endl;
    //       std::stringstream res_id;
    //       std::string res_resource;
    //       std::string res_uri;
    //       res_id << "1";
    //       res_resource = CreateURIResource(gmml::OntSequenceResidue, root_oligo_id, id_prefix, res_id.str());
    //       res_uri = CreateURI(res_resource);
    //       std::string mono_short_name = oligo->root_->sugar_name_.monosaccharide_short_name_;
    //       std::string monoSNFG = oligo->root_->SNFG_name_;
    //       // gmml::log(__LINE__, __FILE__,  gmml::INF, "About to check derivative map");
    //       CheckDerivativesAndPopulate(oligo_stream, residue_stream, mono_short_name, oligo_uri, res_uri, monoSNFG);
    //       gmml::AddTriple(res_uri, Ontology::isConnectedTo, term_uri, oligo_stream);
    //   }
    //   else if(oligo->child_oligos_.size() == 0 && o_name.compare("") != 0)
    //   {
    //     // gmml::log(__LINE__, __FILE__,  gmml::INF, " ");
    //       std::stringstream res_id;
    //       std::string res_resource;
    //       std::string res_uri;
    //       res_id << "1";
    //       res_resource = CreateURIResource(gmml::OntSequenceResidue, root_oligo_id, id_prefix, res_id.str());
    //       res_uri = CreateURI(res_resource);
    //       std::string mono_short_name = oligo->root_->sugar_name_.monosaccharide_short_name_;
    //       std::string monoSNFG = oligo->root_->SNFG_name_;
    //       // gmml::log(__LINE__, __FILE__,  gmml::INF, "About to check derivative map");
    //       CheckDerivativesAndPopulate(oligo_stream, residue_stream, mono_short_name, oligo_uri, res_uri, monoSNFG);
    //   }
    //   // gmml::log(__LINE__, __FILE__,  gmml::INF, " ");
    //   Glycan::Monosaccharide* mono = oligo->root_;
    //   PopulateMonosaccharide(mono_stream, oligo_stream, oligo_uri, id_prefix, mono, side_or_ring_atoms, pdb_uri);
    //
    //   std::vector<Glycan::Oligosaccharide*> child_oligos = oligo->child_oligos_;
    //   PopulateOligosaccharide(pdb_stream, oligo_stream, oligo_stream, mono_stream, linkage_stream, pdb_uri,
    //   id_prefix, link_id, child_oligos, side_or_ring_atoms, visited_oligos, mono_to_short_name_map,
    //   oligo_to_res_uri_map, root_oligo_id);
    //
    //   */
    // }
}

void Assembly::PopulateNotes(std::stringstream& pdb_stream, std::stringstream& note_stream, std::string pdb_uri,
                             NoteVector notes, std::string id_prefix, int note_id)
{
    std::string note_resource = "";
    std::string note_uri      = "";
    for (NoteVector::iterator it = notes.begin(); it != notes.end(); it++)
    {
        int local_debug       = -1;
        Glycan::Note* note    = (*it);
        std::string newPrefix = pdb_uri + "_";
        note_resource         = CreateURIResource(gmml::OntNote, note_id, newPrefix, "");

        // combining these got weird and i need a quick fix
        if (note_resource.find(':') == std::string::npos)
        {
            note_uri = CreateURI(note_resource);
        }
        else
        {
            note_uri = note_resource;
        }

        // gmml::AddTriple(pdb_uri, Ontology::hasNote, note_uri, pdb_stream);

        // //        note_stream << Ontology::ENTITY_COMMENT << note_resource << std::endl;
        // gmml::AddTriple(note_uri, Ontology::TYPE, Ontology::Note, note_stream);
        // //        gmml::AddLiteral(note_uri, Ontology::LABEL, note_resource, note_stream);
        // gmml::AddLiteral(note_uri, Ontology::note_type, note->ConvertGlycanNoteType2String(note->type_),
        // note_stream); gmml::AddLiteral(note_uri, Ontology::note_category,
        // note->ConvertGlycanNoteCat2String(note->category_), note_stream); gmml::AddLiteral(note_uri,
        // Ontology::note_description, note->description_, note_stream);

        if (local_debug > 0)
        {
            gmml::log(__LINE__, __FILE__, gmml::INF, "note_resource: " + note_resource);
            gmml::log(__LINE__, __FILE__, gmml::INF, "note_uri: " + note_uri);
        }

        pdb_stream << "\t" << Ontology::hasNote << "\t\t" << note_uri << ";\n";

        CreateSubtitle(note_resource, note_stream);

        note_stream << "  " << note_uri << "\n";

        std::string noteType        = note->ConvertGlycanNoteType2String(note->type_);
        std::string noteCategory    = note->ConvertGlycanNoteCat2String(note->category_);
        std::string noteDescription = note->description_;
        gmml::CreateLiteral(noteType);
        gmml::CreateLiteral(noteCategory);
        gmml::CreateLiteral(noteDescription);
        note_stream << "\t\t" << Ontology::note_type << "\t" << noteType << ";\n";
        note_stream << "\t\t" << Ontology::note_category << "\t" << noteCategory << ";\n";
        note_stream << "\t\t" << Ontology::note_description << "\t" << noteDescription << ";\n";
        note_stream << "\t\t" << Ontology::TYPE << "\t" << Ontology::Note << ".\n\n";

        note_id++;
    }
}

void Assembly::CheckDerivativesAndPopulate(std::stringstream& oligo_stream, std::stringstream& residue_stream,
                                           std::string oligo_uri, std::string res_uri, Glycan::Monosaccharide* mono)
{
    // gmml::AddTriple(oligo_uri, Ontology::hasSequenceResidue, res_uri, oligo_stream);
    // gmml::AddTriple(res_uri, Ontology::TYPE, Ontology::SequenceResidue, oligo_stream);
    // gmml::AddLiteral(res_uri, Ontology::id, mono->cycle_atoms_[0]->GetResidue()->GetId(), oligo_stream);
    // gmml::AddLiteral(res_uri, "gmmo:residueName", mono->cycle_atoms_[0]->GetResidue()->GetName(), oligo_stream);
    // gmml::AddLiteral(res_uri, Ontology::mono_short_name, mono_short_name, oligo_stream);
    // gmml::AddLiteral(res_uri, Ontology::hasSNFGName, monoSNFG, oligo_stream);
    std::string mono_short_name = mono->sugar_name_.monosaccharide_short_name_;
    std::string monoSNFG        = mono->SNFG_name_;

    oligo_stream << "\t" << Ontology::hasSequenceResidue << "\t" << res_uri << ";\n";

    if (hasDerivative(mono_short_name))
    {
        std::vector<std::string> derivatives;
        getDerivatives(mono_short_name, derivatives);
        for (std::vector<std::string>::iterator t = derivatives.begin(); t != derivatives.end(); ++t)
        {
            // gmml::AddLiteral(res_uri, Ontology::seq_derivative, *t, oligo_stream);
            std::string derivative = *t;
            gmml::CreateLiteral(derivative);
            residue_stream << "\t\t" << Ontology::seq_derivative << "\t" << derivative << ";\n";
        }
    }
}

bool Assembly::hasDerivative(std::string mono_short_name)
{
    std::string::size_type loc1 = mono_short_name.find("[", 0);
    std::string::size_type loc2 = mono_short_name.find("]", 0);
    if (loc1 != std::string::npos && loc2 != std::string::npos && loc1 < loc2)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void Assembly::getDerivatives(std::string& mono_short_name, std::vector<std::string>& derivatives)
{
    std::string::size_type loc1 = mono_short_name.find("[", 0);
    std::string::size_type loc2 = mono_short_name.find("]", 0);
    std::string deriv           = "";
    if (loc1 != std::string::npos && loc2 != std::string::npos && loc1 < loc2)
    {
        deriv = mono_short_name.substr(loc1, loc2 - loc1 + 1);
        gmml::FindReplaceString(mono_short_name, deriv, "");
        deriv       = deriv.substr(1, deriv.length() - 2);
        derivatives = gmml::Split(deriv, ",");
    }
}

void Assembly::PopulateLinkage(std::stringstream& linkage_stream, std::stringstream& oligo_stream,
                               std::string oligo_uri, std::string parent_mono_uri, std::string child_mono_uri,
                               int& linkNum, Glycan::GlycosidicLinkage* thisLinkage, Glycan::Monosaccharide* thisMono,
                               Glycan::Monosaccharide* thisMonoNeighbor)
{
    int local_debug                      = -1;
    std::string linkage_resource         = "";
    std::string linkage_uri              = "";
    std::string child_oligo_resource     = "";
    std::string child_oligo_uri          = "";
    std::string child_atom_resource      = "";
    std::string child_atom_uri           = "";
    std::string glycosidic_atom_resource = "";
    std::string glycosidic_atom_uri      = "";
    std::string parent_atom_resource     = "";
    std::string parent_atom_uri          = "";
    std::stringstream linkage_str;
    std::stringstream glycosidic_linkage_str;

    if ((thisMono == thisLinkage->non_reducing_mono_) || (thisMono == thisLinkage->non_reducing_mono_2_))
    {
        linkage_str << oligo_uri << "_link_" << linkNum;
        linkage_uri = linkage_str.str();

        oligo_stream << "\t" << Ontology::hasGlycosidicLinkage << "\t" << linkage_uri << ";\n";

        // add subtitle to the linkage
        CreateSubtitle(linkage_uri, linkage_stream);

        linkage_stream << linkage_uri << "\n";
        linkage_stream << "\t" << Ontology::hasParentMono << "\t\t" << parent_mono_uri << ";\n";
        linkage_stream << "\t" << Ontology::hasChildMono << "\t\t" << child_mono_uri << ";\n";
        std::string linkage_type           = thisLinkage->linkage_type_;
        std::string hydroxyl_configuration = thisLinkage->hydroxyl_configuration_;
        std::string anomeric_configuration = thisLinkage->anomeric_configuration_;
        std::string linkage_name           = thisLinkage->linkage_name_;
        std::string residue_linkage_name   = thisLinkage->residue_linkage_name_;
        gmml::CreateLiteral(linkage_type);
        gmml::CreateLiteral(hydroxyl_configuration);
        gmml::CreateLiteral(anomeric_configuration);
        gmml::CreateLiteral(linkage_name);
        gmml::CreateLiteral(residue_linkage_name);

        linkage_stream << "\t" << Ontology::linkageType << "\t\t" << linkage_type << ";\n";
        linkage_stream << "\t" << Ontology::orientation << "\t\t" << hydroxyl_configuration << ";\n";
        linkage_stream << "\t" << Ontology::configuration << "\t\t" << anomeric_configuration << ";\n";
        linkage_stream << "\t" << Ontology::glycosidic_linkage << "\t" << linkage_name << ";\n";
        linkage_stream << "\t" << Ontology::residue_linkage << "\t\t" << residue_linkage_name << ";\n";

        std::string totalCHIEnergy = gmml::ConvertT(thisLinkage->total_CHI_Energy_);
        std::string phiCHIEnergy   = gmml::ConvertT(thisLinkage->phi_CHI_Energy_);
        std::string psiCHIEnergy   = gmml::ConvertT(thisLinkage->psi_CHI_Energy_);
        std::string psiAngle       = gmml::ConvertT(thisLinkage->psi_angle_);
        std::string phiAngle       = gmml::ConvertT(thisLinkage->phi_angle_);
        gmml::CreateDecimal(totalCHIEnergy);
        gmml::CreateDecimal(phiCHIEnergy);
        gmml::CreateDecimal(psiCHIEnergy);
        gmml::CreateDecimal(psiAngle);
        gmml::CreateDecimal(phiAngle);

        linkage_stream << "\t" << Ontology::totalCHIEnergy << "\t\t" << totalCHIEnergy << ";\n";
        linkage_stream << "\t" << Ontology::hasPhiAngle << "\t\t" << phiAngle << ";\n";
        linkage_stream << "\t" << Ontology::phiCHIEnergy << "\t\t" << phiCHIEnergy << ";\n";
        linkage_stream << "\t" << Ontology::phiCHIFunction << "\t\t" << thisLinkage->phi_CHI_function_ << ";\n";

        if (thisLinkage->anomeric_anomeric_linkage_)
        {
            std::string phiPrimeAngle = gmml::ConvertT(thisLinkage->phi_prime_angle_);
            gmml::CreateDecimal(phiPrimeAngle);
            // gmml::AddDecimal(linkage_uri, Ontology::hasGlycosidicPhiPrimeAngle, thisLinkage->phi_prime_angle_,
            // linkage_stream);
            linkage_stream << "\t" << Ontology::hasPhiPrimeAngle << "\t" << phiPrimeAngle << ";\n";
        }
        else
        {
            linkage_stream << "\t" << Ontology::hasPsiAngle << "\t\t" << psiAngle << ";\n";
            linkage_stream << "\t" << Ontology::psiCHIEnergy << "\t\t" << psiCHIEnergy << ";\n";
            linkage_stream << "\t" << Ontology::psiCHIFunction << "\t\t" << thisLinkage->psi_CHI_function_ << ";\n";
        }
        if ((thisLinkage->linkage_type_ == "1-6") || (thisLinkage->inverse_linkage_type_ == "1-6") ||
            (thisLinkage->linkage_type_ == "2-6") || (thisLinkage->inverse_linkage_type_ == "2-6"))
        {
            std::string omegaAngle     = gmml::ConvertT(thisLinkage->omega_angle_);
            std::string omegaCHIEnergy = gmml::ConvertT(thisLinkage->omega_CHI_Energy_);
            gmml::CreateDecimal(omegaAngle);
            gmml::CreateDecimal(omegaCHIEnergy);

            linkage_stream << "\t" << Ontology::hasOmegaAngle << "\t\t" << omegaAngle << ";\n";
            linkage_stream << "\t" << Ontology::omegaCHIEnergy << "\t\t" << omegaCHIEnergy << ";\n";
            linkage_stream << "\t" << Ontology::omegaCHIFunction << "\t" << thisLinkage->omega_CHI_function_ << ";\n";
        }

        linkage_stream << "\t" << Ontology::TYPE << "\t\t" << Ontology::Linkage << ".\n\n";
        linkNum++;
    }
    /* Keeping for reference
    for(OligosaccharideVector::iterator it = oligo->child_oligos_.begin(); it != oligo->child_oligos_.end(); it++)
    {
        int index = distance(oligo->child_oligos_.begin(), it);

        Glycan::Oligosaccharide* child_oligo = (*it);
        // OligosaccharideVector::iterator it2 = (std::next(it,1));
        // Glycan::Oligosaccharide* parent_oligo = (*it2);
        //        visited_oligos.push_back(child_oligo->root_->mono_id_);

        linkage_resource = CreateURIResource(gmml::OntLinkage, link_id, id_prefix, "");
        linkage_uri = CreateURI(linkage_resource);
        //        linkage_stream << Ontology::ENTITY_COMMENT << linkage_resource << std::endl;
        gmml::AddTriple(linkage_uri, Ontology::TYPE, Ontology::Linkage, linkage_stream);
        //        gmml::AddLiteral(linkage_uri, Ontology::LABEL, linkage_resource, linkage_stream);
        link_id++;

        gmml::AddTriple(linkage_uri, Ontology::hasParent, oligo_uri, linkage_stream);
        child_oligo_resource = CreateURIResource(gmml::OntOligosaccharide, child_oligo->root_->mono_id_, id_prefix, "");
        child_oligo_uri = CreateURI(child_oligo_resource);
        gmml::AddTriple(linkage_uri, Ontology::hasChild, child_oligo_uri, linkage_stream);

        std::vector<std::string> linkage_tokens = gmml::Split(oligo->child_oligos_linkages_.at(index), "-");
        std::string parent_atom_id = linkage_tokens.at(0);
        std::string glycosidic_atom_id = linkage_tokens.at(1);
        std::string child_atom_id = linkage_tokens.at(2);

        int parent_c_index = ExtractLinkageCarbonIndex(oligo, parent_atom_id);
        int child_c_index = ExtractLinkageCarbonIndex(child_oligo, child_atom_id);

        if(child_c_index != 0 && parent_c_index != 0)
        {
            std::stringstream link_indeces_str;
            link_indeces_str << child_c_index << "-" << parent_c_index;
            gmml::AddLiteral(linkage_uri, Ontology::linkageIndeces, link_indeces_str.str(), linkage_stream);
        }

        child_atom_resource = CreateURIResource(gmml::OntAtom, 0, id_prefix, child_atom_id);
        child_atom_uri = CreateURI(child_atom_resource);
        gmml::AddTriple(linkage_uri, Ontology::hasChildAtomLinkage, child_atom_uri, linkage_stream);

        glycosidic_atom_resource = CreateURIResource(gmml::OntAtom, 0, id_prefix, glycosidic_atom_id);
        glycosidic_atom_uri = CreateURI(glycosidic_atom_resource);
        gmml::AddTriple(linkage_uri, Ontology::hasGlycosidicLinkage, glycosidic_atom_uri, linkage_stream);


        gmml::AddLiteral(linkage_uri, Ontology::hasChildMono,
    child_oligo->root_->sugar_name_.monosaccharide_short_name_, linkage_stream); gmml::AddLiteral(linkage_uri,
    Ontology::hasParentMono, oligo->root_->sugar_name_.monosaccharide_short_name_, linkage_stream); double
    glycosidic_phi_angle = CalculatePhiAngle(child_oligo, parent_atom_id, child_atom_id, glycosidic_atom_id);
        glycosidic_phi_angle = gmml::ConvertRadian2Degree(glycosidic_phi_angle);
        if(local_debug > 0)
        {
          gmml::log(__LINE__, __FILE__, gmml::INF, "Phi");
          gmml::log(__LINE__, __FILE__, gmml::INF, std::to_string(glycosidic_phi_angle));
        }
        gmml::AddTriple(linkage_uri, Ontology::hasGlycosidicPhiAngle, std::to_string(glycosidic_phi_angle),
    linkage_stream);


        double glycosidic_psi_angle = CalculatePsiAngle(child_oligo, parent_atom_id, child_atom_id, glycosidic_atom_id);
        glycosidic_psi_angle = gmml::ConvertRadian2Degree(glycosidic_psi_angle);
        if(local_debug > 0)
        {
          gmml::log(__LINE__, __FILE__, gmml::INF, "Psi");
          gmml::log(__LINE__, __FILE__, gmml::INF, std::to_string(glycosidic_psi_angle));
        }
        gmml::AddTriple(linkage_uri, Ontology::hasGlycosidicPsiAngle, std::to_string(glycosidic_psi_angle),
    linkage_stream);

        if (parent_c_index == 6)
        {
          gmml::log(__LINE__, __FILE__,  gmml::INF, "About to calculate Omega Angle");
          double glycosidic_omega_angle = CalculateOmegaAngle(oligo, parent_atom_id, glycosidic_atom_id);
          glycosidic_omega_angle = gmml::ConvertRadian2Degree(glycosidic_omega_angle);
          if(local_debug > 0)
          {
            gmml::log(__LINE__, __FILE__, gmml::INF, "Omega");
            gmml::log(__LINE__, __FILE__, gmml::INF, std::to_string(glycosidic_omega_angle));
          }
          gmml::AddTriple(linkage_uri, Ontology::hasGlycosidicOmegaAngle, std::to_string(glycosidic_omega_angle),
    linkage_stream);

        }
        // gmml::log(__LINE__, __FILE__,  gmml::INF, "Done with angles");
        parent_atom_resource = CreateURIResource(gmml::OntAtom, 0, id_prefix, parent_atom_id);
        parent_atom_uri = CreateURI(parent_atom_resource);
        // gmml::log(__LINE__, __FILE__,  gmml::INF, "About to add parent linkage");
        gmml::AddTriple(linkage_uri, Ontology::hasParentAtomLinkage, parent_atom_uri, linkage_stream);

        std::vector<std::string> child_atom_id_tokens = gmml::Split(child_atom_id, "_");
        if(child_atom_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
            linkage_str << child_atom_id_tokens.at(2) << "(" << child_atom_id_tokens.at(4) << ")" <<
    child_atom_id_tokens.at(0); else linkage_str << child_atom_id_tokens.at(2) << "(" << child_atom_id_tokens.at(4) <<
    "_" << child_atom_id_tokens.at(3) << ")" << child_atom_id_tokens.at(0);

        std::vector<std::string> parent_atom_id_tokens = gmml::Split(parent_atom_id, "_");
        if(parent_atom_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
            linkage_str << "-" << parent_atom_id_tokens.at(2) << "(" << parent_atom_id_tokens.at(4) << ")"  <<
    parent_atom_id_tokens.at(0); else linkage_str << "-" << parent_atom_id_tokens.at(2) << "(" <<
    parent_atom_id_tokens.at(4) <<  "_" << parent_atom_id_tokens.at(3) << ")"  << parent_atom_id_tokens.at(0);

        //        gmml::AddLiteral(linkage_uri, Ontology::linkage_str, linkage_str.str(), linkage_stream);

        std::vector<std::string> glycosidic_atom_id_tokens = gmml::Split(glycosidic_atom_id, "_");
        if(glycosidic_atom_id_tokens.at(3).at(0) == gmml::BLANK_SPACE)
            glycosidic_linkage_str << glycosidic_atom_id_tokens.at(2) << "(" << glycosidic_atom_id_tokens.at(4) << ")"
    << glycosidic_atom_id_tokens.at(0); else glycosidic_linkage_str << glycosidic_atom_id_tokens.at(2) << "(" <<
    glycosidic_atom_id_tokens.at(4) << "_" << glycosidic_atom_id_tokens.at(3)
                                   << ")"  << glycosidic_atom_id_tokens.at(0);
        // gmml::log(__LINE__, __FILE__,  gmml::INF, "About to add linkage info");
        gmml::AddLiteral(linkage_uri, Ontology::glycosidic_linkage, glycosidic_linkage_str.str(), linkage_stream);
    }
    gmml::log(__LINE__, __FILE__,  gmml::INF, "Done populating linkages");
    */
}

void Assembly::PopulateSequenceLinkage(std::stringstream& oligo_stream, Glycan::Oligosaccharide* oligo,
                                       std::string oligo_uri, std::string id_prefix, std::vector<int>& visited_oligos,
                                       std::map<std::string, std::string>& mono_to_short_name_map,
                                       std::map<std::string, std::string>& oligo_to_res_uri_map, int& root_oligo_id)
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
    // for(std::vector<Glycan::Monosaccharide*>::reverse_iterator rit = oligo->mono_nodes_.rbegin(); rit !=
    // oligo->mono_nodes.rend(); rit++)
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
    //   CheckDerivativesAndPopulate(oligo_stream, residue_stream, monoShortName, root_oligo_uri, parent_res_uri,
    //   monoSNFG); for(std::vector<std::pair<Glycan::GlycosidicLinkage*, Glycan::Monosaccharide*> >::iterator it =
    //   thisMono->mono_neighbors_.begin(); it!=thisMono->mono_neighbors_.end(); it++)
    //   {
    //     thisMonoNeighbor = *it.second;
    //     int MonoNeighborNum = oligo->mono_nodes_.size() - thisMonoNeighbor->oligosaccharide_index_;
    //     std::string resID = std::to_string(MonoNeighborNum);
    //     std::string monoSNFG = thisMonoNeighbor->SNFG_name_;
    //     std::string monoShortName = thisMonoNeighbor->sugar_name_.monosaccharide_short_name_;
    //     parent_res_resource = CreateURIResource(gmml::OntSequenceResidue, root_oligo_id, id_prefix, resID);
    //     parent_res_uri = CreateURI(parent_res_resource);
    //     gmml::AddTriple(child_res_uri, Ontology::isConnectedTo, parent_res_uri, oligo_stream);
    //   }
    //
    // }
    // for(OligosaccharideVector::iterator it = oligo->child_oligos_.begin(); it != oligo->child_oligos_.end(); it++)
    // {
    //     int index = distance(oligo->child_oligos_.begin(), it);
    //
    //     Glycan::Oligosaccharide* child_oligo = (*it);
    //
    //     child_oligo_resource = CreateURIResource(gmml::OntOligosaccharide, child_oligo->root_->mono_id_, id_prefix,
    //     ""); child_oligo_uri = CreateURI(child_oligo_resource);
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
    //         CheckDerivativesAndPopulate(oligo_stream, residue_stream, mono_to_short_name_map[child_oligo_uri],
    //         root_oligo_uri, child_res_uri, monoSNFG);
    //
    //         //gmml::AddTriple(root_oligo_uri, Ontology::hasSequenceResidue, child_res_uri, oligo_stream);
    //         //gmml::AddTriple(child_res_uri, Ontology::TYPE, Ontology::SequenceResidue, oligo_stream);
    //         //gmml::AddLiteral(child_res_uri, Ontology::id, mono_to_short_name_map[child_oligo_uri], oligo_stream);
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
    //         CheckDerivativesAndPopulate(oligo_stream, mono_to_short_name_map[oligo_uri], root_oligo_uri,
    //         child_res_uri, monoSNFG);
    //
    //         //gmml::AddTriple(root_oligo_uri, Ontology::hasSequenceResidue, child_res_uri, oligo_stream);
    //         //gmml::AddTriple(child_res_uri, Ontology::TYPE, Ontology::SequenceResidue, oligo_stream);
    //         //gmml::AddLiteral(child_res_uri, Ontology::id, mono_to_short_name_map[oligo_uri], oligo_stream);
    //         oligo_to_res_uri_map[oligo_uri] = child_res_uri;
    //     }
    //
    //     std::string child_uri = oligo_to_res_uri_map[child_oligo_uri];
    //     std::string parent_uri = oligo_to_res_uri_map[oligo_uri];
    //     gmml::AddTriple(child_uri, Ontology::isConnectedTo, parent_uri, oligo_stream);
    //
    //     if(o_name.compare("") != 0 && o_term_name.compare("") != 0 && !terminalAdded)
    //     {
    //         gmml::AddTriple(parent_uri, Ontology::isConnectedTo, term_uri, oligo_stream);
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
    //         gmml::AddLiteral(oligo_to_res_uri_map[child_oligo_uri], Ontology::sequence_linkage,
    //         link_indeces_str.str(), oligo_stream);
    //     }
    // }
    // if(o_name.compare("") != 0 && o_term_name.compare("") != 0)
    // {
    //     term_resource = CreateURIResource(gmml::OntTerminal, root_oligo_id, id_prefix, "");
    //     term_uri = CreateURI(term_resource);
    //     gmml::AddTriple(root_oligo_uri, Ontology::hasTerminal, term_uri, oligo_stream);
    //     gmml::AddTriple(term_uri, Ontology::TYPE, Ontology::Terminal, oligo_stream);
    //     gmml::AddLiteral(term_uri, Ontology::id, o_term_name, oligo_stream);
    //     //std::cout << "Terminalll " << o_term_name << std::endl;
    // }
}

int Assembly::ExtractLinkageCarbonIndex(Glycan::Oligosaccharide* oligo, std::string linkage_carbon_id)
{
    int c_index                                = 0;
    std::vector<std::string> cycle_atom_tokens = gmml::Split(oligo->root_->cycle_atoms_str_, "-");

    if (oligo->root_->side_atoms_.at(0).at(0) != NULL)
    {
        c_index++;
        Atom* anomeric_side_carbon = oligo->root_->side_atoms_.at(0).at(0);
        if (anomeric_side_carbon->GetId().compare(linkage_carbon_id) == 0)
        {
            return c_index;
        }
    }
    if (!cycle_atom_tokens.empty())
    {
        for (unsigned int i = 0; i < cycle_atom_tokens.size() - 1;
             i++) /// cycle_atom_tokens.size() - 1 > because the ring oxygen is not considered
        {
            c_index++;
            if (cycle_atom_tokens.at(i).compare(linkage_carbon_id) == 0)
            {
                return c_index;
            }
        }
    }
    AtomVector side_atoms_of_last_ring_carbon = oligo->root_->side_atoms_.at(oligo->root_->side_atoms_.size() - 1);
    if (!side_atoms_of_last_ring_carbon.empty())
    {
        for (AtomVector::iterator it1 = side_atoms_of_last_ring_carbon.begin();
             it1 != side_atoms_of_last_ring_carbon.end(); it1++)
        {
            Atom* side_atom = (*it1);
            c_index++;
            if (side_atom != NULL)
            {
                if (side_atom->GetId().compare(linkage_carbon_id) == 0)
                {
                    return c_index;
                }
            }
        }
    }
    return c_index;
}

void Assembly::PopulateMonosaccharide(std::stringstream& mono_stream, std::stringstream& oligo_stream,
                                      std::stringstream& pdb_stream, std::string oligo_uri, std::string id_prefix,
                                      Glycan::Monosaccharide* mono, std::vector<std::string>& side_or_ring_atoms,
                                      std::string pdb_uri, int numR)
{
    int local_debug = -1;

    std::string mono_resource = CreateURIResource(gmml::OntMonosaccharide, mono->mono_id_, id_prefix, "");
    std::string mono_uri      = CreateURI(mono_resource);

    pdb_stream << "\t" << Ontology::hasMono << "\t" << mono_uri << ";" << std::endl;

    oligo_stream << "\t" << Ontology::hasMono << "\t\t" << mono_uri << ";"
                 << "\n";

    CreateSubtitle(mono_resource, mono_stream);
    mono_stream << mono_uri << "\n";
    mono_stream << "\t" << Ontology::hasOligoParent << "\t\t" << oligo_uri << ";"
                << "\n"; // might be redundant

    std::string mono_index            = gmml::ConvertT<int>(mono->mono_id_);
    std::string IUPAC_index           = gmml::ConvertT<int>(mono->IUPAC_index_);
    std::string oligosaccharide_index = gmml::ConvertT<int>(mono->oligosaccharide_index_);
    gmml::CreateInteger(mono_index);
    gmml::CreateInteger(IUPAC_index);
    gmml::CreateInteger(oligosaccharide_index);
    mono_stream << "\t" << Ontology::hasIndex << "\t\t" << mono_index << ";"
                << "\n";
    mono_stream << "\t" << Ontology::hasIUPACIndex << "\t\t" << IUPAC_index << ";"
                << "\n";
    mono_stream << "\t" << Ontology::hasNameIndex << "\t\t" << oligosaccharide_index << ";"
                << "\n";

    // Chemical modification info
    std::stringstream mod_stream;
    std::string subtitle = mono_resource + " R Groups";
    CreateSubtitle(subtitle, mod_stream);
    if (mono->on_R_ > 0)
    {
        mono_stream << "\t" << Ontology::isChemicallyModified << "\t"
                    << "true"
                    << ";\n";

        for (std::vector<std::pair<std::string, std::string>>::iterator derivative = mono->unknown_derivatives_.begin();
             derivative != mono->unknown_derivatives_.end(); derivative++)
        {

            std::string Rgroup = (*derivative).second;
            if (Rgroup != "" &&
                mono->sugar_name_.monosaccharide_short_name_ != mono->cycle_atoms_[0]->GetResidue()->GetName())
            {
                numR++;
                // std::stringstream RgroupStream, RnumStream;
                // RgroupStream << oligo_uri << "_R" << numR;
                // RnumStream << ":has" << "R" << numR;
                // gmml::AddTriple(oligo_uri, RnumStream.str(), RgroupStream.str(), oligo_stream);
                // gmml::AddLiteral(RgroupStream.str(), Ontology::hasFormula,(*derivative).second, oligo_stream);

                // todo - I'm sure there is a better way for the database to see these so the
                // SPARQL queries can be more efficient

                oligo_stream << "\t"
                             << ":hasR" << numR << "\t\t\t";
                oligo_stream << oligo_uri << "_R" << numR << ";\n";

                mono_stream << "\t"
                            << ":hasR" << numR << "\t\t\t";
                mono_stream << oligo_uri << "_R" << numR << ";\n";

                gmml::CreateLiteral(Rgroup);
                mod_stream << "  " << oligo_uri << "_R" << numR << "\n";
                mod_stream << "\t\t" << Ontology::hasFormula << "\t\t" << Rgroup << ";\n";
                mod_stream << "\t\t" << Ontology::hasIndex << "\t\t" << numR << ";\n";
                mod_stream << "\t\t" << Ontology::TYPE << "\t\t" << Ontology::ChemicalModification << ".\n\n";
            }
        }
    }

    std::string residueName = mono->residue_name_;
    gmml::CreateLiteral(residueName);
    mono_stream << "\t" << Ontology::hasResidueName << "\t\t" << residueName << ";"
                << "\n";

    std::string residueId = mono->cycle_atoms_[0]->GetResidue()->GetId();
    gmml::CreateLiteral(residueId);
    mono_stream << "\t" << Ontology::id << "\t\t" << residueId << ";"
                << "\n";

    if (mono->cycle_atoms_[0]->GetResidue()->CheckIfNucleicAcid())
    {
        mono_stream << "\t" << Ontology::isNucleotide << "\t\t"
                    << "true"
                    << ";"
                    << "\n";
    }
    else
    {
        mono_stream << "\t" << Ontology::isNucleotide << "\t\t"
                    << "false"
                    << ";"
                    << "\n";
    }

    if (mono->cycle_atoms_[0]->GetResidue()->CheckIfSaccharide())
    {
        mono_stream << "\t" << Ontology::isSaccharide << "\t\t"
                    << "true"
                    << ";"
                    << "\n";
    }
    else
    {
        mono_stream << "\t" << Ontology::isSaccharide << "\t\t"
                    << "false"
                    << ";"
                    << "\n";
    }

    // Quick patch to make it easier to clean the data; we get lots of false positives still.
    // These should be checked and assigned in the monosaccharide class.
    if (mono->anomeric_carbon_pointer_ != NULL)
    {
        MolecularModeling::AtomVector neighbors = mono->anomeric_carbon_pointer_->GetNode()->GetNodeNeighbors();
        for (MolecularModeling::AtomVector::iterator it = neighbors.begin(); it != neighbors.end(); it++)
        {
            MolecularModeling::Atom* neighbor = (*it);
            if ((neighbor->GetIsRing()) && (neighbor->GetElementSymbol().compare("C") != 0))
            {
                if (local_debug > 0)
                {
                    gmml::log(__LINE__, __FILE__, gmml::INF,
                              "Found non-carbon ring atom: " + neighbor->GetElementSymbol());
                }

                if (neighbor->GetElementSymbol().compare("O") == 0)
                {
                    mono_stream << "\t" << Ontology::hasRingO << "\t\t"
                                << "true"
                                << ";"
                                << "\n";
                }
                else if (neighbor->GetElementSymbol().compare("N") == 0)
                {
                    mono_stream << "\t" << Ontology::hasRingN << "\t\t"
                                << "true"
                                << ";"
                                << "\n";
                }
                else
                {
                    mono_stream << "\t" << Ontology::hasRingO << "\t\t"
                                << "false"
                                << ";"
                                << "\n";
                    mono_stream << "\t" << Ontology::hasRingN << "\t\t"
                                << "false"
                                << ";"
                                << "\n";
                }
            }
            else if (!(neighbor->GetIsRing()))
            {
                if (local_debug > 0)
                {
                    gmml::log(__LINE__, __FILE__, gmml::INF, "Found non-ring atom: " + neighbor->GetElementSymbol());
                }

                std::string neighborElement = neighbor->GetElementSymbol();
                gmml::CreateLiteral(neighborElement);
                mono_stream << "\t" << Ontology::anomericNeighborElement << "\t" << neighborElement << ";"
                            << "\n";

                if (neighborElement.compare("O") == 0)
                {
                    mono_stream << "\t" << Ontology::hasNonRingO << "\t\t"
                                << "true"
                                << ";"
                                << "\n";
                }
                else if (neighborElement.compare("N") == 0)
                {
                    mono_stream << "\t" << Ontology::hasNonRingN << "\t\t"
                                << "true"
                                << ";"
                                << "\n";
                }
                else
                {
                    mono_stream << "\t" << Ontology::hasNonRingN << "\t\t"
                                << "false"
                                << ";"
                                << "\n";
                    mono_stream << "\t" << Ontology::hasNonRingO << "\t\t"
                                << "false"
                                << ";"
                                << "\n";
                }
            }
        }
    }
    else
    {
        mono_stream << "\t" << Ontology::anomericProperlyAssigned << "\t\t"
                    << "false"
                    << ";"
                    << "\n";
    }

    if (mono->bfmp_ring_conformation_.compare("") != 0)
    {
        std::string bfmp = mono->bfmp_ring_conformation_;

        std::size_t index;
        // This is to only output two decimal places for BFMP values which have a decimal value. It also only returns
        // the 1st bfmp value if there are multiple.
        if ((index = mono->bfmp_ring_conformation_.find(".")) != std::string::npos)
        {
            bfmp = bfmp.substr(0, index + 3) + ")";
        }
        gmml::FindReplaceString(bfmp, "\t", " ");
        gmml::CreateLiteral(bfmp);
        mono_stream << "\t" << Ontology::BFMP << "\t\t\t" << bfmp << ";"
                    << "\n";

        bfmp = mono->bfmp_ring_conformation_;
        gmml::FindReplaceString(bfmp, "\t\t", "\t");
        gmml::FindReplaceString(bfmp, "\t", ", ");
        if (bfmp[0] == ',')
        {
            // bfmp from 1 to the last ')'
            bfmp = bfmp.substr(1, bfmp.find_last_of(")"));
        }
        else
        {
            // bfmp from 0 to the last ')'
            bfmp = bfmp.substr(0, bfmp.find_last_of(")"));
        }
        if (bfmp[0] == ' ')
        {
            bfmp = bfmp.substr(1);
        }
        gmml::CreateLiteral(bfmp);
        mono_stream << "\t" << Ontology::fullBFMP << "\t\t" << bfmp << ";"
                    << "\n";
    }

    if (mono->b_factor_ > 0)
    {
        std::string b_factor = std::to_string(mono->b_factor_);
        gmml::CreateDecimal(b_factor);
        mono_stream << "\t" << Ontology::hasBFactor << "\t\t" << b_factor << ";"
                    << "\n";
    }

    std::string chemical_code = mono->sugar_name_.chemical_code_string_;
    gmml::CreateLiteral(chemical_code);
    mono_stream << "\t" << Ontology::stereochemistry_chemical_code << "\t\t" << chemical_code << ";"
                << "\n";

    std::string stereoName = mono->sugar_name_.monosaccharide_stereochemistry_name_;
    gmml::CreateLiteral(stereoName);
    mono_stream << "\t" << Ontology::mono_stereo_name << "\t\t" << stereoName << ";"
                << "\n";

    std::string stereoShortName = mono->sugar_name_.monosaccharide_stereochemistry_short_name_;
    gmml::CreateLiteral(stereoShortName);
    mono_stream << "\t" << Ontology::mono_stereo_short_name << "\t" << stereoShortName << ";"
                << "\n";

    std::string monoName = mono->sugar_name_.monosaccharide_name_;
    gmml::CreateLiteral(monoName);
    mono_stream << "\t" << Ontology::mono_name << "\t\t" << monoName << ";"
                << "\n";

    std::string monoShortName = mono->sugar_name_.monosaccharide_short_name_;
    gmml::CreateLiteral(monoShortName);
    mono_stream << "\t" << Ontology::mono_short_name << "\t\t" << monoShortName << ";"
                << "\n";

    std::string isomer = mono->sugar_name_.isomer_;
    gmml::CreateLiteral(isomer);
    mono_stream << "\t" << Ontology::isomer << "\t\t\t" << isomer << ";"
                << "\n";

    if (mono->sugar_name_.configuration_.compare("a") == 0)
    {
        std::string configuration = "alpha";
        gmml::CreateLiteral(configuration);
        mono_stream << "\t" << Ontology::configuration << "\t\t" << configuration << ";"
                    << "\n";
    }
    else if (mono->sugar_name_.configuration_.compare("b") == 0)
    {
        // gmml::AddLiteral(sugar_name_uri, Ontology::configuration, "beta", sugar_name_stream);
        std::string configuration = "beta";
        gmml::CreateLiteral(configuration);
        mono_stream << "\t" << Ontology::configuration << "\t\t" << configuration << ";"
                    << "\n";
    }

    if (mono->sugar_name_.ring_type_.compare("") != 0)
    {
        std::string ringType = mono->sugar_name_.ring_type_;
        gmml::CreateLiteral(ringType);
        mono_stream << "\t" << Ontology::ring_type << "\t\t" << ringType << ";"
                    << "\n";
    }

    if (mono->SNFG_name_ != "")
    {
        std::string SNFGName = mono->SNFG_name_;
        gmml::CreateLiteral(SNFGName);
        mono_stream << "\t" << Ontology::hasSNFGName << "\t\t" << SNFGName << ";"
                    << "\n";
    }

    if (mono->author_SNFG_name_ != "")
    {
        std::string authorSNFGName = mono->author_SNFG_name_;
        gmml::CreateLiteral(authorSNFGName);
        mono_stream << "\t" << Ontology::hasAuthorSNFGName << "\t\t" << authorSNFGName << ";"
                    << "\n";
    }

    NoteVector notes = mono->mono_notes_;
    std::stringstream monoNotes;
    if (notes.size() != 0)
    {
        int note_id = 1;
        // gmml::log(__LINE__, __FILE__,  gmml::INF, "Poulating notes");
        // std::string id_prefix = mono_uri.substr(5,mono_uri.length()) + "_";
        PopulateNotes(mono_stream, monoNotes, mono_uri, notes, id_prefix, note_id);
    }

    // Always end with rdf:type
    mono_stream << "\t" << Ontology::TYPE << "\t\t" << Ontology::Monosaccharide << "."
                << "\n\n";

    // Adding modifications and notes if they exist
    if (mono->on_R_ > 0)
    {
        mono_stream << mod_stream.str();
    }
    if (notes.size() != 0)
    {
        mono_stream << monoNotes.str();
    }

    // The below code is being replaced by the code above

    // gmml::AddLiteral(mono_uri, Ontology::id, mono->cycle_atoms_[0]->GetResidue()->GetId(), mono_stream);
    // gmml::AddLiteral(mono_uri, "gmmo:residueName", mono->cycle_atoms_[0]->GetResidue()->GetName(), mono_stream);
    // gmml::AddTriple(mono_uri, Ontology::hasOligoParent, oligo_uri, mono_stream);
    // int Index = mono->IUPAC_index_; //change to IUPAC_index_
    // gmml::AddLiteral(mono_uri, Ontology::hasIndex, std::to_string(Index), mono_stream);
    // gmml::AddLiteral(mono_uri, Ontology::hasNameIndex, std::to_string(mono->oligosaccharide_index_), mono_stream);
    //
    // ////////////////////////////////////////////////////////////////////
    // // Write bool for is Nucleotide, Saccharide, etc. to the ontology //
    // ////////////////////////////////////////////////////////////////////
    //
    // // if(checkIfNucleotide(mono))
    // if(mono->cycle_atoms_[0]->GetResidue()->CheckIfNucleicAcid())
    // {
    //   gmml::AddTriple(mono_uri, Ontology::isNucleotide, "true", mono_stream);
    // }
    // else
    // {
    //   gmml::AddTriple(mono_uri, Ontology::isNucleotide, "false", mono_stream);
    // }
    //
    // if(mono->cycle_atoms_[0]->GetResidue()->CheckIfSaccharide())
    // {
    //   gmml::AddTriple(mono_uri, Ontology::isSaccharide, "true", mono_stream);
    // }
    // else
    // {
    //   gmml::AddTriple(mono_uri, Ontology::isSaccharide, "false", mono_stream);
    // }
    // // Quick patch to make it easier to clean the data; we get lots of false positives still
    // if(mono->anomeric_carbon_pointer_ != NULL)
    // {
    //   MolecularModeling::AtomVector neighbors =      mono->anomeric_carbon_pointer_->GetNode()->GetNodeNeighbors();
    //   for(MolecularModeling::AtomVector::iterator it = neighbors.begin(); it != neighbors.end(); it++)
    //   {
    //     MolecularModeling::Atom* neighbor = (*it);
    //     if((neighbor->GetIsRing()) && (neighbor->GetElementSymbol().compare("C") != 0))
    //     {
    //       if(local_debug > 0)
    //       {
    //         gmml::log(__LINE__, __FILE__,  gmml::INF, "Found non-carbon ring atom: " + neighbor->GetElementSymbol());
    //       }
    //       if(neighbor->GetElementSymbol().compare("O") == 0)
    //       {
    //         gmml::AddTriple(mono_uri, "gmmo:hasRingO", "true", mono_stream);
    //       }
    //       else
    //       {
    //         gmml::AddTriple(mono_uri, "gmmo:hasRingO", "false", mono_stream);
    //       }
    //       if(neighbor->GetElementSymbol().compare("N") == 0)
    //       {
    //         gmml::AddTriple(mono_uri, "gmmo:hasRingN", "true", mono_stream);
    //       }
    //       else
    //       {
    //         gmml::AddTriple(mono_uri, "gmmo:hasRingN", "false", mono_stream);
    //       }
    //     }
    //     else if(!(neighbor->GetIsRing()))
    //     {
    //       if(local_debug > 0)
    //       {
    //         gmml::log(__LINE__, __FILE__,  gmml::INF, "Found non-ring atom: " + neighbor->GetElementSymbol());
    //       }
    //       if(neighbor->GetElementSymbol().compare("O") == 0)
    //       {
    //         gmml::AddTriple(mono_uri, "gmmo:hasNonRingO", "true", mono_stream);
    //       }
    //       else if(neighbor->GetElementSymbol().compare("N") == 0)
    //       {
    //         gmml::AddTriple(mono_uri, "gmmo:hasNonRingN", "true", mono_stream);
    //       }
    //       else
    //       {
    //         gmml::AddLiteral(mono_uri, "gmmo:anomericNeighborElement", neighbor->GetElementSymbol(), mono_stream);
    //       }
    //     }
    //   }
    // }
    // else
    // {
    //   gmml::AddTriple(mono_uri, "gmmo:anomericProperlyAssigned", "false", mono_stream);
    // }
    //
    // // gmml::AddLiteral(mono_uri, Ontology::hasSNFGName, mono->SNFG_name_, mono_stream);
    // //    gmml::AddLiteral(mono_uri, Ontology::LABEL, mono_resource, mono_stream);
    //
    // AtomVector ring_atoms = mono->cycle_atoms_;
    // object.str(std::string());
    // int ring_index = 1;
    // std::stringstream ring_atom_stream;
    // // if(ring_atoms.size() > 0)
    // // {
    // //   for(AtomVector::iterator it = ring_atoms.begin(); it != ring_atoms.end(); it++)
    // //   {
    // //       Atom* ring_atom = (*it);
    // //       if (ring_atom != NULL)
    // //       {
    // //         ring_resource = CreateURIResource(gmml::OntAtom, 0, id_prefix, ring_atom->GetId());
    // //         ring_uri = CreateURI(ring_resource);
    // //         gmml::AddTriple(mono_uri, Ontology::hasRingAtom, ring_uri, mono_stream);
    // //         //None of this is used and it clutters the file (and takes of a lot of disk space)
    // //         //PopulateRingAtom(ring_atom_stream, id_prefix, ring_uri, ring_resource, ring_index, ring_atom, mono,
    // side_or_ring_atoms);
    // //         ring_index++;
    // //         if(it == ring_atoms.end() - 1)
    // //             object << ring_resource;
    // //         else
    // //             object << ring_resource << "-";
    // //       }
    // //   }
    // // }
    // // gmml::AddLiteral(mono_uri, Ontology::ring_atoms, object.str(), mono_stream);
    // // object.str(std::string());
    // // object << mono->anomeric_status_ << " " << CreateURIResource(gmml::OntAtom, 0, id_prefix,
    // mono->cycle_atoms_.at(0)->GetId());
    // // gmml::AddLiteral(mono_uri, Ontology::anomeric_status, object.str(), mono_stream);
    // gmml::AddLiteral(mono_uri, Ontology::stereochemistry_chemical_code, mono->sugar_name_.chemical_code_string_,
    // mono_stream); if(mono->bfmp_ring_conformation_.compare("") != 0) {
    //   std::string bfmp = mono->bfmp_ring_conformation_;
    //   std::size_t index;
    //   // This is to only output two decimal places for BFMP values which have a decimal value.
    //   if( ( index = mono->bfmp_ring_conformation_.find( "." ) ) != std::string::npos ) {
    //     bfmp = bfmp.substr( 0, index + 3 ) + ")";
    //   }
    //   gmml::AddLiteral(mono_uri, Ontology::bfmp_ring_conformation, bfmp, mono_stream);
    //   gmml::AddLiteral(mono_uri, "gmmo:fullBFMP", mono->bfmp_ring_conformation_, mono_stream);
    // }
    // if(mono->b_factor_ > 0)
    // {
    //   gmml::AddLiteral(mono_uri, "gmmo:monoBFactor", std::to_string(mono->b_factor_), mono_stream);
    // }
    // gmml::AddLiteral(mono_uri, Ontology::hasSNFGName, mono->SNFG_name_, mono_stream);
    // gmml::AddLiteral(mono_uri, Ontology::hasAuthorSNFGName, mono->author_SNFG_name_, mono_stream);
    // std::size_t offset = 0;
    // int numR;
    // if(mono->on_R_ > 0)
    // {
    //   for(std::vector<std::pair<std::string, std::string> >::iterator derivative =
    //   mono->unknown_derivatives_.begin(); derivative != mono->unknown_derivatives_.end(); derivative++)
    //   {
    //     std::size_t found = mono->sugar_name_.monosaccharide_short_name_.find("<", offset);
    //     if(found != std::string::npos)
    //     {
    //       found+=2;
    //       if(mono->on_R_ < 10)
    //         numR = std::stoi(mono->sugar_name_.monosaccharide_short_name_.substr(found, 1));
    //       else
    //         numR = std::stoi(mono->sugar_name_.monosaccharide_short_name_.substr(found,2));
    //       if((*derivative).second != "")
    //       {
    //         std::stringstream RgroupStream, RnumStream;
    //         RgroupStream << mono_uri << "_R" << numR;
    //         RnumStream << "gmmo:hasR" << numR;
    //         gmml::AddTriple(mono_uri, RnumStream.str(), RgroupStream.str(), mono_stream);
    //         gmml::AddLiteral(RgroupStream.str(), Ontology::hasFormula,(*derivative).second, mono_stream);
    //       }
    //       offset = found;
    //     }
    //     else
    //       break;
    //   }
    // }
    // gmml::AddLiteral(mono_uri, Ontology::author_mono_name, mono->author_sugar_name_.monosaccharide_name_,
    // mono_stream); NoteVector notes = mono->mono_notes_; if(notes.size() != 0)
    // {
    //     int note_id = 1;
    //     // gmml::log(__LINE__, __FILE__,  gmml::INF, "Poulating notes");
    //     std::string id_prefix = mono_uri.substr(5,mono_uri.length()) + "_";
    //     PopulateNotes(mono_stream, mono_stream, mono_uri, notes, id_prefix, note_id);
    // }
    // Glycan::SugarName sugar_name = mono->sugar_name_;
    // PopulateSugarName(mono_stream, id_prefix, mono_uri, mono->mono_id_, sugar_name);
    // mono_stream << ring_atom_stream.str();
}

void Assembly::PopulateRingAtom(std::stringstream& ring_atom_stream, std::string id_prefix, std::string ring_uri,
                                std::string ring_resource, unsigned int ring_index, Atom* ring_atom,
                                Glycan::Monosaccharide* mono, std::vector<std::string>& side_or_ring_atoms)
{
    if (ring_index <= mono->side_atoms_.size())
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
        coord_stream << gmml::ConvertT<double>(coords->GetX()) << ", " << gmml::ConvertT<double>(coords->GetY()) << ", "
                     << gmml::ConvertT<double>(coords->GetZ());
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
        if ((ring_atom->GetName().substr(0, 1).compare("O") !=
             0)) /// side atoms for the oxygen of the ring are not saved
        {
            std::vector<AtomVector> all_sides = mono->side_atoms_;
            AtomVector sides                  = all_sides.at(ring_index - 1);
            std::string side_resource         = "";
            std::string side_uri              = "";
            int side_index                    = 1;
            for (AtomVector::iterator it = sides.begin(); it != sides.end(); it++)
            {
                Atom* side_atom = (*it);
                if (side_atom != NULL)
                {
                    side_resource = CreateURIResource(gmml::OntAtom, 0, id_prefix, side_atom->GetId());
                    side_uri      = CreateURI(side_resource);
                    gmml::AddTriple(ring_uri, Ontology::hasSideAtom, side_uri, ring_atom_stream);

                    PopulateSideAtom(side_atom_stream, id_prefix, side_uri, side_resource, ring_index, side_index,
                                     side_atom, mono, side_or_ring_atoms);
                    side_index++;
                }
            }
        }
        ring_atom_stream << side_atom_stream.str();
    }
}

void Assembly::PopulateSideAtom(std::stringstream& side_atom_stream, std::string id_prefix, std::string side_uri,
                                std::string side_resource, int ring_index, int side_index, Atom* side_atom,
                                Glycan::Monosaccharide* mono, std::vector<std::string>& side_or_ring_atoms)
{
    std::stringstream object;

    std::string chemical_code_str = mono->sugar_name_.chemical_code_string_;
    if ((mono->sugar_name_.ring_type_.compare("P") == 0 && ring_index == 5) ||
        (mono->sugar_name_.ring_type_.compare("F") == 0 && ring_index == 4))
    {
        object << "+" << side_index;
    }
    else if (ring_index == 1 && (side_atom->GetName().substr(0, 1).compare("C") != 0))
    {
        object << "1";
    }
    else if (ring_index == 1 && (side_atom->GetName().substr(0, 1).compare("C") == 0))
    {
        object << "-1";
    }
    else
    {
        object << ring_index;
    }

    if (find(side_or_ring_atoms.begin(), side_or_ring_atoms.end(), side_atom->GetId()) ==
        side_or_ring_atoms.end()) /// if this side atom has not been added to the ontology as side atom of another mono
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
        coord_stream << gmml::ConvertT<double>(coords->GetX()) << ", " << gmml::ConvertT<double>(coords->GetY()) << ", "
                     << gmml::ConvertT<double>(coords->GetZ());
        gmml::AddLiteral(side_uri, Ontology::coordinate, coord_stream.str(), side_atom_stream);

        if (object.str().compare("1") == 0)
        {
            // if(mono->derivatives_map_.find("a") != mono->derivatives_map_.end())
            std::vector<std::pair<std::string, std::string>>::iterator thisPosition =
                std::find_if(mono->derivatives_map_.begin(), mono->derivatives_map_.end(),
                             [&](const std::pair<std::string, std::string>& element)
                             {
                                 return element.first == "a";
                             });
            if (thisPosition != mono->derivatives_map_.end())
            {
                gmml::AddLiteral(side_uri, Ontology::derivative, (*thisPosition).second, side_atom_stream);
            }
        }
        else
        {
            std::vector<std::pair<std::string, std::string>>::iterator thisPosition =
                std::find_if(mono->derivatives_map_.begin(), mono->derivatives_map_.end(),
                             [&](const std::pair<std::string, std::string>& element)
                             {
                                 return element.first == object.str();
                             });
            // if(mono->derivatives_map_.find(object.str()) != mono->derivatives_map_.end())
            if (thisPosition != mono->derivatives_map_.end())
            {
                gmml::AddLiteral(side_uri, Ontology::derivative, (*thisPosition).second, side_atom_stream);
            }
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
    if (object.str().compare("1") == 0)
    {
        index = chemical_code_str.find("a");
    }
    else
    {
        index = chemical_code_str.find(object.str().c_str());
    }
    if (index != std::string::npos)
    {
        index--;
        if (chemical_code_str.at(index) == '^')
        {
            gmml::AddLiteral(side_uri, Ontology::orientation, "Up", side_atom_stream);
        }
        else if (chemical_code_str.at(index) == '_')
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

void Assembly::PopulateSugarName(std::stringstream& mono_stream, std::string id_prefix, std::string mono_uri,
                                 int mono_id_, Glycan::SugarName sugar_name)
{
    std::stringstream sugar_name_stream;

    std::string sugar_name_resource = "";
    std::string sugar_name_uri      = "";

    sugar_name_resource = CreateURIResource(gmml::OntSugarName, mono_id_, id_prefix, "");
    sugar_name_uri      = CreateURI(sugar_name_resource);

    gmml::AddTriple(mono_uri, Ontology::hasSugarName, sugar_name_uri, mono_stream);

    //    sugar_name_stream << Ontology::ENTITY_COMMENT << sugar_name_resource << std::endl;
    gmml::AddTriple(sugar_name_uri, Ontology::TYPE, Ontology::SugarName, sugar_name_stream);
    gmml::AddLiteral(sugar_name_uri, Ontology::id, sugar_name_resource, sugar_name_stream);
    //    gmml::AddLiteral(sugar_name_uri, Ontology::LABEL, sugar_name_resource, sugar_name_stream);
    gmml::AddLiteral(sugar_name_uri, Ontology::mono_stereo_name, sugar_name.monosaccharide_stereochemistry_name_,
                     sugar_name_stream);
    gmml::AddLiteral(sugar_name_uri, Ontology::mono_stereo_short_name,
                     sugar_name.monosaccharide_stereochemistry_short_name_, sugar_name_stream);
    gmml::AddLiteral(sugar_name_uri, Ontology::mono_name, sugar_name.monosaccharide_name_, sugar_name_stream);

    gmml::AddLiteral(sugar_name_uri, Ontology::mono_short_name, sugar_name.monosaccharide_short_name_,
                     sugar_name_stream);
    gmml::AddLiteral(sugar_name_uri, Ontology::isomer, sugar_name.isomer_, sugar_name_stream);
    if (sugar_name.configuration_.compare("a") == 0)
    {
        gmml::AddLiteral(sugar_name_uri, Ontology::configuration, "alpha", sugar_name_stream);
    }
    else if (sugar_name.configuration_.compare("b") == 0)
    {
        gmml::AddLiteral(sugar_name_uri, Ontology::configuration, "beta", sugar_name_stream);
    }
    if (sugar_name.ring_type_.compare("") != 0)
    {
        gmml::AddLiteral(sugar_name_uri, Ontology::ring_type, sugar_name.ring_type_, sugar_name_stream);
    }

    mono_stream << sugar_name_stream.str();
}

void Assembly::PopulateResidue(std::stringstream& pdb_stream, std::stringstream& residue_stream, std::string pdb_uri,
                               std::string id_prefix, ResidueVector residues,
                               std::vector<std::string> side_or_ring_atoms)
{
    // All this seems to do is write out the atomic coordinates and
    // a TON of unnecessary triples which waste a lot of time and space.
    std::string res_resource  = "";
    std::string res_uri       = "";
    std::string atom_resource = "";
    std::string atom_uri      = "";

    std::stringstream atom_stream;
    for (ResidueVector::iterator it = residues.begin(); it != residues.end(); it++)
    {
        Residue* residue = (*it);
        res_resource     = CreateURIResource(gmml::OntResidue, 0, id_prefix, residue->GetId());
        res_uri          = CreateURI(res_resource);

        // gmml::AddTriple(pdb_uri, Ontology::hasResidue, res_uri, pdb_stream);
        pdb_stream << "\t" << Ontology::hasResidue << "\t" << res_uri << ";\n";

        gmml::AddTriple(res_uri, Ontology::TYPE, Ontology::Residue, residue_stream);
        gmml::AddLiteral(res_uri, Ontology::id, res_resource, residue_stream);

        // Way too granular for our needs.  If they want this they can look at the pdb file

        // AtomVector res_atoms = residue->GetAtoms();
        // for(AtomVector::iterator it1 = res_atoms.begin(); it1 != res_atoms.end(); it1++)
        // {
        //     Atom* atom = (*it1);
        //     atom_resource = CreateURIResource(gmml::OntAtom, 0, id_prefix, atom->GetId());
        //     atom_uri = CreateURI(atom_resource);

        //     gmml::AddTriple(res_uri, Ontology::hasAtom, atom_uri, residue_stream);

        //     if(find(side_or_ring_atoms.begin(), side_or_ring_atoms.end(), atom->GetId()) == side_or_ring_atoms.end())
        //     {
        //         PopulateAtom(atom_stream, atom_uri, atom_resource, id_prefix, atom);
        //     }
        // }
    }
    residue_stream << atom_stream.str();
}

void Assembly::PopulateAtom(std::stringstream& atom_stream, std::string atom_uri, std::string atom_resource,
                            std::string id_prefix, Atom* atom)
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
    coord_stream << gmml::ConvertT<double>(coords->GetX()) << ", " << gmml::ConvertT<double>(coords->GetY()) << ", "
                 << gmml::ConvertT<double>(coords->GetZ());
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

    int title_length  = pdb_resource.length() + 2;
    int title_padding = (48 - title_length) / 2;

    pdb_stream << std::setw(4) << std::setfill(' ') << " ";
    pdb_stream << std::setw(48) << std::setfill('#') << "";
    pdb_stream << std::setw(4) << std::setfill(' ') << " "
               << "\n";

    pdb_stream << std::setw(4) << std::setfill(' ') << " ";
    pdb_stream << "#" << std::setw(title_padding) << " ";
    pdb_stream << pdb_resource;
    pdb_stream << std::setw(title_padding) << " "
               << "#";
    pdb_stream << std::setw(4) << std::setfill(' ') << " "
               << "\n";

    pdb_stream << std::setw(4) << std::setfill(' ') << " ";
    pdb_stream << std::setw(48) << std::setfill('#') << "";
    pdb_stream << std::setw(4) << std::setfill(' ') << " "
               << "\n\n";
}

void Assembly::CreateSubtitle(std::string pdb_resource, std::stringstream& pdb_stream)
{
    // adding 12 for information and 2 for #
    int title_length  = pdb_resource.length() + 12 + 2;
    int title_padding = (48 - title_length) / 2;

    pdb_stream << std::setw(8) << std::setfill(' ') << " ";
    pdb_stream << std::setw(48) << std::setfill('#') << "";
    pdb_stream << std::setw(8) << std::setfill(' ') << " "
               << "\n";

    pdb_stream << std::setw(8) << std::setfill(' ') << " ";
    pdb_stream << "#" << std::setw(title_padding) << " ";
    pdb_stream << pdb_resource << " Information";
    pdb_stream << std::setw(title_padding) << " "
               << "#";
    pdb_stream << std::setw(8) << std::setfill(' ') << " "
               << "\n";

    pdb_stream << std::setw(8) << std::setfill(' ') << " ";
    pdb_stream << std::setw(48) << std::setfill('#') << "";
    pdb_stream << std::setw(8) << std::setfill(' ') << " "
               << "\n";
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

std::string Assembly::CreateURIResource(gmml::URIType resource, int number, std::string id_prefix, std::string id)
{
    std::stringstream uri_resource;
    switch (resource)
    {
        case gmml::OntPDB:
            uri_resource
                << (gmml::Split(this->GetSourceFile().substr(this->GetSourceFile().find_last_of('/') + 1), ".").at(0));
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
            replace(id.begin(), id.end(), '?', 'n'); // replace all '?' with 'n'
            gmml::FindReplaceString(id, "\'", "q");
            gmml::FindReplaceString(id, ",", "c");
            // gmml::FindReplaceString(id, "*", "s");
            replace(id.begin(), id.end(), '*', 's'); // replace all '*' with 's'
            uri_resource << id_prefix << id;
            break;
        case gmml::OntResidue:
            replace(id.begin(), id.end(), '?', 'n'); // replace all '?' with 'n'
            replace(id.begin(), id.end(), '*', 's'); // replace all '*' with 's'
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

    if (output_file_type.compare("csv") == 0)
    {
        curl << Ontology::CSV_OUTPUT_FORMAT;
    }
    else if (output_file_type.compare("xml") == 0)
    {
        curl << Ontology::XML_OUTPUT_FORMAT;
    }
    else if (output_file_type.compare("json") == 0)
    {
        curl << Ontology::JSON_OUTPUT_FORMAT;
    }

    curl << Ontology::DATA_STORE_ADDRESS << Ontology::QUERY_PREFIX << query << Ontology::QUERY_POSTFIX;
    std::string tmp  = curl.str();
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

    if (output_file_type.compare("csv") == 0)
    {
        curl << Ontology::CSV_OUTPUT_FORMAT;
    }
    else if (output_file_type.compare("xml") == 0)
    {
        curl << Ontology::XML_OUTPUT_FORMAT;
    }
    else if (output_file_type.compare("json") == 0)
    {
        curl << Ontology::JSON_OUTPUT_FORMAT;
    }

    // Had to change this from DATA_STORE_ADDRESS_GF to the url parameter to allow for Docker to work.
    // TODO Functions like this should not be in GMML C++. They should be at the GEMS scripting level.
    curl << url << Ontology::QUERY_PREFIX << query << Ontology::QUERY_POSTFIX;
    std::string tmp = curl.str();
    // std::cout << std::endl << "RESULTS: " << std::endl;
    //  const char* cstr = tmp.c_str();
    // system(cstr);
    // std::cout << std::endl;
    return tmp;
}
