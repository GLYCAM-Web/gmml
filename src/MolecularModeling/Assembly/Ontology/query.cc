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

#include <unistd.h>
#include <errno.h>
#include <string.h>

using MolecularModeling::Assembly;

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////
void Assembly::ExtractOntologyInfoByNameOfGlycan(std::string stereo_name, std::string stereo_condensed_name, std::string name, std::string condensed_name, std::string output_file_type)
{
    if(stereo_name.compare("") == 0 && stereo_condensed_name.compare("") == 0 && name.compare("") == 0 && condensed_name.compare("") == 0)
    {
//        std::cout << "Please specify at least one of the arguments and set the others as \"\" " << std::endl;
        return;
    }
    std::stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE << " ?pdb ?residue_id ?stereo_name ?stereo_condensed_name ?name ?condensed_name ?shape " << Ontology::WHERE_CLAUSE;
    query << "?pdb_file     :hasOligo   ?oligo.\n";
    query << "?oligo        :hasCore    ?mono.\n";
    query << "?mono         :hasSugarName   ?sugarName.\n";
    if(stereo_name.compare("") != 0)
        query << "?sugarName    :monosaccharideStereochemName   \"" << stereo_name << "\".\n";
    if(stereo_condensed_name.compare("") != 0)
        query << "?sugarName    :monosaccharideStereochemShortName   \"" << stereo_condensed_name << "\".\n";
    if(name.compare("") != 0)
        query << "?sugarName    :monosaccharideName   \"" << name << "\".\n";
    if(condensed_name.compare("") != 0)
        query << "?sugarName    :monosaccharideShortName   \"" << condensed_name << "\".\n";
    query << "?pdb_file     :identifier   ?pdb.\n";

    query << "?mono         :hasRingAtom   ?ring_atoms.\n";
    query << "?residue      :hasAtom       ?ring_atoms.\n";
    query << "?residue      :identifier    ?residue_id.\n";

    query << "?sugarName    :monosaccharideStereochemName   ?stereo_name.\n";
    query << "?sugarName    :monosaccharideStereochemShortName   ?stereo_condensed_name.\n";
    query << "?sugarName    :monosaccharideName   ?name.\n";
    query << "?sugarName    :monosaccharideShortName   ?condensed_name.\n";
    query << "OPTIONAL { ?mono         :BFMPRingConformation   ?shape.}\n";
    query << Ontology::END_WHERE_CLAUSE;

    FormulateCURL(output_file_type, query.str());
}

void Assembly::ExtractOntologyInfoByNamePartsOfGlycan(std::string isomer, std::string ring_type, std::string configuration, std::string output_file_type)
{
    if(isomer.compare("") == 0 && ring_type.compare("") == 0 && configuration.compare("") == 0)
    {
//        std::cout << "Please specify at least one of the arguments and set the others as \"\" " << std::endl;
        return;
    }
    std::stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE << " ?pdb ?stereo_name ?stereo_condensed_name ?name ?condensed_name ?shape " << Ontology::WHERE_CLAUSE;
    query << "?pdb_file     :hasOligo   ?oligo.\n";
    query << "?oligo        :hasCore    ?mono.\n";
    query << "?mono         :hasSugarName   ?sugarName.\n";
    if(isomer.compare("") != 0)
        query << "?sugarName    :isomer   \"" << isomer << "\".\n";
    if(ring_type.compare("") != 0)
        query << "?sugarName    :ringType   \"" << ring_type << "\".\n";
    if(configuration.compare("") != 0)
        query << "?sugarName    :configuration   \"" << configuration << "\".\n";
    query << "?pdb_file     :identifier   ?pdb.\n";

//    query << "?mono         :hasRingAtom   ?ring_atoms.\n";
//    query << "?residue      :hasAtom       ?ring_atoms.\n";
//    query << "?residue      :identifier    ?residue_id.\n";

    query << "?sugarName    :monosaccharideStereochemName   ?stereo_name.\n";
    query << "?sugarName    :monosaccharideStereochemShortName   ?stereo_condensed_name.\n";
    query << "?sugarName    :monosaccharideName   ?name.\n";
    query << "?sugarName    :monosaccharideShortName   ?condensed_name.\n";
    query << "OPTIONAL { ?mono         :BFMPRingConformation   ?shape.}\n";
    query << Ontology::END_WHERE_CLAUSE;

    FormulateCURL(output_file_type, query.str());
}

void Assembly::ExtractOntologyInfoByPDBID(std::string pdb_id, std::string output_file_type)
{
    if(pdb_id.compare("") == 0)
    {
//        std::cout << "Please specify the input argument." << std::endl;
        return;
    }
    std::stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE << " ?oligo_sequence ?residue_links ?glycosidic_linkage " << Ontology::WHERE_CLAUSE;
    query <<  ":" << pdb_id << "    :hasOligo   ?oligo.\n";
    query << "?oligo    :oligoName 	?oligo_sequence.\n";
    query << "OPTIONAL { ?oligo	:oligoResidueLinks	?residue_links.\n";
    query << "?linkage 	:hasParent 	?oligo.\n";
    query << "?linkage	:glycosidicLinkage    ?glycosidic_linkage.}\n";

    //query << "?oligo	:hasCore	?mono.\n";
    //query << "?mono     :anomericStatus    ?anomeric_status.\n";

    query << Ontology::END_WHERE_CLAUSE;

    FormulateCURL(output_file_type, query.str());
}

void Assembly::ExtractOntologyInfoByStringChemicalCode(std::string chemical_code, std::string output_file_type)
{
    if(chemical_code.compare("") == 0)
    {
//        std::cout << "Please specify the input argument." << std::endl;
        return;
    }
    std::stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE << " ?pdb ?name ?short_name ?stereo_name ?stereo_short_name " << Ontology::WHERE_CLAUSE;
    query << "?mono         :stereochemistryChemicalCode	   \"" << chemical_code << "\".\n";
    query << "?pdb_file     :hasOligo	?oligo.\n";
    query << "?oligo        :hasCore	?mono.\n";
    query << "?mono         :hasSugarName	?sn.\n";
    query << "?sn           :monosaccharideName 	?name.\n";
    query << "?sn           :monosaccharideShortName 	?short_name.\n";
    query << "?sn           :monosaccharideStereochemName 	?stereo_name.\n";
    query << "?sn           :monosaccharideStereochemShortName 	?stereo_short_name.\n";
    query << "?pdb_file     :identifier   ?pdb.\n";
    query << Ontology::END_WHERE_CLAUSE;

    FormulateCURL(output_file_type, query.str());
}

void Assembly::ExtractOntologyInfoByOligosaccharideNameSequence(std::string oligo_name, std::string output_file_type)
{
    std::stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE << " ?pdb ?residue_links ?glycosidic_linkage " << Ontology::WHERE_CLAUSE;

    query << "?pdb_file     :hasOligo	?oligo.\n";
    query << "?oligo        :oligoName	\"" << oligo_name << "\"\n";
    query << "OPTIONAL { ?oligo        :oligoResidueLinks	?residue_links.\n";
    query << "?linkage      :hasParent 	?oligo.\n";
    query << "?linkage      :glycosidicLinkage    ?glycosidic_linkage.}\n";
    query << "?pdb_file     :identifier   ?pdb.\n";

    ///To DO: std::string manipulation: split the oligo_name by _ and for each of them write the following to represent the names of the monos:
    //    query << "?oligo	:hasCore	?mono.\n";
    //    query << "?mono     :hasSugarName	?sn.\n";
    //    query << "?sn       :monosaccharideName 	?name.\n";
    //    query << "?sn       :monosaccharideShortName 	?short_name.\n";
    //    query << "?sn       :monosaccharideStereochemName 	?stereo_name.\n";
    //    query << "?sn       :monosaccharideStereochemShortName 	?stereo_short_name.\n";
    query << Ontology::END_WHERE_CLAUSE;

    FormulateCURL(output_file_type, query.str());
}

void Assembly::ExtractOntologyInfoByOligosaccharideNameSequenceByRegex(std::string oligo_name_pattern, std::string output_file_type)
{
    gmml::FindReplaceString(oligo_name_pattern, "[", "\\\\[");
    gmml::FindReplaceString(oligo_name_pattern, "]", "\\\\]");
    if(oligo_name_pattern.compare("") == 0)
    {
//        std::cout << "Please specify the input argument. (you can use up to two * in the name pattern)" << std::endl;
        return;
    }
    if(count(oligo_name_pattern.begin(), oligo_name_pattern.end(), '*') > 3)
    {
//        std::cout << "Wrong name pattern format. Please use only up tp three * in the input argument." << std::endl;
        return;
    }

    std::stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE << " ?pdb ?oligo_sequence ?residue_links ?glycosidic_linkage " << Ontology::WHERE_CLAUSE;
    query << "?oligo        :oligoName	?oligo_sequence.\n";

    size_t first = oligo_name_pattern.find_first_of("*");
    size_t last = oligo_name_pattern.find_last_of("*");

    std::string filter1 = oligo_name_pattern.substr(0, first);
    std::string filter2 = oligo_name_pattern.substr(first + 1, last - 1);
    std::string filter3 = oligo_name_pattern.substr(last + 1, oligo_name_pattern.size() - 1);
    if(count(oligo_name_pattern.begin(), oligo_name_pattern.end(), '*') == 0) ///No *
        query << "?oligo	:oligoName	\"" << oligo_name_pattern << "\".\n";
    else if(count(oligo_name_pattern.begin(), oligo_name_pattern.end(), '*') == 1) ///Only one *
    {
        if(first == 0) ///* at the beginning
            query << "FILTER regex(?oligo_sequence, \"" << filter2 << "$\", \"i\")\n";
        else if(first == oligo_name_pattern.size()-1) ///* at the end
            query << "FILTER regex(?oligo_sequence, \"^" << filter1 << "\", \"i\")\n";
        else if(first < oligo_name_pattern.size()-1 && first > 0) ///* in the middle
            query << "FILTER regex(?oligo_sequence, \"^" << filter1 << ".+" << filter3 << "$\", \"i\")\n";
    }
    else if (count(oligo_name_pattern.begin(), oligo_name_pattern.end(), '*') == 2)
    {
        if(first == 0 && last == oligo_name_pattern.size() - 1) ///* at the beginning and end
            query << "FILTER regex(?oligo_sequence, \"" << filter2 << "\", \"i\")\n";
        else if(first == 0 && last < oligo_name_pattern.size() - 1)///one * at the beginning another in the middle
            query << "FILTER regex(?oligo_sequence, \"" << filter2 << ".+" << filter3 << "$\", \"i\")\n";
        else if(first > 0 && last == oligo_name_pattern.size() - 1)///one * in the middle another at the end
            query << "FILTER regex(?oligo_sequence, \"^" << filter1 << ".+" << filter2 << "\", \"i\")\n";
    }
    else
    {
        std::vector<std::string> pattern_tokens = gmml::Split(filter2, "*");
        query << "FILTER regex(?oligo_sequence, \"" << pattern_tokens.at(0) << ".+" << pattern_tokens.at(1) << "\", \"i\")\n";
    }

    query << "?pdb_file     :hasOligo	?oligo.\n";
    query << "OPTIONAL { ?oligo        :oligoResidueLinks	?residue_links.\n";
    query << "?linkage      :hasParent 	?oligo.\n";
    query << "?linkage      :glycosidicLinkage    ?glycosidic_linkage.}\n";
    query << "?pdb_file     :identifier	?pdb.\n";

    query << Ontology::END_WHERE_CLAUSE;

    FormulateCURL(output_file_type, query.str());
}

void Assembly::ExtractOntologyInfoByGlycanStructure(std::string ring_type, std::string anomeric_orientation, std::string minus_one_orientation, std::string index_two_orientation, std::string index_three_orientation,
                                                    std::string index_four_orientation, std::string plus_one_orientation, std::string output_file_type)
{
    if(ring_type.compare("") == 0)
    {
//        std::cout << "Please specify the ring type which is the first argument of the function as either \"P\" or \"F\" " << std::endl;
        return;
    }
    std::stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE << " ?stereo_name ?stereo_condensed_name ?condensed_name ?name ?oligo_sequence ?residue_links ?glycosidic_linkage " << Ontology::WHERE_CLAUSE;
    query << "?pdb_file     :hasOligo       ?oligo.\n";
    query << "?oligo	    :hasCore        ?mono.\n";

    if(anomeric_orientation.compare("") != 0)
    {
        query << "?mono         :hasRingAtom	?anomeric.\n";
        query << "?anomeric	    :ringIndex  	\"1\".\n";
        query << "?anomeric	    :hasSideAtom    ?a_side.\n";
        query << "?a_side	    :sideIndex      \"1\".\n";
        query << "?a_side	    :orientation	\"" << anomeric_orientation << "\".\n";
    }
    if(anomeric_orientation.compare("") != 0 && minus_one_orientation.compare("") != 0)
    {
        query << "?anomeric     :hasSideAtom    ?a_minus_side.\n";
        query << "?a_minus_side	:sideIndex      \"-1\".\n";
        query << "?a_side	    :orientation	\"" << anomeric_orientation << "\".\n";
    }
    else if(minus_one_orientation.compare("") != 0)
    {
        query << "?mono         :hasRingAtom	?anomeric.\n";
        query << "?anomeric     :ringIndex  	\"1\".\n";
        query << "?anomeric     :hasSideAtom    ?a_minus_side.\n";
        query << "?a_minus_side	:sideIndex      \"-1\".\n";
        query << "?a_minus_side	:orientation	\"" << minus_one_orientation << "\".\n";
    }
    if(index_two_orientation.compare("") != 0)
    {
        query << "?mono         :hasRingAtom	?two.\n";
        query << "?two  	    :ringIndex  	\"2\".\n";
        query << "?two          :hasSideAtom    ?two_side.\n";
        query << "?two_side	    :sideIndex      \"2\".\n";
        query << "?two_side	    :orientation	\"" << index_two_orientation << "\".\n";
    }
    if(index_three_orientation.compare("") != 0)
    {
        query << "?mono         :hasRingAtom	?three.\n";
        query << "?three        :ringIndex  	\"3\".\n";
        query << "?three        :hasSideAtom    ?three_side.\n";
        query << "?three_side	:sideIndex      \"3\".\n";
        query << "?three_side	:orientation	\"" << index_three_orientation << "\".\n";
    }
    if(index_four_orientation.compare("") != 0)
    {
        query << "?mono         :hasRingAtom	?four.\n";
        query << "?four         :ringIndex  	\"4\".\n";
        query << "?four         :hasSideAtom    ?four_side.\n";
        query << "?four_side    :sideIndex      \"4\".\n";
        query << "?four_side	:orientation	\"" << index_four_orientation << "\".\n";
    }
    if(plus_one_orientation.compare("") != 0)
    {
        query << "?mono         :hasRingAtom	?last_c.\n";
        if(ring_type.compare("P") == 0 )
            query << "?last_c       :ringIndex  	\"5\".\n";
        else
            query << "?last_c       :ringIndex  	\"4\".\n";
        query << "?last_c        :hasSideAtom    ?plus_one.\n";
        query << "?plus_one      :sideIndex      \"+1\".\n";
        query << "?plus_one    	 :orientation	\"" << plus_one_orientation << "\".\n";
    }

    query << "?oligo        :oligoName  ?oligo_sequence.\n";
    query << "?pdb_file     :identifier	?pdb.\n";
    query << "OPTIONAL { ?oligo        :oligoResidueLinks	?residue_links.\n";
    query << "?linkage      :hasParent 	?oligo.\n";
    query << "?linkage      :glycosidicLinkage    ?glycosidic_linkage.}\n";
    query << "?mono         :hasSugarName	?sn.\n";
    query << "?sn           :monosaccharideName 	?name.\n";
    query << "?sn           :monosaccharideShortName 	?condensed_name.\n";
    query << "?sn           :monosaccharideStereochemName 	?stereo_name.\n";
    query << "?sn           :monosaccharideStereochemShortName 	?stereo_condensed_name.\n";
    query << Ontology::END_WHERE_CLAUSE;

    FormulateCURL(output_file_type, query.str());

}

void Assembly::ExtractOntologyInfoByDerivativeModificationMap(std::string ring_type, DerivativeModificationMap derivative_modification_map, std::string output_file_type)
{
    if(ring_type.compare("") == 0)
    {
//        std::cout << "Please specify the ring type as the first argument of the function as either \"P\" or \"F\" " << std::endl;
        return;
    }
    std::stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE << "?pdb ?stereo_name ?stereo_condensed_name ?condensed_name " << Ontology::WHERE_CLAUSE;
    for(DerivativeModificationMap::iterator it = derivative_modification_map.begin(); it != derivative_modification_map.end(); it++)
    {
        std::string index = (*it).first;
        std::string pattern = (*it).second;
        query << "?mono	       :hasRingAtom    ?ring_atom.\n";
        if(index.compare("-1") != 0 && index.compare("+1") != 0)
        {
            query << "?ring_atom    :ringIndex      \"" << index << "\".\n";
            query << "?ring_atom    :hasSideAtom    ?side.\n";
        }
        else if(index.compare("-1") == 0 )
        {
            query << "?ring_atom       :ringIndex      \"1\".\n";
            query << "?ring_atom       :hasSideAtom    ?side.\n";
        }
        else if(index.compare("+1") == 0 )
        {
            query << "?mono         :hasRingAtom	?ring_atom.\n";
            if(ring_type.compare("P") == 0 )
                query << "?ring_atom       :ringIndex  	\"5\".\n";
            else
                query << "?ring_atom       :ringIndex  	\"4\".\n";
            query << "?ring_atom        :hasSideAtom    ?side.\n";
        }
        query << "?side         :sideIndex      \"" << index << "\".\n";
        query << "?side         :derivative     \"" << pattern << "\".\n";
        query << "?mono         :hasSugarName    ?sn.\n";
        query << "?sn           :monosaccharideName 	?name.\n";
        query << "?sn           :monosaccharideShortName 	?condensed_name.\n";
        query << "?sn           :monosaccharideStereochemName 	?stereo_name.\n";
        query << "?sn           :monosaccharideStereochemShortName 	?stereo_condensed_name.\n";
        query << "?pdb_file     :hasOligo 	?oligo.\n";
        query << "?oligo        :hasCore 	?mono.\n";
        query << "?pdb_file     :identifier	?pdb.\n";
    }
    query << Ontology::END_WHERE_CLAUSE;

    FormulateCURL(output_file_type, query.str());
}

void Assembly::ExtractOntologyInfoByAttachedGlycanStructures(AttachedGlycanStructuresVector attached_structures, std::string output_file_type)
{
    std::stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE << " ?pdb ?linkage_indices ?stereo_short_name0 ?short_name0 ?stereo_short_name1 ?short_name1 ?oligo_sequence ?residue_links "
          << Ontology::WHERE_CLAUSE;
    unsigned int i = 0;
    std::vector<std::string> oligos = std::vector<std::string>();
    for(AttachedGlycanStructuresVector::iterator it = attached_structures.begin(); it != attached_structures.end(); it++)
    {
        std::vector<std::string> structure = (*it);
        if(structure.size() < 7)
        {
//            std::cout << "Missing arguments! All should be set even as an empty value (for empty values set \"\")" << std::endl;
            return;
        }
        if(structure.at(0).compare("") == 0)
        {
//            std::cout << "Please specify the ring type which is the first argument of the function as either \"P\" or \"F\" " << std::endl;
            return;
        }
        std::stringstream oligo;
        oligo << "?oligo" << i;
        std::stringstream mono;
        mono << "?mono" << i;
        query << oligo.str() << "		:hasCore	" << mono.str() << ".\n";

        for(int j = 1; j < 7; j++)
        {
            if(structure.at(j).compare("") != 0)
            {
                std::stringstream ring_atom;
                std::stringstream side_atom;
                switch (j)///anomeric, -1 and +1 are special cases
                {
                    case 1:///anomeric
                        ring_atom << mono.str() << "_anomeric";
                        side_atom << mono.str() << "_anomeric_side_atom";
                        query << mono.str() << "     :hasRingAtom	" << ring_atom.str() << ".\n";
                        query << ring_atom.str() << "  	:ringIndex  	\"" << j << "\".\n";
                        query << ring_atom.str() << "   :hasSideAtom    " << side_atom.str() << " .\n";
                        query << side_atom.str() << "	:sideIndex      \"" << j << "\".\n";
                        query << side_atom.str() << "	:orientation	\"" << structure.at(j) << "\".\n";
                        break;
                    case 2:///minus one
                        ring_atom << mono.str() << "_anomeric";
                        side_atom << mono.str() << "_anomeric_minus_one_side_atom";
                        if(structure.at(1).compare("") == 0)///anomeric has not been set
                        {
                            query << mono.str() << "     :hasRingAtom	" << ring_atom.str() << ".\n";
                            query << ring_atom.str() << "  	:ringIndex  	\""<< j - 1 << "\".\n";
                        }
                        query << ring_atom.str() << "   :hasSideAtom    " << side_atom.str() << " .\n";
                        query << ring_atom.str() << "	:sideIndex     \"-1\".\n";
                        query << side_atom.str() << "	:orientation	\"" << structure.at(j) << "\".\n";
                        break;
                    case 6:///plus one
                        ring_atom << mono.str() << "_last_c";
                        side_atom << mono.str() << "_last_c_side_atom";
                        query << mono.str() << "         :hasRingAtom	" << ring_atom.str() << ".\n";
                        if(structure.at(0).compare("P") == 0 )
                            query << ring_atom.str() << "       :ringIndex  	\"" << j - 1 << "\".\n";
                        else
                            query << ring_atom.str() << "       :ringIndex  	\"" << j - 2 << "\".\n";
                        query << ring_atom.str() << "         :hasSideAtom    " << side_atom.str() << ".\n";
                        query << side_atom.str() << "       :sideIndex      \"+1\".\n";
                        query << side_atom.str() << "    	  :orientation	\"" << structure.at(j) << "\".\n";
                        break;
                    default:
                        ring_atom << mono.str() << "_ring_atom" << j - 1;
                        side_atom << mono.str() << "_side_atom" << j - 1;
                        query << mono.str() << "     :hasRingAtom	" << ring_atom.str() << ".\n";
                        query << ring_atom.str() << "  	:ringIndex  	\"" << j - 1 << "\".\n";
                        query << ring_atom.str() << "   :hasSideAtom    " << side_atom.str() << " .\n";
                        query << side_atom.str() << "	:sideIndex      \"" << j - 1 << "\".\n";
                        query << side_atom.str() << "	:orientation	\"" << structure.at(j) << "\".\n";
                }
            }
        }
        i++;
        oligos.push_back(oligo.str());
    }
    for(i = 0; i < oligos.size(); i++)
    {
        if(i + 1 < oligos.size())
        {
            query << "{\n";
            query << "?linkage" << i << " :hasParent " << oligos.at(i) << ".\n";
            query << "?linkage" << i << " :hasChild " << oligos.at(i + 1) << ".\n";
            query << oligos.at(i) << " :oligoResidueLinks ?residue_links.\n";
            query << "OPTIONAL {" << oligos.at(i) << "  :oligoName ?oligo_sequence}\n";
            query << "} UNION {\n";
            query << "?linkage" << i << " :hasParent " << oligos.at(i + 1) << ".\n";
            query << "?linkage" << i << " :hasChild " << oligos.at(i) << ".\n";
            query << oligos.at(i + 1) << " :oligoResidueLinks ?residue_links.\n";
            query << "OPTIONAL {" << oligos.at(i + 1) << "  :oligoName ?oligo_sequence}\n";
            query << "}\n";

            query << "?linkage" << i << ":linkageIndeces ?linkage_indices.\n";

            query << "?mono0    :hasSugarName ?sn0.\n";
            query << "?sn0      :monosaccharideShortName ?short_name0.\n";
            query << "?sn0      :monosaccharideStereochemShortName ?stereo_short_name0.\n";

            query << "?mono1    :hasSugarName ?sn1.\n";
            query << "?sn1      :monosaccharideShortName ?short_name1.\n";
            query << "?sn1      :monosaccharideStereochemShortName ?stereo_short_name1.\n";
        }
    }
    for(i = 0; i < oligos.size(); i++)
    {
        query << "?pdb            :hasOligo       " << oligos.at(i) << ".\n";
        ///optional {?oligo0    :oligoName ?oligoName0.}???
        ///optional {?oligo1    :oligoName ?oligoName1.}??? only one of these should be matched
    }
    query << Ontology::END_WHERE_CLAUSE;

    FormulateCURL(output_file_type, query.str());
}

void Assembly::ExtractOntologyInfoByNote(std::string pdb_id, std::string note_type, std::string note_category, std::string output_file_type)
{
    std::stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE << " ?pdb ?note_type ?note_category ?description "<< Ontology::WHERE_CLAUSE;

    if(pdb_id.compare("") != 0)
    {
        query <<  ":" << pdb_id << "    :hasNote   ?note.\n";
        query <<  ":" << pdb_id << "    :identifier    ?pdb.\n";
    }
    else
    {
        query << "?pdb_file      :hasNote    ?note.\n";
        query << "?pdb_file      :identifier    ?pdb.\n";
    }
    if(note_type.compare("") != 0)
        query << "?note	         :NoteType      \"" << note_type << "\".\n";
    query << "?note	       :NoteType    ?note_type.\n";

    if(note_category.compare("") != 0)
        query << "?note	       :NoteCategory      \"" << note_category << "\".\n";
    query << "?note	       :NoteCategory    ?note_category.\n";
    query << "?note	       :description    ?description.\n";

    query << Ontology::END_WHERE_CLAUSE;
    FormulateCURL(output_file_type, query.str());
}

void Assembly::ExtractOntologyInfoByCustomQuery(std::string query_file, std::string output_file_type)
{
    std::string line;
    std::stringstream query;
    std::ifstream in(query_file.c_str());

    if (!in.is_open())
    {
//        std::cout << "Error in reading the query file" << std::endl;
        return;
    }

    while (getline (in, line))
    {
        query << line << std::endl;
    }
    in.close();

    FormulateCURL(output_file_type, query.str());
}

void Assembly::ExtractAtomCoordinatesForTorsionAnglesFromOntologySlow(std::string disaccharide_pattern/*, std::string output_file_type*/)
{

    int link_index = disaccharide_pattern.find_first_of("-");
    std::string child_mono = disaccharide_pattern.substr(0, link_index - 1); /// e.g DNeupNAca in DNeupNAca2-3DGalpb
    std::string parent_mono = disaccharide_pattern.substr(link_index + 2, disaccharide_pattern.size()); /// e.g DGalpb in DNeupNAca2-3DGalpb
    std::string linkage_indeces = disaccharide_pattern.substr(link_index - 1, 3); /// e.g 2-3 in DNeupNAca2-3DGalpb
    bool omega = false;
    if(linkage_indeces.find("6") != std::string::npos || linkage_indeces.find("7") != std::string::npos
             || linkage_indeces.find("8") != std::string::npos || linkage_indeces.find("9") != std::string::npos) /// Disaccharides with any of 1-6, 1-7, 2-6, and 2-7 linkages might have omega torsion angles
        omega = true;

    std::stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE << " ?pdb ?O5_crd ?C1_crd ?Ox_crd ?Cx ?Cx_crd ?Cx_neighbor ?Cx_neighbor_crd ?O5_prime_crd "<< Ontology::WHERE_CLAUSE;

    if(child_mono.compare("*") == 0)
        query <<  "?oligo1        :hasCore    ?mono1. \n";
    else
    {
        if(child_mono.find("*") != std::string::npos)
        {
            query <<  "?sn1           :monosaccharideShortName    ?short_name1. \n";
            query << "FILTER regex(?short_name1, \"" << gmml::Split(child_mono, "*").at(0) << "\", \"i\")\n";
        }
        else
            query <<  "?sn1           :monosaccharideShortName    \"" << child_mono << "\". \n";
        query <<  "?mono1         :hasSugarName    ?sn1.\n";
        query <<  "?oligo1        :hasCore    ?mono1. \n";
    }

    if(parent_mono.compare("*") == 0)
        query <<  "?oligo2        :hasCore    ?mono2. \n";
    else
    {
        if(parent_mono.find("*") != std::string::npos)
        {
            query <<  "?sn2           :monosaccharideShortName    ?short_name2. \n";
            query << "FILTER regex(?short_name2, \"" << gmml::Split(parent_mono, "*").at(0) << "\", \"i\")\n";
        }
        else
            query <<  "?sn2           :monosaccharideShortName    \"" << parent_mono << "\". \n";
        query <<  "?mono2         :hasSugarName    ?sn2. \n";
        query <<  "?oligo2        :hasCore    ?mono2. \n";
    }

    query <<  "?link          :hasParent   ?oligo2. \n";
    query <<  "?link          :hasChild    ?oligo1. \n";

    if(linkage_indeces.find("?") == std::string::npos)
        query <<  "?link          :linkageIndeces   \"" << linkage_indeces << "\". \n";
    else if(linkage_indeces.compare("?-?") != 0 && linkage_indeces.find("?") != std::string::npos)
    {
        query <<  "?link          :linkageIndeces   ?linkage_indeces. \n";
        query << "FILTER regex(?linkage_indeces, \"" << gmml::Split(linkage_indeces, "?").at(0) << "\", \"i\")\n";
    }

    query <<  "?pdb           :hasOligo    ?oligo1. \n";
    query <<  "?pdb           :hasOligo    ?oligo2. \n";

    query <<  "?mono1         :hasRingAtom    ?O5. \n";
    query <<  "?O5            :ringIndex    \"6\". \n";
    query <<  "?O5            :coordinate    ?O5_crd. \n";

    query <<  "?link          :hasChildAtomLinkage    ?C1. \n";
    query <<  "?C1            :coordinate    ?C1_crd. \n";

    query <<  "?link          :hasGlycosidicLinkage    ?Ox. \n";
    query <<  "?Ox            :coordinate    ?Ox_crd. \n";

    query <<  "?link          :hasParentAtomLinkage    ?Cx. \n";
    query <<  "?Cx            :coordinate     ?Cx_crd. \n";

    query <<  "?Cx            :hasNeighbor    ?Cx_neighbor. \n";
    query <<  "?Cx_neighbor   :coordinate     ?Cx_neighbor_crd. \n";

    if(omega)
    {
        query <<  "?mono2         :hasRingAtom    ?O5_prime.\n";
        query <<  "?O5_prime      :ringIndex    \"6\".\n";
        query <<  "?O5_prime      :coordinate    ?O5_prime_crd.\n";
    }

    query << Ontology::END_WHERE_CLAUSE;
    std::stringstream curl;
    curl << Ontology::CURL_PREFIX;
    curl << Ontology::CSV_OUTPUT_FORMAT;

    curl << Ontology::DATA_STORE_ADDRESS << Ontology::QUERY_PREFIX << query.str() << Ontology::QUERY_POSTFIX << " \\>\\> result.txt";   //Has ewarning: unknown escape sequence: '\>' [enabled by default] so changed \>\> to \\>\\> by Ayush on 06/22/2017
    std::string tmp = curl.str();
    const char* cstr = tmp.c_str();
    system(cstr);
}

void Assembly::ExtractAtomCoordinatesForTorsionAnglesFromOntologyFast(std::string disaccharide_pattern/*, std::string output_file_type*/)
{
    int link_index = disaccharide_pattern.find_first_of("-");

    std::string child_mono = disaccharide_pattern.substr(0, link_index - 1); /// e.g DNeupNAca in DNeupNAca2-3DGalpb
    std::string parent_mono = disaccharide_pattern.substr(link_index + 2, disaccharide_pattern.size()); /// e.g DGalpb in DNeupNAca2-3DGalpb
    std::string linkage_indeces = disaccharide_pattern.substr(link_index - 1, 3);
    bool omega = false;
    if(linkage_indeces.find("-6") != std::string::npos || linkage_indeces.find("-7") != std::string::npos)
        omega = true;

    std::stringstream query;
    query << "sparql PREFIX : <http://gmmo.uga.edu/#> " <<
             "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> " <<
             "PREFIX owl: <http://www.w3.org/2002/07/owl#> " <<
             "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> " <<
             "PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>  " <<
             "SELECT ?pdb ?O5_crd ?C1_crd ?Ox_crd ?Cx ?Cx_crd ?Cx_neighbor ?Cx_neighbor_crd ?O5_prime_crd WHERE { " ;

    query <<  "?sn1           :monosaccharideShortName    \"" << child_mono << "\". ";
    query <<  "?mono1         :hasSugarName    ?sn1. ";
    query <<  "?oligo1        :hasCore    ?mono1. ";

    query <<  "?sn2           :monosaccharideShortName    \"" << parent_mono << "\". ";
    query <<  "?mono2         :hasSugarName    ?sn2. ";
    query <<  "?oligo2        :hasCore    ?mono2. ";

    query <<  "?link          :hasParent   ?oligo2. ";
    query <<  "?link          :hasChild    ?oligo1. ";
    query <<  "?link          :linkageIndeces   \"" << linkage_indeces << "\". ";

    query <<  "?pdb           :hasOligo    ?oligo1. ";
    query <<  "?pdb           :hasOligo    ?oligo2. ";

    query <<  "?mono1         :hasRingAtom    ?O5. ";
    query <<  "?O5            :ringIndex    \"6\". ";
    query <<  "?O5            :coordinate    ?O5_crd. ";

    query <<  "?link          :hasChildAtomLinkage    ?C1. ";
    query <<  "?C1            :coordinate    ?C1_crd. ";

    query <<  "?link          :hasGlycosidicLinkage    ?Ox. ";
    query <<  "?Ox            :coordinate    ?Ox_crd. ";

    query <<  "?link          :hasParentAtomLinkage    ?Cx. ";
    query <<  "?Cx            :coordinate     ?Cx_crd. ";

    query <<  "?Cx            :hasNeighbor    ?Cx_neighbor. ";
    query <<  "?Cx_neighbor   :coordinate     ?Cx_neighbor_crd. ";

    if(omega)
    {
        query <<  "?mono2         :hasRingAtom    ?O5_prime. ";
        query <<  "?O5_prime      :ringIndex    \"6\". ";
        query <<  "?O5_prime      :coordinate    ?O5_prime_crd. ";
    }

    query << "};";

    std::ofstream sparql;
    sparql.open("sparql.sparql", std::fstream::app);
    sparql << query.str() ;
    sparql.close();

//    std::stringstream ss;
//    ss << "/home/delaram/virtuoso-7.2.4/bin/isql 1111 dba dba \< sparql.sparql \>  result.txt";
//    std::cout << ss.str() << std::endl;
//    std::string tmp = ss.str();
//    const char* cstr = tmp.c_str();
 //   system("/home/delaram/virtuoso-7.2.4/bin/isql 1111 dba dba \< sparql.sparql \>  result.txt");
   // remove("sparql.sparql");

//    std::stringstream curl;
//    curl << Ontology::CURL_PREFIX;
//    curl << Ontology::CSV_OUTPUT_FORMAT;

//    curl << Ontology::DATA_STORE_ADDRESS << Ontology::QUERY_PREFIX << query.str() << Ontology::QUERY_POSTFIX << " >> result.txt";
//    std::string tmp = curl.str();
//    const char* cstr = tmp.c_str();
//    system(cstr);
//    std::cout << std::endl;



 /*   int link_index = disaccharide_pattern.find_first_of("-");

    std::string child_mono = disaccharide_pattern.substr(0, link_index - 1); /// e.g DNeupNAca in DNeupNAca2-3DGalpb
    std::string parent_mono = disaccharide_pattern.substr(link_index + 2, disaccharide_pattern.size()); /// e.g DGalpb in DNeupNAca2-3DGalpb
    std::string linkage_indeces = disaccharide_pattern.substr(link_index - 1, 3);
    bool omega = false;
    if(linkage_indeces.find("6") != std::string::npos || linkage_indeces.find("7") != std::string::npos
             || linkage_indeces.find("8") != std::string::npos || linkage_indeces.find("9") != std::string::npos)
        omega = true;

    std::stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE << " ?pdb ?RA1 ?RA1_crd ?SA1 ?SA1_crd ?RA2 ?RA2_crd ?SA2 ?SA2_crd "<< Ontology::WHERE_CLAUSE;
    query << "?sn1           :monosaccharideShortName    \"" << child_mono << "\".\n";
    query << "?mono1         :hasSugarName    ?sn1.\n";
    query << "?oligo1        :hasCore    ?mono1.\n";
    query << "?sn2           :monosaccharideShortName    \"" << parent_mono << "\".\n";
    query << "?mono2         :hasSugarName    ?sn2.\n";
    query << "?oligo2        :hasCore    ?mono2.\n";
    query << "?link          :hasParent   ?oligo2.\n";
    query << "?link          :hasChild    ?oligo1.\n";
    query << "?link          :linkageIndeces   \"" << linkage_indeces << "\".\n";
    query << "?pdb           :hasOligo    ?oligo1.\n";
    query << "?pdb           :hasOligo    ?oligo2.\n";

    query << "?mono1      :hasRingAtom  ?RA1.\n";
    query << "?RA1        :coordinate    ?RA1_crd.\n";

    query << "optional {?RA1        :hasSideAtom    ?SA1.\n";
    query << "?SA1        :coordinate     ?SA1_crd.}\n";


    query << "?mono2      :hasRingAtom  ?RA2.\n";
    query << "?RA2        :coordinate    ?RA2_crd.\n";

    query << "optional {?RA2        :hasSideAtom    ?SA2.\n";
    query << "?SA2        :coordinate     ?SA2_crd.}\n";


    query << Ontology::END_WHERE_CLAUSE;
    FormulateCURL(output_file_type, query.str());
    */
}

void Assembly::ExtractTorsionAnglesFromSlowQueryResult()
{

    ///Uncomment the following section and substitute std::cout with out_file in order to write the results into a file
    /*
    ///Open file to append
    std::ofstream out_file;
    out_file.open("torsions.txt", std::fstream::app);

    ///If the file is empty, write the titles
    std::ifstream check_file("torsions.txt");
    size_t out_file_size = 0;
    check_file.seekg(0,ios_base::end);
    out_file_size = check_file.tellg();
    check_file.close();
    if(out_file_size == 0)
    {
        out_file << "ϕ (O5′-C1′-Ox-Cx)" << std::endl;
        out_file << "ψ (C1′-Ox-Cx-Cx−1)" << std::endl;
        out_file << "Ω (O1-C6′-C5′-O5′)" << std::endl;
        out_file << std::left << std::setw(15) << "PDB" << std::setw(15) << "Phi Angle" << std::setw(15) << "Psi Angle" << std::setw(15) << "Omega Angle" << std::endl;
    } */


    ///Outputting the result in std out. In order to write the results into a file comment the following section
//    std::cout << "ϕ (O5′-C1′-Ox-Cx)" << std::endl;
//    std::cout << "ψ (C1′-Ox-Cx-Cx−1)" << std::endl;
//    std::cout << "Ω (O1-C6′-C5′-O5′)" << std::endl;
//    std::cout << std::left << std::setw(15) << "PDB" << std::setw(15) << "Phi Angle" << std::setw(15) << "Psi Angle" << std::setw(15) << "Omega Angle" << std::endl;

    ///Read query result file
    std::string line;
    std::ifstream in("result.txt");

    while (getline (in, line))
    {
        if(line.find("\"pdb\",\"O5_crd\",\"C1_crd\",\"Ox_crd\",\"Cx\",\"Cx_crd\",\"Cx_neighbor\",\"Cx_neighbor_crd\",\"O5_prime_crd\"") != std::string::npos)
            break;
    }
    std::string last_Cx = "";

    GeometryTopology::Coordinate* O5_crd = new GeometryTopology::Coordinate();
    GeometryTopology::Coordinate* C1_crd = new GeometryTopology::Coordinate();
    GeometryTopology::Coordinate* Ox_crd = new GeometryTopology::Coordinate();
    GeometryTopology::Coordinate* Cx_crd = new GeometryTopology::Coordinate();
    GeometryTopology::Coordinate* Cx_neighbor_crd = new GeometryTopology::Coordinate();
    GeometryTopology::Coordinate* O5_prime_crd = NULL;
    std::string Cx_id = "";
    std::string Cx_neighbor_id = "";
    double phi_angle = 0.0;
    double psi_angle = 0.0;
    double omega_angle = 0.0;

    while (getline (in, line))
    {
        std::vector<std::string> line_tokens = gmml::Split(line, ",");

        if(last_Cx.compare("") == 0 || last_Cx.compare(line_tokens.at(10)) != 0)//If it is the first line or the line with info about a new oligosaccharide
        {
            ///e.g. "0.686, -15.194, 26.371" --after splits--> x=0.686 y=-15.194 z=26.371
            O5_crd->SetX(gmml::ConvertString<double>(gmml::Split(line_tokens.at(1), "\"").at(0)));
            O5_crd->SetY(gmml::ConvertString<double>(gmml::Split(line_tokens.at(2), " ").at(0)));
            O5_crd->SetZ(gmml::ConvertString<double>(gmml::Split(line_tokens.at(3), " \"").at(0)));

            C1_crd->SetX(gmml::ConvertString<double>(gmml::Split(line_tokens.at(4), "\"").at(0)));
            C1_crd->SetY(gmml::ConvertString<double>(gmml::Split(line_tokens.at(5), " ").at(0)));
            C1_crd->SetZ(gmml::ConvertString<double>(gmml::Split(line_tokens.at(6), " \"").at(0)));

            Ox_crd->SetX(gmml::ConvertString<double>(gmml::Split(line_tokens.at(7), "\"").at(0)));
            Ox_crd->SetY(gmml::ConvertString<double>(gmml::Split(line_tokens.at(8), " ").at(0)));
            Ox_crd->SetZ(gmml::ConvertString<double>(gmml::Split(line_tokens.at(9), " \"").at(0)));

            Cx_crd->SetX(gmml::ConvertString<double>(gmml::Split(line_tokens.at(11), "\"").at(0)));
            Cx_crd->SetY(gmml::ConvertString<double>(gmml::Split(line_tokens.at(12), " ").at(0)));
            Cx_crd->SetZ(gmml::ConvertString<double>(gmml::Split(line_tokens.at(13), " \"").at(0)));

            if(line_tokens.size() > 18)
            {
                O5_prime_crd = new GeometryTopology::Coordinate();
                O5_prime_crd->SetX(gmml::ConvertString<double>(gmml::Split(line_tokens.at(18), "\"").at(0)));
                O5_prime_crd->SetY(gmml::ConvertString<double>(gmml::Split(line_tokens.at(19), " ").at(0)));
                O5_prime_crd->SetZ(gmml::ConvertString<double>(gmml::Split(line_tokens.at(20), " \"").at(0)));
            }

            Cx_id = gmml::Split(line_tokens.at(10), "#").at(1); ///e.g. spliting 5BO9_C3_4773_GAL_A_410_n_n_1 from http://gmmo.uga.edu/#5BO9_C3_4773_GAL_A_410_n_n_1
        }

        last_Cx = line_tokens.at(10);
        Cx_neighbor_id = gmml::Split(line_tokens.at(14), "#").at(1);

        if(gmml::Split(Cx_neighbor_id, "_").at(1).find("C") != std::string::npos)///If the neighbor is a carbon
        {
            /// e.g. removing C * , and ' from the atom name to get the index
            int Cx_index = gmml::ConvertString<int>(gmml::Split(gmml::Split(Cx_id, "_").at(1), "C*,\'").at(0)); ///e.g. 5BO9_C3_4773_GAL_A_410_n_n_1 -split-> C3 -split-> 3
            int Cx_neighbor_index = gmml::ConvertString<int>(gmml::Split(gmml::Split(Cx_neighbor_id, "_").at(1), "C*,\'").at(0));

            if(Cx_index > Cx_neighbor_index) ///if the neighbor is Cx-1
            {
                Cx_neighbor_crd->SetX(gmml::ConvertString<double>(gmml::Split(line_tokens.at(15), "\"").at(0)));
                Cx_neighbor_crd->SetY(gmml::ConvertString<double>(gmml::Split(line_tokens.at(16), " ").at(0)));
                Cx_neighbor_crd->SetZ(gmml::ConvertString<double>(gmml::Split(line_tokens.at(17), "\"").at(0)));

                phi_angle = CalculateTorsionAngleByCoordinates(O5_crd, C1_crd, Ox_crd, Cx_crd); /// ϕ (O5′-C1′-Ox-Cx)
                psi_angle = CalculateTorsionAngleByCoordinates(C1_crd, Ox_crd, Cx_crd,Cx_neighbor_crd); /// ψ (C1′-Ox-Cx-Cx−1)

//                std::cout << std::left << std::setw(15) << gmml::Split(gmml::Split(line_tokens.at(0), "#").at(1), "\"").at(0)
//                         << std::setw(15) << gmml::ConvertRadian2Degree(phi_angle) << std::setw(15)
//                         << gmml::ConvertRadian2Degree(psi_angle);

                if(O5_prime_crd != NULL)
                {
                    omega_angle = CalculateTorsionAngleByCoordinates(Ox_crd, Cx_crd,Cx_neighbor_crd, O5_prime_crd); /// Ω (O1-C6′-C5′-O5′)
//                    std::cout << std::setw(15) << gmml::ConvertRadian2Degree(omega_angle);
                }
//                std::cout << std::endl;
            }
        }
    }
    in.close();
}

void Assembly::ExtractTorsionAnglesFromFastQueryResult()
{

    /* Sample query result

     */

    ///Uncomment the following section and substitute std::cout with out_file in order to write the results into a file
    /*
    ///Open file to append
    std::ofstream out_file;
    out_file.open("torsions.txt", std::fstream::app);

    ///If the file is empty, write the titles
    std::ifstream check_file("torsions.txt");
    size_t out_file_size = 0;
    check_file.seekg(0,ios_base::end);
    out_file_size = check_file.tellg();
    check_file.close();
    if(out_file_size == 0)
    {
        out_file << "ϕ (O5′-C1′-Ox-Cx)" << std::endl;
        out_file << "ψ (C1′-Ox-Cx-Cx−1)" << std::endl;
        out_file << "Ω (O1-C6′-C5′-O5′)" << std::endl;
        out_file << std::left << std::setw(15) << "PDB" << std::setw(15) << "Phi Angle" << std::setw(15) << "Psi Angle" << std::setw(15) << "Omega Angle" << std::endl;
    } */


    ///Outputting the result in std out. In order to write the results into a file comment the following section
//    std::cout << "ϕ (O5′-C1′-Ox-Cx)" << std::endl;
//    std::cout << "ψ (C1′-Ox-Cx-Cx−1)" << std::endl;
//    std::cout << "Ω (O1-C6′-C5′-O5′)" << std::endl;
//    std::cout << std::left << std::setw(15) << "PDB" << std::setw(15) << "Phi Angle" << std::setw(15) << "Psi Angle" << std::setw(15) << "Omega Angle" << std::endl;

    ///Read query result file
    std::string line;
    std::ifstream in("result.txt");

    while (getline (in, line))
    {
        if(line.find("http://gmmo.uga.edu/") != std::string::npos)
            break;
    }
    std::string last_Cx = "";

    GeometryTopology::Coordinate* O5_crd = new GeometryTopology::Coordinate();
    GeometryTopology::Coordinate* C1_crd = new GeometryTopology::Coordinate();
    GeometryTopology::Coordinate* Ox_crd = new GeometryTopology::Coordinate();
    GeometryTopology::Coordinate* Cx_crd = new GeometryTopology::Coordinate();
    GeometryTopology::Coordinate* Cx_neighbor_crd = new GeometryTopology::Coordinate();
    GeometryTopology::Coordinate* O5_prime_crd = NULL;
    std::string Cx_id = "";
    std::string Cx_neighbor_id = "";
    double phi_angle = 0.0;
    double psi_angle = 0.0;
    double omega_angle = 0.0;

    do
    {
        std::vector<std::string> line_tokens = gmml::Split(line, " ");
        if(last_Cx.compare("") == 0 || last_Cx.compare(line_tokens.at(10)) != 0)//If it is the first line or the line with info about a new oligosaccharide
        {
            ///e.g. "0.686, -15.194, 26.371" --after splits--> x=0.686 y=-15.194 z=26.371

            O5_crd->SetX(gmml::ConvertString<double>(gmml::Split(line_tokens.at(1), ",").at(0)));
            O5_crd->SetY(gmml::ConvertString<double>(gmml::Split(line_tokens.at(2), ",").at(0)));
            O5_crd->SetZ(gmml::ConvertString<double>(line_tokens.at(3)));

            C1_crd->SetX(gmml::ConvertString<double>(gmml::Split(line_tokens.at(4), ",").at(0)));
            C1_crd->SetY(gmml::ConvertString<double>(gmml::Split(line_tokens.at(5), ",").at(0)));
            C1_crd->SetZ(gmml::ConvertString<double>(line_tokens.at(6)));

            Ox_crd->SetX(gmml::ConvertString<double>(gmml::Split(line_tokens.at(7), ",").at(0)));
            Ox_crd->SetY(gmml::ConvertString<double>(gmml::Split(line_tokens.at(8), ",").at(0)));
            Ox_crd->SetZ(gmml::ConvertString<double>(line_tokens.at(9)));

            Cx_crd->SetX(gmml::ConvertString<double>(gmml::Split(line_tokens.at(11), ",").at(0)));
            Cx_crd->SetY(gmml::ConvertString<double>(gmml::Split(line_tokens.at(12), ",").at(0)));
            Cx_crd->SetZ(gmml::ConvertString<double>(line_tokens.at(13)));

            if(line_tokens.at(18).compare("") != 0)
            {
                O5_prime_crd->SetX(gmml::ConvertString<double>(gmml::Split(line_tokens.at(18), ",").at(0)));
                O5_prime_crd->SetY(gmml::ConvertString<double>(gmml::Split(line_tokens.at(19), ",").at(0)));
                O5_prime_crd->SetZ(gmml::ConvertString<double>(line_tokens.at(20)));
            }

            Cx_id = gmml::Split(line_tokens.at(10), "#").at(1); ///e.g. spliting 5BO9_C3_4773_GAL_A_410_n_n_1 from http://gmmo.uga.edu/#5BO9_C3_4773_GAL_A_410_n_n_1
        }

        last_Cx = line_tokens.at(10);
        Cx_neighbor_id = gmml::Split(line_tokens.at(14), "#").at(1);

        if(gmml::Split(Cx_neighbor_id, "_").at(1).find("C") != std::string::npos)///If the neighbor is a carbon
        {
            /// e.g. removing C * , and ' from the atom name to get the index
            int Cx_index = gmml::ConvertString<int>(gmml::Split(gmml::Split(Cx_id, "_").at(1), "C*,\'").at(0)); ///e.g. 5BO9_C3_4773_GAL_A_410_n_n_1 -split-> C3 -split-> 3
            int Cx_neighbor_index = gmml::ConvertString<int>(gmml::Split(gmml::Split(Cx_neighbor_id, "_").at(1), "C*,\'").at(0));

            if(Cx_index > Cx_neighbor_index) ///if the neighbor is Cx-1
            {
                Cx_neighbor_crd->SetX(gmml::ConvertString<double>(gmml::Split(line_tokens.at(15), ",").at(0)));
                Cx_neighbor_crd->SetY(gmml::ConvertString<double>(gmml::Split(line_tokens.at(16), ",").at(0)));
                Cx_neighbor_crd->SetZ(gmml::ConvertString<double>(line_tokens.at(17)));

                phi_angle = CalculateTorsionAngleByCoordinates(O5_crd, C1_crd, Ox_crd, Cx_crd); /// ϕ (O5′-C1′-Ox-Cx)
                psi_angle = CalculateTorsionAngleByCoordinates(C1_crd, Ox_crd, Cx_crd,Cx_neighbor_crd); /// ψ (C1′-Ox-Cx-Cx−1)

//                std::cout << std::left << std::setw(15) << gmml::Split(gmml::Split(line_tokens.at(0), "#").at(1), "\"").at(0)
//                         << std::setw(15) << gmml::ConvertRadian2Degree(phi_angle) << std::setw(15)
//                         << gmml::ConvertRadian2Degree(psi_angle);

                if(O5_prime_crd != NULL)
                {
                    omega_angle = CalculateTorsionAngleByCoordinates(Ox_crd, Cx_crd,Cx_neighbor_crd, O5_prime_crd); /// Ω (O1-C6′-C5′-O5′)
//                    std::cout << std::setw(15) << gmml::ConvertRadian2Degree(omega_angle);
                }
//                std::cout << std::endl;
            }
        }
    }
    while (getline (in, line) && !line.empty());
    in.close();
}

std::string Assembly::ExtractOntologyInfoByOligosaccharideNameSequenceGF(std::string oligo_name, std::string output_file_type)
{
    std::stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE << " ?pdb ?residue_links " << Ontology::WHERE_CLAUSE;

    query << "?pdb_file     :hasOligo	?oligo.\n";
    query << "?oligo        :oligoName	\"" << oligo_name << "\"\n";
    query << "OPTIONAL { ?oligo        :oligoResidueLinks	?residue_links.}\n";
    query << "?pdb_file     :identifier   ?pdb.\n";

    ///To DO: std::string manipulation: split the oligo_name by _ and for each of them write the following to represent the names of the monos:
    //    query << "?oligo	:hasCore	?mono.\n";
    //    query << "?mono     :hasSugarName	?sn.\n";
    //    query << "?sn       :monosaccharideName 	?name.\n";
    //    query << "?sn       :monosaccharideShortName 	?short_name.\n";
    //    query << "?sn       :monosaccharideStereochemName 	?stereo_name.\n";
    //    query << "?sn       :monosaccharideStereochemShortName 	?stereo_short_name.\n";
    query << Ontology::END_WHERE_CLAUSE;
    return FormulateCURLGF(output_file_type, query.str(), Ontology::DATA_STORE_ADDRESS_GF);
}

std::string Assembly::ExtractOntologyInfoByOligosaccharideNameSequenceByRegexGF(std::string oligo_name_pattern, std::string url, std::string output_file_type)
{
    gmml::FindReplaceString(oligo_name_pattern, "[", "\\\\[");
    gmml::FindReplaceString(oligo_name_pattern, "]", "\\\\]");
    gmml::FindReplaceString(oligo_name_pattern, "-OH", "-ROH");
    if(oligo_name_pattern.compare("") == 0)
    {
//        std::cout << "Please specify the input argument. (you can use up to two * in the name pattern)" << std::endl;
        return "Please specify the input argument. (you can use up to two * in the name pattern)";
    }
    if(count(oligo_name_pattern.begin(), oligo_name_pattern.end(), '*') > 3)
    {
//        std::cout << "Wrong name pattern format. Please use only up to three * in the input argument." << std::endl;
        return "Wrong name pattern format. Please use only up to three * in the input argument.";
    }

    std::stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE << " ?pdb ?oligo_sequence ?residue_links " << Ontology::WHERE_CLAUSE;
    query << "?oligo        :oligoName	?oligo_sequence.\n";

    size_t first = oligo_name_pattern.find_first_of("*");
    size_t last = oligo_name_pattern.find_last_of("*");

    std::string filter1 = oligo_name_pattern.substr(0, first);
    std::string filter2 = oligo_name_pattern.substr(first + 1, last - 1);
    std::string filter3 = oligo_name_pattern.substr(last + 1, oligo_name_pattern.size() - 1);
    if(count(oligo_name_pattern.begin(), oligo_name_pattern.end(), '*') == 0) ///No *
        query << "?oligo	:oligoName	\"" << oligo_name_pattern << "\".\n";
    else if(count(oligo_name_pattern.begin(), oligo_name_pattern.end(), '*') == 1) ///Only one *
    {
        if(first == 0) ///* at the beginning
            query << "FILTER regex(?oligo_sequence, \"" << filter2 << "$\")\n";
        else if(first == oligo_name_pattern.size()-1) ///* at the end
            query << "FILTER regex(?oligo_sequence, \"^" << filter1 << "\")\n";
        else if(first < oligo_name_pattern.size()-1 && first > 0) ///* in the middle
            query << "FILTER regex(?oligo_sequence, \"^" << filter1 << ".+" << filter3 << "$\")\n";
    }
    else if (count(oligo_name_pattern.begin(), oligo_name_pattern.end(), '*') == 2)
    {
        if(first == 0 && last == oligo_name_pattern.size() - 1) ///* at the beginning and end
            query << "FILTER regex(?oligo_sequence, \"" << filter2 << "\")\n";
        else if(first == 0 && last < oligo_name_pattern.size() - 1)///one * at the beginning another in the middle
            query << "FILTER regex(?oligo_sequence, \"" << filter2 << ".+" << filter3 << "$\")\n";
        else if(first > 0 && last == oligo_name_pattern.size() - 1)///one * in the middle another at the end
            query << "FILTER regex(?oligo_sequence, \"^" << filter1 << ".+" << filter2 << "\")\n";
    }
    else
    {
        std::vector<std::string> pattern_tokens = gmml::Split(filter2, "*");
        query << "FILTER regex(?oligo_sequence, \"" << pattern_tokens.at(0) << ".+" << pattern_tokens.at(1) << "\")\n";
    }

    query << "?pdb_file     :hasOligo	?oligo.\n";
    query << "OPTIONAL { ?oligo        :oligoResidueLinks	?residue_links.}\n";
    query << "?pdb_file     :identifier	?pdb.\n";

    query << Ontology::END_WHERE_CLAUSE;

    return FormulateCURLGF(output_file_type, query.str(), url);
}

std::string Assembly::ExtractOntologyInfoByPDBIDGF(std::string pdb_id, std::string url, std::string output_file_type)
{
    if(pdb_id.compare("") == 0)
    {
//        std::cout << "Please specify the input argument." << std::endl;
        return "Please specify the input argument.";
    }
    std::stringstream query;
    query << Ontology::PREFIX << Ontology::SELECT_CLAUSE << " ?oligo_sequence ?residue_links " << Ontology::WHERE_CLAUSE;
    query <<  ":" << pdb_id << "    :hasOligo   ?oligo.\n";
    query << "?oligo    :oligoName 	?oligo_sequence.\n";
    query << "OPTIONAL { ?oligo	:oligoResidueLinks	?residue_links.}\n";

    //query << "?oligo	:hasCore	?mono.\n";
    //query << "?mono     :anomericStatus    ?anomeric_status.\n";

    query << Ontology::END_WHERE_CLAUSE;

    return FormulateCURLGF(output_file_type, query.str(), url);
}
