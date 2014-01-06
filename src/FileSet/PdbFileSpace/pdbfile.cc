// Author: Alireza Khatamian

#include <fstream>
#include <iostream>
#include <iomanip>

#include "../../../includes/FileSet/PdbFileSpace/pdbfile.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbfileprocessingexception.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbFile::PdbFile(const std::string &pdb_file)
{
    path_ = pdb_file;
    std::ifstream in_file;
    try
    {
        in_file.open(pdb_file.c_str());
    }
    catch(...)
    {
        throw PdbFileProcessingException(__LINE__,"File not found");
    }
    Read(in_file);
    in_file.close();            /// Close the parameter files
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
string PdbFile::GetPath()
{
    return path_;
}

PdbHeaderCard* PdbFile::GetHeader()
{
    return header_;
}

PdbTitleCard* PdbFile::GetTitle()
{
    return title_;
}

PdbCompoundCard* PdbFile::GetCompound()
{
    return compound_;
}

PdbNumModelCard* PdbFile::GetNumberOfModels()
{
    return number_of_models_;
}

PdbModelTypeCard* PdbFile::GetModelType()
{
    return model_type_;
}

PdbResidueSeqenceCard* PdbFile::GetResidueSequence()
{
    return residues_sequence_;
}

PdbResidueModificationCard* PdbFile::GetResidueModification()
{
    return residue_modification_;
}

PdbHeterogenCard* PdbFile::GetHeterogens()
{
    return heterogens_;
}

PdbHeterogenSynonymCard* PdbFile::GetHeterogenSynonyms()
{
    return heterogen_synonyms_;
}

PdbFormulaCard* PdbFile::GetFormulas()
{
    return formulas_;
}

PdbHelixCard* PdbFile::GetHelixes()
{
    return helixes_;
}

PdbSheetCard* PdbFile::GetSheets()
{
    return sheets_;
}

PdbDisulfideBondCard* PdbFile::GetDisulfideBonds()
{
    return disulfide_bonds_;
}

PdbLinkCard* PdbFile::GetLinks()
{
    return links_;
}

PdbSiteCard* PdbFile::GetSites()
{
    return sites_;
}

PdbCrystallographicCard* PdbFile::GetCrystallography()
{
    return crystallography_;
}

PdbOriginXnCard* PdbFile::GetOrigins()
{
    return origins_;
}

PdbScaleNCard* PdbFile::GetScales()
{
    return scales_;
}

PdbMatrixNCard* PdbFile::GetMatrices()
{
    return matrices_;
}

PdbModelCard* PdbFile::GetModels()
{
    return models_;
}

PdbConnectCard* PdbFile::GetConnectivities()
{
    return connectivities_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////
void PdbFile::Read(ifstream &in_file)
{
    this->ParseSections(in_file);
}

void PdbFile::ParseCards(ifstream &in_stream)
{
    string line;

    /// Unable to read file
    if (!getline(in_stream, line))
    {
        throw PdbFileProcessingException("Error reading file");
    }

    string record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "HEADER")
    {
        ParseHeaderCard(in_stream, line);
    }
    if(record_name == "OBSLTE")
    {
        ParseObsoleteCard(in_stream, line);
    }
    if(record_name == "TITLE")
    {
        ParseTitleCard(in_stream, line);
    }
    if(record_name == "SPLIT")
    {
        ParseSplitCard(in_stream, line);
    }
    if(record_name == "CAVEAT")
    {
        ParseCaveatCard(in_stream, line);
    }
    if(record_name == "COMPND")
    {
        ParseCompoundCard(in_stream, line);
    }
    if(record_name == "SOURCE")
    {
        ParseSourceCard(in_stream, line);
    }
    if(record_name == "KEYWDS")
    {
        ParseKeywordCard(in_stream, line);
    }
    if(record_name == "EXPDTA")
    {
        ParseExpirationDateCard(in_stream, line);
    }
    if(record_name == "NUMMDL")
    {
        ParseNumModelCard(in_stream, line);
    }
    if(record_name == "MDLTYP")
    {
        ParseModelTypeCard(in_stream, line);
    }
    if(record_name == "AUTHOR")
    {
        ParseAuthorCard(in_stream, line);
    }
    if(record_name == "REVDAT")
    {
        ParseRevisionDateCard(in_stream, line);
    }
    if(record_name == "SPRSDE")
    {
        ParseSupersededEntriesCard(in_stream, line);
    }
    if(record_name == "JRNL")
    {
        ParseJournalCard(in_stream, line);
    }
    if(record_name == "REMARK")
    {
        ParseRemarkCard(in_stream, line);
    }
    if(record_name.find("DBREF") != string::npos)
    {
        ParseDatabaseReferenceCard(in_stream, line);
    }
    if(record_name == "SEQADV")
    {
        ParseSequenceAdvancedCard(in_stream, line);
    }
    if(record_name == "SEQRES")
    {
        ParseSequenceResidueCard(in_stream, line);
    }
    if(record_name == "MODRES")
    {
        ParseModificationResidueCard(in_stream, line);
    }
    if(record_name == "HET")
    {
        ParseHeterogenCard(in_stream, line);
    }
    if(record_name == "HETNAM")
    {
        ParseHeterogenNameCard(in_stream, line);
    }
    if(record_name == "HETSYN")
    {
        ParseHeterogenSynonymCard(in_stream, line);
    }
    if(record_name == "FORMUL")
    {
        ParseFormulaCard(in_stream, line);
    }
    if(record_name == "HELIX")
    {
        ParseHelixCard(in_stream, line);
    }
    if(record_name == "SHEET")
    {
        ParseSheetCard(in_stream, line);
    }
    if(record_name == "SSBOND")
    {
        ParseDisulfideBondCard(in_stream, line);
    }
    if(record_name == "LINK")
    {
        ParseLinkCard(in_stream, line);
    }
    if(record_name == "CISPEP")
    {
        ParseCISPeptideCard(in_stream, line);
    }
    if(record_name == "SITE")
    {
        ParseSiteCard(in_stream, line);
    }
    if(record_name == "CRYST1")
    {
        ParseCrystallographyCard(in_stream, line);
    }
    if(record_name.find("ORIGX") != string::npos)
    {
        ParseOriginCard(in_stream, line);
    }
    if(record_name.find("SCALE") != string::npos)
    {
        ParseScaleCard(in_stream, line);
    }
    if(record_name.find("MTRIX") != string::npos)
    {
        ParseMatrixCard(in_stream, line);
    }
    if(record_name == "MODEL")
    {
        ParseModelCard(in_stream, line);
    }
    if(record_name == "CONECT")
    {
        ParseConnectivityCard(in_stream, line);
    }
    if(record_name == "MASTER")
    {
        ParseMasterCard(in_stream, line);
    }
    if(record_name == "END")
    {
        ParseEndCard(in_stream, line);
    }
}

void PdbFile::ParseHeaderCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseObsoleteCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseTitleCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseSplitCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseCaveatCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseCompoundCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseSourceCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseKeywordCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseExpirationDateCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseNumModelCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseModelTypeCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseAuthorCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseRevisionDateCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseSupersededEntriesCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseJournalCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseRemarkCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseDatabaseReferenceCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseSequenceAdvancedCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseSequenceResidueCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseModificationResidueCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseHeterogenCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseHeterogenNameCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseHeterogenSynonymCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseFormulaCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseHelixCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseSheetCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseDisulfideBondCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseLinkCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseCISPeptideCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseSiteCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseCrystallographyCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseOriginCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseScaleCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseMatrixCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseModelCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseConnectivityCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseMasterCard(std::ifstream& stream, string& line)
{

}

void PdbFile::ParseEndCard(std::ifstream& stream, string& line)
{

}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////

