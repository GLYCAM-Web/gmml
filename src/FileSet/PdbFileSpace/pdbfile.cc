// Author: Alireza Khatamian

#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <exception>

#include "../../../includes/FileSet/PdbFileSpace/pdbfile.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbheadercard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbtitlecard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbcompoundcard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbcompoundspecification.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbnummodelcard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbmodeltypecard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbresiduesequencecard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbresiduesequence.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbresiduemodificationcard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbresiduemodification.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbheterogencard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbheterogen.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbheterogennamecard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbheterogenname.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbheterogensynonymcard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbheterogensynonym.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbformulacard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbformula.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbhelixcard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbhelix.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbhelixresidue.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbsheetcard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbsheet.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbsheetstrand.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbsheetstrandresidue.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbdisulfidebondcard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbdisulfideresiduebond.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbdisulfideresidue.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdblinkcard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdblink.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdblinkresidue.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbsitecard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbsite.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbsiteresidue.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbcrystallographiccard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdboriginxncard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdboriginxn.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbscalencard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbscalen.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbmatrixncard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbmatrixn.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbmodelcard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbmodel.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbconnectcard.hpp"
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
    header_ = NULL;
    title_ = NULL;
    compound_ = NULL;	
    number_of_models_ = NULL;
    model_type_ = NULL;
	residues_sequence_ = NULL;
	residue_modification_ = NULL;
	heterogens_ = NULL;
	heterogens_name_ = NULL;
	heterogen_synonyms_ = NULL;
	formulas_ = NULL;
	helixes_ = NULL;
	sheets_ = NULL;
	disulfide_bonds_ = NULL;
	links_ = NULL;
	sites_ = NULL;
	crystallography_ = NULL;
	origins_ = NULL;
	scales_ = NULL;
	matrices_ = NULL;
	models_ = NULL;
	connectivities_ = NULL;

    std::ifstream in_file;
    try
    {
        in_file.open(pdb_file.c_str());
    }
    catch(exception &ex)
    {
        throw PdbFileProcessingException(__LINE__, "File not found");
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

PdbResidueSequenceCard* PdbFile::GetResidueSequence()
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

PdbHeterogenNameCard* PdbFile::GetHeterogensName()
{
    return heterogens_name_;
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
    this->ParseCards(in_file);
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
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "OBSLTE")
    {
        ParseObsoleteCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "TITLE")
    {
        ParseTitleCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "SPLIT")
    {
        ParseSplitCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "CAVEAT")
    {
        ParseCaveatCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "COMPND")
    {
        ParseCompoundCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "SOURCE")
    {
        ParseSourceCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "KEYWDS")
    {
        ParseKeywordCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "EXPDTA")
    {
        ParseExpirationDateCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "NUMMDL")
    {
        ParseNumModelCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "MDLTYP")
    {
        ParseModelTypeCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "AUTHOR")
    {
        ParseAuthorCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "REVDAT")
    {
        ParseRevisionDateCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "SPRSDE")
    {
        ParseSupersededEntriesCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "JRNL")
    {
        ParseJournalCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "REMARK")
    {
        ParseRemarkCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.find("DBREF") != string::npos)
    {
        ParseDatabaseReferenceCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "SEQADV")
    {
        ParseSequenceAdvancedCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "SEQRES")
    {
        ParseSequenceResidueCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "MODRES")
    {
        ParseModificationResidueCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "HET")
    {
        ParseHeterogenCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "HETNAM")
    {
        ParseHeterogenNameCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "HETSYN")
    {
        ParseHeterogenSynonymCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "FORMUL")
    {
        ParseFormulaCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "HELIX")
    {
        ParseHelixCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "SHEET")
    {
        ParseSheetCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "SSBOND")
    {
        ParseDisulfideBondCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "LINK")
    {
        ParseLinkCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "CISPEP")
    {
        ParseCISPeptideCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "SITE")
    {
        ParseSiteCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "CRYST1")
    {
        ParseCrystallographyCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.find("ORIGX") != string::npos)
    {
        ParseOriginCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.find("SCALE") != string::npos)
    {
        ParseScaleCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name.find("MTRIX") != string::npos)
    {
        ParseMatrixCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "MODEL")
    {
        ParseModelCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "CONECT")
    {
        ParseConnectivityCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "MASTER")
    {
        ParseMasterCard(in_stream, line);
    }
    record_name = line.substr(0,6);
    record_name = Trim(record_name);
    if(record_name == "END")
    {
        ParseEndCard(in_stream, line);
    }
}

void PdbFile::ParseHeaderCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "HEADER")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }

    header_ = new PdbHeaderCard(stream_block);
}

void PdbFile::ParseObsoleteCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "OBSLTE")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }
}

void PdbFile::ParseTitleCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "TITLE")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }

    title_ = new PdbTitleCard(stream_block);
}

void PdbFile::ParseSplitCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "SPLIT")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }
}

void PdbFile::ParseCaveatCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "CAVEAT")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }
}

void PdbFile::ParseCompoundCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "COMPND")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }

    compound_ = new PdbCompoundCard(stream_block);
}

void PdbFile::ParseSourceCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "SOURCE")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }
}

void PdbFile::ParseKeywordCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "KEYWDS")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }
}

void PdbFile::ParseExpirationDateCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "EXPDTA")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }
}

void PdbFile::ParseNumModelCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "NUMMDL")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }

    number_of_models_ = new PdbNumModelCard(stream_block);
}

void PdbFile::ParseModelTypeCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "MDLTYP")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }

    model_type_ = new PdbModelTypeCard(stream_block);
}

void PdbFile::ParseAuthorCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "AUTHOR")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }
}

void PdbFile::ParseRevisionDateCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "REVDAT")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }
}

void PdbFile::ParseSupersededEntriesCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "SPRSDE")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }
}

void PdbFile::ParseJournalCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "JRNL")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }
}

void PdbFile::ParseRemarkCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "REMARK")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }
}

void PdbFile::ParseDatabaseReferenceCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "DBREF")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }
}

void PdbFile::ParseSequenceAdvancedCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "SEQADV")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }
}

void PdbFile::ParseSequenceResidueCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "SEQRES")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }

    residues_sequence_ = new PdbResidueSequenceCard(stream_block);
}

void PdbFile::ParseModificationResidueCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "MODRES")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }

    residue_modification_ = new PdbResidueModificationCard(stream_block);
}

void PdbFile::ParseHeterogenCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "HET")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }

    heterogens_ = new PdbHeterogenCard(stream_block);
}

void PdbFile::ParseHeterogenNameCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "HETNAM")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }

    heterogens_name_ = new PdbHeterogenNameCard(stream_block);
}

void PdbFile::ParseHeterogenSynonymCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "HETSYN")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }

    heterogen_synonyms_ = new PdbHeterogenSynonymCard(stream_block);
}

void PdbFile::ParseFormulaCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "FORMUL")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }

    formulas_ = new PdbFormulaCard(stream_block);
}

void PdbFile::ParseHelixCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "HELIX")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }

    helixes_ = new PdbHelixCard(stream_block);
}

void PdbFile::ParseSheetCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "SHEET")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }

    sheets_ = new PdbSheetCard(stream_block);
}

void PdbFile::ParseDisulfideBondCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "SSBOND")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }

    disulfide_bonds_ = new PdbDisulfideBondCard(stream_block);
}

void PdbFile::ParseLinkCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "LINK")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }

    links_ = new PdbLinkCard(stream_block);
}

void PdbFile::ParseCISPeptideCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "CISPEP")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }
}

void PdbFile::ParseSiteCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "SITE")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }

    sites_ = new PdbSiteCard(stream_block);
}

void PdbFile::ParseCrystallographyCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "CRYST1")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }

    crystallography_ = new PdbCrystallographicCard(stream_block);
}

void PdbFile::ParseOriginCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,5);
    record_name = Trim(record_name);

    while(record_name == "ORIGX")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,5);
        record_name = Trim(record_name);
    }

    origins_ = new PdbOriginXnCard(stream_block);
}

void PdbFile::ParseScaleCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,5);
    record_name = Trim(record_name);

    while(record_name == "SCALE")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,5);
        record_name = Trim(record_name);
    }

    scales_ = new PdbScaleNCard(stream_block);
}

void PdbFile::ParseMatrixCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,5);
    record_name = Trim(record_name);

    while(record_name == "MTRIX")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,5);
        record_name = Trim(record_name);
    }

    matrices_ = new PdbMatrixNCard(stream_block);
}

void PdbFile::ParseModelCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "MODEL" || record_name == "ATOM" || record_name == "ANISOU"
          || record_name == "TER" || record_name == "HETATM" || record_name == "ENDMDL" )
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }

    models_ = new PdbModelCard(stream_block);
}

void PdbFile::ParseConnectivityCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "CONECT")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }

    connectivities_ = new PdbConnectCard(stream_block);
}

void PdbFile::ParseMasterCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "MASTER")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }
}

void PdbFile::ParseEndCard(std::ifstream& stream, string& line)
{
    stringstream stream_block;
    stream_block << line << endl;
    getline(stream, line);
    string record_name = line.substr(0,6);
    record_name = Trim(record_name);

    while(record_name == "END")
    {
        stream_block << line << endl;
        getline(stream, line);
        record_name = line.substr(0,6);
        record_name = Trim(record_name);
    }
}


void PdbFile::Write(const std::string& pdb_file)
{
    std::ofstream out_file;
    try
    {
        out_file.open(pdb_file.c_str());
    }
    catch(...)
    {
        throw PdbFileProcessingException(__LINE__,"File could not be created");
    }
    try
    {
        this->ResolveCards(out_file);
    }
    catch(...)
    {
        out_file.close();            /// Close the parameter files
    }
}

void PdbFile::ResolveCards(std::ofstream& out_stream)
{
    if(this->header_ != NULL)
    {
        this->ResolveHeaderCard(out_stream);
    }
    if(this->title_ != NULL)
    {
        this->ResolveTitleCard(out_stream);
    }
    if(this->compound_ != NULL)
    {
        this->ResolveCompoundCard(out_stream);
    }
    if(this->number_of_models_ != NULL)
    {
        this->ResolveCompoundCard(out_stream);
    }
    if(this->model_type_ != NULL)
    {
        this->ResolveNumModelCard(out_stream);
    }
    if(this->residues_sequence_ != NULL)
    {
        this->ResolveSequenceResidueCard(out_stream);
    }
    if(this->residue_modification_ != NULL)
    {
        this->ResolveModificationResidueCard(out_stream);
    }
    if(this->heterogens_ != NULL)
    {
        this->ResolveHeterogenCard(out_stream);
    }
    if(this->heterogens_name_ != NULL)
    {
        this->ResolveHeterogenNameCard(out_stream);
    }
    if(this->heterogen_synonyms_ != NULL)
    {
        this->ResolveHeterogenSynonymCard(out_stream);
    }
    if(this->formulas_ != NULL)
    {
        this->ResolveFormulaCard(out_stream);
    }
    if(this->helixes_ != NULL)
    {
        this->ResolveHelixCard(out_stream);
    }
    if(this->sheets_ != NULL)
    {
        this->ResolveSheetCard(out_stream);
    }
    if(this->disulfide_bonds_ != NULL)
    {
        this->ResolveDisulfideBondCard(out_stream);
    }
    if(this->links_ != NULL)
    {
        this->ResolveLinkCard(out_stream);
    }
    if(this->sites_ != NULL)
    {
        this->ResolveSiteCard(out_stream);
    }
    if(this->crystallography_ != NULL)
    {
        this->ResolveCrystallographyCard(out_stream);
    }
    if(this->origins_ != NULL)
    {
        this->ResolveOriginCard(out_stream);
    }
    if(this->scales_ != NULL)
    {
        this->ResolveScaleCard(out_stream);
    }
    if(this->matrices_ != NULL)
    {
        this->ResolveMatrixCard(out_stream);
    }
    if(this->models_ != NULL)
    {
        this->ResolveModelCard(out_stream);
    }
    if(this->connectivities_ != NULL)
    {
        this->ResolveConnectivityCard(out_stream);
    }
    this->ResolveEndCard(out_stream);
}

void PdbFile::ResolveHeaderCard(std::ofstream& stream)
{
    stream << left << setw(6) << header_->GetRecordName()
           << left << setw(4) << " "
           << left << setw(40) << header_->GetClassification()
           << left << setw(9) << header_->GetDepositionDate()
           << left << setw(3) << " "
           << right << setw(4) << header_->GetIdentifierCode()
           << left << setw(14) << " "
           << endl;
}

void PdbFile::ResolveObsoleteCard(std::ofstream& stream)
{
}

void PdbFile::ResolveTitleCard(std::ofstream& stream)
{
    const int MAX_TITLE_LENGTH_IN_LINE = 70;
    stream << left << setw(6) << title_->GetRecordName()
           << left << setw(2) << " ";
    if(title_->GetTitle().length() > MAX_TITLE_LENGTH_IN_LINE)
    {
        stream << right << setw(2) << " "
               << left << setw(70) << title_->GetTitle().substr(0,MAX_TITLE_LENGTH_IN_LINE)
               << endl;

        int counter = ceil((double)(title_->GetTitle().length()) / MAX_TITLE_LENGTH_IN_LINE);
        for(unsigned int i = 2; i <= counter; i++)
        {
            if(i != counter)
            {
                stream << left << setw(6) << title_->GetRecordName()
                       << left << setw(2) << " "
                       << right << setw(2) << i
                       << left << setw(70) << title_->GetTitle().substr(MAX_TITLE_LENGTH_IN_LINE*(i-1), MAX_TITLE_LENGTH_IN_LINE)
                       << endl;
            }
            else
            {
                stream << left << setw(6) << title_->GetRecordName()
                       << left << setw(2) << " "
                       << right << setw(2) << i
                       << left << setw(70) << title_->GetTitle().substr(MAX_TITLE_LENGTH_IN_LINE*(i-1), title_->GetTitle().length())
                       << endl;
            }
        }
    }
    else
    {
        stream << right << setw(2) << " "
               << left << setw(70) << title_->GetTitle()
               << endl;
    }
}

void PdbFile::ResolveSplitCard(std::ofstream& stream)
{

}

void PdbFile::CaveatCard(std::ofstream& stream)
{

}

void PdbFile::ResolveCompoundCard(std::ofstream& stream)
{
    stream << left << setw(6) << compound_->GetRecordName()
           << left << setw(1) << " "
           << right << setw(3) << " ";

    PdbCompoundCard::PdbCompoundSpecificationMap compound_specification_map = compound_->GetCompoundSpecifications();
    if((*(compound_specification_map.begin())).second->GetMoleculeId() != "")
    {
        stringstream ss;
        ss << "MOL_ID: " << (*(compound_specification_map.begin())).second->GetMoleculeId() << ";";
        stream << left << setw(70) << ss.str() << endl;
    }
    bool first = true;
    int counter = 2;
    for(PdbCompoundCard::PdbCompoundSpecificationMap::iterator it = compound_specification_map.begin(); it != compound_specification_map.end(); it++)
    {
        PdbCompoundSpecification* compound_specification = (*it).second;
        if(!first)
        {
            if(compound_specification->GetMoleculeId() != "")
            {
                stringstream ss;
                ss << " MOL_ID: " << compound_specification->GetMoleculeId() << ";";
                stream << left << setw(6) << compound_->GetRecordName()
                       << left << setw(1) << " "
                       << right << setw(3) << counter
                       << left << setw(70) << ss.str() << endl;
                counter++;
            }
        }

        /// Molecule name specification
        if(compound_specification->GetMoleculeName() != "")
        {
            stringstream ss;
            ss << " MOLECULE: " << compound_specification->GetMoleculeName() << ";";
            stream << left << setw(6) << compound_->GetRecordName()
                   << left << setw(1) << " "
                   << right << setw(3) << counter
                   << left << setw(70) << ss.str() << endl;
            counter++;
        }

        /// Molecule chain ids specification
        if(compound_specification->GetChainIds().size() > 0)
        {
            vector<string> chain_ids = compound_specification->GetChainIds();
            stringstream ss;
            ss << " CHAIN: ";
            for(vector<string>::iterator it1 = chain_ids.begin(); it1 != chain_ids.end(); it1++)
            {
                if(it1 < chain_ids.end()-1)
                    ss << (*it1) << ", ";
                else
                    ss << (*it1);
            }
            ss << ";";
            stream << left << setw(6) << compound_->GetRecordName()
                   << left << setw(1) << " "
                   << right << setw(3) << counter
                   << left << setw(70) << ss.str() << endl;
            counter++;
        }

        /// Fragment specification
        if(compound_specification->GetFragment() != "")
        {
            stringstream ss;
            ss << " FRAGMENT: " << compound_specification->GetFragment() << ";";
            stream << left << setw(6) << compound_->GetRecordName()
                   << left << setw(1) << " "
                   << right << setw(3) << counter
                   << left << setw(70) << ss.str() << endl;
            counter++;
        }

        /// Molecule synonyms specification
        if(compound_specification->GetMoleculeSynonyms().size() > 0)
        {
            vector<string> molecule_synonyms = compound_specification->GetMoleculeSynonyms();
            stringstream ss;
            ss << " SYNONYM: ";
            for(vector<string>::iterator it1 = molecule_synonyms.begin(); it1 != molecule_synonyms.end(); it1++)
            {
                if(it1 < molecule_synonyms.end()-1)
                    ss << (*it1) << ", ";
                else
                    ss << (*it1);
            }
            ss << ";";
            stream << left << setw(6) << compound_->GetRecordName()
                   << left << setw(1) << " "
                   << right << setw(3) << counter
                   << left << setw(70) << ss.str() << endl;
            counter++;
        }

        /// Enzyme commission numbers specification
        if(compound_specification->GetEnzymeCommissionNumbers().size() > 0)
        {
            vector<int> enzyme_commission_numbers = compound_specification->GetEnzymeCommissionNumbers();
            stringstream ss;
            ss << " EC: ";
            for(vector<int>::iterator it1 = enzyme_commission_numbers.begin(); it1 != enzyme_commission_numbers.end(); it1++)
            {
                if(it1 < enzyme_commission_numbers.end()-1)
                    ss << (*it1) << ", ";
                else
                    ss << (*it1);
            }
            ss << ";";
            stream << left << setw(6) << compound_->GetRecordName()
                   << left << setw(1) << " "
                   << right << setw(3) << counter
                   << left << setw(70) << ss.str() << endl;
            counter++;
        }

        /// Engineered specification
        if(compound_specification->GetIsEngineered())
        {
            stringstream ss;
            ss << " ENGINEERED: YES;";
            stream << left << setw(6) << compound_->GetRecordName()
                   << left << setw(1) << " "
                   << right << setw(3) << counter
                   << left << setw(70) << ss.str() << endl;
            counter++;
        }

        /// Mutation specification
        if(compound_specification->GetHasMutation())
        {
            stringstream ss;
            ss << " MUTATION: YES;";
            stream << left << setw(6) << compound_->GetRecordName()
                   << left << setw(1) << " "
                   << right << setw(3) << counter
                   << left << setw(70) << ss.str() << endl;
            counter++;
        }

        /// Other comments specification
        if(compound_specification->GetComments() != "")
        {
            stringstream ss;
            ss << " OTHER_DETAILS: " << compound_specification->GetComments() << ";";
            stream << left << setw(6) << compound_->GetRecordName()
                   << left << setw(1) << " "
                   << right << setw(3) << counter
                   << left << setw(70) << ss.str() << endl;
            counter++;
        }

        first = false;

    }

}

void PdbFile::ResolveSourceCard(std::ofstream& stream)
{

}

void PdbFile::ResolveKeywordCard(std::ofstream& stream)
{

}

void PdbFile::ResolveExpirationDateCard(std::ofstream& stream)
{

}

void PdbFile::ResolveNumModelCard(std::ofstream& stream)
{
    stream << left << setw(6) << number_of_models_->GetRecordName()
           << left << setw(4) << " "
           << right << setw(4) << number_of_models_->GetNumberOfModels()
           << left << setw(66) << " "
           << endl;
}

void PdbFile::ResolveModelTypeCard(std::ofstream& stream)
{
    stream << left << setw(6) << model_type_->GetRecordName()
           << left << setw(2) << " ";
    stringstream ss;
    for(vector<string>::iterator it = model_type_->GetComments().begin(); it != model_type_->GetComments().end(); it++)
    {
        if(it != model_type_->GetComments().end() - 1)
        {
            ss << (*it) << "; ";
        }
        else
        {
            ss << (*it);
        }
    }
    if(ss.str().length() > 70)
    {
        stream << right << setw(2) << " "
               << left << setw(70) << ss.str().substr(0,70)
               << endl;

        int counter = ceil((double)(ss.str().length()) / 70);
        for(unsigned int i = 2; i <= counter; i++)
        {
            if(i != counter)
            {
                stream << left << setw(6) << model_type_->GetRecordName()
                       << left << setw(2) << " "
                       << right << setw(2) << i
                       << left << setw(70) << ss.str().substr(70*(i-1), 70)
                       << endl;
            }
            else
            {
                stream << left << setw(6) << model_type_->GetRecordName()
                       << left << setw(2) << " "
                       << right << setw(2) << i
                       << left << setw(70) << ss.str().substr(70*(i-1), ss.str().length())
                       << endl;
            }
        }
    }
    else
    {
        stream << right << setw(2) << " "
               << left << setw(70) << ss.str()
               << endl;
    }
}

void PdbFile::ResolveAuthorCard(std::ofstream& stream)
{

}

void PdbFile::ResolveRevisionDateCard(std::ofstream& stream)
{

}

void PdbFile::ResolveSupersededEntriesCard(std::ofstream& stream)
{

}

void PdbFile::ResolveJournalCard(std::ofstream& stream)
{

}

void PdbFile::ResolveRemarkCard(std::ofstream& stream)
{

}

void PdbFile::ResolveDatabaseReferenceCard(std::ofstream& stream)
{

}

void PdbFile::ResolveSequenceAdvancedCard(std::ofstream& stream)
{

}

void PdbFile::ResolveSequenceResidueCard(std::ofstream& stream)
{
    PdbResidueSequenceCard::ResidueSequenceMap residue_sequence_map = residues_sequence_->GetResidueSequenceChain();
    for(PdbResidueSequenceCard::ResidueSequenceMap::iterator it = residue_sequence_map.begin(); it != residue_sequence_map.end(); it++)
    {

        PdbResidueSequence* residue_sequence = (*it).second;
        int serial_number = 1;
        const int MAX_RESIDUE_IN_SINGLE_LINE = 13;
        vector<string> residue_names = residue_sequence->GetResidueNames();
        if(residue_sequence->GetNumberOfResidues() <= MAX_RESIDUE_IN_SINGLE_LINE)
        {
            stringstream ss;
            for(vector<string>::iterator it1 = residue_names.begin(); it1 != residue_names.end(); it1++)
            {
                ss << right << setw(1) << " "
                   << right << setw(3) << (*it1);
            }
            for(unsigned int i = 0; i < MAX_RESIDUE_IN_SINGLE_LINE - residue_sequence->GetNumberOfResidues(); i++)
            {
                ss << right << setw(1) << " "
                   << right << setw(3) << " ";
            }
            stream << left << setw(6) << residues_sequence_->GetRecordName()
                   << left << setw(1) << " "
                   << right << setw(3) << serial_number
                   << right << setw(1) << " "
                   << right << setw(1) << residue_sequence->GetChainId()
                   << right << setw(1) << " "
                   << right << setw(4) << residue_sequence->GetNumberOfResidues()
                   << right << setw(1) << " "
                   << right << setw(52) << ss.str()
                   << right << setw(10) << " "
                   << endl;
        }
        else
        {
            int number_of_lines = ceil((double)(residue_sequence->GetNumberOfResidues()) / MAX_RESIDUE_IN_SINGLE_LINE);
            for(unsigned int i = 0; i < number_of_lines; i++)
            {
                stringstream ss;
                if(i != number_of_lines - 1)
                {
                    for(vector<string>::iterator it1 = residue_names.begin() + i * MAX_RESIDUE_IN_SINGLE_LINE;
                        it1 != residue_names.begin() + (i+1) * MAX_RESIDUE_IN_SINGLE_LINE; it1++)
                    {
                        ss << right << setw(1) << " "
                           << right << setw(3) << (*it1);
                    }
                }
                else
                {
                    for(vector<string>::iterator it1 = residue_names.begin() + i * MAX_RESIDUE_IN_SINGLE_LINE;
                        it1 != residue_names.end(); it1++)
                    {
                        ss << right << setw(1) << " "
                           << right << setw(3) << (*it1);
                    }
                    for(unsigned int i = 0; i < MAX_RESIDUE_IN_SINGLE_LINE - residue_sequence->GetNumberOfResidues() % MAX_RESIDUE_IN_SINGLE_LINE; i++)
                    {
                        ss << right << setw(1) << " "
                           << right << setw(3) << " ";
                    }
                }
                stream << left << setw(6) << residues_sequence_->GetRecordName()
                       << left << setw(1) << " "
                       << right << setw(3) << serial_number
                       << right << setw(1) << " "
                       << right << setw(1) << residue_sequence->GetChainId()
                       << right << setw(1) << " "
                       << right << setw(4) << residue_sequence->GetNumberOfResidues()
                       << right << setw(1) << " "
                       << right << setw(52) << ss.str()
                       << right << setw(10) << " "
                       << endl;
                serial_number++;
            }
        }

    }
}

void PdbFile::ResolveModificationResidueCard(std::ofstream& stream)
{
    PdbResidueModificationCard::ResidueModificationMap residue_modification_map = residue_modification_->GetResidueModifications();
    for(PdbResidueModificationCard::ResidueModificationMap::iterator it = residue_modification_map.begin(); it != residue_modification_map.end(); it++)
    {
        PdbResidueModification* residue_modification = (*it).second;
        stream << left << setw(6) << residue_modification_->GetRecordName()
               << left << setw(1) << " "
               << right << setw(4) << residue_modification->GetIdCode()
               << left << setw(1) << " "
               << right << setw(3) << residue_modification->GetResidueName()
               << left << setw(1) << " "
               << right << setw(1) << residue_modification->GetChainIdentifier()
               << left << setw(1) << " "
               << right << setw(4) << residue_modification->GetSequenceNumber()
               << right << setw(1) << residue_modification->GetInsertionCode()
               << left << setw(1) << " "
               << right << setw(3) << residue_modification->GetStandardResidueName()
               << left << setw(2) << " "
               << left << setw(41) << residue_modification->GetDscr()
               << left << setw(10) << " "
               << endl;
    }
}

void PdbFile::ResolveHeterogenCard(std::ofstream& stream)
{
    PdbHeterogenCard::HeterogenMap heterogen_map = heterogens_->GetHeterogens();
    for(PdbHeterogenCard::HeterogenMap::iterator it = heterogen_map.begin(); it != heterogen_map.end(); it++)
    {

        PdbHeterogen* heterogen = (*it).second;
        stream << left << setw(6) << heterogens_->GetRecordName()
               << left << setw(1) << " "
               << right << setw(3) << heterogen->GetHeterogenId()
               << left << setw(2) << " "
               << right << setw(1) << heterogen->GetChainIdentifier()
               << right << setw(4) << heterogen->GetSequenceNumber()
               << right << setw(1) << heterogen->GetInsertionCode()
               << left << setw(2) << " "
               << right << setw(5) << heterogen->GetNumberOfHeterogenAtoms()
               << left << setw(5) << " "
               << left << setw(40) << heterogen->GetDscr()
               << left << setw(10) << " "
               << endl;
    }
}

void PdbFile::ResolveHeterogenNameCard(std::ofstream& stream)
{
    PdbHeterogenNameCard::HeterogenNameMap heterogen_name_map = heterogens_name_->GetHeterogenNames();
    for(PdbHeterogenNameCard::HeterogenNameMap::iterator it = heterogen_name_map.begin(); it != heterogen_name_map.end(); it++)
    {
        PdbHeterogenName* heterogen_name = (*it).second;
        const int MAX_NAME_LENGTH_IN_LINE = 55;
        if(heterogen_name->GetHeterogenName().length() > MAX_NAME_LENGTH_IN_LINE)
        {
            stream << left << setw(6) << heterogens_name_->GetRecordName()
                   << left << setw(2) << " "
                   << right << setw(2) << " "
                   << left << setw(1) << " "
                   << right << setw(3) << heterogen_name->GetHeterogenIdentifier()
                   << left << setw(1) << " "
                   << left << setw(55) << heterogen_name->GetHeterogenName().substr(0,MAX_NAME_LENGTH_IN_LINE)
                   << left << setw(10) << " "
                   << endl;
            int counter = ceil((double)(heterogen_name->GetHeterogenName().length()) / MAX_NAME_LENGTH_IN_LINE);
            for(unsigned int i = 2; i <= counter; i++)
            {
                if(i != counter)
                {
                    stream << left << setw(6) << heterogens_name_->GetRecordName()
                           << left << setw(2) << " "
                           << right << setw(2) << i
                           << left << setw(1) << " "
                           << right << setw(3) << heterogen_name->GetHeterogenIdentifier()
                           << left << setw(1) << " "
                           << left << setw(55) << heterogen_name->GetHeterogenName().substr(MAX_NAME_LENGTH_IN_LINE*(i-1),MAX_NAME_LENGTH_IN_LINE)
                           << left << setw(10) << " "
                           << endl;
                }
                else
                {
                    stream << left << setw(6) << heterogens_name_->GetRecordName()
                           << left << setw(2) << " "
                           << right << setw(2) << i
                           << left << setw(1) << " "
                           << right << setw(3) << heterogen_name->GetHeterogenIdentifier()
                           << left << setw(1) << " "
                           << left << setw(55) << heterogen_name->GetHeterogenName().substr(MAX_NAME_LENGTH_IN_LINE*(i-1),heterogen_name->GetHeterogenName().length())
                           << left << setw(10) << " "
                           << endl;
                }
            }
        }
        else
        {
            stream << left << setw(6) << heterogens_name_->GetRecordName()
                   << left << setw(2) << " "
                   << right << setw(2) << " "
                   << left << setw(1) << " "
                   << right << setw(3) << heterogen_name->GetHeterogenIdentifier()
                   << left << setw(1) << " "
                   << left << setw(55) << heterogen_name->GetHeterogenName()
                   << left << setw(10) << " "
                   << endl;
        }
    }
}

void PdbFile::ResolveHeterogenSynonymCard(std::ofstream& stream)
{
    PdbHeterogenSynonymCard::HeterogenSynonymMap heterogen_synonym_map = heterogen_synonyms_->GetHeterogensSynonyms();
    for(PdbHeterogenSynonymCard::HeterogenSynonymMap::iterator it = heterogen_synonym_map.begin(); it != heterogen_synonym_map.end(); it++)
    {
        PdbHeterogenSynonym* heterogen_synonym = (*it).second;
        stringstream ss;
        vector<string> synonyms = heterogen_synonym->GetHeterogenSynonyms();
        for(vector<string>::iterator it = synonyms.begin(); it != synonyms.end(); it++)
        {
            if(it != synonyms.end() - 1)
            {
                ss << (*it) << "; ";
            }
            else
            {
                ss << (*it) << ";";
            }
        }
        const int MAX_SYNONYM_LENGTH_IN_LINE = 55;
        if(ss.str().length() > MAX_SYNONYM_LENGTH_IN_LINE)
        {
            stream << left << setw(6) << heterogen_synonyms_->GetRecordName()
                   << left << setw(2) << " "
                   << right << setw(2) << " "
                   << left << setw(1) << " "
                   << right << setw(3) << heterogen_synonym->GetHeterogenIdentifier()
                   << left << setw(1) << " "
                   << left << setw(55) << ss.str().substr(0,MAX_SYNONYM_LENGTH_IN_LINE)
                   << left << setw(10) << " "
                   << endl;
            int counter = ceil((double)(ss.str().length()) / MAX_SYNONYM_LENGTH_IN_LINE);
            for(unsigned int i = 2; i <= counter; i++)
            {
                if(i != counter)
                {
                    stream << left << setw(6) << heterogen_synonyms_->GetRecordName()
                           << left << setw(2) << " "
                           << right << setw(2) << i
                           << left << setw(1) << " "
                           << right << setw(3) << heterogen_synonym->GetHeterogenIdentifier()
                           << left << setw(1) << " "
                           << left << setw(55) << ss.str().substr(MAX_SYNONYM_LENGTH_IN_LINE*(i-1),MAX_SYNONYM_LENGTH_IN_LINE)
                           << left << setw(10) << " "
                           << endl;
                }
                else
                {
                    stream << left << setw(6) << heterogen_synonyms_->GetRecordName()
                           << left << setw(2) << " "
                           << right << setw(2) << i
                           << left << setw(1) << " "
                           << right << setw(3) << heterogen_synonym->GetHeterogenIdentifier()
                           << left << setw(1) << " "
                           << left << setw(55) << ss.str().substr(MAX_SYNONYM_LENGTH_IN_LINE*(i-1),ss.str().length())
                           << left << setw(10) << " "
                           << endl;
                }
            }
        }
        else
        {
            stream << left << setw(6) << heterogen_synonyms_->GetRecordName()
                   << left << setw(2) << " "
                   << right << setw(2) << " "
                   << left << setw(1) << " "
                   << right << setw(3) << heterogen_synonym->GetHeterogenIdentifier()
                   << left << setw(1) << " "
                   << left << setw(55) << ss.str()
                   << left << setw(10) << " "
                   << endl;
        }
    }
}

void PdbFile::ResolveFormulaCard(std::ofstream& stream)
{
    PdbFormulaCard::FormulaMap formula_map = formulas_->GetFormulas();
    for(PdbFormulaCard::FormulaMap::iterator it = formula_map.begin(); it != formula_map.end(); it++)
    {
        PdbFormula* formula = (*it).second;
        const int MAX_NAME_LENGTH_IN_LINE = 51;
        if(formula->GetChemicalFormula().length() > MAX_NAME_LENGTH_IN_LINE)
        {
            stream << left << setw(6) << formulas_->GetRecordName()
                   << left << setw(2) << " "
                   << right << setw(2) << formula->GetComponentNumber()
                   << left << setw(2) << " "
                   << right << setw(3) << formula->GetHeterogenIdentifier()
                   << left << setw(1) << " "
                   << right << setw(2) << " "
                   << left << setw(1) << " "
                   << left << setw(51) << formula->GetChemicalFormula().substr(0,MAX_NAME_LENGTH_IN_LINE)
                   << left << setw(10) << " "
                   << endl;
            int counter = ceil((double)(formula->GetChemicalFormula().length()) / MAX_NAME_LENGTH_IN_LINE);
            for(unsigned int i = 2; i <= counter; i++)
            {
                if(i != counter)
                {
                    stream << left << setw(6) << formulas_->GetRecordName()
                           << left << setw(2) << " "
                           << right << setw(2) << formula->GetComponentNumber()
                           << left << setw(2) << " "
                           << right << setw(3) << formula->GetHeterogenIdentifier()
                           << left << setw(1) << " "
                           << right << setw(2) << i
                           << left << setw(1) << " "
                           << left << setw(51) << formula->GetChemicalFormula().substr(MAX_NAME_LENGTH_IN_LINE*(i-1),MAX_NAME_LENGTH_IN_LINE)
                           << left << setw(10) << " "
                           << endl;
                }
                else
                {
                    stream << left << setw(6) << formulas_->GetRecordName()
                           << left << setw(2) << " "
                           << right << setw(2) << formula->GetComponentNumber()
                           << left << setw(2) << " "
                           << right << setw(3) << formula->GetHeterogenIdentifier()
                           << left << setw(1) << " "
                           << right << setw(2) << i
                           << left << setw(1) << " "
                           << left << setw(51) << formula->GetChemicalFormula().substr(MAX_NAME_LENGTH_IN_LINE*(i-1),formula->GetChemicalFormula().length())
                           << left << setw(10) << " "
                           << endl;
                }
            }
        }
        else
        {
            stream << left << setw(6) << formulas_->GetRecordName()
                   << left << setw(2) << " "
                   << right << setw(2) << formula->GetComponentNumber()
                   << left << setw(2) << " "
                   << right << setw(3) << formula->GetHeterogenIdentifier()
                   << left << setw(1) << " "
                   << right << setw(2) << " "
                   << left << setw(1) << " "
                   << left << setw(51) << formula->GetChemicalFormula().substr(0,MAX_NAME_LENGTH_IN_LINE)
                   << left << setw(10) << " "
                   << endl;
        }
    }
}

void PdbFile::ResolveHelixCard(std::ofstream& stream)
{
    PdbHelixCard::HelixMap helix_map = helixes_->GetHelixes();
    int counter = helix_map.size();
    int serial_number = 1;
    while(serial_number <= counter)
    {
        for(PdbHelixCard::HelixMap::iterator it = helix_map.begin(); it != helix_map.end(); it++)
        {

            PdbHelix* helix = (*it).second;
            PdbHelix::HelixResidueVector helix_residues = helix->GetHelixResidues();
            if(helix->GetHelixSerialNumber() == serial_number)
            {
                stream << left << setw(6) << helixes_->GetRecordName()
                       << left << setw(1) << " "
                       << right << setw(3) << helix->GetHelixSerialNumber()
                       << left << setw(1) << " "
                       << right << setw(3) << helix->GetHelixId()
                       << left << setw(1) << " "
                       << right << setw(3) << helix_residues.at(0)->GetResidueName()
                       << left << setw(1) << " "
                       << right << setw(1) << helix_residues.at(0)->GetResidueChainId()
                       << left << setw(1) << " "
                       << right << setw(4) << helix_residues.at(0)->GetResidueSequenceNumber()
                       << right << setw(1) << helix_residues.at(0)->GetResidueInsertionCode()
                       << left << setw(1) << " "
                       << right << setw(3) << helix_residues.at(1)->GetResidueName()
                       << left << setw(1) << " "
                       << right << setw(1) << helix_residues.at(1)->GetResidueChainId()
                       << left << setw(1) << " "
                       << right << setw(4) << helix_residues.at(1)->GetResidueSequenceNumber()
                       << right << setw(1) << helix_residues.at(1)->GetResidueInsertionCode()
                       << right << setw(2) << helix->GetHelixClass()
                       << left << setw(30) << helix->GetComment()
                       << left << setw(1) << " "
                       << right << setw(5) << helix->GetHelixLength()
                       << left << setw(4) << " "
                       << endl;
                break;
            }
        }
        serial_number++;
    }
}

void PdbFile::ResolveSheetCard(std::ofstream& stream)
{
    PdbSheetCard::SheetMap sheet_map = sheets_->GetSheets();
    for(PdbSheetCard::SheetMap::iterator it = sheet_map.begin(); it != sheet_map.end(); it++)
    {
        PdbSheet* sheet = (*it).second;
        PdbSheet::SheetStrandVector strands = sheet->GetStrands();
        int serial_number = 1;
        for(PdbSheet::SheetStrandVector::iterator it1 = strands.begin(); it1 != strands.end(); it1++)
        {
            PdbSheetStrand* strand = (*it1);
            PdbSheetStrand::SheetStrandResidueVector strand_residues = strand->GetStrandResidues();
            if(strand->GetSense() != 0)
            {
                stream << left << setw(6) << sheets_->GetRecordName()
                       << left << setw(1) << " "
                       << right << setw(3) << serial_number
                       << left << setw(1) << " "
                       << right << setw(3) << sheet->GetSheetId()
                       << right << setw(2) << sheet->GetNumberOfStrands()
                       << left << setw(1) << " "
                       << right << setw(3) << strand_residues.at(0)->GetResidueName()
                       << left << setw(1) << " "
                       << right << setw(1) << strand_residues.at(0)->GetResidueChainId()
                       << right << setw(4) << strand_residues.at(0)->GetResidueSequenceNumber()
                       << right << setw(1) << strand_residues.at(0)->GetResidueInsertionCode()
                       << left << setw(1) << " "
                       << right << setw(3) << strand_residues.at(1)->GetResidueName()
                       << left << setw(1) << " "
                       << right << setw(1) << strand_residues.at(1)->GetResidueChainId()
                       << right << setw(4) << strand_residues.at(1)->GetResidueSequenceNumber()
                       << right << setw(1) << strand_residues.at(1)->GetResidueInsertionCode()
                       << right << setw(2) << strand->GetSense()
                       << left << setw(1) << " "
                       << left << setw(4) << strand->GetCurrentAtom()
                       << right << setw(3) << strand_residues.at(2)->GetResidueName()
                       << left << setw(1) << " "
                       << right << setw(1) << strand_residues.at(2)->GetResidueChainId()
                       << right << setw(4) << strand_residues.at(2)->GetResidueSequenceNumber()
                       << right << setw(1) << strand_residues.at(2)->GetResidueInsertionCode()
                       << left << setw(1) << " "
                       << left << setw(4) << strand->GetPreviousAtom()
                       << right << setw(3) << strand_residues.at(3)->GetResidueName()
                       << left << setw(1) << " "
                       << right << setw(1) << strand_residues.at(3)->GetResidueChainId()
                       << right << setw(4) << strand_residues.at(3)->GetResidueSequenceNumber()
                       << right << setw(1) << strand_residues.at(3)->GetResidueInsertionCode()
                       << left << setw(10) << " "
                       << endl;
            }
            else
            {
                stream << left << setw(6) << sheets_->GetRecordName()
                       << left << setw(1) << " "
                       << right << setw(3) << serial_number
                       << left << setw(1) << " "
                       << right << setw(3) << sheet->GetSheetId()
                       << right << setw(2) << sheet->GetNumberOfStrands()
                       << left << setw(1) << " "
                       << right << setw(3) << strand_residues.at(0)->GetResidueName()
                       << left << setw(1) << " "
                       << right << setw(1) << strand_residues.at(0)->GetResidueChainId()
                       << right << setw(4) << strand_residues.at(0)->GetResidueSequenceNumber()
                       << right << setw(1) << strand_residues.at(0)->GetResidueInsertionCode()
                       << left << setw(1) << " "
                       << right << setw(3) << strand_residues.at(1)->GetResidueName()
                       << left << setw(1) << " "
                       << right << setw(1) << strand_residues.at(1)->GetResidueChainId()
                       << right << setw(4) << strand_residues.at(1)->GetResidueSequenceNumber()
                       << right << setw(1) << strand_residues.at(1)->GetResidueInsertionCode()
                       << right << setw(2) << strand->GetSense()
                       << left << setw(40) << " "
                       << endl;
            }
            serial_number++;
        }

    }
}

void PdbFile::ResolveDisulfideBondCard(std::ofstream& stream)
{
    PdbDisulfideBondCard::DisulfideResidueBondMap disulfide_bond_map = disulfide_bonds_->GetDisulfideResidueBonds();
    for(PdbDisulfideBondCard::DisulfideResidueBondMap::iterator it = disulfide_bond_map.begin(); it != disulfide_bond_map.end(); it++)
    {
        PdbDisulfideResidueBond* disulfide_bonds = (*it).second;
        PdbDisulfideResidueBond::DisulfideResidueVector disulfide_bonds_residues = disulfide_bonds->GetResidues();
        stream << left << setw(6) << disulfide_bonds_->GetRecordName()
               << left << setw(1) << " "
               << right << setw(3) << disulfide_bonds->GetSerialNumber()
               << left << setw(1) << " "
               << right << setw(3) << disulfide_bonds_residues.at(0)->GetResidueName()
               << left << setw(1) << " "
               << right << setw(1) << disulfide_bonds_residues.at(0)->GetResidueChainIdentifier()
               << left << setw(1) << " "
               << right << setw(4) << disulfide_bonds_residues.at(0)->GetResidueSequenceNumber()
               << right << setw(1) << disulfide_bonds_residues.at(0)->GetResidueInsertionCode()
               << left << setw(3) << " "
               << right << setw(3) << disulfide_bonds_residues.at(1)->GetResidueName()
               << left << setw(1) << " "
               << right << setw(1) << disulfide_bonds_residues.at(1)->GetResidueChainIdentifier()
               << left << setw(1) << " "
               << right << setw(4) << disulfide_bonds_residues.at(1)->GetResidueSequenceNumber()
               << right << setw(1) << disulfide_bonds_residues.at(1)->GetResidueInsertionCode()
               << left << setw(23) << " "
               << right << setw(6) << disulfide_bonds_residues.at(0)->GetSymmetryOperator()
               << left << setw(1) << " "
               << right << setw(6) << disulfide_bonds_residues.at(1)->GetSymmetryOperator()
               << left << setw(1) << " "
               << right << setw(5) << fixed << setprecision(2) << disulfide_bonds->GetBondLength()
               << left << setw(2) << " "
               << endl;
    }
}

void PdbFile::ResolveLinkCard(std::ofstream& stream)
{
    PdbLinkCard::LinkVector links = links_->GetResidueLinks();
    for(PdbLinkCard::LinkVector::iterator it = links.begin(); it != links.end(); it++)
    {
        PdbLink* link = (*it);
        PdbLink::LinkResidueVector link_residues = link->GetResidues();
        stream << left << setw(6) << links_->GetRecordName()
               << left << setw(6) << " "
               << left << setw(4) << link_residues.at(0)->GetAtomName()
               << right << setw(1) << link_residues.at(0)->GetAlternateLocationIndicator()
               << right << setw(3) << link_residues.at(0)->GetResidueName()
               << left << setw(1) << " "
               << right << setw(1) << link_residues.at(0)->GetResidueChainIdentifier()
               << right << setw(4) << link_residues.at(0)->GetResidueSequenceNumber()
               << right << setw(1) << link_residues.at(0)->GetResidueInsertionCode()
               << left << setw(15) << " "
               << left << setw(4) << link_residues.at(1)->GetAtomName()
               << right << setw(1) << link_residues.at(1)->GetAlternateLocationIndicator()
               << right << setw(3) << link_residues.at(1)->GetResidueName()
               << left << setw(1) << " "
               << right << setw(1) << link_residues.at(1)->GetResidueChainIdentifier()
               << right << setw(4) << link_residues.at(1)->GetResidueSequenceNumber()
               << right << setw(1) << link_residues.at(1)->GetResidueInsertionCode()
               << left << setw(2) << " "
               << right << setw(6) << link_residues.at(0)->GetSymmetryOperator()
               << left << setw(1) << " "
               << right << setw(6) << link_residues.at(1)->GetSymmetryOperator()
               << left << setw(1) << " "
               << right << setw(5) << fixed << setprecision(2) << link->GetLinkLength()
               << left << setw(2) << " "
               << endl;
    }
}

void PdbFile::ResolveCISPeptideCard(std::ofstream& stream)
{

}

void PdbFile::ResolveSiteCard(std::ofstream& stream)
{
    PdbSiteCard::PdbSiteMap site_map = sites_->GetResidueSites();
    for(PdbSiteCard::PdbSiteMap::iterator it = site_map.begin(); it != site_map.end(); it++)
    {
        PdbSite* site = (*it).second;
        PdbSite::SiteResidueVector site_residues = site->GetResidues();
        const int MAX_RESIDUE_IN_LINE = 4;
        const int RESIDUE_LENGHT_IN_LINE = 11;
        int number_of_residues = site_residues.size();

        int sequence_number = 1;
        if(number_of_residues > MAX_RESIDUE_IN_LINE)
        {
            int number_of_lines = ceil((double)(number_of_residues) / MAX_RESIDUE_IN_LINE);
            while(sequence_number <= number_of_lines)
            {
                if(sequence_number != number_of_lines)
                {
                    stringstream ss;
                    for(PdbSite::SiteResidueVector::iterator it1 = site_residues.begin()+(sequence_number - 1)*MAX_RESIDUE_IN_LINE;
                        it1 != site_residues.begin()+(sequence_number)*MAX_RESIDUE_IN_LINE; it1++)
                    {
                        PdbSiteResidue* residue = (*it1);
                        ss << left << setw(1) << " "
                           << right << setw(3) << residue->GetResidueName()
                           << left << setw(1) << " "
                           << right << setw(1) << residue->GetResidueChainId()
                           << right << setw(4) << residue->GetresidueSequenceNumber()
                           << right << setw(1) << residue->GetResidueInsertionCode();
                    }
                    ss << left << setw(19) << " ";
                    stream << left << setw(6) << sites_->GetRecordName()
                           << left << setw(1) << " "
                           << right << setw(3) << sequence_number
                           << left << setw(1) << " "
                           << right << setw(3) << site->GetSiteId()
                           << left << setw(1) << " "
                           << right << setw(2) << site->GetNumberOfResidues()
                           << left << setw(63) << ss.str()
                           << endl;
                }
                else
                {
                    stringstream ss;
                    for(PdbSite::SiteResidueVector::iterator it1 = site_residues.begin()+(sequence_number - 1)*MAX_RESIDUE_IN_LINE;
                        it1 != site_residues.end(); it1++)
                    {
                        PdbSiteResidue* residue = (*it1);
                        ss << left << setw(1) << " "
                           << right << setw(3) << residue->GetResidueName()
                           << left << setw(1) << " "
                           << right << setw(1) << residue->GetResidueChainId()
                           << right << setw(4) << residue->GetresidueSequenceNumber()
                           << right << setw(1) << residue->GetResidueInsertionCode();
                    }
                    ss << left << setw((sequence_number*MAX_RESIDUE_IN_LINE-number_of_residues)*RESIDUE_LENGHT_IN_LINE) << " ";
                    ss << left << setw(19) << " ";
                    stream << left << setw(6) << sites_->GetRecordName()
                           << left << setw(1) << " "
                           << right << setw(3) << sequence_number
                           << left << setw(1) << " "
                           << right << setw(3) << site->GetSiteId()
                           << left << setw(1) << " "
                           << right << setw(2) << site->GetNumberOfResidues()
                           << left << setw(63) << ss.str()
                           << endl;
                }
                sequence_number++;
            }

        }
        else
        {
            stringstream ss;
            for(PdbSite::SiteResidueVector::iterator it1 = site_residues.begin(); it1 != site_residues.end(); it1++)
            {
                PdbSiteResidue* residue = (*it1);
                ss << left << setw(1) << " "
                   << right << setw(3) << residue->GetResidueName()
                   << left << setw(1) << " "
                   << right << setw(1) << residue->GetResidueChainId()
                   << right << setw(4) << residue->GetresidueSequenceNumber()
                   << right << setw(1) << residue->GetResidueInsertionCode();
            }
            ss << left << setw((MAX_RESIDUE_IN_LINE-number_of_residues)*RESIDUE_LENGHT_IN_LINE) << " ";
            ss << left << setw(19) << " ";

            stream << left << setw(6) << sites_->GetRecordName()
                   << left << setw(1) << " "
                   << right << setw(3) << sequence_number
                   << left << setw(1) << " "
                   << right << setw(3) << site->GetSiteId()
                   << left << setw(1) << " "
                   << right << setw(2) << site->GetNumberOfResidues()
                   << left << setw(63) << ss.str()
                   << endl;
            sequence_number++;
        }
    }
}

void PdbFile::ResolveCrystallographyCard(std::ofstream& stream)
{
    stream << left << setw(6) << crystallography_->GetRecordName()
           << right << setw(9) << fixed << setprecision(3) << crystallography_->GetA()
           << right << setw(9) << fixed << setprecision(3) << crystallography_->GetB()
           << right << setw(9) << fixed << setprecision(3) << crystallography_->GetC()
           << right << setw(7) << fixed << setprecision(2) << crystallography_->GetAlpha()
           << right << setw(7) << fixed << setprecision(2) << crystallography_->GetBeta()
           << right << setw(7) << fixed << setprecision(2) << crystallography_->GetGamma()
           << left << setw(1) << " "
           << left << setw(11) << crystallography_->GetSpaceGroup()
           << right << setw(4) << crystallography_->GetZValue()
           << left << setw(10) << " "
           << endl;
}

void PdbFile::ResolveOriginCard(std::ofstream& stream)
{
    PdbOriginXnCard::OriginXnVector origins = origins_->GetOriginXN();
    for(PdbOriginXnCard::OriginXnVector::iterator it = origins.begin(); it != origins.end(); it++)
    {
        PdbOriginXn* origin = (*it);
        stringstream ss;
        ss << origin->GetRecordName() << origin->GetN();
        stream << left << setw(6) << ss.str()
               << left << setw(4) << " "
               << right << setw(10) << fixed << setprecision(6) << origin->GetOrigin().GetX()
               << right << setw(10) << fixed << setprecision(6) << origin->GetOrigin().GetY()
               << right << setw(10) << fixed << setprecision(6) << origin->GetOrigin().GetZ()
               << left << setw(5) << " "
               << right << setw(10) << fixed << setprecision(5) << origin->GetT()
               << left << setw(25) << " "
               << endl;
    }
}

void PdbFile::ResolveScaleCard(std::ofstream& stream)
{
    PdbScaleNCard::ScaleNVector scales = scales_->GetScaleN();
    for(PdbScaleNCard::ScaleNVector::iterator it = scales.begin(); it != scales.end(); it++)
    {
        PdbScaleN* scale = (*it);
        stringstream ss;
        ss << scale->GetRecordName() << scale->GetN();
        stream << left << setw(6) << ss.str()
               << left << setw(4) << " "
               << right << setw(10) << fixed << setprecision(6) << scale->GetScaleVector().GetX()
               << right << setw(10) << fixed << setprecision(6) << scale->GetScaleVector().GetY()
               << right << setw(10) << fixed << setprecision(6) << scale->GetScaleVector().GetZ()
               << left << setw(5) << " "
               << right << setw(10) << fixed << setprecision(5) << scale->GetU()
               << left << setw(25) << " "
               << endl;
    }
}

void PdbFile::ResolveMatrixCard(std::ofstream& stream)
{
    PdbMatrixNCard::MatrixNVectorVector matrices = matrices_->GetMatrixN();    
    int number_of_matrix_entries = matrices.at(0).size();
    for(unsigned int i = 0; i < number_of_matrix_entries; i++)
    {
        for(unsigned int j = 0; j < 3; j++)
        {
	    PdbMatrixNCard::MatrixNVector matrix_vector = matrices.at(j);
            PdbMatrixN* matrix = matrix_vector.at(i);
            stringstream ss;
            ss << matrix->GetRecordName() << matrix->GetN();
            stream << left << setw(6) << ss.str()
                   << left << setw(1) << " "
                   << right << setw(3) << matrix->GetSerialNumber()
                   << right << setw(10) << fixed << setprecision(6) << matrix->GetTransformationVector().GetX()
                   << right << setw(10) << fixed << setprecision(6) << matrix->GetTransformationVector().GetY()
                   << right << setw(10) << fixed << setprecision(6) << matrix->GetTransformationVector().GetZ()
                   << left << setw(5) << " "
                   << right << setw(10) << fixed << setprecision(5) << matrix->GetV()
                   << left << setw(4) << " "
                   << right << setw(1) << matrix->GetIGiven()
                   << left << setw(20) << " "
                   << endl;
        }
    }

}

void PdbFile::ResolveModelCard(std::ofstream& stream)
{
//    PdbModelCard::PdbModelMap models = models_->GetModels();
//    int number_of_models = models.size();
//    if(number_of_models == 1)
//    {
//        for(PdbModelCard::PdbModelMap::iterator it = models.begin(); it != models.end(); it++)
//        {
//            PdbModel* model = (*it).second;
//            PdbModelResidueSet* residue_set = model->GetModelResidueSet();
//            PdbModelResidueSet::AtomCardVector atoms = residue_set->GetAtoms();
//            for(PdbModelResidueSet::AtomCardVector::iterator it1 = atoms.begin(); it1 != atoms.end(); it1++)
//            {
//                PdbAtomCard* atom = (*it);
//            }
//            PdbModelResidueSet::HeterogenAtomCardVector heterogen_atoms = residue_set->GetHeterogenAtoms();


//        }
//    }
//    else
//    {

//    }
}

void PdbFile::ResolveConnectivityCard(std::ofstream& stream)
{

}

void PdbFile::ResolveMasterCard(std::ofstream& stream)
{

}

void PdbFile::ResolveEndCard(std::ofstream& stream)
{
    stream << left << setw(6) << "END" << left << setw(74) << " " << endl;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////

void PdbFile::Print(ostream &out)
{
    if(header_ != NULL)
    {
        out << "******************************* HEADER *******************************" << endl;
        header_->Print(out);
    }
    if(title_ != NULL)
    {
        out << "******************************** TITLE *******************************" << endl;
        title_->Print(out);
    }
    if(compound_ != NULL)
    {
        out << "****************************** COMPOUND ******************************" << endl;
        compound_->Print(out);
    }
    if(number_of_models_ != NULL)
    {
        out << "************************** NUMBER OF MODELS **************************" << endl;
        number_of_models_->Print(out);
    }
    if(model_type_ != NULL)
    {
        out << "***************************** MODEL TYPE *****************************" << endl;
        model_type_->Print(out);
    }
    if(residues_sequence_ != NULL)
    {
        out << "************************** RESIDUE SEQUENCE **************************" << endl;
        residues_sequence_->Print(out);
    }
    if(residue_modification_ != NULL)
    {
        out << "************************ RESIDUE MODIFICATION ************************" << endl;
        residue_modification_->Print(out);
    }
    if(heterogens_ != NULL)
    {
        out << "***************************** HETEROGEN ******************************" << endl;
        heterogens_->Print(out);
    }
    if(heterogens_name_ != NULL)
    {
        out << "*************************** HETEROGEN NAME ***************************" << endl;
        heterogens_name_->Print(out);
    }
    if(heterogen_synonyms_ != NULL)
    {
        out << "************************** HETEROGEN SYNONYM *************************" << endl;
        heterogen_synonyms_->Print(out);
    }
    if(formulas_ != NULL)
    {
        out << "******************************* FORMULA ******************************" << endl;
        formulas_->Print(out);
    }
    if(helixes_ != NULL)
    {
        out << "******************************** HELIX *******************************" << endl;
        helixes_->Print(out);
    }
    if(sheets_ != NULL)
    {
        out << "******************************** SHEET *******************************" << endl;
        sheets_->Print(out);
    }
    if(disulfide_bonds_ != NULL)
    {
        out << "*************************** DISULFIDE BOND ***************************" << endl;
        disulfide_bonds_->Print(out);
    }
    if(links_ != NULL)
    {
        out << "******************************** LINK ********************************" << endl;
        links_->Print(out);
    }
    if(sites_ != NULL)
    {
        out << "******************************** SITE ********************************" << endl;
        sites_->Print(out);
    }
    if(crystallography_ != NULL)
    {
        out << "************************** CRYSTALLOGRAPHIC **************************" << endl;
        crystallography_->Print(out);
    }
    if(origins_ != NULL)
    {
        out << "******************************* ORIGIN *******************************" << endl;
        origins_->Print(out);
    }
    if(scales_ != NULL)
    {
        out << "******************************** SCALE *******************************" << endl;
        scales_->Print(out);
    }
    if(matrices_ != NULL)
    {
        out << "******************************* MATRIX *******************************" << endl;
        matrices_->Print(out);
    }
    if(models_ != NULL)
    {
        out << "******************************* MODEL ********************************" << endl;
        models_->Print(out);
    }
    if(connectivities_ != NULL)
    {
        out << "******************************* CONNECT ******************************" << endl;
        connectivities_->Print(out);
    }
}
