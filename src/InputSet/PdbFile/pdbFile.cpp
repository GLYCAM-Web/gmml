#include "includes/InputSet/PdbFile/pdbFile.hpp"
#include "includes/CodeUtils/files.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/InputSet/PdbFile/databaseReferenceRecord.hpp"

using pdb::PdbFile;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbFile::PdbFile(const std::string &pdbFile)
{
    codeutils::ensureFileExists(pdbFile);
}

std::string PdbFile::GetUniprotIDs() const
{
    std::string UniprotIDs = "";
    for (auto &databaseReference : this->GetDatabaseReferences())
    {
        UniprotIDs += databaseReference.GetUniprotID();
    }
    return UniprotIDs;
}
