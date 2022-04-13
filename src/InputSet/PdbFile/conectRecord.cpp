#include <iomanip> // setw
#include "includes/InputSet/PdbFile/conectRecord.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/CodeUtils/logging.hpp"

using pdb::ConectRecord;

//////////////////////////////////////////////////////////
//                    CONSTRUCTOR                       //
//////////////////////////////////////////////////////////
ConectRecord::ConectRecord(std::string &line, pdb::PdbModel& pdbModel)
{
    std::vector<std::string> possibleSerialNumberStrings = {line.substr(6,5), line.substr(11,5), line.substr(16,5), line.substr(21,5), line.substr(26,5) };
    std::vector<AtomRecord*> atomRecords;
    for (auto &serialNumberString : possibleSerialNumberStrings)
    {
        int serialNumber = 0;
        try
        {
            serialNumber = std::stoi(codeUtils::RemoveWhiteSpace(serialNumberString));
        }
        catch (...) {} // this is fine, they might not all be present.
        if (serialNumber != 0)
        {
            AtomRecord* foundAtom = pdbModel.FindAtom(serialNumber);
            if (foundAtom != nullptr)
            {
                atomRecords.push_back(foundAtom);
            }
            else
            {
                gmml::log(__LINE__, __FILE__, gmml::WAR, "Could not find atom with serial number: " + serialNumberString);
            }
        }
    }
}
ConectRecord::ConectRecord(std::vector<AtomRecord*> atoms)
{
    atomRecordPtrs_ = atoms;
}
//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
void ConectRecord::Print(std::ostream& out) const
{
    out << "CONECT ";
    out << std::setw(5);
    for (auto &atomRecordPtr : atomRecordPtrs_)
    {
        out << std::right << atomRecordPtr->GetSerialNumber();
    }
    out << "\n";
}

void ConectRecord::Write(std::ostream& stream) const
{
    stream << "CONECT";
    stream << std::setw(5);
    for (auto &atomRecordPtr : atomRecordPtrs_)
    {
        stream << " " << std::right << atomRecordPtr->GetSerialNumber();
    }
    stream << "\n";
}

