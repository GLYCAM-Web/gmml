#include "includes/CentralDataStructure/Readers/Pdb/pdbFunctions.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"

int pdb::checkShiftFromSerialNumberOverrun(std::string line)
{
    int shift = 0;
    if (isdigit(line[12]) && line[22] != ' ')
    {
        shift = codeUtils::GetSizeOfIntInString(line.substr(12));
        std::stringstream ss;
        ss << "Shift detected as position 12 is a digit: " << line[12] << " and this isn't blank: " << line[22] << "\n"
           << line << "\n";
        gmml::log(__LINE__, __FILE__, gmml::ERR, ss.str());
    }
    return shift;
}
