#include "includes/CentralDataStructure/Readers/Pdb/pdbFunctions.hpp"
#include "includes/CodeUtils/logging.hpp"
#include "includes/CodeUtils/strings.hpp"

int pdb::checkShiftFromSerialNumberOverrun(std::string line)
{
    int shift = 0;
    if (isdigit(line[11]) && line[20] != ' ')
    {
        shift =
            codeUtils::GetSizeOfIntInString(line.substr(12)); // The shift starts when this gets overrun, not 11. Right?
        std::stringstream ss;
        ss << "Shift of size " << shift << " detected as position 12 is a digit: " << line[11] << " and position 21 >>>"
           << line[20] << "<<< isn't blank: " << line[20] << ":\n           v        v\n"
           << line << "\n";
        gmml::log(__LINE__, __FILE__, gmml::WAR, ss.str());
    }
    return shift;
}
