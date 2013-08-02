#include <sstream>
#include <stdexcept>

#include "../../../includes/common.hpp"
#include "../../../includes/ParameterSet/LibraryFileSpace/libraryfileprocessingexception.hpp"

using namespace gmml;
using namespace LibraryFileSpace;

///////////////////////////////////////// CONSTRUCTOR //////////////////////////////////////////
LibraryFileProcessingException::LibraryFileProcessingException(const std::string &message)
    : line_number_(kNotSet), message_(message) {}

LibraryFileProcessingException::LibraryFileProcessingException(int line_number, const std::string &message)
    : line_number_(line_number), message_(message) {}

/////////////////////////////////////// FUNCTIONS /////////////////////////////////////////////

/// Exception handler for parameter file exceptions
const char* LibraryFileProcessingException::what() const throw()
{
    what_ = "ParameterFile: " + message_;
    if (line_number_ != kNotSet)
    {
        std::stringstream ss;
        if(ss << line_number_)
        {
            what_ += " (line " + ss.str() + ")";
            return what_.c_str();
        }
        throw std::invalid_argument("to_string: invalid conversion");       // Invalid conversion from int to string
    }
    return "";
}

LibraryFileProcessingException::~LibraryFileProcessingException() throw() {}

