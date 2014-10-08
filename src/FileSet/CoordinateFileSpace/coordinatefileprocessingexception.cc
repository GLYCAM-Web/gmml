#include <sstream>
#include <stdexcept>

#include "../../../includes/common.hpp"
#include "../../../includes/FileSet/CoordinateFileSpace/coordinatefileprocessingexception.hpp"

using namespace gmml;
using namespace CoordinateFileSpace;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
CoordinateFileProcessingException::CoordinateFileProcessingException(const std::string &message)
    : line_number_(dNotSet), message_(message) {}

CoordinateFileProcessingException::CoordinateFileProcessingException(int line_number, const std::string &message)
    : line_number_(line_number), message_(message) {}

//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
/// Exception handler for parameter file exceptions
const char* CoordinateFileProcessingException::what() const throw()
{
    what_ = "CoordinateFile: " + message_;
    if (line_number_ != dNotSet)
    {
        std::stringstream ss;
        if(ss << line_number_)
        {
            what_ += " (line " + ss.str() + ")";
            return what_.c_str();
        }
        else
            throw std::invalid_argument("to_string: invalid conversion");       /// Invalid conversion from int to string
    }
    return what_.c_str();
}

CoordinateFileProcessingException::~CoordinateFileProcessingException() throw() {}


