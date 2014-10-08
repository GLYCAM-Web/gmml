// Author: Alireza Khatamian

#include <sstream>
#include <stdexcept>

#include "../../../includes/common.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbfileprocessingexception.hpp"

using namespace gmml;
using namespace PdbFileSpace;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
PdbFileProcessingException::PdbFileProcessingException(const std::string &message)
    : line_number_(dNotSet), message_(message) {}

PdbFileProcessingException::PdbFileProcessingException(int line_number, const std::string &message)
    : line_number_(line_number), message_(message) {}

//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
/// Exception handler for parameter file exceptions
const char* PdbFileProcessingException::what() const throw()
{
    what_ = "PdbFile: " + message_;
    if (line_number_ != dNotSet)
    {
        std::stringstream ss;
        if(ss << line_number_)
        {
            what_ += " (line " + ss.str() + ")";
            return what_.c_str();
        }
        else
        {
            throw std::invalid_argument(__LINE__ + "to_string: invalid conversion");       /// Invalid conversion from int to string
        }
    }
    return what_.c_str();
}

PdbFileProcessingException::~PdbFileProcessingException() throw() {}


