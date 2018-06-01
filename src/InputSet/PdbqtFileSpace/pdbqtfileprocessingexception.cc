// Author: Alireza Khatamian

#include <sstream>
#include <stdexcept>

#include "../../../includes/common.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtfileprocessingexception.hpp"

using PdbqtFileSpace::PdbqtFileProcessingException;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
PdbqtFileProcessingException::PdbqtFileProcessingException(const std::string &message)
    : line_number_(gmml::dNotSet), message_(message) {}

PdbqtFileProcessingException::PdbqtFileProcessingException(int line_number, const std::string &message)
    : line_number_(line_number), message_(message) {}

//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
/// Exception handler for parameter file exceptions
const char* PdbqtFileProcessingException::what() const throw()
{
    what_ = "PdbqtFile: " + message_;
    if (line_number_ != gmml::dNotSet)
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

PdbqtFileProcessingException::~PdbqtFileProcessingException() throw() {}
