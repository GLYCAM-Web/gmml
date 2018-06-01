// Author: Alireza Khatamian

#include <sstream>
#include <stdexcept>

#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbfileprocessingexception.hpp"

using PdbFileSpace::PdbFileProcessingException;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
PdbFileProcessingException::PdbFileProcessingException(const std::string &message)
    : line_number_(gmml::dNotSet), message_(message) {}

PdbFileProcessingException::PdbFileProcessingException(int line_number, const std::string &message)
    : line_number_(line_number), message_(message) {}

//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
/// Exception handler for parameter file exceptions
const char* PdbFileProcessingException::what() const throw()
{
    what_ = "PdbFile: " + message_;
    if (line_number_ != gmml::dNotSet)
    {
        std::stringstream ss;
        if(ss << line_number_)
        {
            what_ += " (line " + ss.str() + ")";
            gmml::log(__LINE__, __FILE__, gmml::ERR, what_.c_str());
            return what_.c_str();
        }
        else
        {
            gmml::log(__LINE__, __FILE__, gmml::ERR, __LINE__ + "to_string: invalid conversion");
            throw std::invalid_argument(__LINE__ + "to_string: invalid conversion");       /// Invalid conversion from int to string
        }
    }
    gmml::log(__LINE__, __FILE__, gmml::ERR, what_.c_str());
    return what_.c_str();
}

PdbFileProcessingException::~PdbFileProcessingException() throw() {}
