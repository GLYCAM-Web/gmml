#include <sstream>
#include <stdexcept>

#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceprocessingexception.hpp"

using CondensedSequenceSpace::CondensedSequenceProcessingException;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
CondensedSequenceProcessingException::CondensedSequenceProcessingException(const std::string &message)
    : line_number_(gmml::dNotSet), message_(message) {}

CondensedSequenceProcessingException::CondensedSequenceProcessingException(int line_number, const std::string &message)
    : line_number_(line_number), message_(message) {}

//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
/// Exception handler for parameter file exceptions
const char* CondensedSequenceProcessingException::what() const throw()
{
    what_ = "CondensedSequence: " + message_;
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
            gmml::log(__LINE__, __FILE__, gmml::ERR, "to_string: invalid conversion");
            throw std::invalid_argument("to_string: invalid conversion");       /// Invalid conversion from int to string
        }
    }
    gmml::log(__LINE__, __FILE__, gmml::ERR, what_.c_str());
    return what_.c_str();
}

CondensedSequenceProcessingException::~CondensedSequenceProcessingException() throw() {}
