#include <sstream>
#include <stdexcept>

#include "../../../includes/common.hpp"
#include "../../../includes/FileSet/TopologyFileSpace/topologyfileprocessingexception.hpp"

using namespace gmml;
using namespace TopologyFileSpace;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
TopologyFileProcessingException::TopologyFileProcessingException(const std::string &message)
    : line_number_(dNotSet), message_(message) {}

TopologyFileProcessingException::TopologyFileProcessingException(int line_number, const std::string &message)
    : line_number_(line_number), message_(message) {}

//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
/// Exception handler for parameter file exceptions
const char* TopologyFileProcessingException::what() const throw()
{
    what_ = "TopologyFile: " + message_;
    if (line_number_ != dNotSet)
    {
        std::stringstream ss;
        if(ss << line_number_)
        {
            what_ += " (line " + ss.str() + ")";
            return what_.c_str();
        }
        else
            throw std::invalid_argument(__LINE__ + "to_string: invalid conversion");       /// Invalid conversion from int to string
    }
    return what_.c_str();
}

TopologyFileProcessingException::~TopologyFileProcessingException() throw() {}
