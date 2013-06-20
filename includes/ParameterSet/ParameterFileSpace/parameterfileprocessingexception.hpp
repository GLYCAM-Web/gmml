#ifndef PARAMETERFILEPROCESSINGEXCEPTION_HPP
#define PARAMETERFILEPROCESSINGEXCEPTION_HPP

#include <exception>
#include <string>

namespace ParameterFileSpace
{
    class ParameterFileProcessingException : public std::exception
    {
        public:
            /////////////////////////////////// CONSTRUCTOR ////////////////////////////////////
            ParameterFileProcessingException(const std::string& message);
            ParameterFileProcessingException(int line_number, const std::string& message);

            ////////////////////////////////// FUNCTIONS ///////////////////////////////////////
            virtual const char *what() const throw();
            virtual ~ParameterFileProcessingException() throw();

            /////////////////////////////////// ATTRIBUTES /////////////////////////////////////
            int line_number_;               // Line number
            std::string message_;           // message
            mutable std::string what_;      // Explanation
    };
}

#endif // PARAMETERFILEPROCESSINGEXCEPTION_HPP
