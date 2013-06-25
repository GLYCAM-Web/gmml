#ifndef PREPFILEPROCESSINGEXCEPTION_HPP
#define PREPFILEPROCESSINGEXCEPTION_HPP

#include <exception>
#include <string>

namespace PrepFileSpace
{
    class PrepFileProcessingException : public std::exception
    {
        public:
            /////////////////////////////////// CONSTRUCTOR ////////////////////////////////////
            PrepFileProcessingException(const std::string& message);
            PrepFileProcessingException(int line_number, const std::string& message);

            ////////////////////////////////// FUNCTIONS ///////////////////////////////////////
            virtual const char *what() const throw();
            virtual ~PrepFileProcessingException() throw();

            /////////////////////////////////// ATTRIBUTES /////////////////////////////////////
            int line_number_;               // Line number
            std::string message_;           // message
            mutable std::string what_;      // Explanation
    };
}

#endif // PREPFILEPROCESSINGEXCEPTION_HPP
