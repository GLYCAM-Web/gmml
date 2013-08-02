#ifndef COORDINATEFILEPROCESSINGEXCEPTION_HPP
#define COORDINATEFILEPROCESSINGEXCEPTION_HPP

#include <exception>
#include <string>

namespace CoordinateFileSpace
{
    class CoordinateFileProcessingException : public std::exception
    {
        public:
            /////////////////////////////////// CONSTRUCTOR ////////////////////////////////////
            CoordinateFileProcessingException(const std::string& message);
            CoordinateFileProcessingException(int line_number, const std::string& message);

            ////////////////////////////////// FUNCTIONS ///////////////////////////////////////
            virtual const char *what() const throw();
            virtual ~CoordinateFileProcessingException() throw();

            /////////////////////////////////// ATTRIBUTES /////////////////////////////////////
            int line_number_;               // Line number
            std::string message_;           // message
            mutable std::string what_;      // Explanation
    };
}

#endif // COORDINATEFILEPROCESSINGEXCEPTION_HPP
