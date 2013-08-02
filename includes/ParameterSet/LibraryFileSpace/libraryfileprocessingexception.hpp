#ifndef LIBRARYFILEPROCESSINGEXCEPTION_HPP
#define LIBRARYFILEPROCESSINGEXCEPTION_HPP

#include <exception>
#include <string>

namespace LibraryFileSpace
{
    class LibraryFileProcessingException : public std::exception
    {
        public:
            /////////////////////////////////// CONSTRUCTOR ////////////////////////////////////
            LibraryFileProcessingException(const std::string& message);
            LibraryFileProcessingException(int line_number, const std::string& message);

            ////////////////////////////////// FUNCTIONS ///////////////////////////////////////
            virtual const char *what() const throw();
            virtual ~LibraryFileProcessingException() throw();

            /////////////////////////////////// ATTRIBUTES /////////////////////////////////////
            int line_number_;               // Line number
            std::string message_;           // message
            mutable std::string what_;      // Explanation
    };
}

#endif // LIBRARYFILEPROCESSINGEXCEPTION_HPP
