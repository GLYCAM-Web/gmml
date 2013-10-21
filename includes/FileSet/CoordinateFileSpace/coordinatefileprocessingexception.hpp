#ifndef COORDINATEFILEPROCESSINGEXCEPTION_HPP
#define COORDINATEFILEPROCESSINGEXCEPTION_HPP

#include <exception>
#include <string>

namespace CoordinateFileSpace
{
    class CoordinateFileProcessingException : public std::exception
    {
        public:
            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Constructor of exception handler of coordinate file class
              * @param message An appropriate message corresponding to an exception
              */
            CoordinateFileProcessingException(const std::string& message);
            /*! \fn
              * Constructor of exception handler of coordinate file class feeded by the line number in which exception has been fired
              * @param line_number The line number in which the exception has been occured
              * @param message An appropriate message corresponding to an exception
              */
            CoordinateFileProcessingException(int line_number, const std::string& message);

            //////////////////////////////////////////////////////////
            //                         FUNCTIONS                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Virtual function to explain what is the reason of the occured exception
              */
            virtual const char *what() const throw();
            /*! \fn
              * Destructor
              */
            virtual ~CoordinateFileProcessingException() throw();

            //////////////////////////////////////////////////////////
            //                         ATTRIBUTES                   //
            //////////////////////////////////////////////////////////
            int line_number_;               /*!< Line number */
            std::string message_;           /*!< message */
            mutable std::string what_;      /*!< Explanation */
    };
}

#endif // COORDINATEFILEPROCESSINGEXCEPTION_HPP
