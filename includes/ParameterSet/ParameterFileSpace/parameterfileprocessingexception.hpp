#ifndef PARAMETERFILEPROCESSINGEXCEPTION_HPP
#define PARAMETERFILEPROCESSINGEXCEPTION_HPP

#include <exception>
#include <string>

namespace ParameterFileSpace
{
    class ParameterFileProcessingException : public std::exception
    {
        public:
            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Constructor of exception handler of parameter file class
              * @param message An appropriate message corresponding to an exception
              */
            ParameterFileProcessingException(const std::string& message);
            /*! \fn
              * Constructor of exception handler of parameter file class feeded by the line number in which exception has been fired
              * @param line_number The line number in which the exception has been occured
              * @param message An appropriate message corresponding to an exception
              */
            ParameterFileProcessingException(int line_number, const std::string& message);

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
            virtual ~ParameterFileProcessingException() throw();

            //////////////////////////////////////////////////////////
            //                         ATTRIBUTES                   //
            //////////////////////////////////////////////////////////
            int line_number_;               /*!< Line number*/
            std::string message_;           /*!< message*/
            mutable std::string what_;      /*!< Explanation*/
    };
}

#endif // PARAMETERFILEPROCESSINGEXCEPTION_HPP
