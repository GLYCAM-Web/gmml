#ifndef PREPFILEPROCESSINGEXCEPTION_HPP
#define PREPFILEPROCESSINGEXCEPTION_HPP

#include <exception>
#include <string>

namespace PrepFileSpace
{
    class PrepFileProcessingException : public std::exception
    {
        public:
            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Constructor of exception handler of prep file class
              * @param message An appropriate message corresponding to an exception
              */
            PrepFileProcessingException(const std::string& message);
            /*! \fn
              * Constructor of exception handler of prep file class feeded by the line number in which exception has been fired
              * @param line_number The line number in which the exception has been occured
              * @param message An appropriate message corresponding to an exception
              */
            PrepFileProcessingException(int line_number, const std::string& message);

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
            virtual ~PrepFileProcessingException() throw();

            //////////////////////////////////////////////////////////
            //                         ATTRIBUTES                   //
            //////////////////////////////////////////////////////////
            int line_number_;               /*!< Line number */
            std::string message_;           /*!< message */
            mutable std::string what_;      /*!< Explanation */
    };
}

#endif // PREPFILEPROCESSINGEXCEPTION_HPP
