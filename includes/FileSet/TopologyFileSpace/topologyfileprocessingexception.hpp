#ifndef TOPOLOGYFILEPROCESSINGEXCEPTION_HPP
#define TOPOLOGYFILEPROCESSINGEXCEPTION_HPP

#include <exception>
#include <string>

namespace TopologyFileSpace
{
    class TopologyFileProcessingException : public std::exception
    {
        public:
            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Constructor of exception handler of topology file class
              * @param message An appropriate message corresponding to an exception
              */
            TopologyFileProcessingException(const std::string& message);
            /*! \fn
              * Constructor of exception handler of topology file class feeded by the line number in which exception has been fired
              * @param line_number The line number in which the exception has been occured
              * @param message An appropriate message corresponding to an exception
              */
            TopologyFileProcessingException(int line_number, const std::string& message);

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
            virtual ~TopologyFileProcessingException() throw();

            //////////////////////////////////////////////////////////
            //                         ATTRIBUTES                   //
            //////////////////////////////////////////////////////////
            int line_number_;               /*!< Line number */
            std::string message_;           /*!< message */
            mutable std::string what_;      /*!< Explanation */
    };
}

#endif // TOPOLOGYFILEPROCESSINGEXCEPTION_HPP
