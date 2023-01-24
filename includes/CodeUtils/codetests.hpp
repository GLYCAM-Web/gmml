#ifndef CODE_UTILS_HPP
#define CODE_UTILS_HPP
// ToDo This belongs in gmml/tests, doesn't need to include common and can be a function instead of a class.
#include "../common.hpp"

namespace CodeUtils
{
    class CodeTests
    {
        public:
            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            CodeTests();

            //////////////////////////////////////////////////////////
            //                           ACCESSOR                   //
            //////////////////////////////////////////////////////////
/** \addtogroup Code_Utils
               * @{
               */
            /*! \fn
              * Return a list of available tests
              */
            StringVector ListCodeTests();
/** @}*/
            //////////////////////////////////////////////////////////
            //                           MUTATOR                    //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                         FUNCTIONS                    //
            //////////////////////////////////////////////////////////

	    void ProduceSegmentationFault();

            //////////////////////////////////////////////////////////
            //                     DISPLAY FUNCTIONS                //
            //////////////////////////////////////////////////////////

    private:
            //////////////////////////////////////////////////////////
            //                         ATTRIBUTES                   //
            //////////////////////////////////////////////////////////

    };
}

#endif // CODE_UTILS_HPP
