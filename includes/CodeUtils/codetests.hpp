#ifndef GMML_INCLUDES_CODEUTILS_CODETESTS_HPP
#define GMML_INCLUDES_CODEUTILS_CODETESTS_HPP

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

#endif // GMML_INCLUDES_CODEUTILS_CODETESTS_HPP
