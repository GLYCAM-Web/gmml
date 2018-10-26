#ifndef SUGAR_GLYCODE_META_HPP
#define SUGAR_GLYCODE_META_HPP

/** File glycode.hpp begun on 16 June 2018 by BLFoley */

/**
    The 'glycode' is GMML's internal representation for the molecular
    structures of monosaccharides.

    A detailed description of the meaning of the code can be found at:

        glycam.org/glycode

    The code is ring-centric.  That is, the numbering is relative to 
    the ring.  At the time of this writing, the code does not handle
    linear saccharide structures, but it might be made to do so later.

    GLYCODE, briefly:

      Ring forms supported:   pyranose (P) and furanose (F)

      Regarding numbers in the code:

          * Unsigned integers refer to positions in the ring, starting
            with 1 (one) being the anomeric position.  Typically, the
            integer 1 is not used.  Instead the lower-case 'a' is used
            to indicate the orientation of the anomeric oxygen.
          * Signede integers (+1, +2, etc. and -1, 02, etc.) refer to 
            extensions of the carbon chain beyond the ring.
              *  If positive, the extension is off the ring on the 
                 non-anomeric side of the ring oxygen.
                 Examples:
                  -  For a glucopyranose, +1 designates the C6 position
                  -  For NeuNAc, +1, +2 and +3 designate the C7, C8 and
                     C9 positions. 
              *  If negative, the extension is off the ring on the 
                 anomeric side of the ring oxygen. This will apply to 
                 ketoses and similar (e.g., uloses).
                 Examples:
                  -  For fructopyranose and NeuNAc, -1 designates the 
                     C1 position.  Note that in these, the C2 position
                     is the anomeric carbon and the first ring carbon. 

      Sample Glycode:

          alpha-D-glucopyranose

            Human-readable form:

                 3   +1
                   P
               2,4   a

            Linearized form:

               _2^3_4P_a^+1

               Variants of this are acceptable, per the convenience of 
               the coder, e.g., _2,4^3P_1^+1


    This file:

      Defines a nomenclature for tagging atoms with respect to their 
      position in the glycode.
 */

#include <string>
#include <map>

namespace gmml
{
        namespace MolecularMetadata
        {
                namespace Sugars
                {
		// The ring atoms
                { "r1"  ;  "ring position 1"           }  //< typically C1 (aldose) or C2 (furnaose, ulose)
                { "r1"  ;  "anomeric carbon position"  }  //< so noted because it might not be a carbon
                { "r2"  ;  "ring position 2"           }  //< the count goes away from the ring oxygen
                { "r3"  ;  "ring position 3"           }  //< the count goes away from the ring oxygen
                { "r4"  ;  "ring position 4"           }  //< the count goes away from the ring oxygen
                { "r5"  ;  "ring position 5"           }  //< the count goes away from the ring oxygen
                { "r5"  ;  "ring oxygen F position"    }  //< so noted because it might not be an oxygen
                { "r6"  ;  "ring position 6"           }  //< the count goes away from the ring oxygen
                { "r6"  ;  "ring oxygen P"             }  //< the count goes away from the ring oxygen
                { "r6"  ;  "ring oxygen P position"    }  //< so noted because it might not be an oxygen
		// The immediate side atoms
                { "r1.O"  ;  "anomeric oxygen position"  }  //< so noted because it might not be an oxygen
                { "r1.R"  ;  "anomeric other position"  }  //< typically hydrogen (aldo) or carbon (keto, ulo)
                { "r1.O.R"  ;  "whatever is after the anomeric oxygen position"  }  //< typically hydrogen (aldo) or carbon (keto, ulo)
                {   ;   }
		} // close namespace Sugars
	} // close namespace MolecularMetadata
} // close namespace gmml

#endif // SUGAR_GLYCODE_META_HPP
