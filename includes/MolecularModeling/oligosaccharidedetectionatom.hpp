#ifndef	OLIGOSACCHARIDEDETECTIONATOM_HPP
#define OLIGOSACCHARIDEDETECTIONATOM_HPP

#include <iostream>
#include <string>
#include <map>

namespace MolecularModeling
{
    class OligoSaccharideDetectionAtom
    {
	typedef std::multimap <std::string, std::string> OligosaccharidePropertyTags;
        public:
            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            OligoSaccharideDetectionAtom();

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup OligoSaccharideDetectionAtom
               * @{
               */
            /*! \fn
              * An accessor function in order to determind if this atom belongs to a cycle(ring). Default is false.
              * @return IsCycle_ attribute of the current object of this class
              */
	    bool GetIsCycle();
            /*! \fn
              * An accessor function in order to determind if this atom belongs to a sidechain. Default is false.
              * @return IsSideChain_ attribute of the current object of this class
              */
            bool GetIsSideChain();
            /*! \fn
              * An accessor function in order to determind if this atom is an anomeric carbon. Default is false.
              * @return IsAnomericCarbon_ attribute of the current object of this class
              */
            bool GetIsAnomericCarbon();

	    std::string GetNaming();
	
	    OligosaccharidePropertyTags GetAllOligosaccharidePropertyTags();
/** @}*/

             /////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup OligoSaccharideDetectionAtom
               * @{
               */
	    /*! \fn
              * An mutator function in order to change the value of IsCycle boolean. Default is false.
	      * Set the IsCycle_attribute of the current object of this class
	      * @param is_cycle A boolean indicating if this atom is part of a cycle(ring).
              */
            void SetIsCycle(bool is_cycle);
	    /*! \fn
              * An mutator function in order to change the value of IsSideChain boolean. Default is false.
	      * Set the IsSideChain_attribute of the current object of this class
	      * @param is_side_chain A boolean indicating if this atom is part of a sidechain.
              */
            void SetIsSideChain(bool is_side_chain);
            /*! \fn
              * An accessor function in order to set if this atom is an anomeric carbon. Default is false.
	      * Set the IsAnomericCarbon_ attribute of the current object of this class
	      * @param is_anomeric_carbon A boolean indicating if this atom is anomeric.
              */
            void SetIsAnomericCarbon(bool is_anomeric_carbon);

	    void SetNaming(std::string naming);

	    void AddOligosaccharidePropertyTag (std::pair<std::string, std::string> role_tag_pair);

	    void RemoveOligosaccharidePropertyTag (std::pair<std::string, std::string> role_tag_pair);

/** @}*/

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the oligosaccharidedetectionatom contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

	private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
	    bool IsCycle_;						// A boolean indicating if this atom is part of a cycle(ring).
	    bool IsSideChain_;						// A boolean indicating if this atom is part of a sidechain.
	    bool IsAnomericCarbon_;
	    std::string naming_;
	    bool is_aglycon_;
	    OligosaccharidePropertyTags oligosaccharide_properties_;

    };//class

}//namespace
#endif
