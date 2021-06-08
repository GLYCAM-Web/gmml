#ifndef OFFFILERESIDUE_HPP
#define OFFFILERESIDUE_HPP

#include <string>
#include <vector>
#include <map>
#include <iostream>

namespace OffFileSpace
{
    class OffFileAtom;
    class OffFileResidue
    {
        public:
            //////////////////////////////////////////////////////////
            //                     TYPE DEFINITION                  //
            //////////////////////////////////////////////////////////
            /*! \def
              * A mapping between the atom order and its atom object
              */
          //  typedef std::map<int, OffFileAtom*> AtomMap;
            typedef std::vector<OffFileSpace::OffFileAtom*> OffFileAtomVector;

            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            OffFileResidue();

            //////////////////////////////////////////////////////////
            //                           ACCESSOR                   //
            //////////////////////////////////////////////////////////
            /*! \fn
              * An accessor function in order to access to residue name of the current object
              * @return name_ attribute of the current residue
              */
            std::string GetName();
            /*! \fn
              * An accessor function in order to access to the atoms belonging to the current object
              * @return atoms_ attribute of the current residue
              */
            OffFileAtomVector GetAtoms();
            /*! \fn
              * An accessor function in order to access to a atom belonging to the current object by its atom index
              * @param index Index of the target atom in the current residue
              * @return A pointer to the atom belonging to the current residue having the given atom index
              */
            OffFileAtom* GetAtomByIndex(int index);
            /*! \fn
              * An accessor function in order to access to a atom belonging to the current object by its atom order
              * @param order Order of the target atom in the current residue
              * @return A pointer to the atom belonging to the current residue having the given atom order
              */
            OffFileAtom* GetAtomByOrder(int order);
            /*! \fn
              * An accessor function in order to access to the box angle of the current object
              * @return box_angle_ attribute of the current residue
              */
            double GetBoxAngle();
            /*! \fn
              * An accessor function in order to access to the box length of the current object
              * @return box_length_ attribute of the current residue
              */
            double GetBoxLength();
            /*! \fn
              * An accessor function in order to access to the box width of the current object
              * @return box_width_ attribute of the current residue
              */
            double GetBoxWidth();
            /*! \fn
              * An accessor function in order to access to the box height of the current object
              * @return box_height_ attribute of the current residue
              */
            double GetBoxHeight();
            /*! \fn
              * An accessor function in order to access to the index of the head atom of the current object
              * @return head_atom_index_ attribute of the current residue
              */
            int GetHeadAtomIndex();
            /*! \fn
              * An accessor function in order to access to the tail atom of the current object
              * @return tail_atom_index attribute of the current residue
              */
            int GetTailAtomIndex();
            /*! \fn
              * A function in order to access to Off file atom by a atom name
              * @param atom_name The name of the atom
              * @return Off_file_atom
              */
            OffFileAtom* GetOffAtomByAtomName(std::string atom_name);
            /*! \fn
              * An accessor function in order to access to listing index of the current object
              * @return listing_index attribute of the current residue
              */
            int GetListingIndex();

            //////////////////////////////////////////////////////////
            //                           MUTATOR                    //
            //////////////////////////////////////////////////////////

            /*! \fn
              * A mutator function in order to set the residue name of the current object
              * Set the name_ attribute of the current atom
              * @param name Residue name
              */
            void SetName(std::string name);
            /*! \fn
              * A mutator function in order to set the atoms belonging to the current object
              * Set the atoms_ attribute of the current atom
              * @param atoms A set of atoms belonging to the current residue
              */
            void SetAtoms(OffFileAtomVector atoms);
            /*! \fn
              * A mutator function in order to add an atom belonging to the current object
              * Add a new entry to the atoms_ attribute of the current atom
              * @param atom A new atom belonging to the current residue
              */
            void AddAtom(OffFileAtom* atom);
            /*! \fn
              * A mutator function in order to set the box angle of the current object
              * Set the box_angle_ of the current atom
              * @param box_angle A double value of the box angle of the current residue
              */
            void SetBoxAngle(double box_angle);
            /*! \fn
              * A mutator function in order to set the box length of the current object
              * Set the box_length_ of the current atom
              * @param box_length A double value of the box length of the current residue
              */
            void SetBoxLength(double box_length);
            /*! \fn
              * A mutator function in order to set the box width of the current object
              * Set the box_width_ of the current atom
              * @param box_width A double value of the box width of the current residue
              */
            void SetBoxWidth(double box_width);
            /*! \fn
              * A mutator function in order to set the box height of the current object
              * Set the box_height_ of the current atom
              * @param box_height A double value of the box_height of the current residue
              */
            void SetBoxHeight(double box_height);
            /*! \fn
              * A mutator function in order to set the index of the head atom of the current object
              * Set the head_atom_index_ of the current atom
              * @param head_atom_index An integer that indicates the head atom index of the current residue
              */
            void SetHeadAtomIndex(int head_atom_index);
            /*! \fn
              * A mutator function in order to set the index of the tail atom of the current object
              * Set the tail_atom_index_ of the current atom
              * @param tail_atom_index An integer that indicates the tail atom index of the current residue
              */
            void SetTailAtomIndex(int tail_atom_index);
            /*! \fn
              * A mutator function in order to set the listing index of the current object
              * Set the listing_index_ of the current atom
              * @param listing_index An integer that indicates the tail atom index of the current residue
              */
            void SetListingIndex(int listing_index);

            //////////////////////////////////////////////////////////
            //                     DISPLAY FUNCTIONS                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the current object information in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                         ATTRIBUTES                   //
            //////////////////////////////////////////////////////////
            std::string name_;          /*!< Residue name.Set by each line of the index section of a Off file */
            OffFileAtomVector atoms_;   /*!< List of atoms belonging to the residue mapped to their order numbers; set by various sections (atom, connectivity, positions) of a Off file*/
            double box_angle_;          /*!< Box angle; set by the 2nd line of the boundbox section of a Off file*/
            double box_length_;         /*!< Box length; set by the 3rd line of the boundbox section of a Off file*/
            double box_width_;          /*!< Box width; set by the 4th line of the boundbox section of a Off file*/
            double box_height_;         /*!< Box height; set by the 5th line of the boundbox section of a Off file*/
            int head_atom_index_;       /*!< Head atom index; set by the 1st line of the connect section of a Off file*/
            int tail_atom_index_;       /*!< Tail atom index; set by the 2nd line of the connect section of a Off file*/
            int listing_index_;
    };
}

#endif // OFFFILERESIDUE_HPP
