// Created by: Alireza Khatamian
// Modified by: Alireza Khatamian, Delaram Rahbarinia

#ifndef PDBHELIX_HPP
#define PDBHELIX_HPP

#include <string>
#include <vector>

namespace PdbFileSpace
{
    enum PdbHelixClass
    {
        RIGHT_HANDED_ALPHA = 1,
        RIGHT_HANDED_OMEGA = 2,
        RIGHT_HANDED_PI = 3,
        RIGHT_HANDED_GAMMA = 4,
        RIGHT_HANDED_310 = 5,
        LEFT_HANDED_ALPHA = 6,
        LEFT_HANDED_OMEGA_ = 7,
        LEFT_HANDED_GAMMA_ = 8,
        RIBBON_27 = 9,
        POLYPROLINE = 10
    };

    class PdbHelixResidue;
    class PdbHelix
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            typedef std::vector<PdbHelixResidue*> HelixResidueVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            PdbHelix();
            PdbHelix(const std::string& helix_id, int helix_serial_number, HelixResidueVector helix_residues,
                     PdbHelixClass helix_class, const std::string& comment, double helix_length);
            PdbHelix(std::stringstream& specification_block);

            //////////////////////////////////////////////////////////
            //                         ACCESSOR                     //
            //////////////////////////////////////////////////////////
            std::string GetHelixId();
            int GetHelixSerialNumber();
            HelixResidueVector GetHelixResidues();
            PdbHelixClass GetHelixClass();
            std::string GetComment();
            double GetHelixLength();

            //////////////////////////////////////////////////////////
            //                          MUTATOR                     //
            //////////////////////////////////////////////////////////
            void SetHelixId(const std::string helix_id);
            void SetHelixSerialNumber(int helix_serial_number);
            void SetHelixResidues(const HelixResidueVector helix_residues);
            void AddHelixResidue(PdbHelixResidue* helix_residue);
            void SetHelixClass(PdbHelixClass helix_class);
            void SetComment(const std::string& comment);
            void SetHelixLength(double helix_length);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                      DISPLAY FUNCTION                //
            //////////////////////////////////////////////////////////

        private:
            //////////////////////////////////////////////////////////
            //                        ATTRIBUTES                    //
            //////////////////////////////////////////////////////////
            std::string helix_id_;
            int helix_serial_number_;
            HelixResidueVector helix_residues_;
            PdbHelixClass helix_class_;
            std::string comment_;
            double helix_length_;
    };
}

#endif // PDBHELIX_HPP
