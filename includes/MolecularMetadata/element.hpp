#ifndef GMML_ELEMENT_META_HPP
#define GMML_ELEMENT_META_HPP

#include <string>
#include <iostream>

namespace MolecularMetadata
{
    class Element
    {
        public:

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            Element();
            ~Element();
            Element(Element& element);
            Element(Element *element);


            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
               * @{
               */
            /*! \fn
              * An accessor function in order to access to the symbol
              * @return symbol_ attribute of the current object of this class
              */
            std::string GetSymbol();
            /*! \fn
              * An accessor function in order to access to the name
              * @return name_ attribute of the current object of this class
              */
            std::string GetName();
            /*! \fn
              * An accessor function in order to access to the mass
              * @return mass_ attribute of the current object of this class
              */
            double GetMass();
            /*! \fn
              * An accessor function in order to access to the exact mass
              * @return exact_mass_ attribute of the current object of this class
              */
            double GetExactMass();
            /*! \fn
              * An accessor function in order to access to the ioinization energy
              * @return ionization_energy_ attribute of the current object of this class
              */
            double GetIonizationEnergy();
            /*! \fn
              * An accessor function in order to access to the election affinity
              * @return election_affinity_ attribute of the current object of this class
              */
            double GetElectionAffinity();
            /*! \fn
              * An accessor function in order to access to the election negativity
              * @return election_negativity_ attribute of the current object of this class
              */
            double GetElectionNegativity();
            /*! \fn
              * An accessor function in order to access to the covalent radius
              * @return covalent_radius_ attribute of the current object of this class
              */
            double GetCovalentRadius();
            /*! \fn
              * An accessor function in order to access to the van der waals radius
              * @return van_der_waals_radius_ attribute of the current object of this class
              */
            double GetVanDerWaalsRadius();
            /*! \fn
              * An accessor function in order to access to the boiling point
              * @return boiling_point_ attribute of the current object of this class
              */
            double GetBoilingPoint();
            /*! \fn
              * An accessor function in order to access to the melting point
              * @return melting_point_ attribute of the current object of this class
              */
            double GetMeltingPoint();
/** @}*/
            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
               * @{
               */
            /*! \fn
              * A mutator function in order to set the symbol of the current object
              * Set the symbol_ attribute of the current element
              * @param symbol The symbol attribute of the current object
              */
            void SetSymbol(std::string symbol);
            /*! \fn
              * A mutator function in order to set the name of the current object
              * Set the name_ attribute of the current element
              * @param name The name attribute of the current object
              */
            void SetName(std::string name);
            /*! \fn
              * A mutator function in order to set the mass of the current object
              * Set the mass_ attribute of the current element
              * @param mass The mass attribute of the current object
              */
            void SetMass(double mass);
            /*! \fn
              * A mutator function in order to set the exact mass of the current object
              * Set the exact_mass_ attribute of the current element
              * @param exact_mass The exact_mass attribute of the current object
              */
            void SetExactMass(double exact_mass);
            /*! \fn
              * A mutator function in order to set the ionization energy of the current object
              * Set the ionization_energy_ attribute of the current element
              * @param mass The ionization_energy attribute of the current object
              */
            void SetIonizationEnergy(double ionization_energy);
            /*! \fn
              * A mutator function in order to set the election affinity of the current object
              * Set the election_affinity_ attribute of the current element
              * @param election_affinity The election affinity attribute of the current object
              */
            void SetElectionAffinity(double election_affinity);
            /*! \fn
              * A mutator function in order to set the election negativity of the current object
              * Set the election_negativity_ attribute of the current element
              * @param election_negativity The election negativity attribute of the current object
              */
            void SetElectionNegativity(double election_negativity);
            /*! \fn
              * A mutator function in order to set the covalent radius of the current object
              * Set the covalent_radius_ attribute of the current element
              * @param covalent_radius The covalent radius attribute of the current object
              */
            void SetCovalentRadius(double covalent_radius);
            /*! \fn
              * A mutator function in order to set the van der waals radius of the current object
              * Set the van_der_waals_radius_ attribute of the current element
              * @param van_der_waals_radius The van der waals radius attribute of the current object
              */
            void SetVanDerWaalsRadius(double van_der_waals_radius);
            /*! \fn
              * A mutator function in order to set the boiling point of the current object
              * Set the boiling_point_ attribute of the current element
              * @param boiling_point The boiling point attribute of the current object
              */
            void SetBoilingPoint(double boiling_point);
            /*! \fn
              * A mutator function in order to set the melting point of the current object
              * Set the melting_point_ attribute of the current element
              * @param melting_point The melting point attribute of the current object
              */
            void SetMeltingPoint(double melting_point);
/** @}*/
            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the element contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            std::string symbol_;                    /*!< Chemical symbol for an element >*/
            std::string name_;                      /*!< Name of an element >*/
            double mass_;                           /*!< Mass of an element >*/
            double exact_mass_;                     /*!< Exact mass of an element >*/
            double ionization_energy_;              /*!< Ionizing energy of an element >*/
            double election_affinity_;              /*!< Election affinity of an element >*/
            double election_negativity_;            /*!< Election negativity of an element >*/
            double covalent_radius_;                /*!< Covalent radius of an element >*/
            double van_der_waals_radius_;           /*!< Van-der-waals radius of an element >*/
            double boiling_point_;                  /*!< Boiling point of an element >*/
            double melting_point_;                  /*!< Melting point of an element >*/

    };
}

#endif // GMML_ELEMENT_META_HPP
