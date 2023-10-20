#ifndef INCLUDES_CENTRALDATASTRUCTURE_ENSEMBLE_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_ENSEMBLE_HPP

#include "includes/CentralDataStructure/assembly.hpp"
#include "includes/MolecularModeling/TemplateGraph/GraphStructure/include/Node.hpp"
#include <string>
#include <vector>
#include <memory> // unique_ptr

namespace cds
{
    class Ensemble : public glygraph::Node<Ensemble>
    {
      public:
        //////////////////////////////////////////////////////////
        //                    CONSTRUCTOR                       //
        //////////////////////////////////////////////////////////
        Ensemble()
        {}

        Ensemble(Ensemble&& other) noexcept; // Move Ctor
        Ensemble(const Ensemble& other);     // Copy Ctor
        Ensemble& operator=(Ensemble other); // Move and Copy assignment operator
        virtual ~Ensemble() = default;

        friend void swap(Ensemble& lhs,
                         Ensemble& rhs) // ToDo figure out how to put this in cpp file once everything is working. Yo
                                        // just define it without the friend keyword you bozo.
        {
            using std::swap;
            swap(lhs.assemblies_, rhs.assemblies_);
        }

        //////////////////////////////////////////////////////////
        //                    ACCESSOR                          //
        //////////////////////////////////////////////////////////
        std::vector<Atom*> getAtoms();
        std::vector<Residue*> getResidues();
        std::vector<Molecule*> getMolecules() const;
        std::vector<Assembly*> getAssemblies() const;
        //////////////////////////////////////////////////////////
        //                    MUTATOR                           //
        //////////////////////////////////////////////////////////
        void addAssembly(std::unique_ptr<Assembly> myAssembly);
        void addAssembly(Assembly& myAssembly);
        //////////////////////////////////////////////////////////
        //                    FUNCTIONS                         //
        //////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////
        //                    DISPLAY                           //
        //////////////////////////////////////////////////////////
      private:
        //////////////////////////////////////////////////////////
        //                    ATTRIBUTES                        //
        //////////////////////////////////////////////////////////
        std::vector<std::unique_ptr<Assembly>> assemblies_;
    };
} // namespace cds
#endif // INCLUDES_CENTRALDATASTRUCTURE_ENSEMBLE_HPP
