#include "gmml.hpp"
#include "../../includes/MolecularModeling/ring_shape_detection.hpp"

int main(void)
{
    MolecularModeling::Assembly assembly( "tests/inputs/ring_shape_detection_input.pdb", gmml::InputFileType::PDB);
    assembly.BuildStructureByDistance(4, 1.91); // 4 threads, 1.91 cutoff to allow C-S in Cys and Met to be bonded.
    GeometryTopology::CoordinateVector ring_coords = assembly.GetAllCoordinates();
    std::string shape = glylib::CalculateRingShapeBFMP(ring_coords, 10);
    std::cout << "Shape is " << shape << "\n";
    return 0;
}
