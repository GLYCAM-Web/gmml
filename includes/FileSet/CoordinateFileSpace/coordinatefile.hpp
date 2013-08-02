#ifndef COORDINATEFILE_HPP
#define COORDINATEFILE_HPP

#include <string>
#include <vector>

#include "../../../includes/Geometry/coordinate.hpp"

namespace CoordinateFileSpace
{
    class CoordinateFile
    {
        public:
            //////////////////////////////////// CONSTRUCTOR ///////////////////////////////////////
            CoordinateFile(const std::string& crd_file);

            ///////////////////////////////////// ACCESSOR /////////////////////////////////////////
            const std::string& GetFilePath() const;
            const std::string& GetTitle() const;
            int GetNumberOfCoodinates();
            const std::vector<Geometry::Coordinate*> GetCoordinates() const;

            ///////////////////////////////////// MUTATOR //////////////////////////////////////////
            void SetPath(std::string& path);
            void SetTitle(std::string& title);
            void SetNumberOfCoordinates(int number_of_coordinates);
            void SetCoordinates(std::vector<Geometry::Coordinate*> coordinates);
            void AddCoordinate(Geometry::Coordinate* coordinate);

            ///////////////////////////////// FUNCTIONS ///////////////////////////////////////////
            void Read(std::ifstream& in_file);

            ///////////////////////////////// DISPLAY FUNCTION ////////////////////////////////////
            void Print(std::ostream& out);

        private:
            std::string path_;                                  // Path of the coordinate file
            std::string title_;                                 // Title; set by the first line of a coordinate file
            int number_of_coordinates_;                         // Number of coordinates containing in the file; set by the second line of a coordinate file
            std::vector<Geometry::Coordinate*> coordinates_;    // List of coordinates in a coordinate file; from the 3rd line of a coordinate file to
                                                                // the end of the file, lines are including coordinates
            // An example of a coordinate file:
            /*
            ROH
                10
               6.0393704   8.3738677  -1.5620029   6.2495186   8.9767698  -2.2789038
               6.1158760  10.3255420  -1.8229947   5.0939355  10.5070215  -1.4901242
               6.4638546  11.2544873  -3.0100838   6.1736819  12.2707215  -2.7392541
               7.9666814  11.2358498  -3.3504444   8.2351252  10.2512600  -3.7361927
               8.8381482  11.5303281  -2.1051706   9.8938827  11.4306846  -2.3616407
            */
    };
}

#endif // COORDINATEFILE_HPP
