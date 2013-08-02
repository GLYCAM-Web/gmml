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
            std::string path_;
            std::string title_;
            int number_of_coordinates_;
            std::vector<Geometry::Coordinate*> coordinates_;

    };
}

#endif // COORDINATEFILE_HPP
