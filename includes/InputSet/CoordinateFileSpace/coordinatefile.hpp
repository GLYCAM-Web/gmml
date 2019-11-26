#ifndef COORDINATEFILE_HPP
#define COORDINATEFILE_HPP

#include <string>
#include <vector>
#include <iostream>

#include "../../GeometryTopology/coordinate.hpp"

/*! \namespace CoordinateFileSpace */
namespace CoordinateFileSpace
{
    /*! \class CoordinateFile */
    class CoordinateFile
    {
        public:            
            //////////////////////////////////////////////////////////
            //                       Constructor                    //
            //////////////////////////////////////////////////////////
            /*!
              * Default constructor
              */
            CoordinateFile();
            /*!
              * Constructor
              * @param crd_file An existing coordinate file path to be read
              */
            CoordinateFile(const std::string& crd_file);
            /*!
              * load coordinate file
              */
            CoordinateFile* LoadCoordinateFile();
            /*!
              * @param crd_file An existing coordinate file path to be loaded
              */
            CoordinateFile* LoadCoordinateFile(const std::string& crd_file);

            //////////////////////////////////////////////////////////
            //                           ACCESSOR                   //
            //////////////////////////////////////////////////////////
/** \addtogroup Molecular_Data_Structure
              * @{
              */
            /*! \fn
              * An accessor function in order to access to coordinate file path of the current object.
              * @return path_ attribute of the current object of this class
              */
            const std::string& GetFilePath() const;
            /*! \fn
              * An accessor function in order to access to title attribute of the current object
              * The attribute is set by the contents of the given file
              * @return title_ attribute of the current object of this class
              */
            const std::string& GetTitle() const;
            /*! \fn
              * An accessor function in order to access to number of coordinate attribute of the current object
              * The attribute is set by the contents of the given file
              * @return number_of_coordinates_ attribute of the current object of this class
              */
            int GetNumberOfCoodinates();
            /*! \fn
              * An accessor function in order to access to coordinates of the current object
              * The attribute is set by the contents of the given file
              * @return coordinates_ attribute of the current object of this class
              */
            const std::vector<GeometryTopology::Coordinate*> GetCoordinates() const;
/** @}*/
            //////////////////////////////////////////////////////////
            //                           MUTATOR                    //
            //////////////////////////////////////////////////////////
/** \addtogroup Manipulators
              * @{
              */
            /*! \fn
              * A mutator function in order to set the path attribute of the current object
              * Set the path_ attribute of the current object
              * @param path A string defines an actual path of a coordinate file
              */
            void SetPath(std::string& path);
            /*! \fn
              * A mutator function in order to set the title attribute of the current object
              * Set the title_ attribute of the current object
              * @param title A string defines the title of the coordinate file which is extracted from the contents of the file
              */
            void SetTitle(std::string& title);
            /*! \fn
              * A mutator function in order to set the number of coordinates attribute of the current object
              * Set the number_of_coordinates_ attribute of the current object
              * @param number_of_coordinates An integer defines the number of coordinates includeing
              *             in the file which is extracted from the contents of the file
              */
            void SetNumberOfCoordinates(int number_of_coordinates);
            /*! \fn
              * A mutator function in order to set the coordinates attribute of the current object
              * Set the coordinates_ attribute of the current object
              * @param coordinates A vector of coordinates including in the file which which each is extracted
              *             from the contents of the file
              */
            void SetCoordinates(std::vector<GeometryTopology::Coordinate*> coordinates);
            /*! \fn
              * A mutator function in order to add a single coordinate to the coordinates attribute of the current object
              * Add a new entry to the coordinates_ attribute of the current object
              * @param title A string defines the title of the coordinate file which is extracted from the contents of the file
              */
            void AddCoordinate(GeometryTopology::Coordinate* coordinate);
/** @}*/
            //////////////////////////////////////////////////////////
            //                         FUNCTIONS                    //
            //////////////////////////////////////////////////////////
/** \addtogroup Input_File_Reader
              * @{
              */
            /*! \fn
              * A function to parse the contents of a given stream of a file
              * Parse the given stream and set the attributes of the current object accordingly
              * @param in_file A stream contains whole contents of a coordinate file
              */
            void Read(std::ifstream& in_file);
/** @}*/
/** \addtogroup Output_File_Builder
              * @{
              */
            /*! \fn
              * A function to write back the coordinate file
              * @param coordinate_file The output file path
              */
            void Write(const std::string& coordinate_file);
            /*! \fn
              * A function to wrtie into a stream that has been created from a file path
              * @param out_stream The output stream
              */
            void BuildCoordinateFile(std::ofstream& out_stream);
/** @}*/
            //////////////////////////////////////////////////////////
            //                     DISPLAY FUNCTIONS                //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the file contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cerr);

        private:
            //////////////////////////////////////////////////////////
            //                         ATTRIBUTES                   //
            //////////////////////////////////////////////////////////
            std::string path_;                                  /*!< Path of the coordinate file */
            std::string title_;                                 /*!< Title; set by the first line of a coordinate file */
            int number_of_coordinates_;                         /*!< Number of coordinates containing in the file; set by the second line of a coordinate file */
            std::vector<GeometryTopology::Coordinate*> coordinates_;    /*!< List of coordinates in a coordinate file; from the 3rd line of a coordinate file to
                                                                        the end of the file, lines are including coordinates */
            /*! \file
              * An example of a coordinate file:
              *     ROH
              *      10
              *     6.0393704   8.3738677  -1.5620029   6.2495186   8.9767698  -2.2789038
              *     6.1158760  10.3255420  -1.8229947   5.0939355  10.5070215  -1.4901242
              *     6.4638546  11.2544873  -3.0100838   6.1736819  12.2707215  -2.7392541
              *     7.9666814  11.2358498  -3.3504444   8.2351252  10.2512600  -3.7361927
              *     8.8381482  11.5303281  -2.1051706   9.8938827  11.4306846  -2.3616407
              */
    };
}

#endif // COORDINATEFILE_HPP
