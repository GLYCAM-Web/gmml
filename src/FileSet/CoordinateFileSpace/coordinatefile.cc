#include <fstream>
#include <sstream>
#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>

#include "../../../includes/utils.hpp"
#include "../../../includes/FileSet/CoordinateFileSpace/coordinatefile.hpp"
#include "../../../includes/FileSet/CoordinateFileSpace/coordinatefileprocessingexception.hpp"

using namespace Geometry;
using namespace CoordinateFileSpace;
using namespace std;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       Constructor                    //
//////////////////////////////////////////////////////////
CoordinateFile::CoordinateFile()
{
    path_ = "GMML-Generated";
}

CoordinateFile::CoordinateFile(const string &crd_file)
{
    path_ = crd_file;
    std::ifstream in_file;        
    if(std::ifstream(crd_file.c_str()))
        in_file.open(crd_file.c_str());
    else
    {
        throw CoordinateFileProcessingException(__LINE__, "Coordinate file not found");
    }
    Read(in_file);
    in_file.close();            /// Close the parameter files
}
CoordinateFile* CoordinateFile::LoadCoordinateFile()
{
    CoordinateFile* crd = new CoordinateFile();
    return crd;
}

CoordinateFile* CoordinateFile::LoadCoordinateFile(const string &crd_file)
{
    CoordinateFile* crd = new CoordinateFile(crd_file);
    return crd;
}
//////////////////////////////////////////////////////////
//                           ACCESSOR                   //
//////////////////////////////////////////////////////////
/// Return path of the coordinate file
const std::string& CoordinateFile::GetFilePath() const
{
    return path_;
}

/// Return title of the coordinate file
const std::string& CoordinateFile::GetTitle() const
{
    return title_;
}

/// Return the number of coordinates contained in the coordinate file
int CoordinateFile::GetNumberOfCoodinates()
{
    return number_of_coordinates_;
}

/// Return list of coordinates contained in the coordinate file
const std::vector<Geometry::Coordinate*> CoordinateFile::GetCoordinates() const
{
    return coordinates_;
}

//////////////////////////////////////////////////////////
//                           MUTATOR                    //
//////////////////////////////////////////////////////////
/// Set the path of the coordinate file
void CoordinateFile::SetPath(std::string& path)
{
    path_ = path;
}

/// Set the title of the coordinate file
void CoordinateFile::SetTitle(std::string& title)
{
    title_ = title;
}

/// Set the number of coordinates contained in the coordinate file
void CoordinateFile::SetNumberOfCoordinates(int number_of_coordinates)
{
    number_of_coordinates_ = number_of_coordinates;
}

/// Clear and set the list of coordinates in the coordinate file
void CoordinateFile::SetCoordinates(std::vector<Coordinate*> coordinates)
{
    coordinates_.clear();
    for(vector<Coordinate*>::iterator it = coordinates.begin(); it != coordinates.end(); it++)
    {
        coordinates_.push_back(*it);
    }
}

/// Add a new coordinate to the list of the coordinates
void CoordinateFile::AddCoordinate(Geometry::Coordinate* coordinate)
{
    coordinates_.push_back(coordinate);
}

//////////////////////////////////////////////////////////
//                         FUNCTIONS                    //
//////////////////////////////////////////////////////////
void CoordinateFile::Read(std::ifstream& in_file)
{
    string line;

    // Unable to read file
    if (!getline(in_file, line))
    {
        throw CoordinateFileProcessingException("Error reading file");
    }

    /// Set the tile by the first read line
    title_ = line;

    /// Extract the number of coordinates in the file
    getline(in_file, line);                         /// Read the next line
    int number_of_coordinates;
    stringstream ss(line);                          /// Create a stream from the read line
    ss >> number_of_coordinates;
    number_of_coordinates_ = number_of_coordinates; /// Set the number of coordinates attribute

    getline(in_file, line);                         /// Read the next line
    while(!Trim(line).empty())                      /// Read until the end of the file
    {
        // Tokenizing the read line
        boost::char_separator<char> separator(" ");
        boost::tokenizer< boost::char_separator<char> > tokens(line, separator);
        vector<string> vectorTokens = vector<string>();
        vectorTokens.assign(tokens.begin(), tokens.end());
        switch(vectorTokens.size())
        {
            /// One coordinate in the read line
            case 3:
                coordinates_.push_back(new Coordinate(ConvertString<double>(vectorTokens.at(0)), ConvertString<double>(vectorTokens.at(1)),
                                                      ConvertString<double>(vectorTokens.at(2))));
                break;
            /// Two coordinates in the read line
            case 6:
                coordinates_.push_back(new Coordinate(ConvertString<double>(vectorTokens.at(0)), ConvertString<double>(vectorTokens.at(1)),
                                                      ConvertString<double>(vectorTokens.at(2))));
                coordinates_.push_back(new Coordinate(ConvertString<double>(vectorTokens.at(3)), ConvertString<double>(vectorTokens.at(4)),
                                                      ConvertString<double>(vectorTokens.at(5))));
                break;
        }
        getline(in_file, line);
    }
    if((int)coordinates_.size() != number_of_coordinates_)
    {
        throw CoordinateFileProcessingException(__LINE__, "Corrupted file");
    }
}
void CoordinateFile::Write(const string &coordinate_file)
{
    std::ofstream out_file;
    try
    {
        out_file.open(coordinate_file.c_str());
    }
    catch(...)
    {
        throw CoordinateFileProcessingException(__LINE__,"File could not be created");
    }
    try
    {
        this->BuildCoordinateFile(out_file);
    }
    catch(...)
    {
        out_file.close();
    }
}
void CoordinateFile::BuildCoordinateFile(ofstream &stream)
{
    stream << left << setw(4) << GetTitle() << endl
           << right << setw(6) << GetNumberOfCoodinates() << endl;
    for(unsigned int i = 0; i < coordinates_.size(); i++)
    {
        stream << right << setw(12) << fixed << setprecision(7) << coordinates_.at(i)->GetX()
               << right << setw(12) << fixed << setprecision(7) << coordinates_.at(i)->GetY()
               << right << setw(12) << fixed << setprecision(7) << coordinates_.at(i)->GetZ();
        if(i%2 == 1)
            stream << endl;
    }
}

//////////////////////////////////////////////////////////
//                     DISPLAY FUNCTIONS                //
//////////////////////////////////////////////////////////
void CoordinateFile::Print(std::ostream& out)
{
    out << "*********** " << title_ << "***********" << endl;
    for(unsigned int i = 0; i < coordinates_.size(); i++)
    {
        out << (i+1) << ". " << coordinates_.at(i)->GetX() << ", " << coordinates_.at(i)->GetY() << ", " << coordinates_.at(i)->GetZ() << endl;
    }
}
