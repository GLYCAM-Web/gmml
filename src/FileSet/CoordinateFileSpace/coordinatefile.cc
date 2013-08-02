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

//////////////////////////////////// CONSTRUCTOR ///////////////////////////////////////
CoordinateFile::CoordinateFile(const string &crd_file)
{
    path_ = crd_file;
    std::ifstream in_file;
    try
    {
        in_file.open(crd_file.c_str());
    }
    catch(...)
    {
        throw CoordinateFileProcessingException(__LINE__,"File not found");
    }
    Read(in_file);
    in_file.close();            // Close the parameter files
}

///////////////////////////////////// ACCESSOR /////////////////////////////////////////
const std::string& CoordinateFile::GetFilePath() const
{
    return path_;
}

const std::string& CoordinateFile::GetTitle() const
{
    return title_;
}

int CoordinateFile::GetNumberOfCoodinates()
{
    return number_of_coordinates_;
}

const std::vector<Geometry::Coordinate*> CoordinateFile::GetCoordinates() const
{
    return coordinates_;
}

///////////////////////////////////// MUTATOR //////////////////////////////////////////
void CoordinateFile::SetPath(std::string& path)
{
    path_ = path;
}

void CoordinateFile::SetTitle(std::string& title)
{
    title_ = title;
}

void CoordinateFile::SetNumberOfCoordinates(int number_of_coordinates)
{
    number_of_coordinates_ = number_of_coordinates;
}

void CoordinateFile::SetCoordinates(std::vector<Coordinate*> coordinates)
{
    coordinates_.clear();
    for(vector<Coordinate*>::iterator it = coordinates.begin(); it != coordinates.end(); it++)
    {
        coordinates_.push_back(*it);
    }
}

void CoordinateFile::AddCoordinate(Geometry::Coordinate* coordinate)
{
    coordinates_.push_back(coordinate);
}

///////////////////////////////// FUNCTIONS ///////////////////////////////////////////
void CoordinateFile::Read(std::ifstream& in_file)
{
    string line;

    // Unable to read file
    if (!getline(in_file, line))
    {
        throw CoordinateFileProcessingException("Error reading file");
    }

    title_ = line;

    getline(in_file, line);
    int number_of_coordinates;
    stringstream ss(line);
    ss >> number_of_coordinates;
    number_of_coordinates_ = number_of_coordinates;

    getline(in_file, line);
    while(!Trim(line).empty())
    {
        boost::char_separator<char> separator(" ");
        boost::tokenizer< boost::char_separator<char> > tokens(line, separator);
        vector<string> vectorTokens = vector<string>();
        vectorTokens.assign(tokens.begin(), tokens.end());
        switch(vectorTokens.size())
        {
            case 3:
                coordinates_.push_back(new Coordinate(ConvertString<double>(vectorTokens.at(0)), ConvertString<double>(vectorTokens.at(1)),
                                                      ConvertString<double>(vectorTokens.at(2))));
                break;
            case 6:
                coordinates_.push_back(new Coordinate(ConvertString<double>(vectorTokens.at(0)), ConvertString<double>(vectorTokens.at(1)),
                                                      ConvertString<double>(vectorTokens.at(2))));
                coordinates_.push_back(new Coordinate(ConvertString<double>(vectorTokens.at(3)), ConvertString<double>(vectorTokens.at(4)),
                                                      ConvertString<double>(vectorTokens.at(5))));
                break;
        }
        getline(in_file, line);
    }
    if(coordinates_.size() != number_of_coordinates_)
    {
        throw CoordinateFileProcessingException(__LINE__, "Corrupted file");
    }
}

///////////////////////////////// DISPLAY FUNCTION ////////////////////////////////////
void CoordinateFile::Print(std::ostream& out)
{
    out << "*********** " << title_ << "***********" << endl;
    for(unsigned int i = 0; i < coordinates_.size(); i++)
    {
        out << (i+1) << ". " << coordinates_.at(i)->x_ << ", " << coordinates_.at(i)->y_ << ", " << coordinates_.at(i)->z_ << endl;
    }
}
