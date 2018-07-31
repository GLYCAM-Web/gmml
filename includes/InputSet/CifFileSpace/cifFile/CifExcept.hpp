/*
FILE:     CifExcept.hpp
*/
/*
VERSION:  7.105
*/
/*
DATE:     10/21/2013
*/
/*
  Comments and Questions to: sw-help@rcsb.rutgers.edu
*/
/*
COPYRIGHT 1999-2013 Rutgers - The State University of New Jersey
This SOFTWARE has been modified from 
the version originally obtained from RUTGERS.  
*/


/*!
** \file CifExcept.hpp
**
** \brief Header file for CifExcept class.
*/


#ifndef CIFEXCEPT_HPP
#define CIFEXCEPT_HPP


#include <string>


/**
**  \class CifExcept
**
**  \brief Static class that represents some exceptions in CIF files
**    related to data values.
**
*/
class CifExcept
{
  public:
    static bool CanBeUnknown(const std::string& itemName);
    static bool CanBeInapplicable(const std::string& itemName);
    static bool IsBadParentRelation(const std::string& itemName);
    static bool IsBadChildRelation(const std::string& itemName);
};


#endif
