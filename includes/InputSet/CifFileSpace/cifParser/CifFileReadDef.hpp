/*
FILE:     CifFileReadDef.hpp
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
** \file CifFileReadDef.hpp
**
** \brief Header file for CifFileReadDef class.
*/


/* 
  PURPOSE:    Definitions for selective parsing/reading cif file
*/

#ifndef CIFFILEREADDEF_HPP
#define CIFFILEREADDEF_HPP

#include <string>
#include <vector>


#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>

enum type {A, D}; //A-list of accepted categorys/datablocks; D-denied

#define  INVALID_NUM_CATS -1

/**
** \class CifFileReadDef
**
** \brief Private class that represents a CIF parser controller.
*/
class CifFileReadDef
{
 private:
  int _numCatsToRead;
  int _numReadCats;
  void SetNumCatsToRead();
 protected:
  std::vector<std::string> _datablocklist;
  std::vector<std::string> _categorylist;
  type _datablocklisttype;
  type _categorylisttype;


 public:
  CifFileReadDef(std::vector<std::string> dblist,std::vector<std::string>clist,type dbtype = A, type ctype = A);
  CifFileReadDef(){_numCatsToRead = INVALID_NUM_CATS; _numReadCats = 0;};
  ~CifFileReadDef(){};
 
  void SetDataBlockList(std::vector<std::string> dblist,type dbtype = A);
  void SetCategoryList(std::vector<std::string>clist, type ctype = A);

  void SetDataBlockListType(type dbtype = A) { _datablocklisttype=dbtype;};
  void SetCategoryListType(type ctype = A){_categorylisttype=ctype;};

  int AreAllCatsRead();
  void IncreaseNumReadCats();

  int Category_OK(const std::string& categoryName);
  int Datablock_OK(const std::string& datablockName);
};
#endif
