/*
FILE:     CifDataInfo.hpp
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


/**
** \file CifDataInfo.hpp
**
** CIF data information class
*/


#ifndef CIFDATAINFO_HPP
#define CIFDATAINFO_HPP


#include <string>
#include <vector>

#include "../common/DataInfo.hpp"
#include "DicFile.hpp"


class CifDataInfo : public DataInfo
{
  public:
    CifDataInfo(DicFile& dictFile);
    ~CifDataInfo();

    void GetVersion(std::string& version);

    const std::vector<std::string>& GetCatNames();

    const std::vector<std::string>& GetItemsNames();

    bool IsCatDefined(const std::string& catName) const;

    bool IsItemDefined(const std::string& itemName);

    const std::vector<std::string>& GetCatKeys(const std::string& catName);

    const std::vector<std::string>& GetCatAttribute(const std::string& catName,
      const std::string& refCatName, const std::string& refAttrName);

    const std::vector<std::string>&
      GetItemAttribute(const std::string& itemName,
      const std::string& refCatName, const std::string& refAttrName);

    virtual void GetCatItemsNames(std::vector<std::string>& itemsNames,
      const std::string& catName);

  protected:
    DicFile& _dictFile;

  private:
    std::string _version;
    std::vector<std::string> _catsNames;
    std::vector<std::string> _itemsNames;
    std::vector<std::string> _catKeyItems;
    std::vector<std::string> _catAttrib;
    std::vector<std::string> _itemAttrib;
    std::vector<std::string> _itemTypeListAttrib;

    void _GetDictVersion(std::string& dictVer);
    bool _isDictCategory(const std::string& category) const;

    const std::vector<std::string>&
      GetItemAttributeForItemTypeListCat(const std::string& itemName,
      const std::string& refCatName,
      const std::string& refAttrName);
};


#endif
