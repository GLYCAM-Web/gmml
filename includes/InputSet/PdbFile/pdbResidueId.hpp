#ifndef INCLUDES_INPUTSET_PDBFILE_PDBRESIDUEID_HPP
#define INCLUDES_INPUTSET_PDBFILE_PDBRESIDUEID_HPP

#include <string>

namespace pdb
{ // This class is a way for me to centralize where res info is extracted from a pbd file line.
  // This is used by PreprocessorInformation and PdbResidue.
class ResidueId
{
public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    ResidueId() {}
    ResidueId(const std::string &line);
    ResidueId(const std::string &name, const std::string &number, const std::string &insertionCode, const std::string &chainId);
    //////////////////////////////////////////////////////////
    //                       ACCESSORS                      //
    //////////////////////////////////////////////////////////
    inline const std::string& getName() const {return residueName_;}
    inline const std::string& getNumber() const {return sequenceNumber_;}
    inline const std::string& getInsertionCode() const {return insertionCode_;}
    inline const std::string& getChainId() const {return chainId_;}
    //////////////////////////////////////////////////////////
    //                       MUTATORS                       //
    //////////////////////////////////////////////////////////
    inline void setName(const std::string& s) {residueName_ = s;}
    inline void setNumber(const std::string& s) {sequenceNumber_ = s;}
    inline void setInsertionCode(const std::string& s) {insertionCode_ = s;}
    inline void setChainId(const std::string& s) {chainId_ = s;}
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
    const std::string getNumberAndInsertionCode() const;
    //////////////////////////////////////////////////////////
    //                       OVERLOADS                      //
    //////////////////////////////////////////////////////////
    bool operator==( ResidueId& rhs)
        { return (
            (residueName_ == rhs.residueName_) &&
            (sequenceNumber_ == rhs.sequenceNumber_) &&
            (insertionCode_ == rhs.insertionCode_) &&
            (chainId_ == rhs.chainId_));
        }
    bool operator!=( ResidueId& rhs)
            { return (
                (residueName_ != rhs.residueName_) ||
                (sequenceNumber_ != rhs.sequenceNumber_) ||
                (insertionCode_ != rhs.insertionCode_) ||
                (chainId_ != rhs.chainId_));
            }
    //////////////////////////////////////////////////////////
    //                       DISPLAY                        //
    //////////////////////////////////////////////////////////
    std::string print() const;
private:
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////
    std::string residueName_;
    std::string sequenceNumber_;
    std::string insertionCode_;
    std::string chainId_;
};
}
#endif
