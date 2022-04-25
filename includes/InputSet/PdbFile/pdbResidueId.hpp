#ifndef INCLUDES_INPUTSET_PDBFILE_PDBRESIDUEID_HPP
#define INCLUDES_INPUTSET_PDBFILE_PDBRESIDUEID_HPP

#include <string>

namespace pdb
{
class ResidueId
{
public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    ResidueId() {}
    ResidueId(const std::string &line);
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
    std::string getId() const;
    const std::string getNumberWithCode() const;
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
