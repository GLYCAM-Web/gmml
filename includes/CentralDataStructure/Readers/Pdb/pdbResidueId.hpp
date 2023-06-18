#ifndef INCLUDES_INPUTSET_PDBFILE_PDBRESIDUEID_HPP
#define INCLUDES_INPUTSET_PDBFILE_PDBRESIDUEID_HPP

#include <string>
#include <vector>

namespace pdb
{ // This class is a way for me to centralize where res info is extracted from a pdb file line.
  // This is used by PreprocessorInformation and PdbResidue.

    class ResidueId
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        ResidueId()
        {}

        ResidueId(const std::string& line);
        ResidueId(const std::string name, const std::string number, const std::string insertionCode,
                  const std::string chainId);
        ResidueId(std::vector<std::string> inputVector);

        //////////////////////////////////////////////////////////
        //                       ACCESSORS                      //
        //////////////////////////////////////////////////////////
        inline std::string getName() const
        {
            return residueName_;
        }

        inline std::string getNumber() const
        {
            return sequenceNumber_;
        }

        inline std::string getInsertionCode() const
        {
            return insertionCode_;
        }

        inline std::string getChainId() const
        {
            return chainId_;
        }

        inline std::string getAlternativeLocation() const
        {
            return alternativeLocation_;
        }

        //////////////////////////////////////////////////////////
        //                       MUTATORS                       //
        //////////////////////////////////////////////////////////
        inline void setName(const std::string s)
        {
            residueName_ = s;
        }

        inline void setNumber(const std::string s)
        {
            sequenceNumber_ = s;
        }

        inline void setInsertionCode(const std::string s)
        {
            insertionCode_ = s;
        }

        inline void setChainId(const std::string s)
        {
            chainId_ = s;
        }

        inline void setAlternativeLocation(const std::string s)
        {
            alternativeLocation_ = s;
        }

        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        const std::string getNumberAndInsertionCode() const;

        //////////////////////////////////////////////////////////
        //                       OVERLOADS                      //
        //////////////////////////////////////////////////////////
        bool operator==(const ResidueId& rhs)
        {
            return ((residueName_ == rhs.residueName_) && (sequenceNumber_ == rhs.sequenceNumber_) &&
                    (insertionCode_ == rhs.insertionCode_) && (chainId_ == rhs.chainId_));
        }

        bool operator!=(const ResidueId& rhs)
        {
            return ((residueName_ != rhs.residueName_) || (sequenceNumber_ != rhs.sequenceNumber_) ||
                    (insertionCode_ != rhs.insertionCode_) || (chainId_ != rhs.chainId_));
        }

        friend std::ostream& operator<<(std::ostream& os, const ResidueId& rId)
        {
            os << rId.print();
            return os;
        }

        //////////////////////////////////////////////////////////
        //                       DISPLAY                        //
        //////////////////////////////////////////////////////////
        std::string print() const;

      private:
        //////////////////////////////////////////////////////////
        //                       ATTRIBUTES                     //
        //////////////////////////////////////////////////////////
        std::string residueName_         = "";
        std::string sequenceNumber_      = "";
        std::string insertionCode_       = "";
        std::string chainId_             = "";
        std::string alternativeLocation_ = "";
    };
} // namespace pdb
#endif
