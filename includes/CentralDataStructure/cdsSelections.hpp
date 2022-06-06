#ifndef INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_HPP
#define INCLUDES_CENTRALDATASTRUCTURE_SELECTIONS_HPP

#include <iostream>
#include <string>
#include <utility>
#include <vector>
namespace cds
{
template< class residueT>
   typename std::vector<residueT*> selectProteinResidues(vectoriterator start, vectoriterator end);

} // namespace

#endif // CDS SELECTIONS

