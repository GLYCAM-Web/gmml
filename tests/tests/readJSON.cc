#include "../../gmml/includes/gmml.hpp"
#include "../../gmml/includes/External_Libraries/json.hpp"
#include <iostream>
#include <string>

int main()
{


struct selectedRotamers
{
    int linkageIndex;
    std::string dihedralName;
    StringVector angles;
};    

    auto j2 = R"(
{
   "entity": {
       "type": "sequence"
   },
   "responses": [
       {
           "Evaluate": {
               "glycosidicLinkages": [
                   {
                       "20": {
                           "likelyRotamers": {
                               "omg": [
                                   "gg",
                                   "gt"
                               ]
                           },
                           "linkageName": "DManpa1-6DManpb",
                           "possibleRotamers": {
                               "omg": [
                                   "gg",
                                   "gt"
                               ]
                           }
                       },
                       "24": {
                           "likelyRotamers": {
                               "omg": [
                                   "gg",
                                   "gt",
                                   "tg"
                               ],
                               "phi": [
                                   "t",
                                   "-g"
                               ]
                           },
                           "linkageName": "DNeu5Acpa2-6DGalpb",
                           "possibleRotamers": {
                               "omg": [
                                   "gg",
                                   "gt",
                                   "tg"
                               ],
                               "phi": [
                                   "t",
                                   "-g"
                               ]
                           }
                       },
                       "25": {
                           "likelyRotamers": {
                               "omg": [
                                   "gg",
                                   "gt"
                               ]
                           },
                           "linkageName": "DManpa1-6DManpa",
                           "possibleRotamers": {
                               "omg": [
                                   "gg",
                                   "gt"
                               ]
                           }
                       }
                   }
               ],
               "inputSequence": "DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeu5Aca2-6DGalpb1-4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH",
               "officialSequence": "DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeu5Aca2-6DGalpb1-4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH"
           }
       }
   ]
}
)"_json;

    //std::string condensed_sequence = "DManpb1-2DManp[6A]b1-2[DGalpb1-6]DManpb1-2DManp[6Me]b1-2[DGlcpb1-6DManp[2S]a1-4]DGalpa1-OME";
    //CondensedSequenceSpace::carbohydrateBuilder steve2("build", condensed_sequence);
    //steve2.WriteFile("PDB", "output2.pdb");
}
//prep file

