#ifndef NOTE_HPP
#define NOTE_HPP

#include <string>

namespace Glycan
{
    /*! \enum
      * Note type enumerator
      */
    enum NoteType
    {
        COMMENT,
        ERROR,
        WARNING
    };
    /*! \enum
      * Note category enumerator
      */
    enum NoteCat
    {
        MONOSACCHARIDE,
        GLYCOSIDIC,
        ANOMERIC,
        DER_MOD,
        RESIDUE_NAME
    };
    struct Note {
            NoteType type_;             /*!< The type of the issue, it can be a comment, warning, error etc. >*/
            NoteCat category_;          /*!< The category of the note>*/
            std::string description_;   /*!< The description of the note for a specific structure>*/
            /*! \fn
              * Convert a NoteCat enumerator to the string version of it
              * @param category The NoteCat enumerator that has to be converted to string
              * @return String format of the given NoteCat enumerator
              */
           std::string ConvertGlycanNoteCat2String(NoteCat category)
            {
                switch(category)
                {
                    case MONOSACCHARIDE:
                        return "Monosaccharide";
                    case GLYCOSIDIC:
                        return "Glycosidic Linkage";
                    case ANOMERIC:
                        return "Anomeric";
                    case DER_MOD:
                        return "Derivative/Modification";
                    case RESIDUE_NAME:
                        return "Residue Name";
                    default:
                        return "Unknown";
                }
            }
            /*! \fn
              * Convert a NoteType enumerator to the string version of it
              * @param type The NoteType enumerator that has to be converted to string
              * @return String format of the given NoteType enumerator
              */
            std::string ConvertGlycanNoteType2String(NoteType type)
            {
                switch(type)
                {
                    case COMMENT:
                        return "Comment";
                    case WARNING:
                        return "Warning";
                    case ERROR:
                        return "Error";
                    default:
                        return "Unknown";
                }
            }
    } ;
}

#endif // NOTE_HPP
