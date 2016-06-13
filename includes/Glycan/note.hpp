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
                        return "monosaccharide";
                    case GLYCOSIDIC:
                        return "glycosidic linkage";
                    case ANOMERIC:
                        return "anomeric";
                    case DER_MOD:
                        return "derivative/modification";
                    case RESIDUE_NAME:
                        return "residue name";
                    default:
                        return "unknown";
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
                        return "comment";
                    case WARNING:
                        return "warning";
                    case ERROR:
                        return "error";
                    default:
                        return "unknown";
                }
            }
    } ;
}

#endif // NOTE_HPP
