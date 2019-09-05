// ---------------------------------------------------------------------------
//
//           The general purpose cross platform C/C++ framework
//
//                       S x A c c e l e r a t e
//
//           Home:       https://www.sxlib.de
//           License:    Apache 2
//           Authors:    see src/AUTHORS
//
// ---------------------------------------------------------------------------

#ifndef _SX_FILE_PARSER_H_
#define _SX_FILE_PARSER_H_

#include <SxIO.h>
#include <SxString.h>
#include <SxArray.h>

/** \brief S/PHI/nX file parser for non-SPHInX files

    \b SxClass = S/PHI/nX file parser

    \author C. Freysoldt freysoldt@mpie.de */
class SX_EXPORT_IO SxFileParser
{
   public:
      /** \brief File pointer

        \note Do not manipulate the file directly via the file pointer.
        Use the file point only for fscanf, when the functionality provided
        by SxFileParser does not allow to perform the scan in an easy way.
        */
      FILE *fp;

   protected:
      /// Line number
      ssize_t line;

      /// Current position within file
      long fpos;

      /// File name
      SxString name;

      /// Update line number
      void updateLine ();

      /// Check that we are not at end of file
      void noEof ();

      /// Current topic
      SxString parserTopic;

   public:
      /// Whether progress should be printed
      bool verbose;

      /// Constructor
      SxFileParser ();

      /// Constructor (open a file)
      SxFileParser (const SxString &fileName);

      /// Destructor
      ~SxFileParser ();

      /// Open a file
      void open (const SxString &fileName);

      /// Close the file
      void close ();

      /// Print out current position
      void where ();

      /// Get current line number
      ssize_t getLineNo ()
      {
         SX_CHECK (fp);
         updateLine ();
         return line; 
      }

      /** \brief Set topic
          Topics allow to structure the file parsing.
          Topics are printed when errors occur and in verbose mode.
          Example:
          \code
   fp.topic ("wave functions");
   // read wave functions
   ...
   fp.topic ("potential");
   // read potential
   ...
          \endcode
        */
      void topic (const SxString &topicName)
      {
         parserTopic = topicName;
         if (verbose)
            cout << "Reading " << topicName << "..." << endl;
      }

      /// Read a double
      double getDouble ();

      /// Read a double
      SxFileParser& operator>> (double &x)
      {
         x = getDouble ();
         return *this;
      }

      /// Read a float
      SxFileParser& operator>> (float &x)
      {
         x = (float)getDouble ();
         return *this;
      }

      /// Read an integer
      int getInt ();

      /// Read an integer 
      SxFileParser& operator>> (int &x)
      {
         x = getInt ();
         return *this;
      }

      /// Read a long
      long getLong ();

      /// Read an integer 
      SxFileParser& operator>> (long &x)
      {
         x = getLong ();
         return *this;
      }

      /** \brief Read rest of line
        */
      SxString getLine ();

      /** \brief Get to next line
          @param nLine   number of lines to skip
        */
      void nextLine (int nLine = 1);

      /// Read next line
      SxString readNextLine ()
      {
         SX_CHECK (fp);
         nextLine ();
         return getLine ();
      }

      /// Read until a specific character occurs
      void readUntil (char what);

      /// Skip white space
      void skipWhite ();

      /** Read a specific text
        */
      void read (const char *what);

      /** \brief Try to read a specific text.
          If this fails, nothing is read from file, and return value is false
          If it succeeds, next read starts after the specific text, return value is true.
        */
      bool reads (const char *what);

      /** \brief Read a vector
          Reads as many real numbers as possible.
        */
      SxArray<double> getVector ();

      /// Read a vector of known size
      SxArray<double> getVector (ssize_t n);

      /** \brief Find text

          Next read starts after this text
        */
      void find (const char *what);

};

#endif /* _SX_FILE_PARSER_H_ */
