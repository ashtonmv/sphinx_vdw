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

#ifndef _SX_EXCEPTION_BASE_H_
#define _SX_EXCEPTION_BASE_H_

#include <SxUtil.h>
#include <SxMacroLib.h>
#include <SxString.h>
#include <atomic>

extern atomic<bool> dumpCoreFile;

/** \brief SxAccelerate Exception class

    \b SxException = SxAccelerate Exception Handling Class

    This is a standard exception class for all SxAccelerate exceptions.
    Exception handling in C++ can cause spurious bus errors when they are
    not caught properly in a try-catch block. If SFHIngX aborts with a bus
    error which could be related to such a missing try-catch block this class
    allows to handle the exception with an SFHIngX SX_EXIT statement instead.
    In the DEBUG mode the resulting core file can be analyzed
    (\ref dbg_analysis).

    \author Sixten Boeck
  */
class SX_EXPORT_UTIL SxException
{
   public:

      // print mode for toString()
      enum PrintMode {
         Simple=0x00,      ///< top-most exception w/o sender info
         Stack=0x01,       ///< exception stack w/o sender info
         DebugSimple=0x02, ///< top-most exception w/ sender info
         DebugStack=0x03,  ///< exception stack w/ sender info
      };

      /// \brief Standard constructor
      SxException ();
      /// \brief Copy constructor
      SxException (const SxException &);
      SxException (const SxString &message,
                   const SxString &senderFile="unknown",
                   int senderLine=0);
      SxException (const SxString &tag,
                   const SxString &message,
                   const SxString &file="unknown", int senderLine=0);
      SxException (const SxException &chain, const SxString &message,
                   const SxString &senderFile="unknown", int senderLine=0);
      SxException (const SxException &chain,
                   const SxString &tag, const SxString &message,
                   const SxString &senderFile="unknown", int senderLine=0);

      /// \brief Destructor
      ~SxException ();

      /// \brief assignment operator
      SxException &operator= (const SxException &);

      /** \brief Print the message body

          if \em printSender is true additionally the file and the line
          from where this exception was trown will be printed. */
      void print (bool printSender=false) const;

      SxString toString (int mode = (int)Simple,
                         const SxString &delimiter="\n") const;

      void printStack (const SxString &delimiter="\n") const;

      bool hasTag (const SxString &tag_) const;   // current instance only
      bool findTag (const SxString &tag_) const;  // search in chain recusively


      /// \brief returns the message body.
      const SxString &getMessage () const  { return message; }

      /// \brief returns the sender file.
      const SxString &getFileName () const  { return senderFile; }

      const SxString &getTag () const { return tag; }

      /// \brief returns the sender line.
      int getLineNumber () const { return senderLine; }

      /// \brief returns next exception in the chain.
      const SxException *getNext () const { return next; }


      /** \brief Treat all exception as SX_EXIT statements.

          Usually an uncaught exception causes an "Abort" signal. So no
          core file will be prompted and backtracing in the debug mode 
          becomes impossible.
          By calling this static function every exception violates the
          memory explicily in order to dump a qualified core file.
          For debugging only
          \sa SxException::throwExceptions
          \sa SxException::dumpCoreFile */
      static void causeSegFault   ();
      /** \brief Treat all exceptions normally (no SX_EXIT emulation)

          For debugging only. Switch off dumping of core files rather
          than thwing exceptions.
          \sa SxException::causeSegFault
          \sa SxException::dumpCodeFile */
      static void throwExceptions ();

      /** Program global variable.

          If SxException::dumpCoreFile is set to true all exception are
          treated with a SFHIngX SX_EXIT statement instead. This allows the
          developer to trace missing try-catch blocks.

          This variable shouldn't be set directly. Use better
          SxException::causeSegFaults or SxException::throwExceptions.

          \code
             SxException::dumpSegFaults();
          \endcode

          \sa \ref page_debug  */
      static atomic<bool> dumpCoreFile;

   protected:

      /// \brief The exception message body
      SxString message;

      /// \brief Optional, exception key, such as Permission, BrokenPipe
      SxString tag;

      /// \brief Contains the file name from where the exception was trown.
      SxString senderFile;

      /// \brief Contains the line from where the exception was trown.
      int senderLine;

      /// \brief Next exception in the chain
      SxException *next;

      /// \brief previous exception in the chain
      SxException *prev;

      /// \brief last exception in the chain
      SxException *tail;

      void init (const SxException &);
      void init (const SxString &message,
                 const SxString &tag,
                 const SxString &senderFile,
                 int            senderLine);
      void destroy ();
      void copy (const SxException &);
};

#define _SxThrow1(msg)                                                  \
   throw SxException(msg, __FILE__, __LINE__)
#define _SxThrow2(var,msg)                                              \
   throw SxException(var, msg, __FILE__, __LINE__)

#define _SxThrow3(e,tag,msg)                                            \
   throw SxException(e, tag, msg, __FILE__, __LINE__)

#define SX_THROW(...) SX_VMACRO(_SxThrow, __VA_ARGS__)


#endif /* _SX_EXCEPTION_BASE_H_ */
