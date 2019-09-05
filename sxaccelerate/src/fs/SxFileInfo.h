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

// Defining the header-guard
#ifndef _SX_FILE_INFO_H_
#define _SX_FILE_INFO_H_

// Including header-files
#include <SxString.h>
#include <SxFS.h>
#include <iostream>
#include <fstream>

// Defining a macro that saves lots of lines
#define TRY_N_PRINT_EX(func) try {\
   (func);\
} catch (SxException ex)  {\
   ex.print ();\
}
#define SX_DBG_TRACE(param) {;}/*{\
                                 ofstream outfile;\
                                 outfile.open ("sxdebug.log",\
                                 ios_base::app);\
                                 outfile << SxString((const char *)(param))\
                                 << std::endl;\
                                 outfile.close ();\
                                 }*/

#define SX_FS_CACHING_ENABLED

#ifndef WIN32
typedef uid_t sxuid_t;
typedef gid_t sxgid_t;
#else
typedef int sxuid_t;
typedef int sxgid_t;
#endif /*not WIN32*/

class SxDir;

/** \brief Information about a file

  \b SxFileInfo = SPHInX File and Directory Information class

  This class is used to inquire information about files and directories

  \ingroup group_os
  \author  Sixten Boeck, boeck@mpie.de
  \author  Thomas Uchdorf, t.uchdorf@mpie.de
 */
class SX_EXPORT_FS SxFileInfo
{
   public:

      enum AccessGroup { Owner=0x01, Group=0x02, Other=0x04, All=0x07 };
      enum AccessType  { Readable=0x01, Writable=0x02, Executable=0x04 };

      static const size_t MAX_PATH_LEN = 10240;
      static const size_t BUFFER_SIZE  = 10240;

      // --- Constructors
      // Constructor
      /** \brief Initialization with a relative or absolute path. */
      SxFileInfo (const SxString &path_ = SxString ()); // #exception

      SxFileInfo (const SxFileInfo &);

      // Destructor
      virtual ~SxFileInfo ();

      /** \brief Returns the separator used for navigation between files,
        directories and symbolic links on the current platform, i.e., '/' on
        unix based system and '\' on windows */
      static char getSeparator ();

      static bool equals (const SxFileInfo &,
                          const SxFileInfo &);
      static bool equals (const SxString &,
                          const SxString &);

      static int compareModificationTime (const SxFileInfo &lhs,
                                          const SxFileInfo &rhs);

      static int compareSize (const SxFileInfo &lhs,
                              const SxFileInfo &rhs);

      // --- Methods
      void setDirty () const;

      /** \brief Check whether a certain file exists.

        This function tries to check, whether a file exist or not.
       */
      virtual bool exists () const;

      /** \brief Check if entry is a file.

        This function returns true if the given entry is a file and if
        it exists. */
      bool isFile () const;

      /** \brief Check if entry is a symbolic link.

        This function returns true if the given entry is a symbolic link. */
      bool isSymLink () const;

      /** \brief Check if entry is a directory.

        This function returns true if the given entry is a directory and if
        it exists. */
      bool isDir () const;

      /** \brief Checks for read permisson.

        This function checks if a certain AccessGroup has read permission
        on a file or directory.

        \par Example 1:
        Check if the current user has permission
        \code
        SxFileInfo info = ...;
        if (info.isReadable())  ...;
        \endcode

        \par Example 2:
        Check if users and group have permission, but not all others
        \code
        SxFileInfo info = ...;
        if (info.isReadable(SxFileInfo::User | SxFileInfo::Group))  ...;
        \endcode */
      bool isReadable (AccessGroup access=Owner) const;

      /** \brief Checks for write permisson.

        This function checks if a certain AccessGroup has write permission
        on a file or directory.

        \par Example 1:
        Check if the current user has permission
        \code
        SxFileInfo info = ...;
        if (info.isWritable())  ...;
        \endcode

        \par Example 2:
        Check if users and group have permission, but not all others
        \code
        SxFileInfo info = ...;
        if (info.isWritable(SxFileInfo::User | SxFileInfo::Group))  ...;
        \endcode */
      bool isWritable (AccessGroup access=Owner) const;

      /** \brief Checks for execution permisson.

        This function checks if a certain AccessGroup has executation
        permission on a file or directory.

        \par Example 1:
        Check if the current user has permission
        \code
        SxFileInfo info = ...;
        if (info.isExecutable())  ...;
        \endcode

        \par Example 2:
        Check if users and group have permission, but not all others
        \code
        SxFileInfo info = ...;
        if (info.isExecutable(SxFileInfo::User | SxFileInfo::Group))  ...;
        \endcode */
      bool isExecutable (AccessGroup access = Owner) const;

      SxString getPerms () const;
      // int version of getPerms. returns octal
      int getMode () const;
      static int getMode (const SxString &);
      static SxString getPerms (int);

      /** \brief Return time of last read access

        Whenever a read access on a file or directory occurs its
        modification time is stamped. This function returns the time of
        last read access.
        It returns the number of seconds since 1.1.1970. */
      unsigned long lastAccessed () const;
      // Same as lastAccessed, but with precision to miliseconds
      unsigned long long lastAccessedMS () const;

      /** \brief Return time of last write access

        Whenever a write access on a file or directory occurs its
        modification time is stamped. This function returns the time of
        last write access.
        It returns the number of seconds since 1.1.1970. */
      unsigned long lastModified () const;
      // Same as lastModified, but with precision to miliseconds
      unsigned long long lastModifiedMS () const;




      static SxString getName (const SxString &path);
      static SxString getPath (const SxString &path);

      /** \brief Return current file or directory name */
      SxString getName () const;

      /** \brief Return the current full absolute path without the file name */
      SxString getPath () const;

      /** \brief Returns the current file containing the full absolute path */
      SxString const & getAbsPath () const;

      /** \brief Returns the relative path in respect to either the current
        working directory, or a specified one. Be careful with using this
        function as it simply truncates the absolute path string. */
      SxString getRelPath (const SxString &path = ".") const;
      SxString getRelPath (const SxDir    &path) const;

      /** \brief Returns the absolute uniform resource identifier */
      SxString getAbsURI () const;

      /** \brief Returns the size of an entry.

        SxFileInfo::getSize returns the total size of a file in bytes. */
      int64_t getSize () const;

      /** \brief Return the UID of a file or directory.

        Every file or directory has a UID (user identification number)
        and a GID (group identification number). This function returns
        the first one.
        \return GID or -1 if it was not possible to retrieve the GID */
      sxuid_t getUID () const;

      /** \brief Return the GID of a file or directory.

        Every file or directory has a UID (user identification number)
        and a GID (group identification number). This function returns
        the latter one.
        \return GID or -1 if it was not possible to retrieve the GID */
      sxgid_t getGID () const;

      //int getMode () const;
      //SxString getModeStr () const;

      SxString const & getOrig () const;

      // --- Class-functions
      static bool isAbs (const SxString &);

      static SxString getHomeStr ();

      static SxString getTmpStr ();

      /** \brief Returns simplified path without // */
      static SxString simplifyPath (const SxString &path);

      /** \brief Returns simplified path without // or .. or .
        The path does not need to exist.
        Symbolic links are not followed.
        ex: /usr/sbin/../lib/ -> /usr/lib
                       ../lib -> lib
                        ./lib -> lib */
      static SxString resolvePath (const SxString &path);

      bool operator== (const SxFileInfo &) const;
      bool operator!= (const SxFileInfo &) const;
   protected:

      // --- Members
      SxString orig;

      class SX_EXPORT_FS Cache {

         public:
            enum Type {Nothing = 0x00, File = 0x01, Dir = 0x02, SymLink = 0x04};

            // --- Constructors
            // Constructor
            Cache ();

            Cache (const Cache &rhs);

            // Destructor
            ~Cache ();

            // --- Methods
            const bool & isDirty () const;

            SxString const & getAbsPath () const;

            const bool & isExisting ();

            const int & getStatus ();

            const Type & getType ();

            const sxuid_t & getUID ();

            const sxgid_t & getGID ();

            const int & getMode ();

            const int64_t &getSize ();

            const unsigned long      &getAccessTime ();
            const unsigned long long &getAccessTimeMS ();

            const unsigned long      &getModificationTime ();
            const unsigned long long &getModificationTimeMS ();

            void setDirty ();

            void setAbsPath (const SxString &);

         protected:

            // --- Methods
            int getStatusInfo (SxStructStat *) const;

            int getStatusInfoCompletelyResolved (SxStructStat *) const;

            void update ();

            void setStatus (int);

            void setTypeAttrFile (int);

            void setTypeAttrDir (int);

            void setTypeAttrSymLink ();

            void setType (Type type_);

            void setExisting (bool);

            void setUID (sxuid_t);

            void setGID (sxgid_t);

            void setMode (int);

            void setSize (int64_t);

            void setAccessTime   (unsigned long);
            void setAccessTimeMS (unsigned long long);

            void setModificationTime   (unsigned long);
            void setModificationTimeMS (unsigned long long);


            // --- Members
            bool dirty;

            Type type;

            bool existing;

            int status;

            sxuid_t userID;

            sxgid_t groupID;

            int mode;

            int64_t size;

            unsigned long      accessTime;
            unsigned long long accessTimeMS;

            unsigned long      modificationTime;
            unsigned long long modificationTimeMS;

            /** \brief Represents the absolute path to a file system node
              (file, directory or symbolic link)*/
            SxString absPath;
      };
      mutable SxPtr<Cache> cache;
      //mutable Cache cache;

      /* SxString newRoot; */
};

class SX_EXPORT_FS SxFISortedByTime : public SxFileInfo
{
   public:
      // --- Constructors
      // Constructor
      SxFISortedByTime (const SxString &path_ = SxString ());

      SxFISortedByTime (const SxFileInfo &);

      // Destructor
      virtual ~SxFISortedByTime ();

      // --- Methods
      bool operator== (const SxFISortedByTime &rhs) const;

      bool operator> (const SxFISortedByTime &rhs) const;

      bool operator< (const SxFISortedByTime &rhs) const;
};

class SX_EXPORT_FS SxFISortedBySize : public SxFileInfo
{
   public:
      // --- Constructors
      // Constructor
      SxFISortedBySize (const SxString &path_ = SxString ());

      SxFISortedBySize (const SxFileInfo &);

      // Destructor
      virtual ~SxFISortedBySize ();

      // --- Methods
      bool operator== (const SxFISortedBySize &rhs) const;

      bool operator> (const SxFISortedBySize &rhs) const;

      bool operator< (const SxFISortedBySize &rhs) const;
};

/** \brief Creates an absolute path to a new file system node, i.e., a file, a
  directory or a symbolic link, by concatenating the absolute path of a file
  system node and a string separated by a slash */
SX_EXPORT_FS SxFileInfo operator/ (const SxFileInfo &a, const SxString &b);

/** \brief Creates an absolute path to a new file system node, i.e., a file, a
  directory or a symbolic link, by concatenating the absolute path of a file
  system node and a string separated by a slash */
SX_EXPORT_FS SxFileInfo operator/ (const SxString &a, const SxFileInfo &b);

/** \brief Concatenation of absolute paths with strings without insertion of
  separating slashes. */
SX_EXPORT_FS SxFileInfo operator+ (const SxFileInfo &a, const SxString &b);

/** \brief Concatenation of absolute paths with strings without insertion of
  separating slashes. */
SX_EXPORT_FS SxFileInfo operator+ (const SxString &a, const SxFileInfo &b);

SX_EXPORT_FS std::ostream &operator<< (std::ostream &s, const SxFileInfo &in);

#ifdef WIN32
   SX_EXPORT_FS std::wostream& operator<< (std::wostream &s, const SxFileInfo &in);
#endif


#endif /* _SX_FILE_INFO_H_ */
