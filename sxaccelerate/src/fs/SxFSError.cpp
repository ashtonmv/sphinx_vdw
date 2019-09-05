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

// Including header files
#include <SxFSError.h>
#include <errno.h>

#ifdef WIN32
#  include <shortcut.h>
#  include <windows.h>
#  include <unknwn.h>
#  include <shlobj.h>
#  include <objbase.h>
#  include <olectl.h>
#  include <initguid.h>
#endif /* WIN32 */

// --- Constructors
// Constructor
SxFSError::SxFSError ()
{
   // empty
}
// Destructor
SxFSError::~SxFSError ()
{
   // empty
}

// --- Class-functions

SxString SxFSError::getUtimeErrMsg ()
{
   SxString errMsg;
   switch (errno)  {
      case(EACCES):
                errMsg = "Either search permission is denied by a component ";
                errMsg += "of the path, or the argument for setting the time ";
                errMsg += "is not set correctly and write access is denied ";
                errMsg += "due to mismatching user IDs.";
                break;
#ifndef WIN32
      case(ELOOP):
                errMsg = "Too many symbolic links were encountered in ";
                errMsg += "resolving path.";
                break;
#endif /* not WIN32 */
      case(ENAMETOOLONG):
                errMsg = "The path exceeds the maximum path length or a ";
                errMsg += "component of it is longer than the maximum name ";
                errMsg += "length.";
                break;
      case(ENOENT):
                errMsg = "Either a component of path does not name an ";
                errMsg += "existing file, or path is an empty string.";
                break;
      case(ENOTDIR):
                errMsg = "A component of the path is not a directory.";
                break;
      case(EPERM):
                errMsg = "The times argument is not a null pointer and the ";
                errMsg += "calling process' effective user ID has write ";
                errMsg += "access to the file but does not match the owner ";
                errMsg += "of the file and the calling process does not have ";
                errMsg += "the appropriate privileges.";
                break;
      case(EROFS):
                errMsg = "The file system containing the file is read-only.";
                break;
      default:
                errMsg = "An unknown error occurred.";
                break;
   }
   return errMsg;
}

SxString SxFSError::getChmodErrMsg ()
{
   SxString errMsg;
   switch (errno)  {
      case (EACCES):
         errMsg = "Search permission is denied on a component of the path ";
         errMsg += "prefix.";
         break;
#ifndef WIN32
      case (ELOOP):
         errMsg = "Too many symbolic links were encountered in resolving ";
         errMsg += "path.";
         break;
#endif /* not WIN32 */
      case (ENAMETOOLONG):
         errMsg = "The length of the path exceeds the maximum one ";
         errMsg += "or a pathname component is longer than the maximum ";
         errMsg += "allowed name length.";
         break;
      case (ENOTDIR):
         errMsg = "A component of the path prefix is not a directory.";
         break;
      case (ENOENT):
         errMsg = "A component of path does not name an existing file or ";
         errMsg += "path is an empty string.";
         break;
      case (EPERM):
         errMsg = "The effective user ID does not match the owner of the ";
         errMsg += "file and the process does not have appropriate ";
         errMsg += "privileges.";
         break;
      case (EROFS):
         errMsg = "The named file resides on a read-only file system.";
         break;
      default:
         errMsg = "An unknown error occurred.";
         break;
   }
   return errMsg;
}

SxString SxFSError::getChownErrMsg ()
{
      SxString errMsg;
      switch (errno)  {
         case(EACCES):
         errMsg = "Search permission is denied on a component of the ";
         errMsg += "path prefix.";
         break;
#ifndef WIN32
         case(ELOOP):
         errMsg = "Too many symbolic links were encountered in ";
         errMsg += "resolving path.";
         break;
#endif /* not WIN32 */
         case(ENAMETOOLONG):
         errMsg = "The length of the path exceeds the maximum one ";
         errMsg += "or a pathname component is longer than the maximum ";
         errMsg += "allowed name length.";
         break;
         case(ENOTDIR):
         errMsg = "A component of the path prefix is not a directory.";
         break;
         case(ENOENT):
         errMsg = "A component of path does not name an existing file ";
         errMsg += "or path is an empty string.";
         break;
         case(EPERM):
         errMsg = "The effective user ID does not match the owner of ";
         errMsg += "the file, or the calling process does not have ";
         errMsg += "appropriate privileges.";
         break;
         case(EROFS):
         errMsg = "The named file resides on a read-only file system.";
         break;
         default:
         errMsg = "An unknown error occurred.";
         break;
      }
     return errMsg;
}

SxString SxFSError::getSymlinkErrMsg ()
{
   SxString errMsg;
   switch (errno)  {
      case(EACCES):
         errMsg = "Write permission is denied in the directory where the ";
         errMsg += "symbolic link is being created, or search permission ";
         errMsg += "is denied for a component of the path prefix of the ";
         errMsg += "desired link.";
         break;
      case(EEXIST):
         errMsg = "The desired link names an existing file or symbolic link.";
         break;
      case(EIO):
         errMsg = "An I/O error occurs while reading from or writing to the ";
         errMsg += "file system.";
         break;
#ifndef WIN32
      case(ELOOP):
         errMsg = "Too many symbolic links were encountered in resolving the ";
         errMsg += "desired link path.";
         break;
#endif /* not WIN32 */
      case(ENAMETOOLONG):
         errMsg = "The length of the path exceeds the maximum one or a ";
         errMsg += "pathname component is longer than the maximum allowed ";
         errMsg += "name length.";
         break;
      case(ENOENT):
         errMsg = "A component of path does not name an existing file or ";
         errMsg += "path is an empty string.";
         break;
      case(ENOSPC):
         errMsg = "The directory in which the entry for the new symbolic ";
         errMsg += "link is being placed cannot be extended because no space ";
         errMsg += "is left on the file system containing the directory, or ";
         errMsg += "the new symbolic link cannot be created because no space ";
         errMsg += "is left on the file system which will contain the link, ";
         errMsg += "or the file system is out of file-allocation resources.";
         break;
      case(ENOTDIR):
         errMsg = "A component of the path prefix of the desired link ";
         errMsg += "is not a directory.";
         break;
      case(EROFS):
         errMsg = "The new symbolic link would reside on a read-only file ";
         errMsg += "system.";
         break;
      default:
         errMsg = "An unknown error occurred.";
         break;
   }
   return errMsg;
}

SxString SxFSError::getRenameErrMsg ()
{
   SxString errMsg;
   switch (errno)  {
      case(EACCES):
                errMsg = "A component of either path prefix denies ";
                errMsg += "search permission; or one of the directories ";
                errMsg += "containing old or new denies write permissions; ";
                errMsg += "or, write permission is required and is denied ";
                errMsg += "for a directory pointed to by the old or new ";
                errMsg += "arguments.";
                break;
      case(EBUSY):
                errMsg = "The directory named by old or new is ";
                errMsg += "currently in use by the system or another ";
                errMsg += "process, and the implementation considers this ";
                errMsg += "an error.";
                break;
      case(EEXIST):
      case(ENOTEMPTY):
                errMsg = "The link named by new is a directory that ";
                errMsg += "is not an empty directory.";
                break;
      case(EINVAL):
                errMsg = "The new directory pathname contains a path ";
                errMsg += "prefix that names the old directory.";
                break;
      case(EIO):
                errMsg = "A physical I/O error has occurred.";
                break;
      case(EISDIR):
                errMsg = "The new argument points to a directory and ";
                errMsg += "the old argument points to a file that is not a ";
                errMsg += "directory.";
                break;
#ifndef WIN32
      case(ELOOP):
                errMsg = "Too many symbolic links were encountered ";
                errMsg += "in resolving either pathname.";
                break;
#endif /* not WIN32 */
      case(EMLINK):
                errMsg = "The file named by old is a directory, and ";
                errMsg += "the link count of the parent directory of new ";
                errMsg += "would exceed {LINK_MAX}.";
                break;
      case(ENAMETOOLONG):
                errMsg = "The length of the old or new argument ";
                errMsg += "exceeds {PATH_MAX} or a pathname component is ";
                errMsg += "longer than {NAME_MAX}.";
                break;
      case(ENOENT):
                errMsg = "The link named by old does not name an ";
                errMsg += "existing file, or either old or new points to ";
                errMsg += "an empty string.";
                break;
      case(ENOSPC):
                errMsg = "The directory that would contain new cannot ";
                errMsg += "be extended.";
                break;
      case(ENOTDIR):
                errMsg = "A component of either path prefix is not a ";
                errMsg += "directory; or the old argument names a directory ";
                errMsg += "and new argument names a non-directory file.";
                break;
      case(EPERM):
                errMsg = "The S_ISVTX flag is set on the directory ";
                errMsg += "containing the file referred to by old and the ";
                errMsg += "caller is not the file owner, nor is the caller ";
                errMsg += "the directory owner, nor does the caller have ";
                errMsg += "appropriate privileges; or new refers to an ";
                errMsg += "existing file, the S_ISVTX flag is set on the ";
                errMsg += "directory containing this file and the ";
                errMsg += "caller is not the file owner, nor is the caller ";
                errMsg += "the directory owner, nor does the caller have ";
                errMsg += "appropriate privileges.";
                break;
      case(EROFS):
                errMsg = "The requested operation requires writing in ";
                errMsg += "a directory on a read-only file system.";
                break;
      case(EXDEV):
                errMsg = "The links named by new and old are on ";
                errMsg += "different file systems and the implementation ";
                errMsg += "does not support links between file systems.";
                break;
      default:
                errMsg = "An unknown error occurred.";
                break;
   }

   return errMsg;
}

SxString SxFSError::getOpenErrMsg ()
{
   SxString errMsg;
   switch (errno)  {
      case(EACCES):
                errMsg = "Search permission is denied on a component of the ";
                errMsg += "path prefix, or the file exists and the ";
                errMsg += "permissions specified by oflag are denied, or the ";
                errMsg += "file does not exist and write permission is denied";
                errMsg += " for the parent directory of the file to be ";
                errMsg += "created, or O_TRUNC is specified and write ";
                errMsg += "permission is denied.";
                break;
      case(EEXIST):
                errMsg = "O_CREAT and O_EXCL are set, and the named file ";
                errMsg += "exists.";
                break;
      case(EINTR):
                errMsg = "A signal was caught during open().";
                break;
      case(EINVAL):
                errMsg = "The implementation does not support synchronised ";
                errMsg += "I/O for this file.";
                break;
      case(EIO):
                errMsg = "The path argument names a STREAMS file and a hangup";
                errMsg += " or error occurred during the open().";
                break;
      case(EISDIR):
                errMsg = "The named file is a directory and oflag includes ";
                errMsg += "O_WRONLY or O_RDWR.";
                break;
#ifndef WIN32
      case(ELOOP):
                errMsg = "Too many symbolic links were encountered in ";
                errMsg += "resolving path.";
                break;
#endif /* not WIN32 */
      case(EMFILE):
                errMsg = "{OPEN_MAX} file descriptors are currently open in ";
                errMsg += "the calling process.";
                break;
      case(ENAMETOOLONG):
                errMsg = "The length of the path argument exceeds {PATH_MAX} ";
                errMsg += "or a pathname component is longer than {NAME_MAX}.";
                break;
      case(ENFILE):
                errMsg = "The maximum allowable number of files is currently ";
                errMsg += "open in the system.";
                break;
      case(ENOENT):
                errMsg = "O_CREAT is not set and the named file does not ";
                errMsg += "exist; or O_CREAT is set and either the path ";
                errMsg += "prefix does not exist or the path argument points ";
                errMsg += "to an empty string.";
                break;
#ifndef WIN32
      case(ENOSR):
                errMsg = "The path argument names a STREAMS-based file and ";
                errMsg += "the system is unable to allocate a STREAM.";
                break;
#endif /* not WIN32 */
      case(ENOSPC):
                errMsg = "The directory or file system that would contain the";
                errMsg += " new file cannot be expanded, the file does not ";
                errMsg += "exist, and O_CREAT is specified.";
                break;
      case(ENOTDIR):
                errMsg = "A component of the path prefix is not a directory.";
                break;
      case(ENXIO):
                errMsg = "O_NONBLOCK is set, the named file is a FIFO, ";
                errMsg += "O_WRONLY is set and no process has the file open ";
                errMsg += "for reading.";
                errMsg = "The named file is a character special or block ";
                errMsg += "special file, and the device associated with this ";
                errMsg += "special file does not exist.";
                break;
#ifndef WIN32
      case(EOVERFLOW):
                errMsg = "The named file is a regular file and the size of ";
                errMsg += "the file cannot be represented correctly in an ";
                errMsg += "object of type off_t.";
                break;
#endif /* not WIN32 */
      case(EROFS):
                errMsg = "The named file resides on a read-only file system ";
                errMsg += "and either O_WRONLY, O_RDWR, O_CREAT (if file ";
                errMsg += "does not exist) or O_TRUNC is set in the oflag ";
                errMsg += "argument.";
                break;
      default:
                errMsg = "An unknown error occurred.";
                break;
   }
   return errMsg;
}

SxString SxFSError::getWriteErrMsg ()
{
   SxString errMsg;
   switch (errno)  {
      case(EAGAIN):
                errMsg = "The O_NONBLOCK flag is set for the file descriptor";
                errMsg += " and the thread would be delayed in the write() ";
                errMsg += "operation.";
                break;
      case(EBADF):
                errMsg = "The fildes argument is not a valid file descriptor";
                errMsg += " open for writing.";
                break;
      case(EFBIG):
                errMsg = "An attempt was made to write a file that exceeds the";
                errMsg += " implementation-dependent maximum file size or the ";
                errMsg += "process' file size limit or ";
                errMsg = "the file is a regular file, nbyte is greater than 0 ";
                errMsg += "and the starting position is greater than or equal ";
                errMsg += "to the offset maximum established in the open file ";
                errMsg += "description associated with fildes.";
                break;
      case(EINTR):
                errMsg = "The write operation was terminated due to the ";
                errMsg += "receipt of a signal, and no data was transferred.";
                break;
      case(EIO):
                errMsg = "A physical I/O error has occurred or ";
                errMsg = "the process is a member of a background process ";
                errMsg += "group attempting to write to its controlling ";
                errMsg += "terminal, TOSTOP is set, the process is neither ";
                errMsg += "ignoring nor blocking SIGTTOU and the process ";
                errMsg += "group of the process is orphaned. This error may ";
                errMsg += "also be returned under implementation-dependent ";
                errMsg += "conditions.";
                break;
      case(ENOSPC):
                errMsg = "There was no free space remaining on the device ";
                errMsg += "containing the file.";
                break;
      case(EPIPE):
                errMsg = "An attempt is made to write to a pipe or FIFO that ";
                errMsg += "is not open for reading by any process, or that ";
                errMsg += "only has one end open. A SIGPIPE signal will also ";
                errMsg += "be sent to the thread.";
                break;
      case(ERANGE):
                errMsg = "The transfer request size was outside the range ";
                errMsg += "supported by the STREAMS file associated with ";
                errMsg += "fildes.";
                break;
      default:
                errMsg = "An unknown error occurred.";
                break;
   }
   return errMsg;
}

SxString SxFSError::getChdirErrMsg ()
{
   SxString errMsg;
   switch (errno)  {
      case(EACCES):
                errMsg = "Search permission is denied for any component of the";
                errMsg += " pathname.";
                break;
#ifndef WIN32
      case(ELOOP):
                errMsg = "Too many symbolic links were encountered in ";
                errMsg += "resolving path.";
                break;
#endif /* not WIN32 */
      case(ENAMETOOLONG):
                errMsg = "The path argument exceeds {PATH_MAX} in length or ";
                errMsg += "a pathname component is longer than {NAME_MAX}.";
                break;
      case(ENOENT):
                errMsg = "A component of path does not name an existing ";
                errMsg += "directory or path is an empty string.";
                break;
      case(ENOTDIR):
                errMsg = "A component of the pathname is not a directory.";
                break;
      default:
                errMsg = "An unknown error occurred.";
                break;
   }
   return errMsg;
}

SxString SxFSError::getGetcwdErrMsg ()
{
   SxString errMsg;
   switch (errno)  {
      case(EINVAL):
                errMsg = "The size argument is 0.";
                break;
      case(ERANGE):
                errMsg = "The size argument is greater than 0, but is smaller";
                errMsg += " than the length of the pathname +1.";
                break;
      default:
                errMsg = "An unknown error occurred.";
                break;
   }
   return errMsg;
}

SxString SxFSError::getOpendirErrMsg ()
{
   SxString errMsg;
   switch (errno)  {
      case(EACCES):
                errMsg = "Search permission is denied for the component of the";
                errMsg += " path prefix of dirname or read permission is ";
                errMsg += "denied for dirname.";
                break;
#ifndef WIN32
      case(ELOOP):
                errMsg = "Too many symbolic links were encountered in ";
                errMsg += "resolving path.";
                break;
#endif /* not WIN32 */
      case(ENAMETOOLONG):
                errMsg = "The length of the dirname argument exceeds ";
                errMsg += "{PATH_MAX}, or a pathname component is longer ";
                errMsg += "than {NAME_MAX}.";
                break;
      case(ENOENT):
                errMsg = "A component of dirname does not name an existing ";
                errMsg += "directory or dirname is an empty string.";
                break;
      case(ENOTDIR):
                errMsg = "A component of dirname is not a directory.";
                break;
      default:
                errMsg = "An unknown error occurred.";
                break;
   }
   return errMsg;
}

SxString SxFSError::getRmdirErrMsg ()
{
   SxString errMsg;
   switch (errno)  {
      case(EACCES):
                errMsg = "Search permission is denied on a component of the";
                errMsg += " path prefix, or write permission is denied on ";
                errMsg += "the parent directory of the directory to be ";
                errMsg += "removed.";
                break;
      case(EBUSY):
                errMsg = "The directory to be removed is currently in use by ";
                errMsg += "the system or another process and the ";
                errMsg += "implementation considers this to be an error.";
                break;
      case(EEXIST):
      case(ENOTEMPTY):
                errMsg = "The path argument names a directory that is not an ";
                errMsg += "empty directory.";
                break;
      case(EIO):
                errMsg = "A physical I/O error has occurred.";
                break;
#ifndef WIN32
      case(ELOOP):
                errMsg = "Too many symbolic links were encountered in ";
                errMsg += "resolving path.";
                break;
#endif /* not WIN32 */
      case(ENAMETOOLONG):
                errMsg = "The length of the path argument exceeds {PATH_MAX} ";
                errMsg += "or a pathname component is longer than {NAME_MAX}.";
                break;
      case(ENOENT):
                errMsg = "A component of path does not name an existing file, ";
                errMsg += "or the path argument names a non-existent directory";
                errMsg += " or points to an empty string.";
                break;
      case(ENOTDIR):
                errMsg = "A component of the path is not a directory.";
                break;
      case(EPERM):
                errMsg = "The S_ISVTX flag is set on the parent directory of ";
                errMsg += "the directory to be removed and the caller is not ";
                errMsg += "the owner of the directory to be removed, nor is ";
                errMsg += "the caller the owner of the parent directory, nor ";
                errMsg += "does the caller have the appropriate privileges.";
                break;
      case(EROFS):
                errMsg = "The directory entry to be removed resides on a ";
                errMsg += "read-only file system.";
                break;

      default:
                errMsg = "An unknown error occurred.";
                break;
   }
   return errMsg;
}

SxString SxFSError::getReadlinkErrMsg ()
{
   SxString errMsg;
   switch (errno)  {
      case(EACCES):
                errMsg = "Search permission is denied for a component of the";
                errMsg += "path prefix of path.";
                break;
      case(EINVAL):
                errMsg = "The path argument names a file that is not a ";
                errMsg += "symbolic link.";
                break;
      case(EIO):
                errMsg = "An I/O error occurred while reading from the file ";
                errMsg += "system.";
                break;
      case(ENOENT):
                errMsg = "A component of path does not name an existing file ";
                errMsg += "or path is an empty string.";
                break;
#ifndef WIN32
      case(ELOOP):
                errMsg = "Too many symbolic links were encountered in ";
                errMsg += "resolving path.";
                break;
#endif /* not WIN32 */
      case(ENAMETOOLONG):
                errMsg = "The length of path exceeds {PATH_MAX}, or a ";
                errMsg += "pathname component is longer than {NAME_MAX}.";
                break;
      case(ENOTDIR):
                errMsg = "A component of the path prefix is not a directory.";
                break;
      default:
                errMsg = "An unknown error occurred.";
                break;

   }
   return errMsg;
}

SxString SxFSError::getUnlinkErrMsg ()
{
   SxString errMsg;
   switch (errno)  {
      case(EACCES):
                errMsg = "Search permission is denied for a component of the ";
                errMsg += "path prefix, or write permission is denied on the ";
                errMsg += "directory containing the directory entry to be ";
                errMsg += "removed.";
                break;
      case(EBUSY):
                errMsg = "The file named by the path argument cannot be ";
                errMsg += "unlinked because it is being used by the system or ";
                errMsg +=  "another process and the implementation considers ";
                errMsg += "this an error.";
                break;
#ifndef WIN32
      case(ELOOP):
                errMsg = "Too many symbolic links were encountered in ";
                errMsg += "resolving path.";
                break;
#endif /* not WIN32 */
      case(ENAMETOOLONG):
                errMsg = "The length of the path argument exceeds {PATH_MAX} ";
                errMsg += "or a pathname component is longer than {NAME_MAX}.";
                break;
      case(ENOENT):
                errMsg = "A component of path does not name an existing file ";
                errMsg += "or path is an empty string.";
                break;
      case(ENOTDIR):
                errMsg = "A component of the path prefix is not a directory.";
                break;
      case(EPERM):
                errMsg = "The file named by path is a directory, and either ";
                errMsg += "the calling process does not have appropriate ";
                errMsg += "privileges, or the implementation prohibits using ";
                errMsg += "unlink() on directories.";
                break;
      case(EROFS):
                errMsg = "The directory entry to be unlinked is part of a ";
                errMsg += "read-only file system.";
                break;
      default:
                errMsg = "An unknown error occurred.";
                break;
   }
   return errMsg;
}

#ifdef WIN32
void SxFSError::throwShortcutException (int err,
                                        const SxString &lnkSrc,
                                        const SxString &lnkDest)
{
	SX_EXIT;
	/*
   switch(err)
   {
      case (Worked):
                std::cout << "Successfully created shortcut from '" << lnkSrc;
                std::cout << "' to '" << lnkDest << "'."<< std::endl;
                break;
      case (ShellLinkCreationFailed):
                SX_THROW (("Error: Can't create a shortcut from '"
                               + lnkSrc + "' to '" + lnkDest
                               + "'. No instance implementing the shell "
                               + "link interface could be fetched."
                              ).getElems (),
                              __FILE__, __LINE__);
                break;
      case (PersistFileCreationFailed):
                SX_THROW (("Error: Can't create a shortcut from '"
                               + lnkSrc + "' to '" + lnkDest
                               + "'. The shortcut couldn't be made "
                               + "persistent."
                              ).getElems (),
                              __FILE__, __LINE__);
                break;
      case (UnicodeConversionFailed):
                SX_THROW ((SxString ("Error: ")
                               + "Can't create a shortcut from '"
                               + lnkSrc + "' to '" + lnkDest
                               + "'. An error occurred during"
                               + " the conversion to an unicode string."
                              ).getElems (),
                              __FILE__, __LINE__);
                break;
      case (TargetDefinitionFailed):
                SX_THROW ("Error: Failed to set the target.",
                              __FILE__, __LINE__);
                break;
      default:
                SX_THROW ("Error: An unknown error occurred.",
                              __FILE__, __LINE__);
                break;
   }
   */
}
#endif /* WIN32 */
