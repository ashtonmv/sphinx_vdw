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

#ifndef _SX_SYSTEM_THREAD_H_
#define _SX_SYSTEM_THREAD_H_

#include <SxConfig.h>
#include <SxPtr.h>
#include <SxFunction.h>
#include <SxException.h>
#include <SxIPC.h>
#include <SxMutex.h>

#ifdef WIN32
#  include <windows.h>
   typedef unsigned int SX_THREAD_ID;
#else
#  include <pthread.h>
   typedef pthread_t SX_THREAD_ID;
#endif


class SxThreadId
{
   public:

      SX_THREAD_ID id;

      SxThreadId () : id (0) { }

      SxThreadId (SX_THREAD_ID id_) : id (id_) { }
      virtual ~SxThreadId () { }

      bool operator== (const SxThreadId &in) const
      {
#         ifdef WIN32
             return id == in.id;
#         else
             return pthread_equal (id, in.id);
#         endif
      }
};


/** \brief Interface for platform-independent operating system threads

    A SxSystemThread is a seperate thread of the program which can run in
    parallel with other threads. Every thread has access to all data of the
    parent process (shared memory). Threading works only if the underlying
    operating system supports multitasking.

    SxSystemThread abstracts whatever threading mechanism the underlying
    operating system provides.
    Its a use-as-best-as-possible approach with essential functionality of
    /threading/ supported, but no attempt to expose all features of each
    implementation.

    \author  Sixten Boeck, boeck@sfhingx.de */
class SX_EXPORT_IPC SxSystemThread
{
   public:

      SxSystemThread ();
      template<class Func>
      SxSystemThread (Func lambda)
         : running (false),
           hasLastException (false),
           hasThreadId (false)
      {
         // SX_TRACE ();
         mainLambda = SxFunction<void>::create (lambda);
      }


      /** \brief Priority of the thread

          This enum value indicates what priority is used when scheduling 
          the thread. The default is to use the same priority as the parent
          process. */
      enum Priority  {
         Idle,
         Lowest,
         Low,
         Normal,
         High,
         Highest,
         SameAsParent
      };


      /** \brief Destructor.

          The destrcutor does not terminate a thread! Waiting for the 
          end of a thread is done using the barrier function SxThread::wait. */
      virtual ~SxSystemThread ();

      /** \brief Start a new thread

          This function sets up the thread environment and execute the
          SxThread::main function of the derived class. 
          \param priority  The priority used for scheduling the new thread */
      virtual void start (Priority priority=SameAsParent);

      /** \brief Block a thread and wait its termination.

           This function creates a barrier and waits until the thread has 
           (been) terminated. 

           \par Example:
\code
   MyThread t1, t2;

   // --- multi-threadding part starts here
   t1.start (); t2.start ();
   t1.end();    t2.end();
   // --- end of multi-threadding part
\endcode 
       */
      virtual void wait ();

      /** \brief Get exception from the last run.

          SxSystemThread::wait() does not automatically rethrow
          the exception because the current code is not
          expecting such behavior.

\code
class MyThread : public SxSystemThread
{
   public:
      virtual void main () {
         SxString data = SxString::read ("nosuchfile");
         // ...
      }
};

int main ()
{
   MyThread t;

   t.start ();
   t.wait ();

   if (t.hasException ())  {
      t.getException().print (true);
      return 1;
   }

   return 0;
}

// Output
// Can't open file 'nosuchfile' for reading: No such file or directory
// This Exception was emitted from /.../src/util/SxString.cpp (line 1401)
\endcode 
      */
      bool hasException () const;

      SxException getException () const;

      /** \brief Terminates a thread

          This function terminates the thread by calling the corresponding
          kill function of the underlying thread library. */
      //void terminate ();

      /** \brief Slot for computation routine of a thread

          When the thread is executed from SxThread::start this virtual
          function is called. When defining a new thread type by deriving
          from SxThread, overload this function. It should contain the
          computational part which runs threadded. */
      virtual void main ();

      /** \brief Sleep for some secondus

        Calling this function makes the current process/thread sleep until
        a certain time has elapsed. This function can be applied to threads
        or to entire processes.

        \par Example 1
        Making a process sleep for 5 seconds
\code
   int main ()
   {
      for (int i=0; i < 100; ++i)  {
         printf ("i=%d\n", i);
         SxThread::sleep(5);
      }
   }
\endcode
 
       \par Example 2 
       Making a single thread sleep for 1 second for observing changes of 
       a variable (which might be controlled from the slave thread.
\code
   class MyMasterThread : public SxThread
   {
      public:
         MyMasterThread () : SxThread (), f(false)  {  }

         void setFlag (bool f)  {  flag = f; }

         virtual void main ()
         {
             for (;;)  {
                if (flag)  doSomeThing ();
                SxThread::sleep(1);
             }
         }
      protected:
         bool flag;
   };

   class MySlaveThread : public SxThread
   {
      public:
         MySlaveThread (MyMasterThread *m) : SxThread (), master(m);

         virtual void main ()
         {
            // do some heavy calculation and communicate with the master
            for (i=0; i<10000; ++i)  {
               // eat CPU time
               if (...)  master->setFlag (true);  // call master process
            }
         }

      protected:
         MyMasterThread *master;
   };

   int main ()
   {
      MyMasterThread tMaster;
      MySlaveThread  tSlave(&tMaster);
      tMaster.start ();
      tSlave.start ();
      tMaster.wait ();
      tSlave.wait ();
      return 0;
   }
\endcode */

      bool isRunning () const;

      //static SxThreadId getCurrentId ();

   protected:
      /**\brief user data for thread
        */
#     ifdef WIN32
         void **userData;
#     endif

      SxFunction<void> mainLambda;

      bool running;

      mutable SxMutex mutex;

      /**\brief Mutex for thread wait

          The wait routine cleans/resets member variables of SxSystemThread
          such as the thread ID or the userData for Windows. Since several
          other threads might wait for this thread to complete, this reset
          must be mutext within the wait routine.
      */
      SxMutex waitMutex;

      bool hasLastException;

      SxException lastException;

      /** \brief Launch a thread

           This hook-in function is called from the static thread callback
           and calls the SxThread::main routine of the derived class. */
      void launcher ();

   private:

      /**\brief thread handle

          The thread handle is the point to the windows thread instance*/
#     ifdef WIN32
         HANDLE       threadHandle;
#     endif

      /** \brief thread identifier

          This is the thread identifier used for internal thread library
          calls. */
      SxThreadId threadId;
      bool hasThreadId;

      /** \brief thread attributes

          The thread attributes specify mainly the scheduling policy of
          the thread. This variable is used for calling the thread library. */
#     ifndef WIN32
         pthread_attr_t  threadAttribs;
#     endif

      /** \brief static function callback

          This static callback function is the main entrance point of the
          thread creation process. */
#     ifdef WIN32
         static unsigned __stdcall loopCallback (void *);
#     else
         static void loopCallback (void *);
#     endif
};

#endif /* _SX_SYSTEM_THREAD_H_ */
