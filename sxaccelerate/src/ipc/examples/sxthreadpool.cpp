#include <SxCLI.h>
#include <SxThreadPool.h>
#include <SxThread.h>

class MyThread : public SxThread
{
   public:
      SxString msg;
      MyThread (const SxString &msg);
      virtual ~MyThread ();
      virtual void main ();
};

MyThread::MyThread (const SxString &msg_)
   : msg(msg_)
{
   // empty
}

MyThread::~MyThread ()
{
   // empty
}

void MyThread::main ()
{
   cout << msg << endl;
   SxThread::sleep (5);
}

int main (int argc, char **argv)
{
   SxCLI cli (argc, argv);
   int nThreads   = cli.option ("--nThreads",
                    "The number of system threads in threadpool."
                    ).toInt (1);
   cli.finalize ();

   // --- create thread pool
   SxThreadPool pool(nThreads);
   
   // --- SxThread will run automatically as a task in the given thread pool
   MyThread a("a");
   a.setThreadPool (pool.getThis ());
   a.start ();
   a.wait ();
   
   // --- submit a task to the thread pool manually
   SxPtr<MyThread> b = SxPtr<MyThread>::create("b");
   pool.submit (b);
   pool.wait (b->getThis ());

   return 0;
}
