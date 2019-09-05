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

#include <SxEventBus.h>
#include <SxSystemThread.h>
#include <SxCLI.h>

struct Publisher : public SxSystemThread
{
   Publisher (int nRuns_ = 10)
      : nRuns(nRuns_)
   {
      // empty
   }
   ~Publisher () { }

   virtual void main () override {
      for (int i = 0; i < nRuns; ++i) {
         SxTime::msleep (150);
         SX_EVENT_BUS_PUB ("StatusUpdate", i, nRuns);
      }
      SxTime::msleep (1'000);
      SxEventBus::terminate ();
   }

   int nRuns;
};

int main (int argc, char **argv)
{
   SxCLI cli (argc, argv);
   cli.finalize ();

   SxEventListener l = {
      { "StatusUpdate", "StatusUpdate"_SX }
   };
   l.setName ("ProgressListener");
   SX_EVENT_BUS_SUB ("StatusUpdate", sx::Int << sx::Int, l);

   Publisher publisher (100);
   publisher.start ();

   bool stop = false;
   while (1) {
      l.waitForEvents<"StatusUpdate"_SX> (&stop, 1, 100);

      if (stop) break;

      SxEvent e = l.getLast ();
      int cur = (int)e.shiftInt ();
      int max = (int)e.shiftInt ();
      sxprintf ("\rprogress: %3d%%", 100 * (cur+1)/max);
      fflush (stdout);
   }
   sxprintf ("\n");
   SxEventBus::unsubscribeAll (&l);

   publisher.wait ();

   return 0;
}
