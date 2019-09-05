#include <SxEventBus.h>
#include <thread>
#include <SxCLI.h>

template<uint32_t event>
SxEvent getEvent1 ()
{
   SxPtr<SxEventData> data = SxPtr<SxEventData>::create ("1");
   data->setEventHash (event);
   SxEvent ev(data);
   return ev;
}

template<uint32_t event>
SxEvent getEvent2 ()
{
   SxString str = "1";
   int64_t integer = 2;
   SxPtr<SxEventData> data = SxPtr<SxEventData>::create (str, integer);
   data->setEventHash (event);
   SxEvent ev(data);
   return ev;
}

template<uint32_t event>
SxEvent getEvent3 ()
{
   SxString str1 = "1";
   int64_t integer = 2;
   SxString str2 = "3";
   SxPtr<SxEventData> data = SxPtr<SxEventData>::create (str1, integer, str2);
   data->setEventHash (event);
   SxEvent ev(data);
   return ev;
}

template<uint32_t event>
SxEvent getEvent4 ()
{
   SxString str1 = "1";
   int64_t int1 = 2;
   SxString str2 = "3";
   int64_t int2 = 4;
   SxPtr<SxEventData> data = SxPtr<SxEventData>::create (str1, int1, str2, int2);
   data->setEventHash (event);
   SxEvent ev(data);
   return ev;
}

template<uint32_t event>
SxEvent getEvent5 ()
{
   SxString str1 = "1";
   SxString str2 = "3";
   SxString str3 = "5";
   SxPtr<SxEventData> data
      = SxPtr<SxEventData>::create (str1, 2, str2, 4, str3);
   data->setEventHash (event);
   SxEvent ev(data);
   return ev;
}

template<uint32_t event>
SxEvent getEvent6 ()
{
   SxString str1 = "1";
   SxString str2 = "3";
   SxString str3 = "5";
   SxPtr<SxEventData> data
      = SxPtr<SxEventData>::create (str1, 2, str2, 4, str3, 6);
   data->setEventHash (event);
   SxEvent ev(data);
   return ev;
}

#define LAYOUT(a) \
   cout << "Layout: " << SxEventData::getLayout(SxList<uint8_t>() << a) << endl;

int main (int argc, char **argv)
{
   SxCLI cli (argc, argv);
   cli.finalize ();

   sxprintf ("Sizeof(SxEventData) = %ju\n", sizeof (SxEventData));
   sxprintf ("sizeof(SxEvent) = %ju\n", sizeof (SxEvent));
   sxprintf ("sizeof(SxPtr) = %ju\n", sizeof (SxPtr<int64_t>));

   cout << SX_SEPARATOR;
   SxEvent ev = getEvent1<"EV1"_SX> ();
   if (ev == "EV1"_SX) {
      cout << "Layout: " << ev.getLayout () << endl;
      LAYOUT(sx::String);
      cout << ev.shiftString () << endl;
   }

   cout << SX_SEPARATOR;
   ev = getEvent2<"EV2"_SX> ();
   if (ev == "EV2"_SX)  {
      cout << "Layout: " << ev.getLayout () << endl;
      LAYOUT(sx::String << sx::Int);
      cout << ev.shiftString () << endl;
      cout << ev.shiftInt () << endl;
   }

   cout << SX_SEPARATOR;
   ev = getEvent3<"EV3"_SX> ();
   if (ev == "EV3"_SX) {
      cout << "Layout: " << ev.getLayout () << endl;
      LAYOUT(sx::String << sx::Int << sx::String);
      cout << ev.shiftString () << endl;
      cout << ev.shiftInt () << endl;
      cout << ev.shiftString () << endl;
      //sxprintf ("%s %jd %s\n", ev.shiftString ().ascii (), ev.shiftInt (),
      //   ev.shiftString ().ascii ());
   }

   cout << SX_SEPARATOR;
   ev = getEvent4<"EV4"_SX> ();
   if (ev == "EV4"_SX) {
      cout << "Layout: " << ev.getLayout () << endl;
      LAYOUT (sx::String << sx::Int << sx::String << sx::Int);
      cout << ev.shiftString () << endl;
      cout << ev.shiftInt () << endl;
      cout << ev.shiftString () << endl;
      cout << ev.shiftInt () << endl;
   }

   cout << SX_SEPARATOR;
   ev = getEvent5<"EV5"_SX> ();
   if (ev == "EV5"_SX) {
      cout << "Layout: " << ev.getLayout () << endl;
      LAYOUT (sx::String << sx::Int << sx::String << sx::Int << sx::String);
      cout << ev.shiftString () << endl;
      cout << ev.shiftInt () << endl;
      cout << ev.shiftString () << endl;
      cout << ev.shiftInt () << endl;
      cout << ev.shiftString () << endl;
   }

   cout << SX_SEPARATOR;
   ev = getEvent6<"EV6"_SX> ();
   if (ev == "EV6"_SX) {
      cout << "Layout: " << ev.getLayout () << endl;
      LAYOUT (sx::String << sx::Int << sx::String << sx::Int << sx::String << sx::Int);
      cout << ev.shiftString () << endl;
      cout << ev.shiftInt () << endl;
      cout << ev.shiftString () << endl;
      cout << ev.shiftInt () << endl;
      cout << ev.shiftString () << endl;
      cout << ev.shiftInt () << endl;
   }

   cout << SX_SEPARATOR;
   SxEventListener l =
   {
      { "EventA", "EventA"_SX }
   };
   SX_EVENT_BUS_SUB ("EventA", sx::String << sx::Int << sx::String, l);

   std::thread publisher = std::thread (
   []() {
      SxString msg = "Hello?!";
      for (int i = 0; i < 100; ++i) {
         SX_EVENT_BUS_PUB("EventA", msg, i * 2, "Hello");
      }
   });

   bool stop = false;
   l.waitForEvents<"EventA"_SX> (&stop, 100, 10);
   if (l.getSize ()) {
      //SxList<SxEvent> events = l.getAll ();
      // method B
      size_t nEvents = l.getSize ();
      //for (size_t i = 0; i < nEvents; ++i) {
      /*LAST*/ {
         //SxEvent e = l.get ();
         SxEvent e = l.getLast ();
         cout << e.shiftString () << " "
              << e.shiftInt () << " "
              << e.shiftString ()
              << endl;
      }
      cout << "Received " << nEvents << " Events\n";
      cout << "Remaining " << l.getSize () << " Events\n";
   }

   publisher.join ();
   cout << SxEventBus::getGlobalObj ();
   SxEventBus::unsubscribeAll (&l);
   SxEventBus::terminate ();


   // sxdoc
   {
      SxEventListener::terminate(false);

      uint32_t eventHash = "MyEventName"_SX;
      sxprintf ("'MyEventName' has the hash: 0x%X\n", eventHash);

      SxEventListener l =
      {
         { "WindowOpened", "WindowOpened"_SX },
         { "WindowClosed", "WindowClosed"_SX }
      };
      l.setName ("WindowEventListener");

      SX_EVENT_BUS_SUB ("WindowOpened", sx::Int << sx::String, l);

      // --- publish
      SX_EVENT_BUS_PUB ("WindowOpened", 1430, "Marc");
      SX_EVENT_BUS_PUB ("WindowClosed", 1630, "Omar");

      cout << SxEventBus::getGlobalObj ();

      // --- waiting for events
      bool stop = false;
      size_t nEvents = 1;   // for how many events of that type do I want to wait
      size_t sleepMS = 100; // sleep interval between checks
      l.waitForEvents<"WindowOpened"_SX> (&stop, nEvents, sleepMS);

      l.waitForEvents_Or ({ "WindowOpened"_SX, "WindowClosed"_SX }, &stop);
      l.waitForEvents_And ({ "WindowOpened"_SX, "WindowClosed"_SX }, &stop);

      // --- obtaining event data
      SxList<SxEvent> events = l.getAll (); // the events in the listener are cleared
      sxprintf ("got %jd new events\n", events.getSize ());
      for (SxEvent &e : events)  {
         if (e == "WindowClosed"_SX)  {
            int64_t  time = e.shiftInt ();
            SxString user = e.shiftString ();
            sxprintf ("%s closed the window at %jd\n", user.getElems (), time);
         } else if (e.getHash () == "WindowOpened"_SX)  {
            int64_t  time = e.shiftInt ();
            SxString user = e.shiftString ();
            sxprintf ("%s closed the window at %jd\n", user.getElems (), time);
         }
      }

      SxEventBus::unsubscribeAll (&l);
   }
   return 0;
}
