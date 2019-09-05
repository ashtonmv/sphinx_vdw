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

#include <SxEvent.h>

SxEvent::SxEvent ()
   : data()
{
   SX_TRACE ();
}

SxEvent::~SxEvent ()
{
   SX_TRACE ();
}

SxEvent::SxEvent (const SxPtr<SxEventData> &data_)
   :  data(data_)
{
   SX_TRACE ();
}

SxEvent::SxEvent (const SxEvent &in)
{
   SX_TRACE ();

   data = in.data;
}

SxEvent &SxEvent::operator= (const SxEvent &in)
{
   if (&in == this) return *this;
   data = in.data;
   return *this;
}

SxEvent &SxEvent::operator>> (int64_t &i)
{
   SX_TRACE ();
   SxEventData *evd = data.getPtr ();
   SX_CHECK (evd != NULL);
   SX_CHECK (evd->iElem < evd->nElem, evd->iElem, evd->nElem);
   SX_CHECK (evd->layout[evd->iElem] == sx::Int, evd->layout[evd->iElem]);
   evd->read (&i);
   ++evd->iElem;
   return *this;
}

SxEvent &SxEvent::operator>> (SxString &str)
{
   SX_TRACE ();
   SxEventData *evd = data.getPtr ();
   SX_CHECK (evd != NULL);
   SX_CHECK (evd->iElem < evd->nElem, evd->iElem, evd->nElem);
   SX_CHECK (evd->layout[evd->iElem] == sx::String,
        (int)evd->layout[evd->iElem]);
   evd->read (&str);
   ++evd->iElem;
   return *this;
}

SxString SxEvent::shiftString ()
{
   SX_TRACE ();
   SxEventData *evd = data.getPtr ();
   SX_CHECK (evd != NULL);
   SX_CHECK (evd->iElem < evd->nElem, evd->iElem, evd->nElem);
   SX_CHECK (evd->layout[evd->iElem] == sx::String, (int)evd->layout[evd->iElem]);
   SxString str;
   evd->read (&str);
   ++evd->iElem;
   return str;
}

int64_t SxEvent::shiftInt ()
{
   SX_TRACE ();
   SxEventData *evd = data.getPtr ();
   SX_CHECK (evd != NULL);
   SX_CHECK (evd->iElem < evd->nElem, evd->iElem, evd->nElem);
   SX_CHECK (evd->layout[evd->iElem] == sx::Int,
        (int)evd->layout[evd->iElem]);
   int64_t i = 0;
   evd->read (&i);
   ++evd->iElem;
   return i;
}
