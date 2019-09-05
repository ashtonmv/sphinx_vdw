// ---------------------------------------------------------------------------
//
//      The ab-initio based multiscale library
//
//                  S / P H I / n X
//
//      Copyright:  Max-Planck-Institute for Iron Research
//                  40237 Duesseldorf, Germany
//
//      Contact:    https://sxlib.mpie.de
//      Authors:    see sphinx/AUTHORS
//      License:    see sphinx/LICENSE
//
// ---------------------------------------------------------------------------

#ifndef _SX_TESTBED_H_
#define _SX_TESTBED_H_

#include <SxParser.h>
#include <SxTimer.h>
#include <SxExt.h>


class SX_EXPORT_EXT SxTest 
{
   public:
      SxTest (const SxSymbolTable *);
      virtual ~SxTest () { }

      void execute ();
      virtual int evaluate ()=0;

      SxString getTitle () const { return title; }
      SxString getFolder () const { return folder; }
      virtual SxString getType ()=0;
      
   protected:
      SxString title, folder, cmd, args, sxHome, cmdln;
      SxTimer  timer;
};

class SX_EXPORT_EXT SxDiffTest : public SxTest
{
   public:
      SxDiffTest (const SxSymbolTable *);
      virtual ~SxDiffTest () { }
      virtual int evaluate ();

      virtual SxString getType () { return "File Diff"; }

   protected:
      SxString file;
};


class SxXYTest : public SxTest
{
   public:
      SxXYTest (const SxSymbolTable *);
      virtual ~SxXYTest () { }
      virtual int evaluate ();
      virtual SxString getType () { return "XY Plot"; }

   protected:
      SxString xLabel, yLabel;
};

class SxNXYTest : public SxXYTest
{
   public:
      SxNXYTest (const SxSymbolTable *);
      virtual ~SxNXYTest () { }
      virtual int evaluate ();
      virtual SxString getType () { return "nXY Plot"; }

   protected:
      SxString xLabel, yLabel;
};



/** \brief ...

    \b SxClass = S/PHI/nX ...

    ....

    \author Sixten Boeck, boeck@mpie.de */
class SxTestbed
{
   public:

      SxTestbed ();
     ~SxTestbed ();


      void read (const SxParser::Table &);
      void print () const;
      void compute ();

   protected:

      SxList<SxPtr<SxTest> > tests;
      SxTimer timer;
            
};

#endif /* _SX_TESTBED_H_ */
