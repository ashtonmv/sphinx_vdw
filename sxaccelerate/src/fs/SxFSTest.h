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
#ifndef _SX_FS_TEST_H_
#define _SX_FS_TEST_H_

// Including header-files
#include <SxFS.h>
#include <SxString.h>

/**
  \brief Base class for filesystem related unit tests.

  \author Thomas Uchdorf, t.uchdorf@mpie.de
 */
class SxTest;
class SX_EXPORT_FS SxTestUnitBase
{
   public:
      virtual ~SxTestUnitBase ()
      {
         // empty
      }
      virtual ssize_t getTestNo (SxTest *) const = 0;
};



/**
  \brief Class for single filesystem related tests that can belong to a unit
  of tests.

  \author Thomas Uchdorf, t.uchdorf@mpie.de
 */
class SX_EXPORT_FS SxTest  {
   public:
      // --- Constructors
      // Constructor
      SxTest (const SxString &descr = SxString (),
              SxTestUnitBase *par_ = NULL) : actualResult(false), correctResult(false),
                                             failed (false), par (par_)
   {
      description = (descr);
   }
      /*void postProcess ()
        {
        if (par)  {
        std::cout << "--------------------";
        std::cout << "--------------------";
        std::cout << "--------------------";
        std::cout << "--------------------" << std::endl;
        std::cout << "Peforming test " << ((par)?par->getTestNo (this) : 0);
        std::cout << "..." << std::endl;
        }
        }*/
      void setPrecondition (bool expr)
      {
         if (!(expr))  { failed = true; }
      }

      void setPostcondition (bool expr)
      {
         if (!(expr))  { failed = true; }
      }
      void setResult (bool expr)
      {
         correctResult = expr;
      }
      void setFailed ()
      {
         failed = true;
      }
      SxString getDescription () const
      {
         return description;
      }
      bool getActualResult () const
      {
         return actualResult;
      }
      bool getCorrectResult () const
      {
         return correctResult;
      }
      void evaluate ()
      {
         actualResult = !failed;
         if (actualResult == correctResult)  {
            std::cout << "OK" << std::endl;
         } else  {
            std::cout << "FAILED" << std::endl;
         }
         std::cout << description << std::endl;
      }

   protected:

      // --- Members
      bool actualResult;
      bool correctResult;
      bool failed;
      SxString description;
      SxTestUnitBase *par;
};



/**
  \brief Class for units of filesystem related tests.

  \author Thomas Uchdorf, t.uchdorf@mpie.de
 */
class SX_EXPORT_FS SxTestUnit : public SxTestUnitBase  {
   public:
      // --- Constructors
      // Constructor
      SxTestUnit () : SxTestUnitBase (), failed (false), testNo (0)
   {
      // empty
   }
      virtual ~SxTestUnit ()
      {
         // empty
      }

      // --- Methods
      void addTest (const SxString &description)
      {
         tests << SxTest (description, this);
         std::cout << "--------------------";
         std::cout << "--------------------";
         std::cout << "--------------------";
         std::cout << "--------------------" << std::endl;
         std::cout << "Peforming test " << tests.getSize ();
         std::cout << "..." << std::endl;
      }
      void setCurTestPrecondition (bool expr)
      {
         getCur ()->setPrecondition (expr);
      }
      void setCurTestPostcondition (bool expr)
      {
         getCur ()->setPostcondition (expr);
      }
      void setCurTestResult (bool expr)
      {
         getCur ()->setResult (expr);
      }
      void setCurTestFailed ()
      {
         getCur ()->setFailed ();
      }
      void evaluateCurTest ()
      {
         getCur ()->evaluate ();
      }
      SxString getCurTestDescription () const
      {
         return getCur ()->getDescription ();
      }
      virtual ssize_t getTestNo (SxTest *test) const
      {
         SxList<SxTest>::ConstIterator it;
         int i;
         for (it = tests.begin (), i = 0; it != tests.end (); ++it, ++i)  {
            if (&(*it) == test)  {
               return i;
            }
         }
         return -1;//tests.findPos (*test);
      }
      ssize_t getSize () const
      {
         return tests.getSize ();
      }

      void printResults (int offset)
      {
         SxList<SxTest>::Iterator itTest;
         int i;
         for (itTest = tests.begin (), i = offset;
              itTest != tests.end ();
              ++itTest, ++i)  {
            std::cout << "Test " << i << "(" << (*itTest).getDescription ();
            std::cout << "): ";
            if ((*itTest).getCorrectResult () !=
                (*itTest).getActualResult ())
            {
               std::cout << "FAILED" << std::endl;
            } else  {
               std::cout << "OK" << std::endl;
            }
         }
      }

   protected:

      // --- Methods
      SxTest * getCur ()
      {
         return &tests(tests.getSize () - 1);
      }

      SxTest const * getCur () const
      {
         return &tests(tests.getSize () - 1);
      }

      // --- Members
      bool failed;
      size_t testNo;
      SxList<SxTest> tests;
};



/**
  \brief Manages units of filesystem related tests.

  \author Thomas Uchdorf, t.uchdorf@mpie.de
 */
class SX_EXPORT_FS SxTestSuite {

   public:
      // --- Constructors
      // Constructor
      SxTestSuite ()
      {
         // empty
      }
      // Destructor
      ~SxTestSuite ()
      {
         // empty
      }
      // --- Methods
      void addTestUnit ()
      {
         testUnits << SxTestUnit ();
      }
      void addTestToCurTestUnit (const SxString &description)
      {
         getCur ()->addTest (description);
      }
      void setCurTestPrecondition (bool expr)
      {
         getCur ()->setCurTestPrecondition (expr);
      }
      void setCurTestPostcondition (bool expr)
      {
         getCur ()->setCurTestPostcondition (expr);
      }
      void setCurTestResult (bool expr)
      {
         getCur ()->setCurTestResult (expr);
      }
      void setCurTestFailed ()
      {
         getCur ()->setCurTestFailed ();
      }
      void evaluateCurTest ()
      {
         getCur ()->evaluateCurTest ();
      }
      void printResults ()
      {
         SxList<SxTestUnit>::Iterator itUnit;
         int offset;
         for (itUnit = testUnits.begin (), offset = 0;
              itUnit != testUnits.end ();
              ++itUnit)  {
            (*itUnit).printResults (offset);
            offset += int((*itUnit).getSize ());
         }
      }

   protected:
      // --- Methods
      SxTestUnit *getCur ()
      {
         return &testUnits(testUnits.getSize () - 1);
      }

      // --- Members
      SxList<SxTestUnit> testUnits;

};
SX_EXPORT_FS SxTestSuite testSuite, *globTestSuite = &testSuite;

#define SX_CHECK_PRE(expr) { globTestSuite->setCurTestPrecondition ((expr)); }
#define SX_CHECK_POST(expr) { globTestSuite->setCurTestPostcondition ((expr)); }
#define SX_INIT_TESTS() globTestSuite->addTestUnit ();
#define SX_CREATE_TEST(descr) { globTestSuite->addTestToCurTestUnit ((descr)); }
#define SX_MAKE_TEST(expr,correct) {\
   globTestSuite->setCurTestResult ((correct));\
   try  {\
      (expr);\
   } catch (SxException ex)  {\
      globTestSuite->setCurTestFailed ();\
      ex.print ();\
   }\
}
#define SX_QUIT_TEST() { globTestSuite->evaluateCurTest (); }

#endif /* _SX_FS_TEST_H_ */
