#include <stdio.h>
#include <SxString.h>
#include <SxSigSlots.h>
#include <SxSlot.h>
#include <SxCLI.h>


class MySlider
{
   public:
      MySlider () { }
     ~MySlider () { }

   public slots:

   signals:

      SxSignal<double,const char*> SX_SIGNAL(valueChanged);
};

class MySurface
{
   public:
      MySurface ();
      ~MySurface ();
 
   public slots:

      SX_SLOT (slotGenerate, double /*threshold*/, const char *);

};


MySurface::MySurface ()
   : SX_LINK_SLOT(this, slotGenerate)
{
    // empty
}


void MySurface::slotGenerate (double threshold, const char *)
{
   printf ("MyIsosurface::slotGenerate received: %g\n", threshold);
}

MySurface::~MySurface ()
{
   std::cout << "MyIsosurface destructor\n";
}


int main (int argc, char **argv)
{
   SxCLI cli (argc, argv);
   cli.finalize ();

   SxPtr<MySlider> slider  = SxPtr<MySlider>::create ();
   SxPtr<MySurface> isoPtr = SxPtr<MySurface>::create ();

   sxconnect (slider->valueChanged, SX_GET_SLOT(isoPtr,slotGenerate));


//  sxdisconnect (slider->valueChanged, SX_GET_SLOT(isoPtr,slotGenerate));
// slider->valueChanged.disconnect (&SX_GET_SLOT(isoPtr,slotGenerate), SX_GET_SLOT(isoPtr,slotGenerate).cbCheck.getCB());

   slider->valueChanged.send (35.0);
   return 0;
}

