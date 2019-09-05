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

#include <SxLinEquation.h>

PrecTauR SxLinEquation::epsZero = 1e-7;
      
SxLinEquation::SxLinEquation (const SxMatrix3<TPrecTauR> &mat, const Coord &res)
   : soluble (true)
{
   for (int i = 0; i < 3; i++)
      addEquation (mat.row (i), res(i));
}

void SxLinEquation::addEquation (const Coord &coeff, PrecTauR res)
{
   if (!soluble) return;
   PrecTauR coeffSqr = coeff.absSqr().sum();
   if (coeffSqr < epsZero)  {
     if (fabs(res) >= epsZero) soluble = false;
     return;
   }
   if (leftSide.getSize () == 0)  {
      PrecTauR norm = sqrt (coeffSqr);
      leftSide.append (coeff / norm);
      rightSide.append (res / norm);
      return;
   }
   int nEq = getSize ();
   PrecTauR ovlp;
   Coord newEq = coeff;
   PrecTauR newRes = res;
   for (int i = 0; i < nEq; i++)  {
      ovlp = newEq ^ leftSide(i);
      newEq -= ovlp * leftSide(i);
      newRes -= ovlp * rightSide(i);
      if ((coeffSqr = newEq.absSqr().sum()) < epsZero)  {
         // equation coefficients depend linearly on others
         if (fabs(newRes) >= epsZero) soluble = false;
         return;
      }
   }
   PrecTauR norm = sqrt (coeffSqr);
   leftSide.append (newEq / norm);
   rightSide.append (newRes / norm);
}
     
void SxLinEquation::operator&= (const SxLinEquation &linEq)
{
   int nEq = linEq.getSize (); 
   for (int i = 0; i < nEq; i++)
      addEquation (linEq.leftSide(i), linEq.rightSide(i));
}

SxLinEquation SxLinEquation::operator&& (const SxLinEquation &eq2) const
{
   SxLinEquation result = *this;
   result &= eq2;
   return result;
}

SxLinEquation& SxLinEquation::operator= (const SxLinEquation &in)
{
   soluble = in.soluble;
   if (!soluble) return *this; // who cares still about the equations
   leftSide = in.leftSide;
   rightSide = in.rightSide;
   return *this;
}

bool SxLinEquation::operator== (const SxLinEquation &eq2)
{
   SX_CHECK (soluble);
   SX_CHECK (eq2.soluble);
   int nEq = this->getSize ();
   if (nEq != eq2.getSize ()) return false;
   SxLinEquation eq12 = (*this && eq2);
   return (eq12.soluble && nEq == eq12.getSize ());
}

SxString toString(const Coord &c)
{
   SxString res = "[";
   res += c(0); res += ",";
   res += c(1); res += ",";
   res += c(2); res += "]";
   return res;
}

Coord SxLinEquation::getPoint () const
{
   if (!soluble)  {
     cout << "Invalid access of solution space offset point: equation "
             "insoluble" << endl;
     SX_EXIT;
   }
   switch (getSize ())  {
      case 0:
         cout << "Invalid access of solution space offset point: there" << endl;
         cout << "is no non-trivial equation to be solved." << endl;
         SX_EXIT;
      case 1:
         // base point vector is parallel to normal vector
         return leftSide (0) * rightSide (0);
      case 2:
         // get line vector
         {
            Coord vec = leftSide (0).x(leftSide(1));
            vec /= sqrt (vec.absSqr().sum());
            // base point lies in planes given by the two equations
            // and its vector is orthogonal to the line vector
            return (SxMatrix3<TPrecTauR> (leftSide (0), leftSide (1), vec)
                    .inverse () 
                    ^ Coord (rightSide (0), rightSide (1), 0.));
         }
      case 3:
         return (SxMatrix3<TPrecTauR> (leftSide(0),leftSide(1),leftSide(2))
                 .inverse ()
                 ^ Coord (rightSide));
      default:
         cout << "Inconsistent SxLinEquation object: too many equations\n";
   }
   // Should never get here.
   SX_EXIT; 
   return Coord(0.,0.,0.);
}

   
SxString SxLinEquation::getName (const SxString &solName) const
{
   if (!soluble) return "no " + solName;
   SxString name;
   Coord vec;
   switch (leftSide.getSize ())  {
      case 0:
         name = "all " + solName + "s";
         break;
      case 1:
         name = "all " + solName + "s on plane ";
         name += toString(leftSide (0));
         name += " ^ x = ";
         name += rightSide (0);
         break;
      case 2:
         name = "all " + solName + "s on line a * ";
         vec = leftSide (0).x(leftSide(1));
         vec /= sqrt(vec.absSqr().sum());
         name += toString (vec) + " + " + toString (getPoint ());
         break;
      case 3:
         name = solName + " ";
         name += toString(getPoint ());
         break;
      default:
         /* empty */
         break;
   }
   return name;
}
