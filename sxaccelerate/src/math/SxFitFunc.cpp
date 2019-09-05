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

#include <SxFitFunc.h>

SxVector<Double> SxFitFunc::getParamGrad (const SxArray<SxVector<Double> > &xList,
                                     const SxArray<double> &yList)
{
   SxVector<Double> g(getNParam ());
   g.set (0.);
   for (int i = 0; i < xList.getSize (); ++i)  {
      double val = 0.;
      SxVector<Double> dfdParam = getParamDeriv(xList(i), &val);
      g.plus_assign_ax (val - yList(i), dfdParam);
   }
   return g;
}

void SxFitFunc::fitFromValsCG (const SxArray<SxVector<Double> > &xList,
                               const SxArray<double> &yList)
{
   SxVector<Double> X, g, Xold, gOld;

   SX_CHECK (xList.getSize () == yList.getSize (),
             xList.getSize (), yList.getSize ());

   double xTrial = paramFitStep ();
   for (int it = 0; it < 1000; ++it)  {
      // get gradient
      g = getParamGrad (xList, yList);

      // conjugate-gradient search direction
      if (it == 0)  {
         X = g.getCopy ();
      } else {
         double gamma = g.normSqr () / gOld.normSqr ();
         X = g + gamma * Xold;
      }
      // --- trial step
      param.plus_assign_ax(-xTrial, X);
      paramChange ();
      SxVector<Double> gTrial = getParamGrad (xList, yList);

      // --- optimal step
      SxVector<Double> dg = (gTrial - g) / xTrial;
      double xOpt = dot (dg, g) / dg.normSqr ();
      param.plus_assign_ax (xTrial - xOpt, X);
      paramChange ();

      // converged ?
      if (fitOK(xOpt * X)) break;

      // --- update xTrial
      if (1.5 * xOpt > xTrial)  {
         if (0.75 * xOpt > xTrial) xTrial *= 2.;
         else xTrial = sqrt(xTrial * xOpt/1.5);
      } else if (xOpt > 0.)  {
         if (3 * xOpt < xTrial) xTrial /= 2.;
         else xTrial = sqrt(xTrial * xOpt/1.5);
      }

      // --- keep gradient and search direction as "old values"
      Xold = X;
      gOld = g;
   }
}

SxFitPolynomial::SxFitPolynomial (int nVar_, int order_)
   : SxFitFunc (nVar_), order(order_), relAcc(1e-6)
{
   if (nVar == 0) return;
   int nParam = 1;
   // nParam is (nVar + order) over (order)
   // reason: for each term in the polynomial sum
   //         make a comma-separated list of powers,
   //         and write an X for each power of a variable
   // Example: x^3 y^0 z^2 -> XXX,,XX
   // The unused powers are written in the front.
   // For instance, if order=7
   // x^0 y^0 z^3 -> XXXX,,,XXX
   // x^1 y^2 z^4 -> ,X,XX,XXXX
   // -> number of commas = number of variables
   // -> number of Xs     = order of polynomial
   // Each term is uniquely defined by the order of X's and ,
   // => distribute nVar commas onto nVar + order places
   // => (nVar + order) over (order)
   for (int i = 0; i < min(nVar,order); i++) {
      nParam *= (nVar + order - i);
      nParam /= (i+1); // works without rest
   }
   param.resize (nParam);
}

double
SxFitPolynomial::evaluate (const SxVector<Double> &var) const
{
   double res;
   evaluate (var, &res, None);
   return res;
}

SxVector<Double>
SxFitPolynomial::getParamDeriv (const SxVector<Double> &var,
                                double *val) const
{
   return evaluate (var, val, ParamGrad);
}

SxVector<Double>
SxFitPolynomial::getDeriv (const SxVector<Double> &var,
                           double *val) const
{
   return evaluate (var, val, VarGrad);
}

SxVector<Double>
SxFitPolynomial::evaluate (const SxVector<Double> &var,
                           double *val,
                           GradType gradientType) const
{
   SxVector<Double> res;

   if (gradientType == ParamGrad) {
      res.resize (getNParam ());
      res(0) = 1.;
   }
   if (gradientType == VarGrad)  {
      res.resize (getNVar ());
      res.set (0.);
   }

   if (val) *val = param(0);

   // variable index. ii(j) >= ii(j+1)
   SxArray<int> ii(order + 1);
   // prod(j) = PI_{k=j+1..n} var(ii(j))
   SxArray<double> prod(order);
   for (int io = 1, ip = 1; io <= order; ++io)  {
      ii.set (0);
      prod(io-1) = 1.;
      for (int j = io-1; j > 0; j--)  {
         prod(j-1) = prod(j) * var(0);
      }
      while (ii(io) < nVar)  {
         //for (int t = 0; t < io; t++) cout << ii(t) << ' ';
         //cout << prod << endl;
         if (val) *val += param(ip) * var(ii(0)) * prod(0);
         if (gradientType == ParamGrad)  {
            res(ip) = var(ii(0)) * prod(0);
            SX_CHECK_NUM(res(ip));
         }
         if (gradientType == VarGrad)  {
            double prod2 =  param(ip); // param(ip) * Prod_{k=0..j-1} var(ii(k))
            for (int j = 0; j < io; ++j)  {
               res(ii(j)) += prod2 * prod(j);
               prod2 *= var(ii(j));
            }
         }
         ip++;
         int k = 0;
         while (k < io) {
            if (++ii(k) == nVar)
               k++;
            else
               break;
         }
         if (k == io) break;
         for (; k > 0; k--) {
            ii(k-1) = ii(k);
            prod(k-1) = prod(k) * var(ii(k));
         }
      }
   }
   /*
   if (gradientType == VarGrad)  {
      SxVector<Double> res2 = evaluate3(var, 1);
      for (int i = 0; i < nVar; i++)
         cout << res(i) << '=' << res2(i+1) << '+' << -(res2(i+1) - res(i)) << endl;
      if (val) cout << *val << '=' << res2(0) << '+' << (*val - res2(0)) << endl;
   }
   */
   //cout << "---" << endl;
   return res;
}

bool SxFitPolynomial::fitOK (const SxVector<Double> &dParam) const
{
   for (int i = 0; i < dParam.getSize (); ++i)  {
      if (fabs(dParam(i)) > relAcc * fabs(param(i))) return false;
   }
   return true;
}

void SxFitPolynomial::rescale (const SxVector<Double> &scale)
{
   // variable index. ii(j) >= ii(j+1)
   SxArray<int> ii(order + 1);
   for (int io = 1, ip = 1; io <= order; ++io)  {
      ii.set (0);
      while (ii(io) < nVar)  {
         for (int t = 0; t < io; t++)
            param(ip) *= scale(ii(t));
         ip++;
         int k = 0;
         while (k < io) {
            if (++ii(k) == nVar)
               k++;
            else
               break;
         }
         if (k == io) break;
         for (; k > 0; k--) {
            ii(k-1) = ii(k);
         }
      }
   }
}

double
SxFitPolynomial::evaluateHorner0 (const SxVector<Double> &var) const
{
   SxArray<int> paramIdx(order + 1), varId(order + 1);
   SxArray<double> partialSum (order + 1);

   paramIdx(0) = 0;
   for (int io = 0, sn = 1; io < order ; io++)
   {
      paramIdx(io+1) = sn;
      sn += sn * nVar / (io + 1);
   }

   int depth = 0;
   varId(0) = 0;
   while (varId(0) < nVar)  {
      // --- move down the tree to the leftmost remaining node
      for ( ; depth < order; ++depth)  {
         varId(depth + 1) = varId(depth);
         partialSum(depth) = 0.;
      }
      partialSum(depth) = 0.;

      // --- now move up again as many levels as necessary
      do {
         // add param of this depth
         int ip = paramIdx(depth)++;
         partialSum(depth) += param(ip);
         // --------- RETURN -------
         if (depth == 0) return partialSum(0);
         // ------------------------
         // contribute to level up by multiplying with var(varId)
         partialSum(depth - 1) += partialSum(depth) * var(varId(depth));

         // set varId for a potential move to the right
         varId(depth)++;
      } while (varId(depth--) == nVar); // move up

      // move down again
      depth++;
   }
   SX_EXIT;
}

/*
int getParamId (int nVar, SxStack<int> &varId)
{
   SX_CHECK (nVar > 0, nVar);
   if (nVar == 1) { int res = (int)varId.getSize (); varId.removeAll (); return res; }
   int res = 0;
   int binom = 1;
   int vOld = 2;
   int m = 0;
   while (!varId.isEmpty ())  {
      int v = nVar - varId.pop ();
      for (int i = vOld; i < v; ++i)
         (binom *= (m + i - 1) ) /= (i - 1);
      m++;
      if (v > 1)  {
         (binom *= (m + v - 2) ) /= m;
         res += binom;
         vOld = v;
      }
      cout << res << ' ' << v << ' ' << binom << ' ' << m << endl;
   }
   // now binom is (m + vOld - 2) over m
   (binom *= (m + vOld - 1) * (m + vOld)) /= ((vOld -1 ) * vOld);
   // now binom is (m + vOld) over m
   for (int i = vOld + 1; i <= nVar; i++)
      (binom *= (m + i)) /= i;
   // now binom is (m + nVar over m)
   return binom - res - 1;
}
*/

SxVector<Double>
SxFitPolynomial::evaluate3 (const SxVector<Double> &var, int levels) const
{
   SX_CHECK (nVar > 0, nVar);
   SxVector<Double> res = param.getCopy ();
   SxArray<int> varId(order);
   SxStack<int> upperId(order);

   for (int iLevel = 0; iLevel <= levels; ++iLevel)  {
      for (int depth = order - 1; depth >= iLevel; depth--)  {
         int ip = 1;
         for (int i = 1; i <= depth ; i++) (ip *= (i + nVar)) /= i;
         int jpEnd = ip - 1;
         for (int i = 0; i <= depth; i++) varId(i) = 0;
         for ( ; varId(0) < nVar; ip++)  {
            /*
            for (int i = 0; i < depth - iLevel; i++)
               upperId << varId(i);
            for (int i = depth - iLevel + 1; i <= depth; i++)
               upperId << varId(i);
            cout << varId << "->" << SxArray<int> (upperId) << endl;
            int jp = getParamId (nVar, upperId);
            */
            int jp = 0;

            // calculate node index of dropping the (depth-iLevel)-th entry in varId
            if (nVar == 1) {
               jp = depth;
            } else {
               int binom = 1;
               int vOld = 2;
               int m = 0;
               int jFromEnd = 0;
               for (int id = depth; id >= 0; id--)  {
                  if (id == depth - iLevel) continue;
                  int v = nVar - varId(id);
                  for (int i = vOld; i < v; ++i)
                     (binom *= (m + i - 1) ) /= (i - 1);
                  m++;
                  if (v > 1)  {
                     (binom *= (m + v - 2) ) /= m;
                     jFromEnd += binom;
                     vOld = v;
                  }
               }
               jp = jpEnd - jFromEnd;
            }

            // add partial result to jp-th node one level up
            res(jp) += res(ip) * var(varId(depth - iLevel));

            // update varId to reflect next parameter ip
            int k = depth;
            // increase last digits until nVar
            while (k >= 0 && ++varId(k) == nVar) k--;
            if (k < 0) break;
            // set all trailing digits to k-th one
            for ( ; k < depth ; k++) varId(k+1) = varId(k);
         }
      }
   }
   return res;
}
