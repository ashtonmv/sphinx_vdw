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

#ifndef _SX_K_FIN_DIFFS_H_
#define _SX_K_FIN_DIFFS_H_


#include <SxDFT.h>
#include <SxCell.h>
#include <SxGkBasis.h>


/** \brief computes finite difference vectors in (periodic) k-space

    \b SxKFinDiffs = S/PHI/nX finite difference vectors in k-space

    Yields for each k-point index \f$ i_{\mathbf k} \f$
    - the difference vectors \f$ \mathbf{b}^n \;(n\in\{1, \ldots, N\} \f$ with
      \f$ N \f$ denoting the number of difference vectors needed\f$)\f$
      connecting the k-point \f$ \mathbf{k} \f$ with these neighboured
      k-points, that are required in the 1st order finite difference formulas
      for Nabla and Laplace (see Appendix B in
      <a href="http://prola.aps.org/abstract/PRB/v56/i20/p12847_1">
      Phys. Rev. B, 56, p. 12 847 (1997) </a>),
    - the required k-points in the neighbourhood themselves, i.e.,
      \f$ \mathbf{k} + \mathbf{b}^n \;(n\in\{1, \ldots, N\}) \f$,
    - the weights \f$ w_{\mathbf{b}^n} \f$ associated with each finite
      difference vector \f$ \mathbf{b}^n \f$.

    While the difference vectors \f$ \mathbf{b}^n \f$ and their associated
    weights \f$ w_{\mathbf{b}^n} \f$ remain the same for all k-points
    \f$ \mathbf{k} \f$ provided by the specified Monkhorst-Pack mesh, the
    set of neighboured k-points changes. A difficulty occurs, if
    \f$ \mathbf{k} \f$ is near to the edge of the Brioullin zone (or, to be
    more precise, the primitive zone of parallelepiped shape). In this case
    the finite difference formulas require k-points outside the primitive zone.
    By adding a suitable reciprocal lattice vector \f$ \mathbf{G} \f$, these
    points can be mapped to "equivalent" k-points
    \f$ \mathbf{q}'=\mathbf{q}+\mathbf{G} \f$ inside the primitive zone.

    Therefore, the neighboured k-points \f$ {\mathbf q} \f$ are provided by
    two parts of information:
    - the index \f$ i_{\mathbf{q}'} \f$ of an equivalent k-point inside the
      primitive zone,
    - the reciprocal lattice vector \f$ \mathbf{G} \f$ conveying
      \f$ {\mathbf q}' \f$ into \f$ {\mathbf q}=\mathbf{q}'-\mathbf{G} \f$.

    \author Matthias Wahn, wahn@fhi-berlin.mpg.de */

class SX_EXPORT_DFT SxKFinDiffs : public SxKPoints
{
   public:
      /** empty contructor */
      SxKFinDiffs ();

      /** constructor */
      SxKFinDiffs (const SxCell &cell_, const SxGkBasis &gkBasis);

      /** destructor */
      ~SxKFinDiffs ()  { /* empty */ }

      /** tries to get/guess the MP (Monkhorst-Pack) folding numbers, e.g,
          2x3x4, from a complete MP set of k-points

          @param  cell_  unit cell
          @param  gamma  if set to 'true', checks, that the\f$\Gamma\f$-point
                         is included

          \todo   may either go to SxKPoints or cell-argument can be removed
       */
      SxVector3<Int> getFoldingMP (bool gamma=false);

      /** gets the index of a k-point given in integer relative coordinates

          \note only works for MP meshes NOT being symmetry-reduced
          \todo may go to SxKPoints
       */
      int getKPointIdx (const RelVec &kRel);

      /** returns the k-point mesh in integer relative coordinates

          \todo may go to SxKPoints
       */
      void computeKVecsRel ();

      /** computes/updates surrounding k-points needed for the 1st order
          finite-difference formulas (see Appendix B in
          <a href="http://prola.aps.org/abstract/PRB/v56/i20/p12847_1">
          Phys. Rev. B, 56, p. 12 847 (1997) </a>)

          The information is then stored in two arrays:
          - SxKFinDiffs::kSurIdx stores the indeces of the surrounding
            k-points
            \f$ \mathbf{k} + \mathbf{b}^n - \mathbf{G}_n\;
            (n\in\{1, \ldots, N\}) \f$, with the \f$ \mathbf{G}_n \f$
            denoting possibly necessary reciprocal lattice vectors, that
            map \f$ \mathbf{k} + \mathbf{b}^n \f$ to equivalent points inside
            the primitive zone.
          - SxKFinDiffs::deltaGCart and SxKFinDiffs::deltaGRel store the
            vectors \f$ \mathbf{G}_n\; (n\in\{1, \ldots, N\}) \f$ described
            in the point above.
      
          \note The order of the indeces corresponds to the order of finite
                difference vectors (\f$\mathbf{b}\f$-vectors) of
                SxKFinDiffs::bVecs.

          \sa   SxKFinDiffs::kSurIdx
          \sa   SxKFinDiffs::deltaGCart
          \sa   SxKFinDiffs::deltaGRel
          \sa   SxKFinDiffs::bVecs
       */
      void computeSurroundingKPoints (int ik);

      /** Stores the indeces of the surrounding k-points
          \f$ \mathbf{k} + \mathbf{b}^n - \mathbf{G}_n\;
          (n\in\{1, \ldots, N\}) \f$, with the \f$ \mathbf{G}_n \f$
          denoting possibly necessary reciprocal lattice vectors, that
          map \f$ \mathbf{k} + \mathbf{b}^n \f$ to equivalent points inside
          the primitive zone.

          \sa SxKFinDiffs::computeSurroundingKPoints
          \sa SxKFinDiffs::deltaGCart
          \sa SxKFinDiffs::deltaGRel
          \sa SxKFinDiffs::bVecs
       */
      SxArray<int>                  kSurIdx;

      /** Stores the reciprocal lattice vectors \f$ \mathbf{G}_n\;
          (n\in\{1, \ldots, N\}) \f$ described under
          SxKFinDiffs::computeSurroundingKPoints and SxKFinDiffs::kSurIdx.
          The vectors are stored in Cartesian coordinates.

          \sa SxKFinDiffs::computeSurroundingKPoints
          \sa SxKFinDiffs::kSurIdx
          \sa SxKFinDiffs::deltaGRel
          \sa SxKFinDiffs::bVecs
       */
      SxArray<Coord>                deltaGCart;

      /** Stores the reciprocal lattice vectors \f$ \mathbf{G}_n\;
          (n\in\{1, \ldots, N\}) \f$ described under
          SxKFinDiffs::computeSurroundingKPoints and SxKFinDiffs::kSurIdx.
          The vectors are stored in relative coordinates.

          \sa SxKFinDiffs::computeSurroundingKPoints
          \sa SxKFinDiffs::kSurIdx
          \sa SxKFinDiffs::deltaGCart
          \sa SxKFinDiffs::bVecs
       */
      SxArray<RelVec>               deltaGRel;

      /** get the number \f$ N \f$ of finite difference vectors
          \f$ \mathbf{b}^n \;(n\in\{1, \ldots, N\} \f$ or, equally, the number
          of the required neighbour points in the \f$\mathbf{k}\f$-point mesh
       */
      int getNb () const;

      /** the array of finite difference vectors (\f$\mathbf{b}\f$-vectors)
          connecting each k-point with its surrounding k-points needed for the
          1st order finite differnce formulas of Nabla and Laplace
          (see Appendix B in
          <a href="http://prola.aps.org/abstract/PRB/v56/i20/p12847_1">
          Phys. Rev. B, 56, p. 12 847 (1997) </a>)
      
          \note The order of the \f$\mathbf{b}\f$-vectors corresponcds to the
                order of the \f$\mathbf{k}\f$-points referred to by the
                indeces of SxKFinDiffs::getSurroundingKPoints.

          \sa   SxKFinDiffs::computeSurroundingKPoints
          \sa   SxKFinDiffs::kSurIdx
          \sa   SxKFinDiffs::deltaG
       */
      SxArray<SxVector3<Double> >   bVecs;

      /** the array of finite difference vectors (\f$\mathbf{b}\f$-vectors,
          see SxKFinDiffs::bVecs) in integer relative coordinates */
      SxArray<RelVec>               bVecsRel;

      /** contains the corresponding weight \f$ w_{\mathbf{b}} \f$
          to each \f$ \mathbf{b} \f$ vector (see Appendix B of
          <a href="http://prola.aps.org/abstract/PRB/v56/i20/p12847_1">
          Phys. Rev. B, 56, p. 12 847 (1997) </a>)
       */
      SxArray<double>               bWeights;

      /** the k-point mesh (MP mesh) given in integer relative coordinates

          /todo may go to SxKPoints, if computeKVecsRel goes there as well
       */
      SxArray<RelVec>               kVecsRel;

   protected:
      /** unit cell */
      SxCell                   cell;

      /** reciprocal unit cell */
      SxCell                   recCell;

      /** type of the lattice structure */
      SxCell::CellType         lattice;

      /** the number \f$ N \f$ of finite difference vectors
          \f$ \mathbf{b}^n \;(n\in\{1, \ldots, N\} \f$ or, equally, the number
          of the required neighbour points in the \f$\mathbf{k}\f$-point mesh
       */
      int                      nbVecs;

      /** the MP folding as SxVector3<Double> */
      SxVector3<Double>        folding;

      /**@name initialisation for different lattice structures */
      //@{
      void initFCC ();
      void initWurtzite ();
      //@}

   public:
      /**@name functions for debugging */
      //@{
      /** \brief checks SxKFinDiffs::getKPointIdx */
      void checkKIdcs ();

      /** \brief checks the surrounding k-point vectors */
      void checkKSur (int ik);

      /** \brief checks the surrounding k-points vectors for all k-points */
      void checkKSur ();
      //@}
};

#endif /* _SX_K_FIN_DIFFS_H_ */
