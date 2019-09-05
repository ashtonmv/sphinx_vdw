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

#ifndef _SX_STRUCT_OPT_H_
#define _SX_STRUCT_OPT_H_

#include <SxSymbolTable.h>
#include <SxArray.h>
#include <SxPotential.h>
#include <SxMatrix.h>
#include <SxTypes.h>
#include <SxSecondaryStructure.h>
#include <SxRandom.h>
#include <SxStruct.h>

/** \brief Structure optimization routines

    \b SxStructOpt = SPHInX Structure Optimization

    In this class the structure optimization, also known as structure 
    relaxation schemes are defined. It basically requires only 
    -# the current atomic structure
    -# a callback to evaluate the forces at a new set of atomic positions

    The callback to the force update routine has been realized by calling
    the virtual function SxPotential::getForces. This allows to use the
    structure optimization class with any kind of potential. To use a
    new potential simple derive from SxPotential and overload the function
    getForces appropriatly.

    \ingroup group_structure
    \author  Sixten Boeck, boeck@mpie.de */
class SX_EXPORT_STRUCT SxStructOpt
{
   public:

      //---------------------------------------------------------------------
      /**@name Constructors and Destructors */
      //---------------------------------------------------------------------
      //@{
      SxStructOpt ();
      SxStructOpt (SxPotential *, SxAtomicStructure &);
      ~SxStructOpt ();
      //@}

      //---------------------------------------------------------------------
      /**@name Interface to input file
         Controlling the class by the \ref tutor_parser. */
      //---------------------------------------------------------------------
      //@{
      void print (const SxSymbolTable *);
      void execute (const SxSymbolTable *, bool calc=true);
      //@}
      //---------------------------------------------------------------------
      static bool isRegistered (const SxSymbolTable *);


      //---------------------------------------------------------------------
      /**@name Structure optmimization algorithms
         Schemes to relax a given atomic structure. */
      //---------------------------------------------------------------------
      //@{
      
      void dampedNewton (const SxSymbolTable *, bool calc);

      /** \brief Quasi-Newton algorithm according to BFGS scheme \ref Schleyer

          The current quasi-Newton implementation is based on the 
          BFGS method \ref Schleyer. It reads as follows
          -# initialize the approximation \b B of the Hessian \b H 
             \f[
                \mathbf{B} \approx \mathbf{H} = -\mathbf{\hat{1}}
             \f] and it's inverse \f[
                \mathbf{\hat{B}} = \mathbf{H}^{-1}
             \f]
          -# calculate the gradient \b g (forces \b f) \f[
                \mathbf{g} \equiv \mathbf{f}
                           = \frac{\delta E}{\delta\mathbf{x} }
             \f]
          -# move the atoms like \f[
                \mathbf{x}' = \mathbf{x} - \mathbf{\hat{B}} \mathbf{g}
             \f]
          -# calculate the forces at new positions \b x' \f[
                \mathbf{g'} \equiv \mathbf{f'}
                            = \frac{\delta E}{\delta\mathbf{x'} }
             \f]
          -# compute the differences \f[
                \mathbf{s}  = \mathbf{x'} - \mathbf{x}
             \f] 
             \f[
                \mathbf{y}  = \mathbf{g'} - \mathbf{g}
             \f]
          -# update the approximation of the Hessian according to \f[
                \mathbf{B'} = \mathbf{B}
                            - \frac{\mathbf{B}   \mathbf{s}
                                    \mathbf{s}^T \mathbf{B}^T}
                                   {\mathbf{s}^T \mathbf{B} \mathbf{s}}
                            + \frac{\mathbf{y}   \mathbf{y}^T}
                                   {\mathbf{y}^T \mathbf{s}}
             \f]
          -# continue at point 2 unless convergence is reached

           \param cmd   the input file group QN {}
           \param calc  if true calculation will be performed otherwise only
                        parsed input parameters will be printed

           \author Sixten Boeck, boeck@fhi-berlin.mpg.de
       */
      void quasiNewton (const SxSymbolTable *cmd, bool calc);
      /** \brief Linear-scaling quasi-Newton algorithm according to 
                 BFGS scheme \ref Schleyer

          The current quasi-Newton implementation is based on the 
          BFGS method \ref Schleyer. It reads as follows
          -# initialize the approximation \f$\hat B\f$ of the inverse Hessian 
             \f$\b H^{-1}\f$
             \f[
                \mathbf{H^{-1}} \approx \mathbf{\hat B} 
                = \lambda \mathbf{\hat{1}}
             \f]
          -# calculate the gradient \b g (forces \b f) \f[
                \mathbf{g} \equiv \mathbf{f}
                           = \frac{\delta E}{\delta\mathbf{x} }
             \f]
          -# move the atoms like \f[
                \mathbf{x}' = \mathbf{x} - \mathbf{\hat{B}} \mathbf{g}
             \f]
          -# calculate the forces at new positions \b x' \f[
                \mathbf{g'} \equiv \mathbf{f'}
                            = \frac{\delta E}{\delta\mathbf{x'} }
             \f]
          -# compute the differences \f[
                \mathbf{s}  = \mathbf{x'} - \mathbf{x}
             \f] 
             \f[
                \mathbf{y}  = \mathbf{g'} - \mathbf{g}
             \f]
          -# update the approximation of the inverse Hessian according to 
          \f[
                \mathbf{B'} = 
                (1 - \frac{\mathbf s \mathbf y^T}{\mathbf y^T \mathbf s})
                \mathbf{B}
                (1 - \frac{\mathbf y \mathbf s^T}{\mathbf y^T \mathbf s})
                + \frac{\mathbf s \mathbf s^T}{\mathbf y^T \mathbf s}
          \f]
          -# continue at point 2 unless convergence is reached

          By expanding the braces in the update step one can see that
          the matrix \f$\hat B\f$ is indeed a sum over projector pairs
          \b s and \b t = \f$\hat B \mathbf y\f$
          \f[
          \hat B = \lambda \hat 1 + 
          \sum_k (\mathbf s_k \mathbf t_k) 
                 \left(\begin{array}{cc} 
                           a_{ss} & a_{st} \\ 
                           a_{ts} & a_{tt} 
                           \end{array}\right)
                 (\mathbf s_k \mathbf t_k)^T
          \f]
          
          If we restrict the number of projector pairs, the scaling
          becomes linear in the number of degrees of freedom, rather
          than quadratic (for the matrix multiplication) to cubic (for the
          inversion) in the normal QN.
          Furthermore, old curvature information is gradually discarded.

          The implementation ensures that the approximate Hessian remains
          positive definite.

           \param cmd   the input file group linQN {}
           \param calc  if true calculation will be performed otherwise only
                        parsed input parameters will be printed

           \author Christoph Freysoldt freyso@fhi-berlin.mpg.de
       */
      void linearQuasiNewton (const SxSymbolTable *cmd, bool calc);
      /** \brief Another BFGS variant, but on top of an approximate
                 Hessian based on interatomic distances.
                 @seealso SxRedundantCoords
        */
      void ricQuasiNewton (const SxSymbolTable *cmd, bool calc);
      /** \brief External control interface
          This is the interface to external minimizers.

          SPHInX reads commands from a file specified via the
          SX_EXT_CTRL environment variable (which could be a named pipe).
          This is the control channel.
          It outputs results to the file specified via the SX_EXT_RES
          environment variable (which could be a named pipe, too).
          This is the "output channel".

          Currently, the following commands are available
          - 'end': end extControl gracefully
          - 'get energy': print the energy on output channel
          - 'get forces': print the forces (n x 3 numbers) on output channel
          - 'get natoms': print number of atoms on output channel
          - 'get structure': print atomic coordinates (n x 3 numbers)
                             on output channel
          - 'set structure': read new atomic coordinates (n x 3 numbers)
                             from control channel
          - 'shift atom <id>': read xyz (3 numbers) from control channel and
                               shift atom <id> (<id> starting at 1) by xyz
          - 'get nspinconstraints': print number of spin-constrained atoms
          - 'get nu': print atomic fields for the PAW spin contraints
          - 'set spinconstraint': read the atomic target spins from control channel (n numbers)
          - 'run [<id>]': run electronic loop.  Different electronic loops
            can be specified in the input file by giving an 'id' to
            each bornOppenheimer {} group.
          - 'onproblem=crash': crash when control input is invalid (core dump)
          - 'onproblem=stop': stop when control input is invalid (no core dump)
          - 'onproblem=ignore': ignore invalid command, if possible

        */
      void extControl (const SxSymbolTable *cmd, bool calc);
      //@}

   protected:
      /** \brief Generate force/step filters from input file
          -# SxDriftFilter
          -# SxStickyFilter
          -# SxLineFilter

          @return the force filter

        */
      SxOperator<SxAtomicStructure> getFilter (const SxSymbolTable *table);

      SxOperator<SxVector<Double> > stepFilter;

      //---------------------------------------------------------------------
      /**@name Service routines */
      //---------------------------------------------------------------------
      //@{
      void set (SxPotential *hamSolver, SxAtomicStructure &str);
      //@}

   protected:

      int maxSteps;

      SxAtomicStructure    structure;
      SxSpeciesData        species;
      SxPotential         *potentialPtr;
      int nDoF;
      SxArray<const SxSymbolTable *> elMinimCmds;

      SxAtomicStructure getForces ();
      SxAtomicStructure getForces (const SxAtomicStructure &);
      SxAtomicStructure getForces (const SxMatrix<TPrecTauR> &);
      

};

#endif /* _SX_STRUCT_OPT_H_ */
