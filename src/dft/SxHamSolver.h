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
#ifndef _SX_HAM_SOLVER_H_
#define _SX_HAM_SOLVER_H_

#include <SxSymbolTable.h>
#include <SxPrecision.h>
#include <SxHamiltonian.h>
#include <SxPW.h>
#include <SxMuPW.h>
#include <SxRho.h>
#include <SxAtomicStructure.h>
#include <SxPotential.h>
#include <SxPseudoPot.h>
#include <SxFermi.h>
#include <SxConfig.h>
#include <SxDFT.h>
#include <SxSpinConstraint.h>

namespace Timer {
   enum HamSolverTimer {
      LcaoTotal,
      LcaoStep,
      ElMinim,
      SubspaceMatrix,SubspaceDiag,
      SearchDir, WaveUpdate,
      RmmDiis,
      parKLoop,   // khr
      khrDebug1,  // khr
      khrDebug2,  // khr
      khrDebug3,  // khr
      khrDebug4,   // khr
      benchmark
   };
}

SX_REGISTER_TIMERS (Timer::HamSolverTimer)
{
   using namespace Timer;
   regTimer (LcaoTotal,  "LCAO Initialization");
   regTimer (LcaoStep,   "LCAO Step");
   regTimer (ElMinim,    "Electronic Loop");
   regTimer (SearchDir,  "search direction");
   regTimer (WaveUpdate, "wave update");
   regTimer (SubspaceMatrix, "Subspace setup");
   regTimer (SubspaceDiag,   "Subspace diag.");
   regTimer (RmmDiis,        "RMM-DIIS solver");
   regTimer (parKLoop,   "inside par k loop");    // khr
   regTimer (khrDebug1,   "khr Debug marker 1");
   regTimer (khrDebug2,   "khr Debug marker 2");
   regTimer (khrDebug3,   "khr Debug marker 3");
   regTimer (khrDebug4,   "khr Debug marker 4");
   regTimer (benchmark,   "khr Benchmark");
}

/**  \brief Electronic minimizations to compute the Born-Oppenheimer surface

  \b SxHamSolver = S/PHI/nX Hamiltonian Solver

  The Hamiltonian solver contains various electronic minimization schemes.
  Most of them are designed to compute the Born-Oppenheimer surface. In
  addition iterative initialization schemes (e.g. tightBinding) are defined
  here. 

  \author  Sixten Boeck
  */
class SX_EXPORT_DFT SxHamSolver : public SxPotential
{
   public:
      double energy;
      SxPtr<SxSpinConstraint> spinConstraint;
      
      SxHamSolver ();
      SxHamSolver (const SxAtomicStructure &, const SxSymbolTable *);
      SxHamSolver (const SxString &wavesFile, const SxString &rhoFile, SxConstPtr<SxSymbolTable> tablePtr, const SxString &tmpDir = "", bool saveMemory = false);
      SxHamSolver (const SxString &rhoFile, SxConstPtr<SxSymbolTable> tablePtr);

      virtual ~SxHamSolver ();

      void init(const SxString &rhoFile, SxConstPtr<SxSymbolTable> tablePtr);
      void setupHam(const SxString &rhoFile, SxConstPtr<SxSymbolTable> tablePtr, SxPtr<SxGkBasis> gkPtr = SxPtr<SxGkBasis>());
      //---------------------------------------------------------------------
      /**@name File parser interface
         Controlling the class by the \ref tutor_parser. */
      //---------------------------------------------------------------------
      //@{
      virtual bool isRegistered (const SxSymbolTable *) const;
      virtual void execute (const SxSymbolTable *, bool calc=true);

      /** \brief Calculate forces for new structure using minimizer
        @param x       new structure (if uninitialized, used current structure)
        @param cmd     minimizer commands (if NULL, don't run minimizer)
        @return        forces (not symmetrized)
        */
      virtual SxAtomicStructure getForces (const SxAtomicStructure &x,
                                           const SxSymbolTable *cmd);
      /** \brief Try to move waves to new atomic structure
        */
      void moveWaves (const SxAtomicStructure &newStr);

      /** \brief Calculate forces for current structure
      */
      inline SxAtomicStructure getForces ()
      {
         return getForces (SxAtomicStructure (), NULL);
      }

      virtual SxSpeciesData getSpeciesData () const;
      virtual PrecEnergy getEnergy () const;
      //@}


      void initialGuess (const SxSymbolTable *, bool calc=true);

      /** \todo treat hamPW, and hamTB has SxHamiltonian, waves as SxPsiSet
          \todo H should be SxSymMatrix */
      void tightBinding (const SxSymbolTable *, bool calc=true);

      /** \brief Tight-binding initialization for PAW
       */
      void tightBindingPAW (const SxSymbolTable *, bool calc=true);

      /** \brief Reallocate memory for new number of states

          With this function reallocates the entities 
          -# SxHamSolver::waves
          -# SxHamSolver::fermi
          consitently. It is mainly used by the direct diagonalization.
          \param nPerK  number of states per k-point 
          \sa    SxPW::setNStates
          \sa    SxFermi::setNStates */
      void setNStates (const SxArray<int> &nPerK);

      /** \brief reallocate memory for new number of states

          This function changes the number of states consistently in
          -# SxHamSolver::waves,
          -# SxHamSolver::fermi,
          When the argument is set to -1, for every k-point the maximal number
          of states is set up. */
      void setNStates (int nStatesNew);


      //---------------------------------------------------------------------
      /**@name Born-Oppenheimer surface 
         These routines compute the Born-Oppenheimer surface. The results of
         the minimizations are the self-consistent charge density, the 
         wavefunctions, and the forces. */
      //---------------------------------------------------------------------
      //@{
      /** \brief Simple downhill minimization

          The Steepest-descent or downhill scheme is the most simple scheme
          implemented here. In order to improve a wavefunction it "walk"
          along the (negative) gradient direction for a given step width
          (timestep \f$ \Delta t \f$ deltaT) until convergence is reached.

          The algorithms is as follows
          -# Consider a single state. To simplify we will write
             \f[
                | \Psi_{i \sigma \mathbf{k}} \rangle 
                   \rightarrow 
                | \Psi_{i} \rangle
             \f]
          -# compute gradient and one-particle energy
             \f[
                \varepsilon_i = \langle \Psi_i | H | \Psi_i \rangle
             \f]
             \f[ 
                | \xi_i \rangle = - H | \Psi_i \rangle 
                                - \varepsilon_i | \Psi_i \rangle
             \f]
          -# get search direction
             \f[
                | X_i > = \Delta t ( |\xi_i \rangle 
                                   - \varepsilon_i | \Psi_i \rangle )
             \f]
          -# improve state
             \f[
                | \Psi_i^{\mathrm{new}} \rangle
                  = 
                | \Psi_i \rangle - | X_i \rangle
             \f]
          This scheme is repeated for all state belonging to all \b k points
          of all spin channels. The result is a new set of wavefunctions
          which is used to reconstruct the charge density. The entire scheme
          is repeated until convergence is reached.
          \author Sixten Boeck */
      void steepestDescent (const SxSymbolTable *, bool calc=true);
      /** \brief All-state conjugate gradient scheme for semiconductors
      
          The all-state conjugate gradient scheme is a special scheme which
          is optimized for insulators and semiconducting systems. It is
          particulary fast because it can take advantage of blocking
          matrix operations. 
          
          In contrast to state-by-state schemes the wavefunction coefficients 
          are treated as matrix
          \f[
             \langle \mathbf{G} | \Psi^{(i)}_{i\mathbf{k}} \rangle 
                \Rightarrow 
             \mathbf{C}^{(i)}_{\mathbf{G}i} (\sigma\mathbf{k})
          \f]
          with \em i being the iteration number.
          The algorithms reads as follows:
          -# compute gradient \b G:  \f[
                \mathbf{G}(\sigma\mathbf{k}) = -H \mathbf{C}{\sigma\mathbf{k}}
             \f]
          -# evaluate preconditioned gradient \b P:  \f[
                \mathbf{P}(\sigma\mathbf{k}) = K \mathbf{G}(\sigma\mathbf{k})
             \f] with \b K being the preconditioner
          -# get search direction of the \em i th iteration \b X: \f[
                \mathbf{X}^{(i)}(\sigma\mathbf{k})
                   = \mathbf{P}(\sigma\mathbf{k})
                   + \kappa \mathbf{X}^{(i-1)}(\sigma\mathbf{k})
             \f] with \f[
                \kappa = \frac{\Re \mathrm{tr} 
                                  \mathbf{P}^{(i)\dagger} (\sigma\mathbf{k})
                                  \mathbf{G}^{(i)}        (\sigma\mathbf{k})
                              }
                              {\Re \mathrm{tr}
                                  \mathbf{P}^{(i-1)\dagger} (\sigma\mathbf{k})
                                  \mathbf{G}^{(i-1)}        (\sigma\mathbf{k})
                              }
             \f]
          -# improve wavefunction according to \f[
                \mathbf{C}^{(i+1)}(\sigma\mathbf{k})
              = \mathbf{C}^{(i)}(\sigma\mathbf{k})
              - \lambda \mathbf{X}^{(i)}(\sigma\mathbf{k})
             \f]

         The parameter \f$ \lambda \f$ is optained from a quadratic line
         minimization of
         \f[
              E_{\mathrm{trial}} 
            = E_{\mathrm{tot}} (\mathbf{C}-\lambda_{\mathrm{trial}} \mathbf{X})
         \f]
         The minimum value is then
         \f[
            \lambda = \frac{D}{2c}
         \f]
         with
         \f[
            c = \frac{1}{\lambda_{\mathrm{trial}}}
                 (E_{\mathrm{trial}} - (E + \lambda_{\mathrm{trial}} D))
         \f]

         The all-state conjugate gradient does not work with the 
         conventional Grahm-Schmidt orthogonalization. Instead it uses
         a matrix-based Löwdin orthogonalization. The non-othogonalized
         wavefunctions are written as \f$ \hat{\mathbf{C}} \f$:
         -# Compute the overlapp matrix \b S: \f[
               \mathbf{S} = \hat{\mathbf{C}}^\dagger \hat{\mathbf{C}}
            \f]
         -# Solve the eigensystem \f[
               \mathbf{S}\mathbf{v} = s \mathbf{v}
            \f]
         -# get the uniform transformation \f[
               \mathbf{U} = \mathbf{v}^\dagger \frac{1}{\sqrt{s}} \mathbf{v}
            \f]
         -# orthogonalize the wavefunctions \f[
               \mathbf{C}_{\bot} = \mathbf{U} \hat{\mathbf{C}}
            \f]
              

       */
      void allStateCG (const SxSymbolTable *cmd, bool calc=true);
      /** \brief 4-th order polynomial fit
        @return vector of polynomial coefficients
        \note This is an auxiliary function for the all-state
              conjugate-gradient minimization. It fits a 4-th order
              polynomial to 5 points:
              - energy f0 and gradient fp0 at x=0
              - energy f1 and gradient fp1 at x=x1
              - energy f2 at x=x2
      */
      static SxVector<Double> lineFit4th (double f0, double fp0,
                                          double x1, double f1, double fp1,
                                          double x2, double f2);
      /** \brief Hybrid scheme of fixed Hamiltonian diagonalization and
                 charge density mixing */
      void scfDiagonalization (const SxSymbolTable *, bool calc=true);
      //@}

      /** \brief Possible modes for minimizer routines */
      enum MinimizerMode {
         /// Don't calculate, just print parameters
         PrintOnly,
         /// Apply convergence criteria for SCF runs
         SCFrun,
         /// Apply convergence criteria for band structure runs
         BandStructureRun,
         /// Apply convergence criteria for band structure runs
         PropagateRun
      };

      //---------------------------------------------------------------------
      /**@name Fixed Hamiltonian diagonalizations
         In contrast to the Born-Oppenheimer minimizers the routines
         in this section keep the Hamiltonian explicitly fixed during
         the diagonalization. */
      //---------------------------------------------------------------------
      //@{
      /** \brief Diagonalizes a fixed Hamiltonian with a conjugate gradient scheme

          The state-by-state or band-by-band conjugate gradent scheme is
          an algorithm to diagonalize a fixed Hamiltonian as described
          in \ref Kresse93.  Fixed means that the electronic charge density 
          is \b not updated during the loop:
          \f[
             H[\varrho|_{\mathrm{fixed}}] | \Psi_{i \sigma \mathbf k} \rangle
              = 
             \varepsilon_{i \sigma \mathbf{k}}
                | \Psi_{i \sigma \mathbf{k}} \rangle 
          \f]
          The algorithm is applied subsequently for all states:

          -# Consider a single state. To simplify we will write
             \f[
                | \Psi_{i \sigma \mathbf{k}} \rangle 
                   \rightarrow 
                | \Psi_{i} \rangle
             \f]
          -# orthonormalize it to all energetically lower lying states
             \f[
                | \Psi_i \rangle
                   = | \Psi_i \rangle
                   - \sum_{j < i}
                        \langle \Psi_j | \Psi_i \rangle | \Psi_j \rangle
             \f]
          -# compute gradient 
             \f[ 
                | \xi_i \rangle = - H | \Psi_i \rangle 
                                - \varepsilon_i | \Psi_i \rangle
             \f]
          -# make (preconditioned) gradient orthonormal to lower lying states as well as to the current state
             \f[
                | K \xi_i \rangle
                   = | K \xi_i \rangle
                   - \sum_{j \leq i}
                        \langle \Psi_j | K \xi_i \rangle | \Psi_j \rangle
             \f]
          -# compute conjugate search direction (previous entities have primes)
             \f[
                \gamma  
                   = \frac{\mathrm{tr} \langle \xi_i  | K \xi_i  \rangle}
                          {\mathrm{tr} \langle \xi'_i | K \xi'_i \rangle}

             \f]
          -# compute search vector
             \f[
                | X_i \rangle = | K \xi_i \rangle + \gamma | X' \rangle
             \f]
          -# make search vector orthonormal to current state
             \f[
                | X_i \rangle
                   = | X_i \rangle
                   - \langle \Psi_i | X_i \rangle | \Psi_i \rangle
             \f]
          -# instead of line minimization use unitform transformation
             \f[
                \mathbf{M} 
                = \left(
                     \begin{matrix}
                       \varepsilon_i 
                     & \langle X_i | H | \Psi_i \rangle  \\
                       \langle \Psi_i | H | X_i \rangle 
                     & \langle X_i | H | X_i \rangle 
                     \end{matrix}
                 \right)
             \f]
             \f[
               \tan \theta 
                 = 
               \frac{1}{2}\frac{M_{10} + M_{01}}{M_{00} - M_{11}}
             \f]
             \f[
               \mathbf{U}
                 = \left(
                     \begin{matrix}
                         \cos \theta & \sin \theta \\
                        -\sin \theta & \cos \theta
                     \end{matrix}
                     \right) 
             \f]
          -# improve state according to
             \f[
                | \Psi_i^{\mathrm{new}} \rangle
                  = 
                     \mathbf{U}_{00} | \Psi_i \rangle 
                  +  \mathbf{U}_{01} | X_i \rangle
             \f]

          This precedure will be repeated for all states belonging to all
          \b k points and spin channels.

          Since the \f$ \varepsilon_i = \langle \Psi_i | H | \Psi_j \rangle \f$
          are not the correct one-particle energies a subspace diagonalization
          has to be performed afterwards. For each \b k point and each spin
          channel 
          \f[
             \xi_{i\mathbf{G}} = H c_{i\mathbf{G}}
          \f]
          is computed with \em c being the wavefunction coeffients in a matrix
          form. 
          \f[
             \tilde{H} = \Psi_{i\mathbf{G}}^\dagger \xi_{i\mathbf{G}}
          \f]
          The eigensolution gives the wanted rotation
          \f[
             \tilde{h} v = e v
          \f]
          and
          \f[
             c_{i\mathbf{G}}' = c_{i\mathbf{G}} \cdot v
          \f]
          
          \param  sloppyEmptyStates   if \em true, only 2 iterations will
                                      be spent in (almost) empty states.
          \author Sixten Boeck */
      void stateByStateCG (const SxSymbolTable *, bool sloppyConv=false,
                           bool printOnly = false);
      /** \brief Block version of state-by-state CG
          @param printOnly if true, print only parameters, do nothing
          @return number of unconverged states
        */
      int blockStateByStateCG (const SxSymbolTable *, MinimizerMode);
      int blockStateByStateCG (int maxItDiag = 5,
                               int defBlockSize = 64,
                               bool numericalLimit = false,
                               MinimizerMode runType = SCFrun,
                               int nSloppy = 0,
                               bool verbose = false);

      /** \brief RMM-DIIS algorithm for all states
        */
      void rmmDiis (const SxSymbolTable *cmd, bool sloppyConv,
                    bool printOnly = false);
      /** \brief RMM-DIIS algorithm (single state)
          \note The resulting state is not normalized. RMM-DIIS needs
          a final orthonormalization run.
       */
      void rmmDiis (int i, int iSpin, int ik,
                    int maxItDiag, bool verbose, bool sloppyConv);
      /** \brief Solve the non-linear DIIS problem

          @param h2In H^2       matrix elements of phi
          @param hssh (HS + SH) matrix elements of phi 
          @param s2   S^2       matrix elements of phi
          @param s    S         matrix elements of phi
          @param epsPtr pointer to eps (with initial guess)
          @param n    number of states
          @return the optimal alpha vector 

          This routine minimizes
          \f[ 
          \frac{\langle \psi | (H-\epsilon S)^2 |\psi \rangle}
               {\langle \psi | S | \psi \rangle}
          \f]
          where
          \f[ \psi = \sum_i \alpha_i \phi_i \f]
          by solving the eigenvalue problem
          \f[
          \langle \phi_i | H^2 - \epsilon (HS + SH) + \varepsilon^2 S^2| \phi_j\rangle
          \alpha_j = \langle \phi_i | S | \phi_j\rangle \alpha_j
          \f]

          The eigenvalue $\epsilon$ is optimized, too.
          When a solution is found, it is taken out of the
          search space and the search continues with the next state.
          
          \note epsPtr must point to double array of n (or more) elements

        */
      static SxDiracVec<TPrecCoeffG> 
      solveRmmDiis (const SxDiracMat<TPrecCoeffG> &h2In,
                    const SxDiracMat<TPrecCoeffG> &hsshIn,
                    const SxDiracMat<TPrecCoeffG> &s2In,
                    const SxDiracMat<TPrecCoeffG> &s,
                    double *epsPtr,
                    int n = 1);
      /** \brief Blocked version of RMM-DIIS algorithm for all k-points
        */
      void blockRmmDiis (const SxSymbolTable *cmd, bool sloppyConv, 
                         bool printOnly = false);
      /** \brief Blocked version of RMM-DIIS algorithm (single block)
        */
      void blockRmmDiis (int i, int iSpin, int ik, int blockSize, 
                         int maxItDiag, bool verbose, bool sloppyConv);

      /** \brief blockRMM-CG algorithm for all k-points
        */
      void blockRmmCG (const SxSymbolTable *cmd, bool sloppyConv, 
                       bool printOnly = false);
      /** \brief blockRMM-CG algorithm (single block)
        */
      void blockRmmCG (int i, int iSpin, int ik, int blockSize, 
                       int maxItDiag, bool verbose, bool sloppyConv);

      /** \brief Diagonalize a fixed Hamiltonian using a matrix eigensolver */
      void directDiagonalization (const SxSymbolTable *cmd, bool calc=true);
      /**
        \brief Diagonalization scheme optimized for bandstructure calculations
       */
      void bandStructureCG (const SxSymbolTable *cmd, bool calc);
      /**
        \brief Diagonalization scheme optimized for (H - eps)^2
       */
      void diffEqtn (const SxSymbolTable *cmd, bool calc);
      /**
        \brief minimization for differential equations
       */
      void hSqrCG (const SxSymbolTable *cmd, bool calc);
      /**  \brief Perform subspace diagonalization for each k-point
        */
      void subspaceDiagonalization (const SxSymbolTable *cmd, bool calc);
        /** \brief Diagonalize subspace Hamiltonian
          @param wavesPtr pointer to waves (input & output)
          @param subDiagSize If greater one, diagonalize only blocks
                 around the diagonal of this size.
          @param subDiagOvlp Overlap between diagonal blocks
          @return The computed eigenvalues
        Set up Hamiltonian matrix for given wavefunctions. Diagonalize it
        and set the waves correspondingly. 
        */
      SxDiracVec<TPrecEps> subspaceDiagonalization (PsiGI *wavesPtr,
                                                      int subDiagSize = -1,
                                                      int subDiagOvlp = -1);
      /** \brief Copy wave coefficients from k-point ik to k-point jk
        */
      void propagateWaves (int ik, int jk);
      //@}


      //---------------------------------------------------------------------
      /**@name File I/O
         These routines are to read and write \ref tutor_sxb which contain
         atomic structures and wavefunctions required for continuation runs. */
      //---------------------------------------------------------------------
      //@{
     void writeData (bool storeWaves, bool storeRho) const;
      //@}

     const SxPW &getWaves () const;

     //---------------------------------------------------------------------
     /**@name Species and structural data */
     //---------------------------------------------------------------------
     //@{
     SxAtomicStructure getStructure () const;
     
     //@}


     // --------------------------------------------------------------------
     /**@name Consistenccy tests */
     // --------------------------------------------------------------------
     /** @{ */
     /** \brief perform a linearize gradient check

         The linearized gradient check varifies the consistency of the
         total energy contribution and its corresponding gradient. For a
         given state a fixed gradient is taken. Along this gradient line
         both the total energy contribution and its tagential are computed.
         If in the Hamiltonian inconsistencies between the total energy
         contribution and its corresponding gradient exist the linearized
         gradient will not show up as a tangent. 
         Note, that the Hamiltonian contributions can be assembled on-the-fly
         using the hContrib input variable. Also note, that at least one
         attractive and a repulsive potential must be taken into account, e.g.,
         the kinetic energy and the Hartree energy. */
     void linGradTest (const SxSymbolTable *, bool calc=true);

      /** \brief linear extrapolation scheme for the wavefunction
           
        This scheme is used by the MD to improve the initial guess
        for the wavefunction after an ionic step has been performed
        */

     void extrapolateWaves (double, const SxPW &);
     /** @} */

   //protected:

      SxMesh3D  mesh;

      Real8 rhoMixing;
      Real8 foccMixing;
      Real8 deltaT;
      Real8 dPsiConv;
      Real8 dEnergy, dEps, dRelEps, dRelR;
      Real8 ekt, ektDefault;
      int   maxSteps;

      /// \name Convergence statistics 
      //@{
          /// Number of converged states due to absolute energy
      int nConvAbs,
          /// Number of converged states due cos numerical limit
          nCosAcc, 
          /// Number of states not converged (iterations exceeded)
          nConvIt, 
          /// Number of converged states due to low occupation
          nConvOcc, 
          /// Number of converged states due relative convergence
          nConvRel;
      /// Number of steps after convergence in a block algorithm
      int nConvCycle;
          /// Number of steps until convergence due to absolute energy
      int nStepsAbs,
          /// Number of steps until convergence due cos numerical limit
          nStepsCos,
          /// Number of steps before max. iterations exceeded
          nStepsIt,
          /// Number of steps until convergence due to low occupation
          nStepsOcc,
          /// Number of steps until convergence due relative convergence
          nStepsRel;

   protected:
      /// Print a single convergence line
      static void printConvStat (const SxString &name, int nConv, int nStep=0);
      /// Print convergence statistics
      void printConvergenceStatistics () const;

      //@}
      /// Print energies to energy.dat
      static void printEnergyDatLine(int it, double eTot, double freeEnergy,
                                     double eBand, double entropy);
   public:

      SxPtr<SxHamiltonian> hamPtr;
      SxPtr<SxPWSet>  wavesPtr;
      SxFermi         fermi;

      /** \brief Pointer to atomic orbital basis
        \note This need not be set.
        */
      SxPtr<SxAOBasis> aoBasisPtr;

      SxConstPtr<SxMuPW> muBasisPtr;

      int nSpin, nStates, nStatesChi, nValStates;
      Real8 nElectrons;
      SxAtomicStructure  structure;

      /// Pointer to species (i.e. potential) data
      SxPtr<SxSpeciesData> potPtr;
      /** \brief Set up an auxiliary ao basis
          \note If file is empty or does not exist, atomic waves are
          extracted from the potential.
        */
      void setupAOBasis (const SxString &TBInFile, const SxGkBasis &gk);
      bool printHartree;
      bool keepRho, keepOcc;

      SxGBasis  G;
      SxRBasis  R;


      void write (const SxString &) const;
      void printParameters (const SxString &prefix,
                            const SxString &title) const;

};

#endif /* _SX_HAM_SOLVER_H_ */
