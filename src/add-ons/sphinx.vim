" ---------------------------------------------------------------------------
" VIM syntax highlighting file for S/PHI/nX
" ---------------------------------------------------------------------------
" Installation:
"
"   (1) create file $HOME/.vim/syntax/syntax.vim containing
"       the following text:
"
"          augroup filetypedetect
"          au! BufRead,BufNewFile *.sx   setfiletype sphinx
"          augroup END
"
"   (2) If you don't have a $HOME/.vim/syntax folder yet create it now
"
"         cd ~
"         mkdir -p .vim/syntax
"
"   (3) Copy or link this file to $HOME/.vim/syntax/sphinx.vim
"
"         cd ~
"         cd .vim/syntax
"
"         cp     <YOUR_SPHINX_PATH>/add-ons/sphinx.vim .
"       OR
"         ln -sf <YOUR_SPHINX_PATH>/add-ons/sphinx.vim .
"
"   (4) Add to $HOME/.vimrc
"
"         so $HOME/.vim/syntax/syntax.vim
"         sy on
"
" ---------------------------------------------------------------------------
" Authors: Sixten Boeck         <boeck@sphinx.de>
"          Christoph Freysoldt  <freysoldt@mpie.de>
" Date:    20/04/03
" ---------------------------------------------------------------------------

syn clear

" syn case ignore

syn keyword sxCommand     format
syn keyword sxCommand     sphinx nextgroup=sxReqSemicolon
syn keyword sxCommand     include
syn keyword sxCommand     needs requires global verboseMode

" --- structure
syn keyword sxGroup       structure species atom symmetry operator
syn keyword sxVariable    element cell
syn keyword sxVariable    operator 
syn keyword sxVariable    potential name element valenceCharge
syn keyword sxVariable    lMax lLoc lcaoOrbitals atomicRhoOcc
syn keyword sxVariable    rGauss reciprocalMass dampingMass ionicMass
syn keyword sxAttrib      movable relative nextgroup=sxReqSemicolon
syn keyword sxAttrib      movableX movableY movableZ nextgroup=sxReqSemicolon

" --- pseudoPot, skData
syn keyword sxGroup       skData pseudoPot 
syn keyword sxVariable    skFilesPath nextgroup=sxReqSemicolon
syn keyword sxVariable    skElementName
syn keyword SxAttrib    realSpaceProjectors nextgroup=sxReqSemicolon

" --- basis
syn keyword sxGroup     basis kPoint kPoints qPoint from to
syn keyword sxGroup     tbBasis
syn keyword sxVariable  kUnits folding eCut mesh meshAccuracy dK eCutChi gCut
syn keyword sxVariable  coords weight nPoints label
syn keyword sxAttrib    saveMemory

" --- Hamiltonian
syn keyword sxGroup       Hamiltonian PWHamiltonian 
syn keyword sxVariable    nEmptyStates nExcessElectrons ekt xc nlBlockSize
syn keyword sxAttrib      LDA PBE PBE_LDA EXX READ_VXX nextgroup=sxReqSemicolon
syn keyword sxAttrib      spinPolarized nextgroup=sxReqSemicolon
syn keyword sxAttrib      dipoleCorrection linearVacMixing nextgroup=sxReqSemicolon
" these are not followed necessarily by a semicolon, so no 
" nextgroup=sxReqSemicolon
syn keyword sxAttrib      CALC_ALL CALC_NONE CALC_KIN CALC_V_HARTREE
syn keyword sxAttrib      CALC_V_X CALC_V_C CALC_V_LOC CALC_RHO
syn keyword sxAttrib      CALC_V_NL CALC_V_XC CALC_V_EFF CALC_V_SCR 
syn keyword sxAttrib      CALC_V_EXT CALC_X CALC_DEFAULT
syn keyword sxAttrib      EXX_WRITE_ALL EXX_WRITE_NONE EXX_WRITE_CHI_MATRIX
syn keyword sxAttrib      EXX_WRITE_CHI_ROWDIAG EXX_WRITE_E_G EXX_WRITE_VXR
syn keyword sxAttrib      EXX_WRITE_VXG

" --- initialGuess
syn keyword sxGroup       initialGuess waves rho lcao occ atomicSpin exchange
syn keyword sxVariable    file maxSteps rhoMixing spinMoment z
syn keyword sxAttrib      fromWaves atomicOrbitals nextgroup=sxReqSemicolon
syn keyword sxAttrib      keepWavesOnDisk random nextgroup=sxReqSemicolon
syn keyword SxAttrib    basisDecomposition changeCharge pawBasis pulayForces nextgroup=sxReqSemicolon

syn keyword SxGroup       bands occupations xcPotential
syn keyword SxVariable    focc foccMixing range values

" --- tbInitialGuess
syn keyword sxGroup       tbInitialGuess
syn keyword sxVariable    withoutERepulsive nextgroup=sxReqSemicolon

" --- main.elecMinim
syn keyword sxGroup     main
syn keyword sxGroup     SD WS DJ CCG DIAG linGrad SCF SH
syn keyword sxVariable  maxSteps dEnergy dCharge dPsi printSteps nMixingSteps
syn keyword sxVariable  deltaT gamma i iSpin ik lambdaMin lambdaMax nSteps
syn keyword sxVariable  maxStepsCCG dEpsCCG dRelEps dRelRes
syn keyword sxVariable  nPulaySteps kerkerDamping
syn keyword sxVariable  update
syn keyword sxAttrib    KERKER LINDHARD CSRB nextgroup=sxReqSemicolon
syn keyword sxAttrib    noDiag avgDensityDiag avgEpsDiag
syn keyword sxVariable  hContrib mixingMethod spinMixing
syn keyword sxAttrib    LINEAR PULAY nextgroup=sxReqSemicolon
syn keyword sxAttrib    s p d f
syn keyword sxAttrib    keepRhoFixed keepOccFixed keepSpinFixed nextgroup=sxReqSemicolon
syn keyword sxAttrib    useFullBasis calcForces nextgroup=sxReqSemicolon
syn keyword sxAttrib    residueProfile nextgroup=sxReqSemicolon
syn keyword sxVariable  eCutDiag nStatesDiag
syn keyword SxAttrib    finalDiag initialDiag testLineMinim nextgroup=sxReqSemicolon
syn keyword SxVariable  xcMeshDensity dEnergyLow kappa testLineStep

syn keyword sxGroup       bandStructure subspaceDiag
syn keyword sxAttrib      autoSteps verbose printResidue propagateWaves
syn keyword sxVariable    maxSize overlap dEps nSloppy

syn keyword SxGroup     blockCCG nonDiagonal preconditioner scfDiag
syn keyword SxAttrib    noRhoStorage noWavesStorage numericalLimit nextgroup=sxReqSemicolon
syn keyword SxVariable  blockSize dielecConstant dumpTime scaling type
syn keyword SxVariable  dRelR dSpinMoment maxResidue spinScaling

" --- main.dampedNewton, quasiNewton, frozenPhonon, molDyn, synchronousTransit

syn keyword sxGroup      initStructure initHessian output 
syn keyword sxGroup      convergence dofRange performance
syn keyword sxGroup      initHistory randomVel devAtoms
syn keyword sxGroup      integrator thermostat
syn keyword sxGroup      initialStructure finalStructure  dofRange  
syn keyword sxGroup      linQN
syn keyword sxVariable   diag maxStepLength nProjectors
syn keyword SxGroup     extControl
syn keyword SxAttrib    noForces driftFilter nextgroup=sxReqSemicolon
syn keyword SxGroup     ricQN
syn keyword SxVariable  softModeDamping

syn keyword sxVariable   initIdentity  file  samplePoints  saveStructure  
syn keyword sxVariable   saveWaves	saveHist    freezeRot   freezeIt 
syn keyword sxVariable   constraints     dAvgForceComponent  dEnergyStruct  
syn keyword sxVariable   dMaxForceComponent  maxStructSteps
syn keyword sxVariable   dAvgStructComponent   dMaxStructComponent 
syn keyword sxVariable   saveHessian  initEkin  gas deltaE  dof 
syn keyword sxVariable   temperature scheme  order  dt   timeSteps   
syn keyword sxVariable   deviation startDof  endDof  

syn keyword sxAttrib     pass optimizeRho expertOutput thirdOrderCor  nextgroup=sxReqSemicolon
 
syn keyword sxAttrib     writeHist   shiftToSticks	rescale   nextgroup=sxReqSemicolon
           
syn keyword sxAttrib     extrapolateWaves nextgroup=sxReqSemicolon
  
syn keyword sxGroup       dampedNewton quasiNewton molDyn frozenPhonon
syn keyword sxGroup       synchronousTransit 

" from qn.std
syn keyword SxGroup     bornOppenheimer QN
syn keyword SxVariable  dF dFavg dX dXavg
syn keyword SxGroup     hessian
syn keyword SxAttrib    Fischer saveRelaxHist Schlegel withAngles nextgroup=sxReqSemicolon
syn keyword SxVariable  angleConstant bondConstant planeCutLimit rmsThreshold typifyThreshold

" --- EXX specific keywords
syn keyword sxGroup     EXX_LOOP relaxRho
syn keyword sxVariable  writeControl
syn keyword sxAttrib    restart nextgroup=sxReqSemicolon


" --- isixServer and user
syn keyword sxGroup       isixServer user
syn keyword sxVariable    host port uid passwd

" from defects.std
syn keyword SxGroup     chemicalPotentials defect state band
syn keyword SxVariable  configurations nAtoms charge energy kBoltzmann
syn keyword SxVariable  siteConcentration energyUnit lengthUnit
syn keyword SxVariable  mass degeneracy value diffusionBarrier
syn keyword SxVariable  diffusionPrefactor diffusionCharge

" from paw.std
syn keyword SxGroup     aoBasis HubbardU PAWHamiltonian pawPot site spinConstraint vExt xcMesh
syn keyword SxAttrib    checkOverlap fromPotential kjxc useProjG nextgroup=sxReqSemicolon
syn keyword SxVariable  angularGrid lMaxRho nRadGrid projectorType coreX alphaHybrid constraint fixedDipoleZ gyromagneticRatio omegaHSE rSoft shift zField coreWaves potType
syn keyword SxGroup     MO
syn keyword SxVariable  mMO nInterpolate sign cutWidth minDist setupBoxSize

" from structure.std
syn keyword SxVariable  epsSym force movableLine spin

" from sphinx.std
syn keyword SxGroup     charged compute externalField leaveOnRemote nlEES peField rsProj StilWeb strain synchronize
syn keyword SxAttrib    actOnNuclei bandstructure electrons nextgroup=sxReqSemicolon

" from parallelHierarchy.std
syn keyword SxGroup     level
syn keyword SxVariable  members siblings workers

" from pse.std
syn keyword SxGroup     chemElements elem
syn keyword SxVariable  symbol number atomicRadius covalentRadius CPKcolor

" from imagecharge.std
syn keyword SxGroup     background isolated Qslab
syn keyword SxAttrib    oldBroadening nextgroup=sxReqSemicolon
syn keyword SxVariable  betaPara betaZ broadening dropV electrodeZ epsilon fromZ posZ toZ
syn keyword SxVariable  cut Q

" --- numbers
syn match  sxNumber "[-+]\=\(\<\d[[:digit:]_]*L\=\>\|0[xX]\x[[:xdigit:]_]*\>\)"
syn match  sxNumber "[-+]\=\<\d[[:digit:]_]*[eE][\-+]\=\d\+"
syn match  sxNumber "[-+]\=\<\d[[:digit:]_]*\.[[:digit:]_]*\([eE][\-+]\=\d\+\)\="
syn match  sxNumber "[-+]\=\<\.[[:digit:]_]\+\([eE][\-+]\=\d\+\)\="

" --- some constants defined in parameters.sx
syn keyword sxNumber ANGSTROEM_TO_BOHR BOHR_TO_ANGSTROEM
syn keyword sxNumber EV_TO_HARTREE HARTREE_TO_EV

" --- logicals
"syn match sxLogical "\(true\|TRUE\|false\|FALSE\|yes\|YES\|no\|NO\)"
syn keyword sxLogical true TRUE false FALSE yes YES no NO

" --- C++, Fortran, and shell like comments
syn region sxString start=/"/ end=/"/
syn region sxString start=/</ end=/>/ 

" --- C like comments
syn region sxComment    matchgroup=sxCommentStart start="/\*" matchgroup=NONE end="\*/" 
syntax match sxError    "\*/"

" --- check that each '=' is followed by ";", or next line starts with + or - 
syn match sxError /=[^[;]*\n *[^ +-;,]/hs=e
syn match sxError /=\zs[^;]*\ze\n *$/hs=e
" --- check that each + or - continuation line is followed by ";"
"     or next line starts with + or - 
syn match sxError /^ *[+-][^;]*\n *[^ +-[;]/hs=e
syn match sxError /^ *[+-]\zs[^;]*\ze\n *$/hs=e

syn match sxError /=[^;[}]*}\{-1,}/hs=e display
syn match sxError "][^;,]*$" display

syn match  sxComment /#.*$/
syn match  sxComment /!.*$/
syn region sxComment start="//" end="$"


" --- check for missing semicolons ("=" is allowed, however)
syn match sxReqSemicolon "\_[ \n]*\_[^ +\n;=]" contained
hi link sxReqSemicolon sxError

" --- check brackets
"syn region sxBracket start="{" end="}" contains=ALLBUT,sxReqSemicolon,sxBrError
"syn match sxBrError "}" containedin=ALLBUT,sxBracket
"hi link sxBrError sxError


" ---- cross linking
hi link sxCommand       Include    
hi link sxGroup         Statement
hi link sxVariable      Type
hi link sxAttrib        SpecialChar
hi link sxChemElem      Identifier
hi link sxComment       Comment
hi link sxCommentStart  Comment
hi link sxNumber        Number
hi link sxLogical       Identifier
hi link sxString        String
hi link sxPath          String
hi link sxError         Error

let b:current_syntax = "sphinx-input"
                                                 

" --- use always C indention for S/PHI/nX input files
set cindent															" CIndent...
set expandtab
set cinkeys=0{,0},:,0#,!,o,O,e,!<Tab>					" ...when...
set cinoptions="n2=10+17(0)100*100							" ...and how
