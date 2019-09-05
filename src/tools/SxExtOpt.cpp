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

#include <SxCLI.h>
#include <SxPotential.h>
#include <SxNeighbors.h>
#include <SxStructOpt.h>
#include <SxSimpleParser.h>


class SxExtPot : public SxPotential
{
   public:
      double energy;
      SxSpeciesData speciesData;

      /** \brief Check whether the provided minimization command is known */
      virtual bool isRegistered (const SxSymbolTable *) const { return true; }

      virtual void execute (const SxSymbolTable *, bool =true) { SX_EXIT; }

      virtual SxAtomicStructure getForces (const SxAtomicStructure &tau,
                                           const SxSymbolTable *cmd=NULL);
      virtual PrecEnergy getEnergy () const { return energy; }
      virtual SxSpeciesData getSpeciesData () const { return speciesData; }

      SxExtPot (const SxString &controlFile,
                const SxString &responseFile);
      ~SxExtPot ();


   protected:
      FILE *control, *response;
      SxArray<int> map; // map external atom id -> internal one

      void command(const char *cmd) const;
      void responseError (const char* what = NULL);
      bool supported ();
      int readInt (const char *what, int no = -1);

   public:
      void readStructureExtern (SxAtomicStructure *structPtr);

      inline void start () { command ("run"); }
      void end ();
};

void SxExtPot::command(const char *cmd) const
{
   SX_CHECK(control);
   fprintf (control, "%s\n", cmd);
   fflush (control);
   sxprintf ("Sent command: %s\n", cmd);
}

void SxExtPot::responseError (const char* what)
{
   if (what)  cout << "Failed to read " << what << endl;
   cout << "Next part of line is '";
   for (int ic = 0; ic < 100; ic++)  {
      if (feof(response)) {
         cout << "<eof>";
         break;
      }
      int c = fgetc(response);
      if (c == '\n') break;
      cout << c;
   }
   cout << "'" << endl;
   command ("end");
   while (!feof(response)) fgetc(response);
   fclose (control);
   fclose (response);
   control = response = NULL;
   SX_QUIT;
}

SxAtomicStructure SxExtPot::getForces (const SxAtomicStructure &tau,
                                      const SxSymbolTable *)
{
   SxAtomicStructure forces = tau.getNewStr ();
   forces.set (Coord(0.,0.,0.));
   command ("set positions");
   SX_LOOP(ia)  {
      Coord pos = tau.constRef(map(ia));
      fprintf(control, "%.8f %.8f %.8f\n", pos(0), pos(1), pos(2));
   }
   command ("run");
   command ("get forces");
   SX_LOOP(ia)  {
      Coord &f = forces.ref(map(ia));
      if (fscanf (response, " %lf %lf %lf", &f(0), &f(1), &f(2)) != 3)  {
         cout << "Error: failed to read force no. " << (ia + 1) << endl;
         responseError ();
      }
   }
   command("get energy");
   if (fscanf (response, " %lf", &energy) != 1)  {
      cout << "Error: failed to read energy" << endl;
      responseError ();
   }

   return forces;
}

SxExtPot::SxExtPot (const SxString &controlFile,
                    const SxString &responseFile)
{
   (cout << "Opening control file '" << controlFile << "'...").flush ();
   control = sxfopen (controlFile, "w");
   cout << endl;
   (cout << "Opening result file '" << responseFile << "'...").flush ();
   response = sxfopen (responseFile, "r");
   cout << endl;
}

void SxExtPot::end ()
{
   command ("end");
   while (!feof(response))  {
      (void)fscanf (response, " ");
      char buffer[10240];
      buffer[0]=0;
      (void)fgets(buffer, 10240, response);
      for (int i = 0; buffer[i]; i++) 
         if (buffer[i] >= 'A' && buffer[i] <= 'Z')
            buffer[i] |= char('a' - 'A'); // lower case
      if (buffer[0] == 'o' && buffer[1] == 'k') break;
      cout << "end was not confirmed; found: " << buffer << endl;
   }
   (cout << "Closing control file." << endl).flush ();
   fclose(control);
   (cout << "Closing result file." << endl).flush ();
   fclose(response);
   control=response=NULL;
}

SxExtPot::~SxExtPot ()
{
   if (response || control) end ();
}

int SxExtPot::readInt (const char *what, int no)
{
   int res;
   if (fscanf (response, " %d", &res) != 1)  {
      cout << "Failed to read " << what;
      if (no >= 0) cout << " no. " << (no + 1);
      cout << endl;
      responseError ();
   }
   return res;
}

/*
bool SxExtPot::supported ()
{
   if (feof(response)) responseError ("answer");
   char buffer[10];
   if (fscanf(" %s", 10, buffer) != 1) responseError("answer");
   SxString answer(buffer);
   if (answer == "no") return false;
   if (answer == "ok") return true;
   cout << "Failed to get confirmation that last command is supported" << endl;
   cout << "Response is '" << answer << "'" << endl;

   command ("end");
   while (!feof(response)) fgetc(response);
   fclose (command);
   fclose (response);
   SX_QUIT;
}
*/

void SxExtPot::readStructureExtern (SxAtomicStructure *structPtr)
{
   SX_CHECK (structPtr);
   SxAtomicStructure &structure = *structPtr;
   SxCell cell;
   command ("get cell");
   for (int a = 0; a < 3; a++)  {
      if (fscanf (response, " %lf %lf %lf", &cell(0,a),
                                            &cell(1,a),
                                            &cell(2,a)) != 3)
      {
         cout << "Failed to read cell basis no. " << (a+1) << endl;
         responseError ();
      }
   }
   cell.setup ();
   if (structure.getNAtoms () > 0)  {
      if ((structure.cell - cell).absSqr ().sum () > 1e-8)  {
         cout << "Cell mismatch: have " << structure.cell << endl;
         cout << "from control        " << cell << endl;
         command ("end");
         while (!feof(response)) fgetc(response);
         fclose (control);
         fclose (response);
         SX_QUIT;
      }
   }

   command ("get natoms");
   int nAtoms = readInt ("number of atoms");

   command ("get positions");
   SxUniqueList<SxString> elementSymbols;
   SxArray<int> iSpecies (nAtoms);
   SxArray<Coord> coords (nAtoms);
   {
      char symbol[3];
      symbol[2] = 0;
      for (int ia = 0; ia < nAtoms; ++ia)  {
         Coord &pos = coords(ia);
         if (fscanf (response, " %lf %lf %lf %2s",
                     &pos(0), &pos(1), &pos(2), symbol) != 4)
         {
            cout << "Failed to atom no. " << (ia+1) << endl;
            responseError ();
         }

         symbol[0] = char(symbol[0] & 0xdf); // capital letter
         if (symbol[1] > 0) symbol[1] |= 0x20; // small letter
         SxString chem(symbol);
         elementSymbols << chem;
         iSpecies(ia) = (int)elementSymbols.findPos (chem);
         // advance to new line
         while (fgetc(response) != '\n' && !feof(response)) ;
      }
   }

   map.resize (nAtoms);
   if (structure.getNAtoms () > 0)  {
      const SxArray<SxString> &elem = structure.getElements ();
      if (nAtoms != structure.getNAtoms ())  {
         cout << "Number of atoms mismatch: " << endl
              << "From SxExtOpt: " << structure.getNAtoms () << endl
              << "From control : " << nAtoms << endl;
         command ("end");
         while (!feof(response)) fgetc(response);
         fclose (control);
         fclose (response);
         SX_QUIT;
      }
      SxGrid grid (structure, 3);
      // map to existing structure
      for (int ia = 0; ia < nAtoms; ++ia)  {
         map(ia) = structure.find (coords(ia), grid);
         // --- consistency checks
         if (map(ia) < 0)  {
            cout << "Structure mismatch: atom " << (ia + 1)
                 << " from control @ " << coords(ia)
                 << "could not be found" << endl;
            command ("end");
            while (!feof(response)) fgetc(response);
            fclose (control);
            fclose (response);
            SX_QUIT;
         }
         int isOld = structure.getISpecies(map(ia));
         if (elementSymbols(iSpecies(ia)) != elem(isOld))  {
            cout << "Element mismatch for atom " << (ia + 1)
                 << " from control @ " << coords(ia) << ":" << endl;
            cout << "From SxExtOpt: " << elem(isOld) << endl;
            cout << "control file: " << elementSymbols(iSpecies(ia)) << endl;
            command ("end");
            while (!feof(response)) fgetc(response);
            fclose (control);
            fclose (response);
            SX_QUIT;
         }
         // copy in precise position
         structure.ref (map(ia)) = coords(ia);
      }
   } else {
      structure.cell = cell;
      int nSpecies = (int)elementSymbols.getSize ();
      structure.resize (nAtoms);
      // --- set speciesData
      speciesData.chemName = elementSymbols;
      speciesData.elementName = elementSymbols;
      // --- set up atom info
      SxPtr<SxAtomInfo> atomInfo = SxPtr<SxAtomInfo>::create (nSpecies);
      atomInfo->nAtoms.set (0);
      SX_LOOP(ia) atomInfo->nAtoms(iSpecies(ia))++;
      atomInfo->setupOffset ();
      atomInfo->meta.attach (SxAtomicStructure::Elements, speciesData.chemName);

      structure.replaceInfo (atomInfo);

      // --- map atoms
      SxArray<int> idAtom(nSpecies);
      idAtom.set (0);
      SX_LOOP(ia)  {
         int is = iSpecies(ia);
         map(ia) = atomInfo->offset (is) + (idAtom(is)++);
         structure.ref(map(ia)) = coords(ia);
      }

      structure.updateSymmetries ();
   }
}

int main (int argc, char **argv)
{
   // --- init S/PHI/nX Utilities
   initSPHInXMath ();

   // --- parse command line
   SxCLI cli (argc, argv);

   SxString inFile = cli.option ("-i|--input", "input file",
                                 "input file").toString ("optim.sx");
   SxString controlFile = cli.option ("-c|--control", "control file",
                                      "control file for the external program")
                          .toString ("");
   SxString responseFile = cli.option ("-r|--response", "response file",
                                      "response file from the external program")
                          .toString ("");

   cli.finalize ();

   // --- read input
   SxParser parser;
   parser.setValidation (false);
   SxConstPtr<SxSymbolTable> tree = parser.read (inFile, "std/optim.std");
   SxAtomicStructure structure;

   SYMBOLPARSE(&*tree)  {
      if (controlFile.getSize () == 0)
         controlFile = SYMBOLGET("controlFile") || "control";
      if (responseFile.getSize () == 0)
         responseFile = SYMBOLGET("responseFile") || "response";
      SYMBOLGROUP("structure")  {
         structure = SxAtomicStructure (SYMBOLGROUP_TABLE);
      }
   }

   SxExtPot extPot(controlFile, responseFile);
   extPot.readStructureExtern (&structure);

   extPot.start ();

   {
      SxStructOpt structOpt(&extPot, structure);
      SxSymbolTable *main = tree->getGroup("main");
      SxSymbolTable *cmd  = NULL;
      for (cmd  = main->begin();
           cmd != NULL;
           cmd  = cmd->nextSibling())
      {
         if (SxStructOpt::isRegistered (cmd))
            structOpt.execute  (cmd);
      }
   }
   extPot.end ();
}

