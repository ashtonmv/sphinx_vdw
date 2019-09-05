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

#include <SxSymbolTable.h>
#include <SxError.h>
#include <SxParser.tab.hpp>
#include <SxParser.h>
#include <stdlib.h>

#ifdef MSVC
double cbrt (double x)
{
   return pow (x, 0.333333333333333333333333333333);
}
#endif


#ifdef WIN32
#  ifndef MSVC
#     include <math.h>
      double cbrt (double x)  {
         return pow (x, 0.3333333333333333333333);
      }
#  endif /* MSVC */
#endif /* WIN32 */


SxSymbol::SxSymbol ()
{
   name = ""; type = 0; status = Defined; initialized = binary = false;
   parserLineNumber = 0;
   val = 0.; str = "";
   valList = new SxList<SxSymbol>;
}


SxSymbol::SxSymbol (const SxString &name_, double val_)
{
   name   = name_.removeWhiteSpace ();
   type   = VAR;
   status = Defined;
   val    = val_;
   str    = "";
   initialized = true;
   binary = false;
   parserLineNumber = 0;
   valList = new SxList<SxSymbol>;
}

SxSymbol::SxSymbol (const SxString &name_, const SxString &str_, 
                    enum AutoCast autoCast)
{
   name   = name_.removeWhiteSpace ();
   val = 0.;
   str = "";
   valList = new SxList<SxSymbol>;

   if (autoCast == NoAutoCast)  {
      type   = STR;
      str    = str_;
   } else {
      bool typeIdentified = false;
      // --- determine type
      if (str_.isDouble ()) {
         val = str_.toDouble ();
         type = VAR;
         typeIdentified = true;
      }

      // SKETCH:
      // if (!typeIdentified)  {          // --- list [1,2,3]
      //    if (   str_.contains("[") == 1 && str_.contains("]" == 1
      //        && str_.contains(","))  
      //    {
      //       line = str_;
      //       while (line.contains (","))  {
      //          term = regexp ("/^\s*(.+?),/");   // minimal matching
      //          line = replace $1 with ""
      //          *varList << SxSymbol ("tmpTerm", term, AutoCast);
      //       }
      //       typeIdentified = true;
      //       type = LIST;
      //    }
      // }

      // if (!typeIdentified)  {          // --- list [[..],[..],..]
      //    if (str_.contains("[")  {
      //       line = str_;
      //       while (line.contains ("[")  {
      //          term = regexp (/(\[.+\])/         // maximal matching
      //          line = replace $1 with ""
      //         *varList << SxSymbol ("tmpTerm", term, AutoCast);
      //       }
      //       check: line should be empty or contain only "['*]"
      //       typeIdentified = true;
      //       type = LIST;
      //    }     


      if (!typeIdentified)  {             // --- string
         str  = str_;
         type = STR;
      }
   }
   status = Defined;
   initialized = true;
   binary = false;
   parserLineNumber = 0;
}

SxSymbol::SxSymbol (const SxString &name_, double (*func_)(double))
{
   name    = name_.removeWhiteSpace ();
   type    = FUNC;
   status  = Defined;
   func    = func_;
   strfunc = NULL;
   val     = 0.;
   str     = "";
   initialized = true;
   binary = false;
   parserLineNumber = 0;
   valList = new SxList<SxSymbol>;
}

SxSymbol::SxSymbol (const SxString &name_,
                    SxString (*strfunc_)(const SxString &))
{
   name    = name_.removeWhiteSpace ();
   type    = STRFUNC;
   status  = Defined;
   func    = NULL;
   strfunc = strfunc_;
   val     = 0.;
   str     = "";
   initialized = true;
   binary = false;
   parserLineNumber = 0;
   valList = new SxList<SxSymbol>;
}


SxSymbol::SxSymbol (const SxSymbol &in)
{
   name     = in.name;
   type     = in.type;
   status   = in.status;
   val      = in.val;
   str      = in.str;
   func     = in.func;
   strfunc  = in.strfunc;
   initialized = in.initialized;
   binary   = in.binary;
   parserFilename = in.parserFilename;
   parserLineNumber = in.parserLineNumber;
   valList  = new SxList<SxSymbol>;
   *valList = *in.valList;
}


SxSymbol::~SxSymbol ()
{
   delete valList;
}

SxSymbol &SxSymbol::operator= (const SxSymbol &in)
{
   if ( this == &in )  return *this;
   name     = in.name;
   type     = in.type;
   status   = in.status;
   val      = in.val;
   str      = in.str;
   func     = in.func;
   strfunc  = in.strfunc;
   initialized = in.initialized;
   binary   = in.binary;
   parserFilename = in.parserFilename;
   parserLineNumber = in.parserLineNumber;
   if (valList) delete valList;
   valList  = new SxList<SxSymbol>;
   *valList = *in.valList;
   return *this;
}



void SxSymbol::append (double val_)
{
   SX_CHECK (type == LIST);
   valList->append (SxSymbol ("", val_));
}


void SxSymbol::append (const SxString &str_)
{
   SX_CHECK (type == LIST);
   initialized = true;
   binary = false;
   valList->append (SxSymbol ("", str_));
}


void SxSymbol::append (SxSymbol &sym_)
{
   initialized = true;
   binary = false;
   valList->append (sym_);
}

void SxSymbol::prepend (SxSymbol &sym_)
{
   initialized = true;
   binary = false;
   valList->prepend (sym_);
}

SxSymbol SxSymbol::flatten ()
{
   if (type == LIST)  {                       // [....]
      SxSymbol sym ("", 0.);
      SxList<SxSymbol>::Iterator it;
      if (valList->getSize() == 0)  {          // [] 
         return *this;
      }  else if (valList->getSize() == 1)  {  // [[...]]
         sym = ( *valList->begin() );
      }  else  {
         bool toFlatten = true;
         for (it=valList->begin(); it!=valList->end(); it++)  {
            if ( (*it).type != LIST )  {
               toFlatten = false;
               break;
            }
         }
         if (toFlatten)  {                    // [[...],[...]]
            sym.type = LIST;
            for (it=valList->begin(); it!=valList->end(); it++)  {
               sym.valList->append ( (*it).flatten () );
            }
         }  else  {
            return *this;
         }
      }
      return sym;
   }  else  {
      return *this;                            // value 123.456;
   } 
}


SxSymbol::SymbolType SxSymbol::getType () const
{
   // --- this has to be done explicitly because SxParser.tab.h is unknown
   //     when dependencies are created. It is built by flex.
   if      (type == NUM)      return Number;
   else if (type == STR)      return String;
   else if (type == VEC)      return Vector;
   else if (type == VAR)      return Variable;
   else if (type == FUNC)     return Function;
   else if (type == STRFUNC)  return StringFunction;
   else if (type == LIST)     return List;
   else                       return Unknown;
}


const SxString &SxSymbol::getName () const
{
   return name;
}

bool SxSymbol::isDefined () const
{
   return status == Defined;
}


void SxSymbol::setDefined (bool defined)
{
   status = defined ? Defined : Undefined;
}


SxList<double> SxSymbol::toList () const
{
   SxList<double> res;
   SxList<SxSymbol>::ConstIterator it;
   if (valList->getSize() > 0)  {
      for (it=valList->begin(); it!=valList->end(); it++)  {
         if ( (*it).type == VAR )  res.append ( (*it).val );
         if ( (*it).type == LIST)  res.append ( (*it).toList() );
      }
   }  else  {
      res.append (val);
   }
   return res;
}


SxList<int> SxSymbol::toIntList () const
{
   SxList<int> res;
   SxList<SxSymbol>::ConstIterator it;
   if (valList->getSize() > 0)  {
      for (it=valList->begin(); it!=valList->end(); it++)  {
         if ( (*it).type == VAR )  res.append ( ::toInt ((*it).val) );
         if ( (*it).type == LIST)  res.append ( (*it).toIntList() );
      }
   } else {
      res.append ((int)val);
   }
   return res;
}


SxList<SxString> SxSymbol::toStringList () const
{
   SxList<SxString> res;
   SxList<SxSymbol>::ConstIterator it;
   if (valList->getSize() > 0)  {
      for (it=valList->begin(); it!=valList->end(); it++)  {
         if ( (*it).type == STR )  res.append ( ((*it).str) );
         if ( (*it).type == LIST)  res.append ( (*it).toStringList() );
      }
   } else {
      res.append (str);
   }
   return res;
}


const SxString &SxSymbol::toString () const
{
   return str;
}

bool SxSymbol::toBool () const
{
   return ( (int)val != 0);
}


bool SxSymbol::toAttribute () const
{
   int v = (int)val;
   // --- val=-1 set in SxParser.l, tag  {ID}
   if      (v == -1)  // attribute only, e.g.       calcForces;
      return true;
   else if (v ==  0)  // attribute switch off, e.g. calcForces = false;
      return false;
   else               // attribute switch on, e.g.  calcForces = true;
      return true;   
}


int SxSymbol::toInt () const
{
   return (int)lround(val);
}


double SxSymbol::toReal () const
{
   return val;
}


int SxSymbol::getRank () const
{
   if (type == LIST)  {                       // [....]
      SxList<SxSymbol>::ConstIterator it;
      if (valList->getSize() == 0)  {          // [] 
         return 1;  // empty vector
      }  else if (valList->getSize() == 1)  {  // [[...]]
         return 1+(*valList->begin()).getRank();
      }  else  {
         int rk, maxRank = 0;
         for (it=valList->begin(); it!=valList->end(); it++)  {
            rk = (*it).getRank();
            if (rk > maxRank)  maxRank = rk;
         }
         return 1+maxRank;
      }
      SX_EXIT;
      return -1;
   }  else  {
      return 0;                               // scalar 123.456;
   } 
   return 0;                                  // scalar
}


SxList<int> SxSymbol::getDimensions () const
{
   SxList<int> res;
   SxList<SxSymbol>::ConstIterator it;
   if (type == LIST)  {
      if (valList->getSize() > 0)  {
         for (it = valList->begin(); it != valList->end(); ++it)  {
            if ( (*it).type == LIST)  {
               res.append (static_cast<int>((*it).valList->getSize()));
               res.append ((*it).getDimensions());
            }
         }
      }
   }  else  {
      res.append (1);
   }
   return res;
}


SxSymbol SxSymbol::operator+ (const SxSymbol &in)
{
   if (status    == Undefined)  undefVarError (*this);
   if (in.status == Undefined)  undefVarError (in);

   SxList<SxSymbol>::Iterator it;
   if         (type == VAR && in.type == VAR)  {
      return SxSymbol ("", val + in.val);

   }  else if (type == STR && in.type == VAR)  {
      SxSymbol sym (*this);
      sym.str = str + in.val;
      return sym;
   }  else if (type == STR && in.type == STR)  {
      SxSymbol sym (*this);
      sym.str = str + in.str;
      return sym;

   }  else if (type == VAR && in.type == LIST)  {
      SxSymbol sym (in);
      for (it=sym.valList->begin(); it != sym.valList->end(); it++)  {
         (*it) = *this + (*it);
      }
      return sym;

   }  else if (type == LIST && in.type == VAR)  {
      SxSymbol sym (*this);
      for (it=sym.valList->begin(); it != sym.valList->end(); it++)  {
         (*it) = (*it) + in;
      }
      return sym;


   }  else if (type == LIST && in.type == LIST)  {
      if ( getRank() != in.getRank() )   {
         sxprintf ("Ranks of addition operands differ.\n");
         SX_EXIT;
      }
      if ( valList->getSize() != in.valList->getSize() )  {
         sxprintf ("Dimensions of addition operands differ.\n");
         SX_EXIT;
      }
      SxSymbol sym (*this);
      SxList<SxSymbol>::ConstIterator itRhs;
      for (it  = sym.valList->begin(), itRhs=in.valList->begin(); 
           it != sym.valList->end(); 
           it++, itRhs++)  
      {
         (*it) = (*it) + (*itRhs);
      }
      return sym;

   }  else  {
      sxprintf ("Corrupt addition statement\n");
      SX_EXIT;
   }
   return SxSymbol ();
}



SxSymbol SxSymbol::operator- (const SxSymbol &in)
{
   if (status    == Undefined)  undefVarError (*this);
   if (in.status == Undefined)  undefVarError (in);

   SxList<SxSymbol>::Iterator it;
   if         (type == VAR && in.type == VAR)  {
      return SxSymbol ("", val - in.val);

   }  else if (type == VAR && in.type == LIST)  {
      SxSymbol sym (in);
      for (it=sym.valList->begin(); it != sym.valList->end(); it++)  {
         (*it) = *this - (*it);
      }
      return sym;

   }  else if (type == LIST && in.type == VAR)  {
      SxSymbol sym (*this);
      for (it=sym.valList->begin(); it != sym.valList->end(); it++)  {
         (*it) = (*it) - in;
      }
      return sym;


   }  else if (type == LIST && in.type == LIST)  {
      if ( getRank() != in.getRank() )   {
         sxprintf ("Ranks of substraction operands differ.\n");
         SX_EXIT;
      }
      if ( valList->getSize() != in.valList->getSize() )  {
         sxprintf ("Dimensions of substraction operands differ.\n");
         SX_EXIT;
      }
      SxSymbol sym (*this);
      SxList<SxSymbol>::ConstIterator itRhs;
      for (it  = sym.valList->begin(), itRhs=in.valList->begin(); 
           it != sym.valList->end(); 
           it++, itRhs++)  
      {
         (*it) = (*it) - (*itRhs);
      }
      return sym;

   }  else  {
      sxprintf ("Corrupt substraction statement\n");
      SX_EXIT;
   }
   return SxSymbol ();
}


SxSymbol SxSymbol::operator* (const SxSymbol &in)
{
   if (status    == Undefined)  undefVarError (*this);
   if (in.status == Undefined)  undefVarError (in);

   SxList<SxSymbol>::Iterator it;
   if         (type == VAR && in.type == VAR)  {
      return SxSymbol ("", val * in.val);

   }  else if (type == VAR && in.type == LIST)  {
      SxSymbol sym (in);
      for (it=sym.valList->begin(); it != sym.valList->end(); it++)  {
         (*it) = *this * (*it);
      }
      return sym;

   }  else if (type == LIST && in.type == VAR)  {
      SxSymbol sym (*this);
      for (it=sym.valList->begin(); it != sym.valList->end(); it++)  {
         (*it) = (*it) * in;
      }
      return sym;


   }  else if (type == LIST && in.type == LIST)  {
      if ( getRank() != in.getRank() )   {
         sxprintf ("Ranks of multiplication operands differ.\n");
         SX_EXIT;
      }
      if ( valList->getSize() != in.valList->getSize() )  {
         sxprintf ("Dimensions of multiplication operands differ.\n");
         SX_EXIT;
      }
      SxSymbol sym (*this);
      SxList<SxSymbol>::ConstIterator itRhs;
      for (it  = sym.valList->begin(), itRhs=in.valList->begin(); 
           it != sym.valList->end(); 
           it++, itRhs++)  
      {
         (*it) = (*it) * (*itRhs);
      }
      return sym;

   }  else  {
      sxprintf ("Corrupt multiplication statement\n");
      SX_EXIT;
   }
   return SxSymbol ();
}



SxSymbol SxSymbol::operator/ (const SxSymbol &in)
{
   if (status    == Undefined)  undefVarError (*this);
   if (in.status == Undefined)  undefVarError (in);

   SxList<SxSymbol>::Iterator it;
   if         (type == VAR && in.type == VAR)  {
      if (fabs (in.val) < 1e-90)  {
         sxprintf ("Division by zero\n");
         SX_EXIT;
      }
      return SxSymbol ("", val / in.val);

   }  else if (type == VAR && in.type == LIST)  {
      SxSymbol sym (in);
      for (it=sym.valList->begin(); it != sym.valList->end(); it++)  {
         (*it) = *this / (*it);
      }
      return sym;

   }  else if (type == LIST && in.type == VAR)  {
      SxSymbol sym (*this);
      for (it=sym.valList->begin(); it != sym.valList->end(); it++)  {
         (*it) = (*it) / in;
      }
      return sym;


   }  else if (type == LIST && in.type == LIST)  {
      if ( getRank() != in.getRank() )   {
         sxprintf ("Ranks of division operands differ.\n");
         SX_EXIT;
      }
      if ( valList->getSize() != in.valList->getSize() )  {
         sxprintf ("Dimensions of multiplication operands differ.\n");
         SX_EXIT;
      }
      SxSymbol sym (*this);
      SxList<SxSymbol>::ConstIterator itRhs;
      for (it  = sym.valList->begin(), itRhs=in.valList->begin(); 
           it != sym.valList->end(); 
           it++, itRhs++)  
      {
         (*it) = (*it) / (*itRhs);
      }
      return sym;

   }  else  {
      sxprintf ("Corrupt division statement\n");
      SX_EXIT;
   }
   return SxSymbol ();
}



SxSymbol SxSymbol::operator^ (const SxSymbol &in)
{
   if (status    == Undefined)  undefVarError (*this);
   if (in.status == Undefined)  undefVarError (in);

   SxList<SxSymbol>::Iterator it;
   if         (type == VAR && in.type == VAR)  {
      return SxSymbol ("", pow(val,in.val));

   }  else if (type == VAR && in.type == LIST)  {
      SxSymbol sym (in);
      for (it=sym.valList->begin(); it != sym.valList->end(); it++)  {
         (*it) = *this ^ *it;
      }
      return sym;

   }  else if (type == LIST && in.type == VAR)  {
      SxSymbol sym (*this);
      for (it=sym.valList->begin(); it != sym.valList->end(); it++)  {
         (*it) = (*it) ^ in;
      }
      return sym;


   }  else if (type == LIST && in.type == LIST)  {
      sxprintf ("Matrix/Vector multiplication not implemented yet.\n");
      SX_EXIT;

   }  else  {
      sxprintf ("Corrupt exponent statement\n");
      SX_EXIT;
   }
   return SxSymbol ();
}



void SxSymbol::print () const
{
   switch (type) {
      case VAR     : cout << "Variable: "   << name << " = "; break;
      case STR     : cout << "String: "     << name << " = "; break;
      case FUNC    : cout << "Function: "   << name << " = "; break;
      case STRFUNC : cout << "StringFunc: " << name << " = "; break;
      case LIST    : cout << "List: "       << name << " = "; break;
   }
   cout << *this << endl;
}


std::ostream &operator<< (std::ostream &s, const SxSymbol &in)
{
   switch (in.type) {
      case VAR     : s << in.val;  break;
      case STR     : s << '"' << in.str << '"';  break;
      case FUNC    : s << &in.func << "(double)"; break;
      case STRFUNC : s << &in.strfunc << "(string)"; break;
      case LIST    : cout << "[ ";
                     ssize_t i, n = in.valList->getSize();
                     for (i=0; i < n; i++)  {
                        cout << (*in.valList)(i);
                        if (i < n-1)  cout << ", ";
                     }
                     cout << " ]";
                     break;
   }
   return s;
}



void SxSymbol::undefVarError (const SxSymbol &sym) const
{
   cout << SX_SEPARATOR;
   cout << "| Error: '" << sym.getName () << "' is undefined.\n";
   cout << SX_SEPARATOR;
   SX_QUIT;
}


SxSymbol operator- (const SxSymbol &in)
{
   if (in.status == SxSymbol::Undefined)  in.undefVarError (in);

   SxList<SxSymbol>::Iterator it;
   if         (in.type == VAR)  {
      return SxSymbol ("", -in.val);

   }  else if (in.type == LIST)  {
      SxSymbol sym (in);
      for (it=sym.valList->begin(); it != sym.valList->end(); it++)  {
         (*it) =  -(*it);
      }
      return sym;

   }  else  {
      sxprintf ("Corrupt substraction statement\n");
      SX_EXIT;
   }
   return SxSymbol ();
}






SxSymbolTable *&SxSymbolTable::getGlobalPtr ()
{
   static SxSymbolTable *tablePtr = NULL;
   return tablePtr;
}





SxSymbolTable::SxSymbolTable () 
   : parent (NULL), 
     name (""), 
     level (0),
     deprecateFormat (false),
     parserLineNumber (0)
{
   table("sin")    = SxSymbol("sin",    sin);
   table("arcsin") = SxSymbol("arcsin", asin);
   table("cos")    = SxSymbol("cos",    cos);
   table("arccos") = SxSymbol("arccos", acos);
   table("tan")    = SxSymbol("tan",    tan);
   table("arctan") = SxSymbol("arctan", atan);
   table("sqrt")   = SxSymbol("sqrt",   sqrt);
   table("cbrt")   = SxSymbol("cbrt",   cbrt);
   table("exp")    = SxSymbol("exp",    exp);
   table("log")    = SxSymbol("log",    log);
   table("read")   = SxSymbol("read",   sxReadFile);
   table("PI")     = SxSymbol("PI",     3.14159265358979323846264338327950);
   table("pi")     = SxSymbol("pi",     3.14159265358979323846264338327950);
   table("Pi")     = SxSymbol("Pi",     3.14159265358979323846264338327950);
   table("true")   = SxSymbol("true",   1.);
   table("TRUE")   = SxSymbol("TRUE",   1.);
   table("false")  = SxSymbol("false",  0.);
   table("FALSE")  = SxSymbol("FALSE",  0.);
   table("yes")    = SxSymbol("yes",    1.);
   table("YES")    = SxSymbol("YES",    1.);
   table("no")     = SxSymbol("no",     0.);
   table("NO")     = SxSymbol("NO",     0.);
}


SxSymbolTable::SxSymbolTable (SxSymbolTable *parentPtr, 
                              const SxString &lvlName)
   : parent (parentPtr),
     name (lvlName),
     level (parentPtr->level+1),
     deprecateFormat (false),
     parserLineNumber (0)
{
   SX_CHECK (parentPtr);
   parentPtr->addChild (this);
}


SxSymbolTable::~SxSymbolTable ()
{
   table.removeAll ();
   // --- remove child nodes
   SxList<SxSymbolTable *>::ConstIterator it;
   for (it=children.begin(); it!=children.end(); it++)  delete *it;
   // --- remove tmpList
   SxList<SxSymbol *>::ConstIterator sym;
   for (sym=tmpList.begin(); sym!=tmpList.end(); sym++)  {
      sxprintf ("WARNING: tmpList was not empty: %p\n", (void *)(*sym));
      delete *sym;
   }
}


const SxSymbol *SxSymbolTable::get (const SxString &name_, bool localOnly,
                                    bool throwException,
                                    const SxString &parentPath) const
{
   return 
      (const_cast<SxSymbolTable *>(this))->get (name_, localOnly,
                                                throwException, parentPath);
}

SxSymbol *SxSymbolTable::get (const SxString &name_, bool localOnly,
                              bool throwException,
                              const SxString &parentPath)
{
   SxString fullPath;

   // --- look in local symbol level
   //SX_DBG_MSG ("get " << name_);
   ssize_t idx = table.findKey (name_);
   if (idx >= 0)  {
      return &table.getValueIdx(idx);
   }

   if ( !localOnly )  {
      // --- look for symbol in upper levels
      if ( parent )  {
         fullPath = getName();
         if (parentPath != "")  fullPath += "." + parentPath;
         return parent->get (name_, localOnly, throwException, fullPath);
      }
   }

   // --- not found
   if (throwException)  {
      fullPath = parentPath + "." + name_;
      SxString mode ("");
      if (useDeprecateFormat())  mode = " (using deprecate file format)";
      SX_THROW ("The symbol '" + fullPath + "' was not found" 
               + mode + ".");
   }
   return NULL;
}


SxString SxSymbolTable::getName () const
{
   return name;
}


bool SxSymbolTable::hasAttribute (const SxString &str) const
{
   if (contains(str))  return get(str)->toAttribute();
   return false;
}


int SxSymbolTable::getNItems (const SxString &name_, bool localOnly) const
{
   int n = 0;
   // --- check node
   if (table.containsKey (name_))  {
      n = 1;
   }

   // --- look in local symbol level (sibling by sibling)
   if (parent)  {
      SxList<SxSymbolTable *>::ConstIterator sibling;
      for (sibling  = parent->children.begin(); 
           sibling != parent->children.end(); 
           sibling++)  
      {
         if ( (*sibling)->name == name_ )  n++;
      }
   } 

   if ( !localOnly )  {
      // --- look for symbol in upper levels
      if ( parent )  n += parent->getNItems (name_);
   }

   return n;
}


bool SxSymbolTable::contains (const SxString &name_, bool localOnly) const
{
   // --- check node
   if (table.containsKey (name_))  {
      return true;
   }

   // --- look in local symbol level (sibling by sibling)
   if ( !localOnly )  {
      if (parent)  {
         SxList<SxSymbolTable *>::ConstIterator sibling;
         for (sibling  = parent->children.begin(); 
              sibling != parent->children.end(); 
              sibling++)  
         {
            if ( (*sibling)->name == name_ )  return true;
         }
      } 
   }

   if ( !localOnly )  {
      // --- look for symbol in upper levels
      if ( parent )  return (parent->getNItems (name_) > 0);
   }

   return false;
}

bool SxSymbolTable::containsGroup (const SxString &name_) const
{
   // --- check node
   if (table.containsKey (name_))  {
      return true;
   }

   // --- look in local symbol level (sibling by sibling)
   SxList<SxSymbolTable *>::ConstIterator sibling;
   for (sibling  = children.begin(); 
        sibling != children.end(); 
        sibling++)  
   {
      if ( (*sibling)->name == name_ )  return true;
   }

   return false;
}





SxSymbolTable *
SxSymbolTable::getGroup (const SxString &name_, bool recursive) const
{
   SxSymbolTable *res = NULL;
   SxList<SxSymbolTable *>::ConstIterator it;
   for (it=children.begin(); it!=children.end();it++)   {
      if ( (*it)->name == name_)  { 
         res = *it; break; 
      }
      if (recursive)  {
         res = (*it)->getGroup (name_, recursive);
         if (res)  break;
      }
   }

   if (!res)  {
      SX_THROW ("Group '" + name_ + "' is missing.");
   }
   
   return res;
}


SxSymbolTable *SxSymbolTable::begin () const
{
   if (children.getSize() == 0)  return NULL;
   return *children.begin();
}



SxSymbolTable *SxSymbolTable::nextSibling (const SxString &name_) const
{
   SxSymbolTable *res = NULL;
   if (parent)  {
      bool found = false;
      SxList<SxSymbolTable *>::ConstIterator it;
      for (it  = parent->children.begin(); 
           it != parent->children.end();
           it++)
      {
         if ( found && (*it)->name == name_)  { 
            res = *it; break; 
         }
         if ( *it == this )  found = true;
      }
   }
   return res;
}


SxSymbolTable *SxSymbolTable::nextSibling () const
{
   SxSymbolTable *res = NULL;
   if (parent)  {
      bool found = false;
      SxList<SxSymbolTable *>::ConstIterator it;
      for (it  = parent->children.begin(); 
           it != parent->children.end();
           it++)
      {
         if ( found )  { 
            res = *it; break; 
         }
         if ( *it == this )  found = true;
      }
   }
   return res;
}




SxSymbol *SxSymbolTable::append (const SxString &name_)
{
   //SX_DBG_MSG ("name '" << name_ << "'");
   return &(table(name_) = SxSymbol (name_, 0.));
}


SxSymbol *SxSymbolTable::append (const SxString &var, const SxString &val)
{
   //SX_DBG_MSG ("var '" << var << "'='" << val <<"'");
   // --- determine group
   SxSymbolTable *group = this;
   SxList<SxString>::ConstIterator partIt;
   SxList<SxString> parts = var.tokenize ('.');
   ssize_t i, groupLvl = parts.getSize() - 1;
   SxString symName;
   if (groupLvl > 0)  {
      for (partIt = parts.begin(), i=0; i < groupLvl; ++partIt, ++i)  {
         cout << "search for group " << *partIt << endl;
         group   = group->getGroup (*partIt);
      }
      symName = *parts.fromLast();
   }  else  {
      symName = var;
   }

   SxSymbol &sym = group->table(symName);
   sym = SxSymbol (symName, val, SxSymbol::DoAutoCast);
   return &sym;
}


void SxSymbolTable::append (const SxMap<SxString,SxString> &vars)
{
   SxMap<SxString,SxString>::ConstIterator var;
   for (var = vars.begin(); var != vars.end(); ++var)  {
      append (var.getKey(), var.getValue());
   }
}


void SxSymbolTable::remove (SxSymbol *child)
{
   if (child)  {
      //SX_DBG_MSG ("remove " << child->name);
      table.removeKey (child->name);
   }
//   SxList<SxSymbol>::ConstIterator it;
//   int i=0;
//   for (it = table.begin(); it != table.end(); it++, i++)  {
//      if ( &(*it) == child )  {
//         table.remove(i);
//         break;
//      }
//   }
}


void SxSymbolTable::addChild (SxSymbolTable *child)
{
   children.append (child);
}


SxSymbolTable *SxSymbolTable::topLevel ()
{
   if (parent)  return parent->topLevel ();
   else         return this;
}

const SxSymbolTable *SxSymbolTable::topLevel () const
{
   if (parent)  return parent->topLevel ();
   else         return this;
}



SxSymbol *SxSymbolTable::pushList ()
{
   SxSymbol *sym = new SxSymbol ("", 0.);
   sym->type = LIST;
   tmpList.prepend (sym);
   return sym;
}


SxSymbol *SxSymbolTable::popList ()
{
   SxSymbol *sym = NULL;
   if ( tmpList.getSize() )  {
      sym = tmpList.last();
      tmpList.removeFirst ();
   }
   return sym;
}


void SxSymbolTable::validateTable (const SxSymbolTable &validator)
{
   SX_TRACE ();
   ssize_t validateTl = validate (validator);
   SX_DBG_MSG (validateTl << " validated symbols");
}


ssize_t SxSymbolTable::validate (const SxSymbolTable &validator)
{
   //SX_DBG_MSG ("validate table '" << name << "'");

   bool topLevelDefs = false;
   bool skipUnexpGroups = false;
   if (validator.contains("topLevelDefs"))  
      topLevelDefs = validator.get("topLevelDefs")->toAttribute();
   if (validator.contains("skipUnexpectedGroups"))  
      skipUnexpGroups = validator.get("skipUnexpectedGroups")->toAttribute();

   SxList<SxSymbolTable *>::ConstIterator valLvl, lvl, cntLvl;
   bool found;
   ssize_t nValidatedSymbols = 0;
   SxString valName, symName;
   SxArray<bool> foundItem;
   SxSymbol *sym=NULL, *symVal=NULL, *needs=NULL;
   SxSymbol *symMin=NULL, *symMax=NULL, *symDim=NULL, *notSym=NULL;;
   SxSymbol *minSym=NULL, *maxSym=NULL, *nSym=NULL, *xorSym=NULL;
   double min, max;
   SxSymbolTable *valTab;

   // look in current level of symbol table
   ssize_t n, v, nVal = validator.children.getSize();
   foundItem.resize (nVal); 
   if (nVal > 0)  {
      foundItem.set (false);
   }

   for (valLvl  = validator.children.begin(), v=0;
        valLvl != validator.children.end();
        ++valLvl, v++)  
   {
         valTab = *valLvl;
         //SX_DBG_MSG ("validate " << valTab->name << " " << validateTl);
         nValidatedSymbols++;
         
         const SxString &type = valTab->get("type", true)->toString ();
         if (   type == "group"
             || valTab->get("optional", true, false)
             || valTab->get("deprecate", true, false)) 
         {
            foundItem(v)=true;  // just skip error prompt
         }
         
         sym = get(valTab->name, false, false);
         if (sym)  {
            foundItem(v)=true;
            
            if (valTab->get("deprecate", true, false))  {
//             cout << SX_SEPARATOR;
//             sxprintf ("| WARNING: Old file format recognized.\n");
//             sxprintf ("|          Item '%s' is deprecate.\n", 
//                     sym->name.ascii());
//             cout << SX_SEPARATOR;
               topLevel()->deprecateFormat = true;
            }

            // --- check for 'needs' conditions for items (see also groups)
            needs = valTab->get("needs", true, false);
            if (needs)  {
               SxString needStr = needs->toString();
               SxList<SxString> needList = needStr.tokenize(',');
               ssize_t needIdx, nNeeds = needList.getSize();
               found = true;
               for (needIdx=0; needIdx < nNeeds; needIdx++)  {
                  if (   !contains      ( needList(needIdx) )
                      && !containsGroup ( needList(needIdx) ) )
                     found = false;
               }
               if ( !found )  {
                  validationError ("Presence of item '"+sym->name
                                   +"' requires also '" +needStr+"'", sym);
               }
            }
            // --- check for 'not' conditions for items (s.a. groups)
            notSym = valTab->get("not", true, false);
            if (notSym)  {
               SxString notStr = notSym->toString();
               SxList<SxString> notList = notStr.tokenize(',');
               ssize_t notIdx, nNot = notList.getSize();
               found = false;
               for (notIdx=0; notIdx < nNot; notIdx++)  {
                  if (   contains      (notList(notIdx))
                      || containsGroup (notList(notIdx)))  {
                     found = true;
                     break;
                  }
               }
               if ( found )  {
                  validationError ("Presence of item '"+sym->name
                                   +"' excludes '"+notStr+"'", sym);
               }
            }
            // --- STRING
            if (type=="string")  {
               if (sym->type != STR)  {
                  validationError ( "Value of item '"+sym->name
                                   +"' is not a string", sym);
               }
            }
            // --- COMBO STRING
            else if (type=="combo")  {
               if (sym->str.getSize() == 0)  {
                  validationError ( "Value of item '"+sym->name
                                   +"' is not a string", sym);
               }
               SxList<SxString> strList;
               SxList<SxString>::ConstIterator strIt;
               strList = valTab->get("val",true,false)->str.tokenize(',');
               found = false;
               for (strIt=strList.begin(); strIt!=strList.end(); strIt++)
                  if ( *strIt == sym->str )  { found=true; break; }
               
               if ( !found )  {
                  validationError ( "Value of '"+sym->name+"' is invalid ('"
                                   +sym->str+"').\n"
                                    "    Choose from the following list:\n    "
                                   + valTab->get("val", true)->str, sym);
               }
            }
            // --- ENUMERATION
            else if (type=="enum")  {
               if (sym->str.getSize() == 0)  {
                  validationError ( "Value of item '"+sym->name
                                   +"' is not a string", sym);
               }
               SxList<SxString> strList, symStrList;
               SxList<SxString>::ConstIterator strIt, symStrIt;
               strList = valTab->get("val",true,false)->str.tokenize(',');
               symStrList = sym->str.tokenize(',');
               for (symStrIt  = symStrList.begin();
                    symStrIt != symStrList.end();
                    symStrIt++)
               {
                  found = false;
                  for (strIt=strList.begin(); strIt!=strList.end(); strIt++)
                     if ( *strIt == *symStrIt )  { found=true; break; }

                  if ( !found )
                     validationError ( "Value of '"+sym->name+"' is invalid ('"
                                      +sym->str+"').\n    Choose from the "
                                       "following list:\n    "
                                      +valTab->get("val",true,false)->str,sym);
               } 
            }
            // --- INT
            else if (type=="int")  {
               if (   sym->type != VAR 
                   || fabs(sym->val - ceil(sym->val)) > 1e-50)
               {
                  validationError ("Item '"+sym->name+"' should be an integer",
                                   sym);
               }
               // --- check range
               symMin = valTab->get("min", true,false);
               symMax = valTab->get("max", true,false);
               min = (symMin) ? symMin->val : -1e+50;
               max = (symMax) ? symMax->val :  1e+50;
               if (sym->val < min || sym->val > max)  {
                  validationError ("Item '"+sym->name+"' is out of range "
                                   "("+min+","+max+")", sym);
               }
            }
            // --- REAL
            else if (type=="real")  {
               if ( sym->type != VAR )
                  validationError ("Item '"+sym->name+"' should be a real",sym);
               // --- check range
               symMin = valTab->get("min", true,false);
               symMax = valTab->get("max", true,false);
               min = (symMin) ? symMin->val : -1e+50;
               max = (symMax) ? symMax->val :  1e+50;
               if (sym->val < min || sym->val > max)  {
                  validationError ("Item '"+sym->name+"' is out of range "
                                   "("+min+","+max+")", sym);
               }
            }
            // --- FLAG
            else if (type=="flag")  {
               // nothing to be tested
            }
            // --- LIST
            else if (type=="list")  {
               if ( sym->type != LIST && sym->type != VAR)
                  validationError ("'"+sym->name+"' is not list", sym);
               // --- check rank
               symMin = valTab->get("minRank", true,false);
               symMax = valTab->get("maxRank", true,false);
               if (symMin)
                  if (sym->getRank()  < symMin->val)
                     validationError ("'"+sym->name+"' has too low rank",sym);
               if (symMax)
                  if (sym->getRank()  > symMax->val)
                     validationError ("'"+sym->name+"' has too high rank",sym);
            }
            // --- VECTOR
            else if (type=="vector")  {
               if ( sym->type != LIST && sym->type != VAR)
                  validationError ("'"+sym->name+"' is not vector",sym);
               // --- check rank
               if (sym->getRank() > 1)
                  validationError ("'"+sym->name+"' should be a vector "
                                   "but has a rank of "+sym->getRank(),sym);
               // --- check dimensions
               ssize_t nElem = sym->valList->getSize();
               symMin = valTab->get("minDim", true,false);
               symMax = valTab->get("maxDim", true,false);
               symDim = valTab->get("dim",    true,false);
               if (symMin)
                  if (nElem  < (ssize_t)symMin->val)
                     validationError ("'"+sym->name+"' has too few elements",
                                      sym);
               if (symMax)
                  if (nElem  > (ssize_t)symMax->val)
                     validationError ("'"+sym->name+"' has too many elements",
                                      sym);
               if (symDim)
                  if (nElem != (ssize_t)symDim->val)
                     validationError ("'"+sym->name+"' has wrong number "
                                      "of elements", sym);
            }
            // --- MATRIX
            else if (type=="matrix")  {
               if ( sym->type != LIST )
                  validationError ("'"+sym->name+"' is not matrix", sym);
               // --- check rank
               if (sym->getRank() != 2)
                  validationError ("'"+sym->name+"' should be a matrix "
                                   "but has a rank of "+sym->getRank(), sym);
               // --- check number of rows
               symDim = valTab->get("dims", true,false);
               if (!symDim || symDim->valList->getSize () != 2)  {
                  // --- corrupted std file
                  SxSymbolTable *valParent = valTab->parent;
                  SxString id = valTab->name;
                  while (valParent)  {
                     if (valParent->parent) id = valParent->name + "." + id;
                     valParent = valParent->parent;
                  }
                  cout << "Error in std file for matrix " << id
                       << ": missing dims = [N,M]" << endl;
                  SX_EXIT;
               }

               ssize_t nRows = sym->valList->getSize();
               if (nRows != (ssize_t)(*symDim->valList)(0).val )
                  validationError ( "'"+sym->name+"' has "+nRows+" rows. "
                                   +(*symDim->valList)(0).val+" are expected",
                                   sym);
               // --- check number of columns
               ssize_t iRow, nCols;
               for (iRow=0; iRow < nRows; iRow++)  {
                  nCols = (*sym->valList)(iRow).valList->getSize();
                  if (nCols != (ssize_t)(*symDim->valList)(1).val)  {
                     validationError ( "'"+sym->name+"' has "+nCols+" columns. "
                                      +(*symDim->valList)(1).val
                                      +" are expected", sym);
                  }
               }
            }
         }
   }


   // --- are there locally defined variables?
   if (topLevelDefs && level > 0)  {
      SxMap<SxString, SxSymbol>::Iterator tableIt = table.begin ();
      for (; tableIt != table.end (); ++tableIt)  {
         found = false;
         for (valLvl  = validator.children.begin(), v=0;
              valLvl != validator.children.end();
              ++valLvl, v++)  
         {
            valTab = *valLvl;
            // skip groups and optional items
            if (valTab->name == tableIt.getKey ()) {
               found = true;
               break;
            }
         }
         if (!found)
            validationError ("Unexpected item '"+tableIt.getKey()+"'", this);
      }
   }

   // --- check missing items
   for (v=0; v < nVal; v++)  {
      if ( !foundItem(v) )  {
         xorSym = validator.children(v)->get("xor", true, false);
         if (xorSym)  {
            SxString xorString = xorSym->toString();
            SxList<SxString> xorList = xorString.tokenize(',');
            ssize_t xorIdx, nXorIdx = xorList.getSize();
            int nXorFound = 0;
            for (xorIdx=0; xorIdx < nXorIdx; xorIdx++)  {
               if (   containsGroup(xorList(xorIdx)) 
                   && xorList(xorIdx) != validator.children(v)->name ) 
                  nXorFound++;
            }
            if ( nXorFound == 0)  {
               validationError ("Missing one of these items: '"
                                +xorString+"'", xorSym);
            }
         } else 
            validationError ("Missing item '"+validator.children(v)->name+"'",
                             this);
      }
   }

   // --- check for unexpected groups
   if ( !skipUnexpGroups )  {
      for (lvl = children.begin(); lvl != children.end(); ++lvl) {
         if ( !(validator.containsGroup((*lvl)->name)) )  {
            validationError ("Unexpected group '"+(*lvl)->name+"'", (*lvl));
         }
      }
   }

   // --- check for missing groups
   for (valLvl  = validator.children.begin();
        valLvl != validator.children.end();
        ++valLvl)
   {
      if ( (symVal = (*valLvl)->get ("type", true,false)) )  {
         if (    symVal->str == "group" 
             && !(    (*valLvl)->get("optional", true,false) 
                   || (*valLvl)->get("deprecate", true, false)) )  
         {
            xorSym = (*valLvl)->get ("xor", true, false);
            if ( !(containsGroup((*valLvl)->name)) )  {
               if (xorSym)  {
                  SxString xorString = xorSym->toString();
                  SxList<SxString> xorList = xorString.tokenize(',');
                  ssize_t xorIdx, nXorIdx = xorList.getSize();
                  int nXorFound = 0;
                  for (xorIdx=0; xorIdx < nXorIdx; xorIdx++)  {
                     if (   containsGroup(xorList(xorIdx)) 
                         && xorList(xorIdx) != (*valLvl)->name )  
                        nXorFound++;
                  }
                  if ( nXorFound == 0)  {
                     validationError ("Missing one of these groups: "
                                      +xorString, xorSym);
                  }
               }  else  {
                  validationError ("Missing group '"+(*valLvl)->name+"'",
                                   this);
               }
            }
         }
      }
   }

   // --- validate all groups
   for (valLvl  = validator.children.begin();
        valLvl != validator.children.end();
        ++valLvl)
   {
      //SX_DBG_MSG ("validator '" << (*valLvl)->name << "'");
      symVal = (*valLvl)->get ("type", true,false);
      if (symVal->str == "group")  {
         // --- loop over children
         n = 0;
         for (lvl = children.begin(); lvl != children.end(); ++lvl)  {
            if ( (*lvl)->name == (*valLvl)->name )  {
               //SX_DBG_MSG ("children '" << (*lvl)->name << "'");
               nValidatedSymbols += (*lvl)->validate ( *(*valLvl) );
               n++;
            }
         }
         //SX_DBG_MSG (n << " children");
         if (n == 0 && (*valLvl)->get ("optional", true,false))  {
            continue;
         }
         
         // --- how many items of that group do we need?
         minSym   = (*valLvl)->get ("minItems", true,false);
         maxSym   = (*valLvl)->get ("maxItems", true,false);
         nSym     = (*valLvl)->get ("nItems",   true,false);

         if (maxSym && n > (int)maxSym->val)
            validationError ("Too many groups '"+(*valLvl)->name+"' defined",
                             this);
         if (nSym   && n > (int)nSym->val)
            validationError ("Too many groups '"+(*valLvl)->name+"' defined",
                             this);

         // --- check for 'needs' conditions for groups 
         //     (see also items)
         needs = (*valLvl)->get("needs", true, false);
         if (needs && n > 0)  {
            SxString needStr = needs->toString();
            SxList<SxString> needList = needStr.tokenize(',');
            ssize_t needIdx, nNeeds = needList.getSize();
            found = true;
            for (needIdx=0; needIdx < nNeeds; needIdx++)  {
               if (     !contains      ( needList(needIdx) )
                     && !containsGroup ( needList(needIdx) ) )
                  found = false;
            }
            if ( !found )  {
               validationError ("Presence of group '"+(*valLvl)->name
                     +"' requires also '" +needStr+"'", this);
            }
         }
         
         // --- check for 'not' conditions for items (s.a. groups)
         notSym = (*valLvl)->get("not", true, false);
         if (notSym && n > 0)  {
            SxString notStr = notSym->toString();
            SxList<SxString> notList = notStr.tokenize(',');
            ssize_t notIdx, nNot = notList.getSize();
            found = false;
            for (notIdx=0; notIdx < nNot; notIdx++)  {
               if (   contains      (notList(notIdx))
                   || containsGroup (notList(notIdx)))  {
                  found = true;
                  break;
               }
            }
            if ( found )  {
               validationError ("Presence of group '"+(*valLvl)->name
                                +"' excludes '"+notStr+"'", this);
            }
         }
         
         // --- has the missing group some XOR references?
         xorSym = (*valLvl)->get ("xor", true, false);
         if (xorSym)  {
            SxString xorString = xorSym->toString();
            SxList<SxString> xorList = xorString.tokenize(',');
            ssize_t xorIdx, nXorIdx = xorList.getSize();
            int nXorFound = 0;
            for (xorIdx=0; xorIdx < nXorIdx; xorIdx++)  {
               if (   containsGroup(xorList(xorIdx)) 
                   && xorList(xorIdx) != (*valLvl)->name) 
                  nXorFound++;
            }
            if ( n == 0 && nXorFound == 0)  {
               validationError ("Missing one of these groups: "+
                                xorString, this);
            }
            if ( n > 0 && nXorFound >= 1)  {
               validationError (
                  "Only one of these groups are allowed: "+xorString, this
               );
            }
         }  else  {
            // --- group not found
            if (n == 0 && !minSym && !maxSym && !nSym 
                  && !(*valLvl)->get("optional", true,false) )  {
               validationError (
                     "Missing group '" +(*valLvl)->name+"'", this
               );
            }
            if (nSym && n != (int)nSym->val)
                  validationError (
                     "Too few groups '"+(*valLvl)->name+"' defined", this
                  );
                  
            if (minSym && n < (int)minSym->val)
                  validationError (
                     "Too few groups '"+(*valLvl)->name+"' found", this
                  );
         }
      }
   }
   
   return nValidatedSymbols;
}


void SxSymbolTable::printList () const
{
   for (int i=0; i < tmpList.getSize(); i++)
      sxprintf ("%d: %p\n", i, (void *)(tmpList(i)));
}


void SxSymbolTable::print () const
{
   cout << "Level: " << level << ": " << name << endl;
   SxMap<SxString, SxSymbol>::ConstIterator it;
   for (it=table.begin(); it!=table.end(); ++it)  {

      if ( it.getValue().type == VAR ||
           it.getValue().type == LIST ||
           it.getValue().type == STR )
      {
         cout << " - variable " << it.getValue().name << " = "
                                << it.getValue() << '\n';
      }
      
      if ( it.getValue().type == FUNC )
         cout << " - function " << it.getValue().name << '\n';
      if ( it.getValue().type == STRFUNC )
         cout << " - strfunc " << it.getValue().name << '\n';
   }
   cout << endl;

   // --- print children
   SxList<SxSymbolTable *>::ConstIterator lvl;
   for (lvl=children.begin(); lvl!=children.end(); lvl++)  
      (*lvl)->print ();
}

void SxSymbolTable::printHash (const SxString &indent1,
                               const SxString &indent2,
                               int idx) const
{
   bool first = true;
   SxMap<SxString, SxSymbol>::ConstIterator it;

   if (name == "")  {
      cout << indent1 << "(\n";         // top level '('
   }  else  {
      cout << indent1                   // otherwise 'name => ('
           << "\"" << name << "~" << idx << "\"" 
           << " => {\n";
   }
   for (it = table.begin(); it != table.end(); ++it)  {
      if ( it.getValue().type == VAR ||
           it.getValue().type == LIST ||
           it.getValue().type == STR )
      {
         if (!first)  cout << ",\n";
         cout << indent1 << indent2 
              << "\"" << it.getValue().name << "\""
              << " => ";
         if (it.getValue().type == STR)
            cout << "q(" 
                 << it.getValue().str.substitute("(","\\(")
                                     .substitute(")","\\)")
                 << ')';
         else
            cout << it.getValue();
         first = false;
      }
   }

   // --- print children
   int grpIdx;
   SxList<SxSymbolTable *>::ConstIterator lvl;
   for (lvl=children.begin(), grpIdx=0; lvl!=children.end(); lvl++, grpIdx++) {
      if (!first)  cout << ",\n";
      (*lvl)->printHash (indent1 + indent2, indent2, grpIdx);
      first = false;
   }
   cout << endl;

   if (level == 0)  cout << ");\n";
   else             cout << indent1 << "}";
}


void SxSymbolTable::validationError (const SxString &str,
                                     const SxSymbol *symbol_)
{
   sxprintf ("\n");
   if (symbol_ && symbol_->parserFilename != "")  {
      std::cout << "    file " << symbol_->parserFilename
                << " at line " << symbol_->parserLineNumber << std::endl;
   }
   
   sxprintf ("    Validation error in ");
   if (level == 0)
      sxprintf ("root node!\n");
   else
      sxprintf ("group '%s'!\n", getPathName().ascii());
   sxprintf (">>> %s.\n", str.ascii());
   if (useDeprecateFormat())  {
      sxprintf ("    Perhaps new and deprecate input tags are used.\n");
   }
   SX_EXIT;
}

void SxSymbolTable::validationError (const SxString      &str,
                                     const SxSymbolTable *table_)
{
   sxprintf ("\n");
   if (table_ && table_->parserFilename != "")  {
      std::cout << "    file " << table_->parserFilename
                << " at line " << table_->parserLineNumber << std::endl;
   }
   
   sxprintf ("    Validation error in ");
   if (level == 0)
      sxprintf ("root node!\n");
   else
      sxprintf ("group '%s'!\n", getPathName().ascii());
   sxprintf (">>> %s.\n", str.ascii());
   if (useDeprecateFormat())  {
      sxprintf ("    Perhaps new and deprecate input tags are used.\n");
   }
   SX_EXIT;
}


SxString SxSymbolTable::getPathName () const
{
   if (parent && level > 1)  return parent->getPathName () + "." + name;
   else         return name;
}

