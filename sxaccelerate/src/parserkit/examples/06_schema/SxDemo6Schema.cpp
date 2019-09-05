#include <SxDemo6Schema.h>

SxDemo6Schema::SxDemo6Schema () { }

SxDemo6Schema::~SxDemo6Schema () { }

ssize_t SxDemo6Schema::validate (const SxGraph<SxDemo6AstNode> &dataG) const
{
   auto dIt = dataG.begin (0);

   // skip dummy root
   ++dIt;

   for (; dIt.isValid (); ++dIt) {
      if (dIt->type  == Type::StringVal) {
         if (dIt->data.getString ().getSize () > 10)
            return false;
      }
   }

   return true;
}

