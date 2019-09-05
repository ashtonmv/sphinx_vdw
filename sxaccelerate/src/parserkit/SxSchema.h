#ifndef _SX_SCHEMA_H_
#define _SX_SCHEMA_H_

#include <SxParserKit.h>
#include <SxGraph.h>
#include <SxGProps.h>
#include <SxVariant.h>
#include <limits>


class SX_EXPORT_PARSER_KIT SxSchema
{
   public:
      typedef typename SxVariantType::DataType Type;

      SxSchema ();
     ~SxSchema ();

      SxSchema (const SxPtr<SxGraph<SxGProps> > &schemaG_);

      bool validate (const SxPtr<SxGraph<SxGProps> > &dataG_);

   protected:

      template<class T>
      inline T max () {
         return std::numeric_limits<T>::max ();
      }

      template<class T>
      inline T min () {
         return std::numeric_limits<T>::min ();
      }

      // -- return the int value for 'key' from child nodes
      inline ssize_t getInt (const SxGraph<SxGProps>::Iterator &sIt,
                             const SxString &key,
                             const ssize_t &default_);

      // -- return the double value for 'key' from child nodes
      inline double getDouble (const SxGraph<SxGProps>::Iterator &sIt,
                               const SxString &key,
                               const double &default_);

      // -- return the string value for 'key' from child nodes
      inline SxString getString (const SxGraph<SxGProps>::Iterator &sIt,
                                 const SxString &key,
                                 const SxString &default_);

      void validationError (const SxString &msg,
                            const SxString &sTag,
                            const SxString &dTag);

      bool validateArray (const SxGraph<SxGProps>::Iterator &sIt,
                          const SxGraph<SxGProps>::Iterator &dIt);


      bool validateObj (const SxGraph<SxGProps>::Iterator &sIt,
                        const SxGraph<SxGProps>::Iterator &dIt);

      // -- find a node matching key,value in child nodes
      inline SxGraph<SxGProps>::Iterator get (const SxGraph<SxGProps>::Iterator &it,
                                              const SxString &key,
                                              const SxVariant &val = SxVariant());

      // -- optimize schema tree
      void simplifyInt    (SxGraph<SxGProps>::Iterator &sIt);
      void simplifyDouble (SxGraph<SxGProps>::Iterator &sIt);
      void simplifyString (SxGraph<SxGProps>::Iterator &sIt);
      void simplifyBool   (SxGraph<SxGProps>::Iterator &sIt);
      void simplifyArray  (SxGraph<SxGProps>::Iterator &sIt);
      void simplifyObj    (SxGraph<SxGProps>::Iterator &sIt);
      void simplifySchemaAst ();

      SxPtr<SxGraph<SxGProps> > schemaG;
      SxPtr<SxGraph<SxGProps> > dataG;
};

#endif /* _SX_SCHEMA_H_ */
