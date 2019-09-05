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
#ifndef _SX_TYPE_CASTER_H_
#define _SX_TYPE_CASTER_H_


//------------------------------------------------------------------------------
// Shortcuts for mappers of types supported by BLAS and CLAPACK.
//------------------------------------------------------------------------------
typedef SxTypeMapper<int, int, int>                     Int;
typedef SxTypeMapper<long int, long int, long int>      Lint;
typedef SxTypeMapper<float, float>                      Float;
typedef SxTypeMapper<double, double>                    Double;
typedef SxTypeMapper<SxComplex8, float,  SxComplex8>    Complex8;
typedef SxTypeMapper<SxComplex16, double, SxComplex16>  Complex16;


//------------------------------------------------------------------------------
// Ranking of datatypes used for typecasting.
//------------------------------------------------------------------------------
// Anonymous 'enum' ranking of template classes as used in generic template
// programming is not supported by all compilers 
// (e.g. HP aC++ treats only one iteration of template class specialization). 
// So what follows is the 'hardcoded' ranking, that should work everywhere!
//------------------------------------------------------------------------------
template<class A, class B> struct SxTypeCaster { 
//#  ifdef WIN32
//      /* Provide missing datatypes to help Developers Studio .NET compiler */
//      typedef Int TRes; typedef Int::Type Type;
//#  endif /* WIN32 */
};

template<> struct SxTypeCaster<Int,      Int>       {typedef Int             TRes;
                                                     typedef Int::Type       Type;};
template<> struct SxTypeCaster<Int,      Lint>      {typedef Lint            TRes;
                                                     typedef Lint::Type      Type;};
template<> struct SxTypeCaster<Int,      Float>     {typedef Float           TRes;
                                                     typedef Float::Type     Type;};
template<> struct SxTypeCaster<Int,      Double>    {typedef Double          TRes;
                                                     typedef Double::Type    Type;};
template<> struct SxTypeCaster<Int,      Complex8>  {typedef Complex8        TRes;
                                                     typedef Complex8::Type  Type;};
template<> struct SxTypeCaster<Int,      Complex16> {typedef Complex16       TRes;
                                                     typedef Complex16::Type Type;};

template<> struct SxTypeCaster<Lint,     Int>       {typedef Lint            TRes;
                                                     typedef Lint::Type      Type; };
template<> struct SxTypeCaster<Lint,     Lint>      {typedef Lint            TRes;
                                                     typedef Lint::Type      Type;};
template<> struct SxTypeCaster<Lint,     Float>     {typedef Float           TRes;
                                                     typedef Float::Type     Type; };
template<> struct SxTypeCaster<Lint,     Double>    {typedef Double          TRes;
                                                     typedef Double::Type    Type; };
template<> struct SxTypeCaster<Lint,     Complex8>  {typedef Complex8        TRes;
                                                     typedef Complex8::Type  Type; };
template<> struct SxTypeCaster<Lint,     Complex16> {typedef Complex16       TRes;
                                                     typedef Complex16::Type Type; };

template<> struct SxTypeCaster<Float,    Int>       {typedef Float           TRes;
                                                     typedef Float::Type     Type; };
template<> struct SxTypeCaster<Float,    Lint>      {typedef Float           TRes;
                                                     typedef Float::Type     Type; };
template<> struct SxTypeCaster<Float,    Float>     {typedef Float           TRes;
                                                     typedef Float::Type     Type;};
template<> struct SxTypeCaster<Float,    Double>    {typedef Double          TRes;
                                                     typedef Double::Type    Type;};
template<> struct SxTypeCaster<Float,    Complex8>  {typedef Complex8        TRes;
                                                     typedef Complex8::Type  Type; };
template<> struct SxTypeCaster<Float,    Complex16> {typedef Complex16       TRes;
                                                     typedef Complex16::Type Type; };

template<> struct SxTypeCaster<Double,   Int>       {typedef Double          TRes;
                                                     typedef Double::Type    Type;};
template<> struct SxTypeCaster<Double,   Lint>      {typedef Double          TRes;
                                                     typedef Double::Type    Type; };
template<> struct SxTypeCaster<Double,   Float>     {typedef Double          TRes;
                                                     typedef Double::Type    Type; };
template<> struct SxTypeCaster<Double,   Double>    {typedef Double          TRes;
                                                     typedef Double::Type    Type;};
template<> struct SxTypeCaster<Double,   Complex8>  {typedef Complex8        TRes;
                                                     typedef Complex8::Type  Type;};
template<> struct SxTypeCaster<Double,   Complex16> {typedef Complex16       TRes;
                                                     typedef Complex16::Type Type;};

template<> struct SxTypeCaster<Complex8, Int>       {typedef Complex8        TRes;
                                                     typedef Complex8::Type  Type;};
template<> struct SxTypeCaster<Complex8, Lint>      {typedef Complex8        TRes;
                                                     typedef Complex8::Type  Type;};
template<> struct SxTypeCaster<Complex8, Float>     {typedef Complex8        TRes;
                                                     typedef Complex8::Type  Type;};
template<> struct SxTypeCaster<Complex8, Double>    {typedef Complex8        TRes;
                                                     typedef Complex8::Type  Type;};
template<> struct SxTypeCaster<Complex8, Complex8>  {typedef Complex8        TRes;
                                                     typedef Complex8::Type  Type;};
template<> struct SxTypeCaster<Complex8, Complex16> {typedef Complex16       TRes;
                                                     typedef Complex16::Type Type;};

template<> struct SxTypeCaster<Complex16,Int>       {typedef Complex16       TRes;
                                                     typedef Complex16::Type Type;};
template<> struct SxTypeCaster<Complex16,Lint>      {typedef Complex16       TRes;
                                                     typedef Complex16::Type Type;};
template<> struct SxTypeCaster<Complex16,Float>     {typedef Complex16       TRes;
                                                     typedef Complex16::Type Type;};
template<> struct SxTypeCaster<Complex16,Double>    {typedef Complex16       TRes;
                                                     typedef Complex16::Type Type;};
template<> struct SxTypeCaster<Complex16,Complex8>  {typedef Complex16       TRes;
                                                     typedef Complex16::Type Type;};
template<> struct SxTypeCaster<Complex16,Complex16> {typedef Complex16       TRes;
                                                     typedef Complex16::Type Type;};

#endif /* _SX_TYPE_CASTER_H_ */
