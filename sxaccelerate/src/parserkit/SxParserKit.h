#ifndef _SX_PARSER_KIT_H_
#define _SX_PARSER_KIT_H_

#ifdef WIN32
#  if defined(_EXPORT_sxparserkit)
#     define SX_EXPORT_PARSER_KIT __declspec(dllexport)
#  else
#     define SX_EXPORT_PARSER_KIT __declspec(dllimport)
#  endif
#else
#  define SX_EXPORT_PARSER_KIT
#endif

#endif /* _SX_PARSER_KIT_H_ */
