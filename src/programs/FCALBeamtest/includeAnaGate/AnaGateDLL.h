// ----------------------------------------------------------------------------
// Projekt       : Analytica AnaGate API
// File          : AnaGateDLL.h
// Author        : Axel Schmidt
// Copyright     : (C) 2004 by Analytica GmbH
// ----------------------------------------------------------------------------
// $Id: AnaGateDLL.h,v 1.1.1.1 2010/07/16 18:12:17 ctarbert Exp $
//
// $Log: AnaGateDLL.h,v $
// Revision 1.1.1.1  2010/07/16 18:12:17  ctarbert
//
//
// Revision 1.10  2009/04/30 14:20:07  stefanwelisch
// Define für LUA_LIB hinzugefügt
//
// Revision 1.9  2007/11/21 09:19:22  axelschmidt
// real c interface for DLL
//
// Revision 1.8  2007/11/08 06:55:10  axelschmidt
// Linux support
//
// Revision 1.7  2007/10/29 08:53:11  axelschmidt
// new define for VB6 support (special DLL)
//
// Revision 1.6  2005/07/22 07:12:54  AxelSchmidt
// DeviceInfo function implemented
//
// Revision 1.5  2004/08/13 12:52:33  AxelSchmidt
// import/export-defines changed
//
// Revision 1.4  2004/08/09 15:49:00  AxelSchmidt
// DLLVersion in DLLInfo umbenannt
//
// Revision 1.3  2004/07/23 15:51:32  AxelSchmidt
// neue Methode DLLVersion
//
// Revision 1.2  2004/07/21 16:12:53  AxelSchmidt
// global OpenDevice and CloseDevice removed
//
// Revision 1.1  2004/07/21 06:53:34  AxelSchmidt
// initial version
//
// ----------------------------------------------------------------------------

#ifndef _ANAGATE_DLL_H
#define _ANAGATE_DLL_H

// Includes -------------------------------------------------------------------
#include "AnaGate.h"

// Defines --------------------------------------------------------------------

// Folgender ifdef-Block ist die Standardmethode zum Erstellen von Makros, die das Exportieren
// aus einer DLL vereinfachen. Alle Dateien in der DLL werden mit dem ANAI2CDLL_EXPORTS-Symbol
// kompiliert, das in der Befehlszeile definiert wurde. Das Symbol darf nicht für ein Projekt definiert werden,
// das diese DLL verwendet. Alle anderen Projekte, deren Quelldateien diese Datei beinhalten, erkennen
// ANAI2CDLL_API-Funktionen als aus einer DLL importiert, während die DLL mit diesem Makro
// definierte Symbole als exportiert ansieht.
#ifndef LUA_LIB
   #ifdef WIN32
      #ifdef ANAGATEDLL_EXPORTS
         #define ANAGATEDLL_API __declspec(dllexport)
      #else
         #define ANAGATEDLL_API __declspec(dllimport)
      #endif
   #else
      #define ANAGATEDLL_API
   #endif
#else
   #define ANAGATEDLL_API
#endif

#ifndef LUA_LIB
   #ifdef WIN32
      // VB6 do not support __cdecl calling convention which is the default c calling convention
      // define the ANAGATEDLL_VB6 if you need DLL for VB6 (then the support __stcall conventions is used )
      #ifdef ANAGATEDLL_VB6
         #define ANAGATEDLL_CALLING __stdcall
      #else
         #define ANAGATEDLL_CALLING  __cdecl
      #endif
   #else
      #define ANAGATEDLL_CALLING
   #endif
#else
   #define ANAGATEDLL_CALLING
#endif

#ifndef WIN32
   // make sure some windows defines are present on non windows platforms
   #define WINAPI
   #define APIENTRY
   #define BOOL   char
   #define TRUE   1
   #define FALSE  0
#endif

// Prototyping ----------------------------------------------------------------
#ifdef __cplusplus
extern "C"
{
#endif // __cplusplus

#ifndef LUA_LIB
   /** Retrieves a textual version description of the AnaGate (dynamic link) library.
      @param pcMessage Pointer to a c-style character buffer, in which the version string is stored.
      @param nMessageLen Length of the supplied charcter buffer. If the version string does not fit into the
             buffer, the string is shortened.
      @return The byte of the returned error string.
   */
   ANAGATEDLL_API int ANAGATEDLL_CALLING DLLInfo( char *pcMessage, int nMessageLen );
#endif

   /** Retrieves version informations from the AnaGate device (includes versions
      of hardware and software).
      @param hHandle Device handle (from a successfull #OpenDevice call).
      @param pDeviceInfo Pointer to the SAnaGateDeviceInfo structure.
      @return. If the function succeeds, the return value is zero. If the function fails, the return value is non zero.
   */
   ANAGATEDLL_API int ANAGATEDLL_CALLING DeviceInfo( int hHandle, struct SAnaGateDeviceInfo *pDeviceInfo );

#ifdef __cplusplus
}
#endif //__cplusplus

#endif // #ifndef _ANAGATE_DLL_H
