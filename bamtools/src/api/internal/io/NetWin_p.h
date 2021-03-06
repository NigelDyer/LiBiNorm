// ***************************************************************************
// NetWin_p.h (c) 2011 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 8 December 2011 (DB)
//                24 July 2017 (ND)
// ---------------------------------------------------------------------------
// Provides common networking-related includes, etc. for Windows systems
//
// Note: requires Windows XP or later
// ***************************************************************************
// Modified for LiBiNorm to allow build using Visual C++ 2015
// Nigel Dyer 24 July 2017
// ***************************************************************************
#ifndef NETWIN_P_H
#define NETWIN_P_H

//  -------------
//  W A R N I N G
//  -------------
//
// This file is not part of the BamTools API.  It exists purely as an
// implementation detail. This header file may change from version to version
// without notice, or even be removed.
//
// We mean it.

#ifdef _WIN32 // <-- source files only include the proper Net*_p.h, but this is a double-check

#include <winsock2.h>  // <-- should bring 'windows.h' along with it
#include <Ws2tcpip.h>

//	The following definitions were added to allow building under Visual C++ for LiBiNorm
#include <io.h>

#define close _close
#define read _read
#define write _write
#define ssize_t int
#define ioctl(a,b,c) ioctlsocket(a,b,(u_long*)c)

/* ipc/network software -- operational errors */
#define	ESHUTDOWN	58		/* Can't send after socket shutdown */
#define	ETOOMANYREFS	59		/* Too many references: can't splice */

#define	EHOSTDOWN	64		/* Host is down */
//	End of Visual C++ additions

#ifndef   BT_SOCKLEN_T
#  define BT_SOCKLEN_T int
#endif

#ifdef _MSC_VER
#  pragma comment(lib, "ws2_32.lib")
#endif

namespace BamTools {
namespace Internal {

// use RAII to ensure WSA is initialized
class WindowsSockInit {
    public:
        WindowsSockInit(void) {
            WSAData wsadata;
            WSAStartup(MAKEWORD(2,2), &wsadata); // catch error ?
        }

        ~WindowsSockInit(void) {
            WSACleanup();
        }
};

} // namespace Internal
} // namespace BamTools

#endif // _WIN32

#endif // NETWIN_P_H

