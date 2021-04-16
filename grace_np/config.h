/* config.h.  Generated from config.h.in by configure.  */
/*
 * Grace - GRaphing, Advanced Computation and Exploration of data
 * 
 * Home page: http://plasma-gate.weizmann.ac.il/Grace/
 * 
 * Copyright (c) 1991-1995 Paul J Turner, Portland, OR
 * Copyright (c) 1996-2000 Grace Development Team
 * 
 * Maintained by Evgeny Stambulchik
 * 
 * 
 *                           All Rights Reserved
 * 
 *    This program is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 * 
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 * 
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef __CONFIG_H
#define __CONFIG_H

/* Define if on AIX 3.
   System headers sometimes define this.
   We just want to avoid a redefinition error message.  */
#ifndef _ALL_SOURCE
#define _ALL_SOURCE 1
#endif

/* Define if you need to in order for stat and other things to work.  */
/* #undef _POSIX_SOURCE */

/* Define if you have the ANSI C header files.  */
#define STDC_HEADERS 1

/* Define if char is unsigned   */
/* #undef __CHAR_UNSIGNED__ */
/* Define sizeof(char)          */
#define SIZEOF_CHAR 1
/* Define sizeof(short)         */
#define SIZEOF_SHORT 2
/* Define sizeof(int)           */
#define SIZEOF_INT 4
/* Define sizeof(long)          */
#define SIZEOF_LONG 8
/* Define sizeof(long long)     */
#define SIZEOF_LONG_LONG 8
/* Define sizeof(float)         */
#define SIZEOF_FLOAT 4
/* Define sizeof(double)        */
#define SIZEOF_DOUBLE 8
/* Define sizeof(long double)   */
#define SIZEOF_LONG_DOUBLE 16
/* Define sizeof(void *)        */
#define SIZEOF_VOID_P 8

/* Define to empty if the keyword does not work.  */
/* #undef const */

/* Define to `int' if <sys/types.h> doesn't define.  */
/* #undef pid_t */

/* Define to `unsigned' if <sys/types.h> doesn't define.  */
/* #undef size_t */

/* Define if you have the <unistd.h> header file.  */
#define HAVE_UNISTD_H 1

/* Define to one of _getb67, GETB67, getb67 for Cray-2 and Cray-YMP systems.
   This function is required for alloca.c support on those systems.  */
/* #undef CRAY_STACKSEG_END */

/* Define if using alloca.c.  */
/* #undef C_ALLOCA */

/* Define if you have alloca, as a function or macro.  */
#define HAVE_ALLOCA 1

/* Define if you have <alloca.h> and it should be used (not on Ultrix).  */
#define HAVE_ALLOCA_H 1

/* Define as the return type of signal handlers (int or void).  */
#define RETSIGTYPE void

/* Define if you have <sys/wait.h> that is POSIX.1 compatible.  */
#define HAVE_SYS_WAIT_H 1

/* Define if you have the <fcntl.h> header file.  */
#define HAVE_FCNTL_H 1

/* Define if you have the <sys/param.h> header file.  */
#define HAVE_SYS_PARAM_H 1

/* Define if you have the <sys/time.h> header file.  */
#define HAVE_SYS_TIME_H 1

/* Define if you have the <sys/select.h> header file.  */
#define HAVE_SYS_SELECT_H 1

/* Define if your <sys/time.h> declares struct tm.  */
/* #undef TM_IN_SYS_TIME */

/* Define if <time.h> and <sys/time.h> can be both included.  */
#define TIME_WITH_SYS_TIME 1

/* Define if you have the gettimeofday function.  */
#define HAVE_GETTIMEOFDAY 1

/* Define if you have the getcwd function.  */
#define HAVE_GETCWD 1

/* Define if you have the gethostname function.  */
#define HAVE_GETHOSTNAME 1

/* Define if you have the memcpy function.  */
#define HAVE_MEMCPY 1

/* Define if you have the memmove function.  */
#define HAVE_MEMMOVE 1

/* Define if you have the unlink function.  */
#define HAVE_UNLINK 1

/* Define if you have the fcntl function.  */
#define HAVE_FCNTL 1

/* Define if you have the popen function.  */
#define HAVE_POPEN 1

/* Define if you have the fnmatch function.  */
#define HAVE_FNMATCH 1

/* Define if you have the on_exit function.  */
#define HAVE_ON_EXIT 1

/* Define if you have the strstr function.  */
#define HAVE_STRSTR 1

/* Define if you have the strerror function.  */
#define HAVE_STRERROR 1

/* Define if stdio.h defines the sys_errlist table */
/* #undef HAVE_SYS_ERRLIST_DECL */

/* Define if you have the vsnprintf function.  */
#define HAVE_VSNPRINTF 1

/* Define if you have the dlopen function.  */
#define HAVE_DLOPEN 1

/* Define if dlfcn.h defines the RTLD_NOW constant */
#define HAVE_RTLD_NOW 1

/* Define if you have the shl_load function.  */
/* #undef HAVE_SHL_LOAD */

/* If words are stored with the most significant byte first (like
                   Motorola and SPARC, but not Intel and VAX, CPUs */
/* #undef WORDS_BIGENDIAN */

/* Define if your FPU arithmetics is of the DEC type.  */
/* #undef HAVE_DEC_FPU */

/* Define if your FPU arithmetics is of the little endian IEEE type.  */
#define HAVE_LIEEE_FPU 1

/* Define if your FPU arithmetics is of the big endian IEEE type.  */
/* #undef HAVE_BIEEE_FPU */

/* Define if realloc() is buggy (SunOS4, for example) */
/* #undef REALLOC_IS_BUGGY */

/* Define if you have the drand48 function.  */
#define HAVE_DRAND48 1

/* Define if your system supports locale.  */
#define HAVE_SETLOCALE 1

/* Define if stdlib.h defines the drand48 function.  */
#define HAVE_DRAND48_DECL 1

/* Define if you have the m library (-lm).  */
#define HAVE_LIBM 1

/* Define if you have the <math.h> header file.  */
#define HAVE_MATH_H 1

/* Define if you have <float.h>.  */
#define HAVE_FLOAT_H 1

/* Define if you have <ieeefp.h>.  */
/* #undef HAVE_IEEEFP_H */

/* Define if you have the hypot function.  */
#define HAVE_HYPOT 1

/* Define if the hypot function is declared in math.h.  */
#define HAVE_HYPOT_DECL 1

/* Define if you have the cbrt function.  */
#define HAVE_CBRT 1

/* Define if the cbrt function is declared in math.h.  */
#define HAVE_CBRT_DECL 1

/* Define if you have the log2 function.  */
#define HAVE_LOG2 1

/* Define if the log2 function is declared in math.h.  */
#define HAVE_LOG2_DECL 1

/* Define if you have the rint function.  */
#define HAVE_RINT 1

/* Define if the rint function is declared in math.h.  */
#define HAVE_RINT_DECL 1

/* Define if you have the lgamma function.  */
#define HAVE_LGAMMA 1

/* Define if the lgamma function is declared in math.h.  */
#define HAVE_LGAMMA_DECL 1

/* Define if math.h defines the signgam variable.  */
#define HAVE_SIGNGAM_DECL 1

/* Define if you have the asinh function.  */
#define HAVE_ASINH 1

/* Define if the asinh function is declared in math.h.  */
#define HAVE_ASINH_DECL 1

/* Define if you have the acosh function.  */
#define HAVE_ACOSH 1

/* Define if the acosh function is declared in math.h.  */
#define HAVE_ACOSH_DECL 1

/* Define if you have the atanh function.  */
#define HAVE_ATANH 1

/* Define if the atanh function is declared in math.h.  */
#define HAVE_ATANH_DECL 1

/* Define if you have the erf function.  */
#define HAVE_ERF 1

/* Define if the erf function is declared in math.h.  */
#define HAVE_ERF_DECL 1

/* Define if you have the erfc function.  */
#define HAVE_ERFC 1

/* Define if the erfc function is declared in math.h.  */
#define HAVE_ERFC_DECL 1

/* Define if you have the finite function.  */
#define HAVE_FINITE 1

/* Define if the finite function is declared in math.h.  */
#define HAVE_FINITE_DECL 1

/* Define if you have the isfinite function.  */
/* #undef HAVE_ISFINITE */

/* Define if the isfinite function is declared in math.h.  */
/* #undef HAVE_ISFINITE_DECL */

/* Define if you have the isnan function.  */
#define HAVE_ISNAN 1

/* Define if the isnan function is declared in math.h.  */
#define HAVE_ISNAN_DECL 1

/* Define if you have the Bessel j0 function.  */
#define HAVE_J0 1

/* Define if the j0 function is declared in math.h.  */
#define HAVE_J0_DECL 1

/* Define if you have the Bessel j1 function.  */
#define HAVE_J1 1

/* Define if the j1 function is declared in math.h.  */
#define HAVE_J1_DECL 1

/* Define if you have the Bessel jn function.  */
#define HAVE_JN 1

/* Define if the jn function is declared in math.h.  */
#define HAVE_JN_DECL 1

/* Define if you have the Bessel y0 function.  */
#define HAVE_Y0 1

/* Define if the y0 function is declared in math.h.  */
#define HAVE_Y0_DECL 1

/* Define if you have the Bessel y1 function.  */
#define HAVE_Y1 1

/* Define if the y1 function is declared in math.h.  */
#define HAVE_Y1_DECL 1

/* Define if you have the Bessel yn function.  */
#define HAVE_YN 1

/* Define if the yn function is declared in math.h.  */
#define HAVE_YN_DECL 1

/* Define if netCDF library is available */
#define HAVE_NETCDF 1

/* Define if FFTW Fourier transform library is available */
/* #undef HAVE_FFTW */

/* Define if PNG library is available */
#define HAVE_LIBPNG 1

/* Define if JPEG library is available */
#define HAVE_LIBJPEG 1

/* Define if PDFlib library is available */
/* #undef HAVE_LIBPDF */

/* Define if you want to compile a Fortran wrapper for grace_np lib.  */
#define WITH_F77_WRAPPER 1

/* Define if the X Window System is missing or not being used.  */
/* #undef X_DISPLAY_MISSING */

/* Define if you have Motif.  */
#define HAVE_MOTIF 1

/* Define if Motif library is actually Lesstif.  */
/* #undef HAVE_LESSTIF */

/* Define if Motif library includes _XmVersionString[] string.  */
#define HAVE__XMVERSIONSTRING 1

/* Define if you have the Xpm library (-lXpm).  */
#define HAVE_XPM 1

/* Define if you have the xpm.h.  */
/* #undef HAVE_XPM_H */

/* Define if you have the xpm.h header among X11 includes.  */
#define HAVE_X11_XPM_H 1

/* Define if you have (and want to use) libhelp  */
/* #undef WITH_LIBHELP */

/* Define if you have (and want to use) XmHTML widget  */
/* #undef WITH_XMHTML */

/* Define if you have (and want to use) editres  */
#define WITH_EDITRES 1

#if (defined(HAVE_MOTIF) && !defined(X_DISPLAY_MISSING))
#  define MOTIF_GUI
#else
#  define NONE_GUI
#endif

/* Define if the print spooling command itself unlinks the temporary file */
/* #undef PRINT_CMD_UNLINKS */

/* Define if you want to compile in (a basic) support for debugging  */
/* #undef WITH_DEBUG */

#if defined(WITH_DEBUG)
#  define DEBUG
#else
#  define NDEBUG
#endif

#endif /* __CONFIG_H */
