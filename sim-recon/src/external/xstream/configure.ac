AC_PREREQ(2.59)
dnl AC_INIT(include/xstream.h)
AC_INIT(XStream, 0.0.2)
AM_CONFIG_HEADER(include/xstream/config.h:config.h.in)
AM_INIT_AUTOMAKE()

dnl Checks for programs

AC_LANG([C++])
AC_PROG_CXX
AC_PROG_CXXCPP
AC_PROG_INSTALL
AC_PROG_RANLIB

dnl compiler amd types stuff

AC_C_CONST
AC_C_INLINE
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T
AC_STRUCT_TM

dnl headers

AC_HEADER_STDC
AC_CHECK_HEADERS([stdint.h inttypes.h])

dnl endianess for md5 optimizations

AC_ARG_ENABLE(
	md5_optimization,
	AC_HELP_STRING(
		[--enable-md5-optimization],
		[if md5 should be compiled with little endian optimizations]
	),
	WITH_MD5_OPTIM="$enableval",
	WITH_MD5_OPTIM="yes"
)

LITTLE_ENDIAN="no"

if test "$WITH_MD5_OPTIM" = "yes";
then

AC_C_BIGENDIAN(
	AC_DEFINE([ARCH_LITTLE_ENDIAN],0,[if machine is little endian]),
	AC_DEFINE([ARCH_LITTLE_ENDIAN],1) LITTLE_ENDIAN="yes",
	AC_DEFINE([ARCH_LITTLE_ENDIAN],0)
)
fi

dnl special functions

dnl check if we can support the dater filter

AC_ARG_ENABLE(
	dater,
	AC_HELP_STRING(
		[--enable-dater],
		[if dater filter is to be installed]
	),
	WITH_DATER="$enableval",
	WITH_DATER="yes"
)

if test "$WITH_DATER" = "yes"; then
	dnl on my machine I have strftime and dater compiles but autoconf keeps saying I don't have it
	dnl so I removed the test
	dnl go figure ...
	dnl AC_CHECK_FUNC([strftime], [AC_DEFINE([HAVE_STRFTIME],1,[if strftime was found])],[WITH_DATER="no"])
	AC_CHECK_FUNC([localtime_r], [AC_DEFINE([HAVE_LOCALTIME_R],1,[if localtime_r was found])], [WITH_DATER="no"])

	AC_MSG_CHECKING([if dater support is possible])
	AC_MSG_RESULT([$WITH_DATER]);
fi
AM_CONDITIONAL(WITH_DATER, test "${WITH_DATER}" = "yes")

AC_ARG_ENABLE(
	examples,
	AC_HELP_STRING(
		[--enable-examples],
		[if examples should be compiled]
	),
	WITH_EXAMPLES="$enableval",
	WITH_EXAMPLES="yes"
)
AM_CONDITIONAL(WITH_EXAMPLES, test "${WITH_EXAMPLES}" = "yes")

AC_ARG_ENABLE(
	tests,
	AC_HELP_STRING(
		[--enable-tests],
		[if tests should be compiled]
	),
	WITH_TEST="$enableval",
	WITH_TEST="no"
)
AM_CONDITIONAL(WITH_TEST, test "${WITH_TEST}" = "yes")

AC_ARG_ENABLE(
	docs,
	AC_HELP_STRING(
		[--enable-docs],
		[if examples should be compiled]
	),
	WITH_DOCS="$enableval",
	WITH_DOCS="yes"
)

dnl if docs need to be built and if doxygen if found (dot looked for and doxygen.cfg is updated accordingly)
if test "${WITH_DOCS}" = yes; then
	AC_PATH_PROG([DOXYGEN],[doxygen])
	AC_PATH_PROG([DOT],[dot])
	if test -z "${DOXYGEN}"; then
		WITH_DOCS="no"
	else
		if test -z "${DOT}"; then
			AC_SUBST(DOC_CLASS_DIAGRAMS,"YES")
			AC_SUBST(DOC_DOT,"NO")
		else
			AC_SUBST(DOC_CLASS_DIAGRAMS,"NO")
			AC_SUBST(DOC_DOT,"YES")
		fi
	fi
fi
AM_CONDITIONAL(WITH_DOCS, test "${WITH_DOCS}" = "yes")

AC_PATH_PROG([PERL],[perl])

dnl check if we can compile fd::streambuf

AC_CHECK_FUNC([read], [WITH_FD="yes"], [WITH_FD="no"])
if test "${WITH_FD}" = yes; then
	AC_CHECK_FUNC([write], [WITH_FD="yes"], [WITH_FD="no"])
	if test "${WITH_FD}" = yes; then
		AC_CHECK_FUNC([fsync], [WITH_FD="yes"], [WITH_FD="no"])
	fi
fi

AM_CONDITIONAL(WITH_FD, test "${WITH_FD}" = "yes")
AC_DEFINE([HAVE_FD],1,[if file descriptor support exists])


dnl library checks

CV_CHECK_LIB(bzlib, [/usr/local /usr],[bzlib.h],BZ2_bzDecompressEnd,bz2,
	[
	AC_DEFINE(HAVE_LIBBZ2, 1, [if bzlib was found])
	HAVE_LIBBZ2="yes"
	],
	[
	HAVE_LIBBZ2="no"
	]
)
AM_CONDITIONAL(HAVE_LIBBZ2, test "${HAVE_LIBBZ2}" = yes)

CV_CHECK_LIB(zlib, [/usr/local /usr],[zlib.h],inflateEnd,z,
	[
	AC_DEFINE(HAVE_LIBZ, 1, [if zlib was found])
	HAVE_LIBZ="yes"
	],
	[
	HAVE_LIBZ="no"
	]
)
AM_CONDITIONAL(HAVE_LIBZ, test "${HAVE_LIBZ}" = yes)

dnl developer stuff

AC_ARG_ENABLE(
	logging,
	AC_HELP_STRING([--enable-logging],[logging messages should be enabeled]), 
	[WITH_LOGGING="$enableval"],
	[WITH_LOGGING="no"]
)
if test "${WITH_LOGGING}" = "yes"; then
	AC_DEFINE(ENABLE_LOGGING, 1, if logging should be written to clog)
fi

AC_OUTPUT(Makefile src/Makefile doc/Makefile doc/doxygen.cfg doc/doxygen_devel.cfg examples/Makefile test/Makefile)

AC_MSG_NOTICE([
*****************************************************

	configure result:

	build documentation:    ${WITH_DOCS}
	build examples:         ${WITH_EXAMPLES}
	build tests:            ${WITH_TEST}
	zlib support:           ${HAVE_LIBZ}
	bzlib support:          ${HAVE_LIBBZ2}
	dater filter support:   ${WITH_DATER}
	fd streambuf support:   ${WITH_FD}
	md5 little endian:      ${LITTLE_ENDIAN}
	logging:                ${WITH_LOGGING}

*****************************************************

	have a look at ./configure --help for options
])
