AC_PREREQ([2.65])

AC_INIT("rtracklayer")

### NOTE: ideally, we would generate a config.h here and use the
### results of these checks. Way too lazy to attempt anything, but at
### least the configure output can help diagnose build issues.

## Checks for libraries
AC_CHECK_LIB([z], [gzeof],,[AC_MSG_ERROR([need libz with gzeof])])

## Checks for headers
AC_CHECK_HEADERS([zlib.h arpa/inet.h fcntl.h netdb.h netinet/in.h stdlib.h string.h strings.h sys/file.h sys/socket.h sys/time.h termios.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_OFF_T
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_UID_T

# Checks for library functions.
AC_FUNC_FORK
AC_FUNC_MALLOC
AC_FUNC_MKTIME
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([bzero dup2 getcwd gethostbyname gettimeofday memmove memset mkdir regcomp select socket sqrt strcasecmp strchr strcspn strdup strerror strrchr strspn strstr uname])

AC_OUTPUT(src/Makevars)
