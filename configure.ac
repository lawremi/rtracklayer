AC_PREREQ([2.65])

AC_INIT

PKG_CHECK_MODULES(OPENSSL, [openssl >= 1.0], [OPENSSL="yes"], [OPENSSL="no"])

if test "${OPENSSL}" != "no"; then
  SSL_CPPFLAGS="-DUSE_SSL"
else
  case "$OSTYPE" in "darwin"*)
  BREWDIR=`brew --prefix`
  OPENSSL_CFLAGS="-I$BREWDIR/opt/openssl/include"
  OPENSSL_LIBS="-L$BREWDIR/opt/openssl/lib"
  SSL_CPPFLAGS="-DUSE_SSL"
  ;;
  esac
fi

AC_SUBST(SSL_CPPFLAGS)

AC_OUTPUT(src/Makevars)
