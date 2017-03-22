AC_PREREQ([2.65])

AC_INIT

case "$OSTYPE" in "darwin"*)
  BREWDIR=`brew --prefix`
  export PKG_CONFIG_PATH=$BREWDIR/opt/openssl/lib/pkgconfig
  ;;
esac

PKG_CHECK_MODULES(OPENSSL, [openssl >= 1.0], [OPENSSL="yes"], [OPENSSL="no"])

if test "${OPENSSL}" != "no"; then
  SSL_CPPFLAGS="-DUSE_SSL"
fi

AC_SUBST(SSL_CPPFLAGS)

AC_OUTPUT(src/Makevars)
