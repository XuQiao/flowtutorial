builddir=`pwd`
installdir="${builddir/build/install}"
configure="../EPAnalyzer/autogen.sh --prefix="$installdir
echo $configure
eval $configure
make
make install
