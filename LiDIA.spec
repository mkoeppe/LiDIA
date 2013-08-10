# RPM spec file for LiDIA

# This spec file makes use of macros that come with RPM version >= 3.04.
# As a result, RPM will automatically configure the package to fit into
# the host's filesystem structure.
# Make sure that you have appropriate optflags defined or set CXXFLAGS
# in the environment before starting to build.

Summary:	C++ library for number theoretical computations
Name:		lidia
Version:	2.3.0
Release:	0
Vendor:		The LiDIA Group
Copyright:	noncommercial
Group:		Development/Libraries
Source:		ftp://ftp.informatik.tu-darmstadt.de/pub/TI/systems/LiDIA/current/%{name}-%{version}.tar.gz
URL:		http://www.informatik.tu-darmstadt.de/TI/LiDIA/
# Actually, you can choose among several multiprecision libraries.
#Requires:	gmp
#BuildRequires:	gmp-devel
BuildRequires:	tetex te_latex texinfo
BuildRoot:	%{_tmppath}/%{name}-install

%description
LiDIA is a C++ library for number-theoretical computations.

LiDIA provides parameterized as well as specialized classes together
with advanced methods for computing in a large variety of mathematical
groups, rings, and fields, ranging from arbitrary-length integers,
fractions, and floating-point approximations, over vectors and matrices,
to high-level constructs of finite fields, lattices, quadratic and
higher-order number fields, polynomial rings, and elliptic curves.
In addition to a rich set of methods for basic structural operations,
I/O, arithmetic, and common number-theoretical primitives, there are
functions for computationally intensive tasks: reducing lattice bases;
factoring polynomials, integers, or algebraic ideals; computing discrete
logarithms; determining group orders on elliptic curves; and more...

LiDIA is free for non-commercial use.  See the copyright notice in the
file COPYING.  Contributors are welcome.

%prep
%setup -q

%build
# Consider setting CPPFLAGS and LIBS in the environment if you want to use
# libs in non-standard locations.  If you do not have a proper optflags
# setting in your /etc/rpmrc or ~/.rpmmacros, you may want to set CXXFLAGS
# too.  Feel free to experiment with some of LiDIA's configure options.
CFLAGS="${CFLAGS-$RPM_OPT_FLAGS}" \
CXXFLAGS="${CXXFLAGS-$RPM_OPT_FLAGS}" \
./configure \
	--prefix=%{_prefix} \
	--exec-prefix=%{_exec_prefix} \
	--libdir=%{_libdir} \
	--includedir=%{_includedir} \
	--datadir=%{_datadir} \
	--enable-shared \
	--disable-static
make
#make examples
make pdf	# some of: dvi ps psgz psbz2 pdf dsc
		# If you change this, adjust the %doc entry below

%install
rm -rf "$RPM_BUILD_ROOT"
DESTDIR="$RPM_BUILD_ROOT" make install	# install-examples

%clean
rm -rf "$RPM_BUILD_ROOT"

%files
%defattr(-,root,root)
%doc COPYING NEWS README TODO
%doc doc/LiDIA.pdf
%{_libdir}/*LiDIA*
%{_includedir}/LiDIA
%{_datadir}/LiDIA

%changelog
* Thu May 22 2003 Christian Cornelssen <ccorn@cs.tu-berlin.de>
- Add spec file template to LiDIA's configure system
* Fri Apr  4 2003 Christian Cornelssen <ccorn@cs.tu-berlin.de> 2.1pre7-4
- Upgrade to 2.1pre7 with patch3
* Mon Oct 28 2002 Christian Cornelssen <ccorn@cs.tu-berlin.de> 2.1pre6-1
- first version for 2.1pre6 using the new build system

