SUBDIRS = turing_problem \
	turing_problem_ccfv	

# don't follow the full GNU-standard
# we need automake 1.9
AUTOMAKE_OPTIONS = foreign 1.9

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-geometry=$(DUNE_GEOMETRY_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-typetree=$(DUNE_TYPETREE_ROOT) --with-dune-pdelab=$(DUNE_PDELAB_ROOT)  CXX="$(CXX)" CC="$(CC)"

EXTRA_DIST = CMakeLists.txt

include $(top_srcdir)/am/global-rules

