# Name of executables: The prefix "noinst" means the executables don't
# have to be installed anywhere.
noinst_PROGRAMS= \
two_d_stick_slip_with_face_element

#one_d_poisson_with_face_element
#three_d_poisson_with_face_element

#---------------------------------------------------------------------------

# Local sources that Joe's Poisson code depends on:
two_d_stick_slip_with_face_element_SOURCES = two_d_stick_slip_with_face_element.cc stokes_sing_face_element.h unstructured_stick_slip_mesh.h

# Required libraries: 
two_d_stick_slip_with_face_element_LDADD = -L@libdir@ -lnavier_stokes -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.
AM_CPPFLAGS +=  -I@includedir@  

EXTRA_DIST = check_condition_number.bash  impose_amplitude_runs.bash  validate_rc.bash condition_number.lay




