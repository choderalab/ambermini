#
#************************************************************************
#                              AMBER                                   **
#                                                                      **
#               Copyright (c) 1986, 1991, 1995, 1997, 1999             **
#                Regents of the University of California               **
#                       All Rights Reserved.                           **
#                                                                      **
#  This software provided pursuant to a license agreement containing   **
#  restrictions on its disclosure, duplication, and use. This software **
#  contains confidential and proprietary information, and may not be   **
#  extracted or distributed, in whole or in part, for any purpose      **
#  whatsoever, without the express written permission of the authors.  **
#  This notice, and the associated author list, must be attached to    **
#  all copies, or extracts, of this software. Any additional           **
#  restrictions set forth in the license agreement also apply to this  **
#  software.                                                           **
#************************************************************************
#
# Makefile for SQM (stand-alone quantum program)
#
include ../config.h

# Sqm Fortran source files are free format
LOCALFLAGS = $(FREEFORMAT_FLAG)

SHELL=/bin/sh
VPATH=

#SQMOBJ = Objects specific to the standalone SQM
SQMOBJ= file_io_dat.o constants.o findmask.o sqm.o \
	xmin.o xminC.o \
	amopen.o mexit.o assert.o \
	timer_dummy.o nmlsrc.o \
	qm_print_info.SQM.o \
	qm2_energy.SQM.o \
	qm2_read_nm_and_alloc.SQM.o \
	qm2_scf.SQM.o \
	qm2_dftb_scf.SQM.o

#LIBQMOBJ = Objects specific to building as a library.
LIBQMOBJ = qm_print_info.o \
	qm2_energy.o \
	qm2_read_nm_and_alloc.o \
	qm2_scf.o \
	qm2_dftb_ewald.o \
	qm2_dftb_scf.o

#QMOBJ = Objects common to both SQM and libsqm
QMOBJ = qm2_allocate_e_repul.o qm2_calc_charges.o qm2_calc_dipole.o \
	qm2_calc_rij_and_eqns.o qm2_dihed.o qm2_fock.o \
	qm2_get_qm_forces.o qm2_get_qmmm_forces.o qm2_h1elec.o qm2_hcore_qmqm.o \
	qm2_hcore_qmmm.o qm2_identify_peptide_links.o \
	qm2_load_params_and_allocate.o qm2_repp.o qm2_rotate_qmqm.o \
	qm2_core_core_repulsion.o qm2_core_core_repulsion_dxyz.o \
	qm2_iterator_mod.o qm2_diagonalizer_module.o \
	qm2_setup_orb_exp.o qm2_smallest_number.o qm2_fock_predict.o \
	qm_gb.o qm_zero_charges.o qm_assign_atom_types.o \
	qm_link_atoms.o qm2_print_charges.o \
	qmmm_qmtheorymodule.o qm2_print_bondorders.o \
	qm2_pm6_hof_module.o qm2_params_module.o\
	qm2_dftb_module.o qm2_dftb_broyden.o \
	qm2_dftb_dispersion_egr.o qm2_dftb_dispersion_params.o \
	qm2_dftb_dispersionread.o qm2_dftb_energy.o \
	qm2_dftb_ewevge.o qm2_dftb_externalshift.o \
	qm2_dftb_fermi.o qm2_dftb_get_qm_forces.o qm2_dftb_gamma.o \
	qm2_dftb_gammamat.o qm2_dftb_get_qmmm_forces.o qm2_dftb_gettab.o \
	qm2_dftb_load_params.o \
	qm2_dftb_mulliken.o \
	qm2_dftb_repulsiv.o qm2_dftb_self.o qm2_dftb_shift.o \
	qm2_dftb_skpar.o qm2_dftb_slkode.o \
	qm2_dftb_slktrafo.o qm2_dftb_read_cm3.o qm2_dftb_cm3.o \
	qm2_dftb_gb.o qm2_dftb_3rd_order.o \
	qmmm_module.o qmmm_nml_module.o qmmm_struct_module.o qmmm_vsolv_module.o \
	qm2_parameters.o dh_correction_module.o \
	utilitiesModule.o elementOrbitalIndex.o parameterReader.o \
	slater_overlap.o qm2_h1elec_d.o rotation.o qm2_repp_d.o \
	opnq_Edisp.o  opnq_Erep.o  opnq_Evdw.o  opnq.o  opnq_SwitchMod.o \
	qm2_fock_d.o MNDOChargeSeparation.o qm2_print_energy.o


install: sqm$(SFX)
	mv sqm$(SFX) $(BINDIR)

sqm$(SFX): $(SQMOBJ) $(QMOBJ) netlib sys
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o sqm$(SFX) $(SQMOBJ) $(QMOBJ) \
		-L$(LIBDIR) $(FLIBSF) ../lib/sys.a $(LDFLAGS) $(AMBERLDFLAGS)

#-----------LIBS

xminC.o:  ../sff/xminC.c
	$(CC) -c -DSQM $(COPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) -o xminC.o ../sff/xminC.c

libsqm: $(LIBDIR)/libsqm.a

$(LIBDIR)/libsqm.a: $(QMOBJ) $(LIBQMOBJ)
	$(AR) libsqm.a $(QMOBJ) $(LIBQMOBJ)
	$(RANLIB) libsqm.a
	mv libsqm.a $(LIBDIR)

sys:
	cd ../lib; $(MAKE) sys.a

netlib:
	cd ../lapack; $(MAKE) $(LAPACK)
	cd ../blas; $(MAKE) $(BLAS)
	cd ../arpack && $(MAKE) install

clean:
	/bin/rm -f *.o  *.mod *.d sqm$(SFX)

uninstall:
	-rm -f $(BINDIR)/sqm$(SFX)
	-rm -f $(LIBDIR)/libsqm.a

depend::
	./makedepend > depend

#----------HEADER DEPENDENCIES (constructed from ./makedepend > depend)
include depend
# DO NOT DELETE
