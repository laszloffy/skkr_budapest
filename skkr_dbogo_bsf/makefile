##################################################################
#                                                                #
# Makefile for the Dirac-Bogo SKKR program                       #
#                                                                #
#                                                                #
#       SPINPOL. REL.                                            #
#         make rsp-scf_ksave.run  => small memory, works for any #
#                                    case, optimized for bulk    #
#                                                                #
##################################################################
#
# fortran compiler:
#
FC = mpifort
#LD = mpif90
FCSER = gfortran
#LDSER = ifort
#
# compiler options:
#
#FFLAGS = -O3 -xHost -w -mkl=sequential -unroll -132 -opt-matmul -mcmodel=medium
#LFLAGS =  -lsvml -opt-matmul -mcmodel=small -i-static 
#FFLAGS = -O2 -C -gen-interfaces -warn interfaces -w -unroll -132 -mcmodel=medium -i-dynamic  -g -traceback  -fp-stack-check
#FFLAGS =   -O2 -w  -unroll  -mcmodel=medium #used on strom
# clustr.phy.bme.hu fflags
#FFLAGS =  -O2  -w -unroll  -mcmodel=medium  #for general use, latest
#
# hpc_debrecen fflags:
FFLAGS = -xSSE4.2 -fpp -traceback -O2 -unroll -i4 -132 -mcmodel=medium -fp-model strict -check bounds #for general use, latest
FFLAGS = -O2 -cpp -fbacktrace -w -unroll -mcmodel=medium 
#FFLAGS = -O2 -cpp -w -unroll  -mcmodel=medium  #for general use, latest
#LFLAGS =  -lpthread 
#FFLAGS = -C -O0 -gen-interfaces -warn interfaces -w -unroll -132 -mcmodel=medium -i-dynamic
LAPACK=~/source/lapack-3.10.1
BLAS=~/source/BLAS-3.10.0
LFLAGS = -pthread -I/opt/nce/packages/global/openmpi/3.1.6/lib -Wl,-rpath -Wl,/opt/nce/packages/global/openmpi/3.1.6/lib -Wl,--enable-new-dtags -L/opt/nce/packages/global/openmpi/3.1.6/lib -lmpi_usempi -lmpi_mpifh -lmpi -Wl,--start-group $(BLAS)/blas_LINUX.a $(LAPACK)/liblapack.a -Wl,--end-group
#LIB=$(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a \
     -Wl,--start-group \
     $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
     $(MKLROOT)/lib/intel64/libmkl_sequential.a \
     $(MKLROOT)/lib/intel64/libmkl_core.a \
     $(MKLROOT)/lib/intel64/libmkl_blacs_sgimpt_lp64.a -Wl,--end-group \
     -lpthread -lm \
     -check bounds

#LFLAGS = -Vaxlib -lsvml -mcmodel=medium -i-dynamic
#LFLAGS=$(LIB)
#
# working directory
#
KPATH = ./
#
# subdirectories 
#
SUBDIRS = clu main rsp-scf 
#
# programs
#
PROGRAMS = rsp-scf_ksave_mpi.run rsp-scf_ksave_serial.run
#
# generate a unique object from the subroutines 
# in different subdirectories
# 
MAINPATH = $(KPATH)main/
#
MAINOBJ = \
$(MAINPATH)alloccheck.o \
$(MAINPATH)alphamat.o   $(MAINPATH)averagebpot.o $(MAINPATH)blm.o \
$(MAINPATH)bulkmad1.o   $(MAINPATH)bulkmad2.o \
$(MAINPATH)broyd.o      $(MAINPATH)cerf.o        $(MAINPATH)chcore.o \
$(MAINPATH)coulint.o    $(MAINPATH)csbf.o        $(MAINPATH)delkam.o \
$(MAINPATH)ibz.o        $(MAINPATH)directg.o     $(MAINPATH)dugo.o \
$(MAINPATH)erf.o        $(MAINPATH)error.o       $(MAINPATH)exchange.o \
$(MAINPATH)fac.o        $(MAINPATH)findangles.o  $(MAINPATH)fullbz.o \
$(MAINPATH)gaunt.o      $(MAINPATH)gauntcalc.o \
$(MAINPATH)genlatt2d.o  $(MAINPATH)genlatt3d.o	 $(MAINPATH)gjinv2.o \
$(MAINPATH)getg2d.o 	$(MAINPATH)getsfg.o 	\
$(MAINPATH)getg2d_io.o 	$(MAINPATH)getsfg_io.o 	\
$(MAINPATH)gstore.o     $(MAINPATH)ilkam.o       $(MAINPATH)incgamma.o \
$(MAINPATH)indekx.o     $(MAINPATH)init.o        $(MAINPATH)inouh.o \
$(MAINPATH)inth.o       $(MAINPATH)jobnam.o      $(MAINPATH)kamstr.o \
$(MAINPATH)kmesh2d.o \
$(MAINPATH)latgen2d.o   $(MAINPATH)latgen3d.o \
$(MAINPATH)legpol.o 	$(MAINPATH)legzero.o	 $(MAINPATH)ll.o \
$(MAINPATH)lmfill.o 	$(MAINPATH)madcheck.o \
$(MAINPATH)madelungb.o  $(MAINPATH)madelungi.o   $(MAINPATH)madelungv.o \
$(MAINPATH)matops.o     $(MAINPATH)newfl.o       $(MAINPATH)packer.o \
$(MAINPATH)phi.o        $(MAINPATH)poisson.o     $(MAINPATH)primkr2d.o \
$(MAINPATH)rcore.o      $(MAINPATH)resld.o       $(MAINPATH)rgroup.o \
$(MAINPATH)rotmat.o     $(MAINPATH)rotmat2.o \
$(MAINPATH)rotvec.o     $(MAINPATH)rsimp.o \
$(MAINPATH)scrinv.o     $(MAINPATH)scrstr.o      $(MAINPATH)setspk.o \
$(MAINPATH)signum.o     $(MAINPATH)speck.o       $(MAINPATH)speck2.o \
$(MAINPATH)spherh.o     $(MAINPATH)sphrm.o       $(MAINPATH)surfgf.o \
$(MAINPATH)struc3d.o    $(MAINPATH)sublattpot.o  $(MAINPATH)sublattmat.o \
$(MAINPATH)tcpaio.o 	$(MAINPATH)tau2d.o	 $(MAINPATH)vecrot.o \
$(MAINPATH)vecsrt2d.o 	$(MAINPATH)vecsrt3d.o \
$(MAINPATH)vector2d.o   $(MAINPATH)vector3d.o \
$(MAINPATH)ylag.o       $(MAINPATH)ylmc.o 	 $(MAINPATH)ytrafo.o \
$(MAINPATH)zero.o       $(MAINPATH)zmesh.o \
$(MAINPATH)clebsch.o    $(MAINPATH)gafill.o      $(MAINPATH)relmtrx.o \
$(MAINPATH)g2d_pack.o
#
MAINSRC = \
$(MAINPATH)alloccheck.f \
$(MAINPATH)alphamat.f   $(MAINPATH)averagebpot.f $(MAINPATH)blm.f \
$(MAINPATH)bulkmad1.f   $(MAINPATH)bulkmad2.f \
$(MAINPATH)broyd.f      $(MAINPATH)cerf.f        $(MAINPATH)chcore.f \
$(MAINPATH)coulint.f    $(MAINPATH)csbf.f        $(MAINPATH)delkam.f \
$(MAINPATH)ibz.f        $(MAINPATH)directg.f     $(MAINPATH)dugo.f \
$(MAINPATH)erf.f        $(MAINPATH)error.f       $(MAINPATH)exchange.f \
$(MAINPATH)fac.f        $(MAINPATH)findangles.f  $(MAINPATH)fullbz.f \
$(MAINPATH)gaunt.f      $(MAINPATH)gauntcalc.f \
$(MAINPATH)genlatt2d.f 	$(MAINPATH)genlatt3d.f   $(MAINPATH)gjinv2.f \
$(MAINPATH)getg2d.f     $(MAINPATH)getsfg.f      \
$(MAINPATH)getg2d_io.f  $(MAINPATH)getsfg_io.f   \
$(MAINPATH)gstore.f     $(MAINPATH)ilkam.f       $(MAINPATH)incgamma.f \
$(MAINPATH)indekx.f     $(MAINPATH)init.f        $(MAINPATH)inouh.f \
$(MAINPATH)inth.f       $(MAINPATH)jobnam.f      $(MAINPATH)kamstr.f \
$(MAINPATH)kmesh2d.f \
$(MAINPATH)latgen2d.f   $(MAINPATH)latgen3d.f \
$(MAINPATH)legpol.f \
$(MAINPATH)legzero.f    $(MAINPATH)ll.f          $(MAINPATH)lmfill.f \
$(MAINPATH)madcheck.f \
$(MAINPATH)madelungb.f  $(MAINPATH)madelungi.f   $(MAINPATH)madelungv.f \
$(MAINPATH)matops.f     $(MAINPATH)newfl.f       $(MAINPATH)packer.f \
$(MAINPATH)phi.f        $(MAINPATH)poisson.f     $(MAINPATH)primkr2d.f \
$(MAINPATH)rcore.f      $(MAINPATH)resld.f       $(MAINPATH)rgroup.f \
$(MAINPATH)rotmat.f     $(MAINPATH)rotmat2.f \
$(MAINPATH)rotvec.f     $(MAINPATH)rsimp.f \
$(MAINPATH)scrinv.f     $(MAINPATH)scrstr.f      $(MAINPATH)setspk.f \
$(MAINPATH)signum.f     $(MAINPATH)speck.f       $(MAINPATH)speck2.f \
$(MAINPATH)spherh.f     $(MAINPATH)sphrm.f       $(MAINPATH)surfgf.f \
$(MAINPATH)struc3d.f	 $(MAINPATH)sublattpot.f $(MAINPATH)sublattmat.f \
$(MAINPATH)tcpaio.f	$(MAINPATH)tau2d.f	 $(MAINPATH)vecrot.f \
$(MAINPATH)vecsrt2d.f 	$(MAINPATH)vecsrt3d.f \
$(MAINPATH)vector2d.f   $(MAINPATH)vector3d.f \
$(MAINPATH)ylag.f       $(MAINPATH)ylmc.f	 $(MAINPATH)ytrafo.f \
$(MAINPATH)zero.f       $(MAINPATH)zmesh.f \
$(MAINPATH)clebsch.f    $(MAINPATH)gafill.f      $(MAINPATH)relmtrx.f \
$(MAINPATH)g2d_pack.f
#
RSPSCFPATH = $(KPATH)rsp-scf/
#
RSPSCF-ksave-OBJ = \
$(RSPSCFPATH)bsf.o \
$(RSPSCFPATH)cpacoord_ksave.o \
$(RSPSCFPATH)dens.o        $(RSPSCFPATH)enprint.o     $(RSPSCFPATH)localrot.o \
$(RSPSCFPATH)initzero.o    $(RSPSCFPATH)localrot2.o   $(RSPSCFPATH)locquant.o \
$(RSPSCFPATH)newdir.o      $(RSPSCFPATH)newpot.o      $(RSPSCFPATH)pothandle.o \
$(RSPSCFPATH)printres.o    $(RSPSCFPATH)printscf.o    $(RSPSCFPATH)readgeom.o \
$(RSPSCFPATH)readini.o     $(RSPSCFPATH)restart.o     $(RSPSCFPATH)sorttmat.o \
$(RSPSCFPATH)tmatini.o     $(RSPSCFPATH)vgen.o 	      $(RSPSCFPATH)upm.o \
$(RSPSCFPATH)bop.o         $(RSPSCFPATH)brmat.o \
$(RSPSCFPATH)enblloyd.o    $(RSPSCFPATH)denblloyd.o   $(RSPSCFPATH)dmag.o \
$(RSPSCFPATH)dirsolv_inw.o $(RSPSCFPATH)dirsolv_out.o $(RSPSCFPATH)fullcpa.o \
$(RSPSCFPATH)g2fill.o      $(RSPSCFPATH)initia.o      $(RSPSCFPATH)matlms.o \
$(RSPSCFPATH)gf2.o          $(RSPSCFPATH)matr.o        $(RSPSCFPATH)matr2.o \
$(RSPSCFPATH)matrep.o      $(RSPSCFPATH)matrot.o      $(RSPSCFPATH)matrot2.o \
$(RSPSCFPATH)moment.o     $(RSPSCFPATH)phystau_d.o \
$(RSPSCFPATH)readpot.o     $(RSPSCFPATH)replms.o      $(RSPSCFPATH)replms_global.o \
$(RSPSCFPATH)replmsf.o     $(RSPSCFPATH)replmsf_global.o \
$(RSPSCFPATH)sphbas.o      $(RSPSCFPATH)matnr.o       \
$(RSPSCFPATH)tinvmat.o     $(RSPSCFPATH)tmat.o        $(RSPSCFPATH)wafu.o \
$(RSPSCFPATH)su2mat.o      $(RSPSCFPATH)pairing.o     $(RSPSCFPATH)delta_bcast.o
#
RSPSCFSRC = \
$(RSPSCFPATH)bsf.f \
$(RSPSCFPATH)cpacoord_ksave.f	$(RSPSCFPATH)cpacoord_save0.f \
$(RSPSCFPATH)cpacoord_save.f    $(RSPSCFPATH)cpacoord_io.f \
$(RSPSCFPATH)dens.f        $(RSPSCFPATH)enprint.f     $(RSPSCFPATH)localrot.f \
$(RSPSCFPATH)initzero.f    $(RSPSCFPATH)localrot2.f   $(RSPSCFPATH)locquant.f \
$(RSPSCFPATH)newdir.f      $(RSPSCFPATH)newpot.f      $(RSPSCFPATH)pothandle.f \
$(RSPSCFPATH)printres.f    $(RSPSCFPATH)printscf.f    $(RSPSCFPATH)readgeom.f \
$(RSPSCFPATH)readini.f     $(RSPSCFPATH)restart.f     $(RSPSCFPATH)sorttmat.f \
$(RSPSCFPATH)tmatini.f     $(RSPSCFPATH)vgen.f 	      $(RSPSCFPATH)upm.f \
$(RSPSCFPATH)bop.f         $(RSPSCFPATH)brmat.f \
$(RSPSCFPATH)enblloyd.f    $(RSPSCFPATH)denblloyd.f   $(RSPSCFPATH)dmag.f \
$(RSPSCFPATH)dirsolv_inw.f $(RSPSCFPATH)dirsolv_out.f $(RSPSCFPATH)fullcpa.f \
$(RSPSCFPATH)g2fill.f      $(RSPSCFPATH)initia.f      $(RSPSCFPATH)matlms.f \
$(RSPSCFPATH)gf2.f          $(RSPSCFPATH)matr.f        $(RSPSCFPATH)matr2.f \
$(RSPSCFPATH)matrep.f      $(RSPSCFPATH)matrot.f      $(RSPSCFPATH)matrot2.f \
$(RSPSCFPATH)moment.f      $(RSPSCFPATH)phystau_d.f \
$(RSPSCFPATH)readpot.f     $(RSPSCFPATH)replms.f      $(RSPSCFPATH)replms_global.f \
$(RSPSCFPATH)replmsf.f     $(RSPSCFPATH)replmsf_global.f \
$(RSPSCFPATH)sphbas.f      $(RSPSCFPATH)matnr.f \
$(RSPSCFPATH)tinvmat.f     $(RSPSCFPATH)tmat.f        $(RSPSCFPATH)wafu.f \
$(RSPSCFPATH)su2mat.f      $(RSPSCFPATH)pairing.f     $(RSPSCFPATH)delta_bcast.f
#
RSPCLUPATH = $(KPATH)clu/
# cpacoord_save.f  ecoreimp.f  gettaucl_mkl.f  readclu.f  tau2d_god.f  tauhin.f  tauhout.f  thin.f  thout.f  vhin.f  vhout.f pothandleimp.f
RSPCLUOBJ = \
$(RSPCLUPATH)readclu.o \
$(RSPCLUPATH)vhin.o \
$(RSPCLUPATH)tau2d_god.o \
$(RSPCLUPATH)pothandleimp.o \
$(RSPCLUPATH)phystau2.o \
$(RSPCLUPATH)unsct.o \
$(RSPCLUPATH)tauhout.o \
$(RSPCLUPATH)thout.o \
$(RSPCLUPATH)tauhin.o \
$(RSPCLUPATH)thin.o \
$(RSPCLUPATH)ecoreimp.o \
$(RSPCLUPATH)tmatiniimp.o \
$(RSPCLUPATH)locquant2.o \
$(RSPCLUPATH)dtat.o \
$(RSPCLUPATH)gettaucl.o \
$(RSPCLUPATH)alphamat1.o \
$(RSPCLUPATH)initzero1.o \
$(RSPCLUPATH)collect_data.o  $(RSPCLUPATH)printscf1.o \
$(RSPCLUPATH)printres1.o  $(RSPCLUPATH)triplet.o
#
RSPCLUSRC = \
$(RSPSCFPATH)readclu.f \
$(RSPSCFPATH)vhin.f \
$(RSPSCFPATH)tau2d_god.f \
$(RSPSCFPATH)pothandleimp.f \
$(RSPCLUPATH)phystau2.f \
$(RSPCLUPATH)unsct.f \
$(RSPCLUPATH)tauhout.f \
$(RSPCLUPATH)thout.f \
$(RSPCLUPATH)tauhin.f \
$(RSPCLUPATH)thin.f \
$(RSPCLUPATH)ecoreimp.f \
$(RSPCLUPATH)tmatiniimp.f \
$(RSPCLUPATH)locquant2.f \
$(RSPCLUPATH)dtat.f \
$(RSPCLUPATH)gettaucl.f \
$(RSPCLUPATH)alphamat1.f \
$(RSPCLUPATH)initzero1.f \
$(RSPCLUPATH)collect_data.f  $(RSPCLUPATH)printscf1.f \
$(RSPCLUPATH)printres1.f  $(RSPCLUPATH)triplet.f
#
     skkr-obj   = $(RSPCLUOBJ) $(MAINOBJ) $(RSPSCF-ksave-OBJ) 
     rsp-scf_ksave-obj = $(RSPCLUOBJ) $(MAINOBJ) $(RSPSCF-ksave-OBJ) 
#
     skkr-src   = $(RSPCLUSRC) $(MAINSRC) $(RSPSCFSRC) 
     rsp-scf-src = $(RSPCLUSRC) $(MAINSRC) $(RSPSCFSRC) 
#
# dependencies
#
# .SUFFIXES: .o .f
.SUFFIXES:  .run .o .f
#
rsp-scf_ksave_serial.run: FC=$(FCSER)
#
rsp-scf_ksave_mpi.run: FFLAGS += -DMPIP
#
.f.o:
	$(FC) -c  $(FFLAGS) $< 
	mv $(*F).o $(*D)
# $<   : source files of an out-of-date module
# $(*F): "file part" 
# $(*D): "directory part" of target
#
# activities: executables
#
all: $(PROGRAMS)
#
rsp-scf_ksave_mpi.run: $(RSPSCFPATH)skkr.o $(rsp-scf_ksave-obj)
	@echo " "
	@echo "Linking rsp-scf_ksave_mpi.run ..."
	@echo "-----------------------------"
	@echo " "
	$(FC) -o $@ $(RSPSCFPATH)skkr.o $(rsp-scf_ksave-obj) $(LFLAGS)
	@echo " "
	@echo "done"
#
rsp-scf_ksave_serial.run: $(RSPSCFPATH)skkr.o $(rsp-scf_ksave-obj)
	@echo " "
	@echo "Linking rsp-scf_ksave.run ..."
	@echo "-----------------------------"
	@echo " "
	$(FC) -o $@ $(RSPSCFPATH)skkr.o $(rsp-scf_ksave-obj) $(LFLAGS)
	@echo " "
	@echo "done"
#
#
# delete object files
#
del-all.o:
	@echo " "
	@echo "Deleting ALL object-files ..."
	@rm -f $(skkr-obj) $(RSPSCFPATH)skkr.o
	@echo "done"
	@echo " "
#
# utils
#
dbdg_skkr.tgz:
	@echo " "
	@echo "Making tar-file: skkr.tar ..."
	@echo "-----------------------------"
	tar cvzf $@ ./makefile ./param* ./compile* ./recompile* ./README.txt \
	$(skkr-src) \
	$(RSPSCFPATH)skkr.f
	@echo " "
	@echo "done"
#
