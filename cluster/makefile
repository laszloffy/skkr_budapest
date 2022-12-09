##################################################################
#                                                                #
# makefile for the SKKR programs                                 #
#                                                                #
# version 22.07.2003 by Laszlo Szunyogh                          #
#         added: FTNCHEK and SPAG                                #
#                                                                #
# usage                                                          #
#         make nr-scf.run         => nonrel. nonspinpol. scf     #
#         make nrsp-scf.run       => nonrel. spinpol. scf        #
#    rel. spin+orb. pol. scf versions                            #
#         make rsp-scf_ksave.run  => small memory, works for any #
#                                    cae, optimized for bulk     #
#         make rsp-scf_save0.run  => bigger memory, fast for     #
#                                    surface calculations, but,  #
#                       ! does not work for spinpol. substrate ! #
#         make rsp-scf_save.run   => biggest memory, fast for    #
#                                    surface calculations,       # 
#                                    works for any case          #
#                                                                #
# remarks                                                        #
#                                                                #
#   (1)  using the SPAG source utility the steps are:            #
#        1. make spag-main                                       # 
#        2. make spag-nonrel                                     # 
#        3. make spag-nr-scf                                     # 
#        4. make spag-rel                                        # 
#        5. make spag-spr-scf                                    # 
#                                                                # 
#        because from the 'spag' standard output is no way       #
#                back to the makefile                            #
#                                                                # 
##################################################################
#
# fortran compiler:
#
# gnu 
# FC = f77 
# ifc
FC = ifort
FCPAR = /opt/ohpc/pub/mpi/openmpi4-gnu9/4.0.4/bin/mpifort
#FC = mpif90
#FC = /opt/intel/compiler81/bin/ifort
#
# compiler options:
#
# gnu
# alpha 
# FFLAGS = -fast -static -extend_source -warn noinformational -non_shared -C
# athlon
# FFLAGS = -O2 -funroll-loops -fno-automatic -malign-double -ffast-math -ffixed-line-length-132 
#
# ifc
# athlon
# FFLAGS = -O2 -w -tpp5 -unroll -132 -xM
# FFLAGS = -O2 -w -tpp5 -unroll -xM
# p4 
FFLAGS = -xSSE4.2 -fpp -O2 -w -traceback -unroll -mcmodel=medium -shared-intel
# ifc70
#FFLAGS = -O2 -w -unroll -132 -traceback
# ifc81
#FFLAGS = -O2 -w -arch pn4 -axW -unroll -132
#
# linker options, if the compiler links:
# 
# gnu
# alpha   
# LFLAGS = -fast -static -extend_source -warn noinformational -non_shared
# athlon 
# LFLAGS = -O2 -funroll-loops -fno-automatic -malign-double -ffast-math -ffixed-line-length-132  -static
#
# ifc
# Budapest
# LFLAGS = -Vaxlib -lsvml
# Vienna
# ifc70
#LFLAGS = -Vaxlib -lsvml
#MKL=/opt/intel/mkl90/lib/32/
#LFLAGS = -Vaxlib -lsvml $(MKL)libmkl_lapack.a $(MKL)libmkl_ia32.a $(MKL)libguide.a -lpthread
# ifc81
#LFLAGS = -Vaxlib -lsvml -i_dynamic -mcmodel=medium -shared-intel
LFLAGS=$(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a \
     -Wl,--start-group \
     $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
     $(MKLROOT)/lib/intel64/libmkl_sequential.a \
     $(MKLROOT)/lib/intel64/libmkl_core.a \
     $(MKLROOT)/lib/intel64/libmkl_blacs_sgimpt_lp64.a -Wl,--end-group \
     -lpthread -lm
#
# LFLAGS = -Vaxlib -static -L/usr/lib/ATLAS -lf77blas -latlas
#
#        -fast -static -extend_source -warn noinformational -non_shared
#
#         -improves run-time performance 
#          with =assume noaccuracy_sensitive, =align dcommons,
#               =math_library fast, and =O4
#               -all local variables are statically allocated
#                        - statements end in column 132
#                                       -suppresses all informational messages 
#                                                            -object files from archives are included
#
# check utility
#
FU = ftnchek
#
# check options
#
FUFLAGS = -calltree=13 
#FUFLAGS = -calltree=9 -crossref -nodeclare -nodivision -noportabil -nolist  
#         -calling tree=sum of below numbersi, e.g. 9 or 13 
#             0:  do not produce any call graph.
#             1:  produce the call graph in tree format.
#             2:  produce the call graph in who-calls-who format (same as -reference switch).
#             3:  produce the call graph in VCG format (same as -vcg switch).
#             4:  do not prune repeated subtrees (applicable only for format 1).
#             8:  do not sort children into alphabetical order.
#                     -cross-reference list
#                                 -undeclared identifiers
#                                            -warn if division occurs
#                                                        -warn about non-portable usages
#                                                                     -listing of program
#
#
FUPRO =  -project 
#        -create project file
# 
# source tools
#
FS = /usr/local/plusFORT/spag
FX = gxchk
#
# tool options in configuration-file
#
PLFCONF = plusFORT.config
#
# working directory
#
KPATH = ./
#
# subdirectories 
#
SUBDIRS = main nonrel nr-scf nrsp-scf nrdlm-scf rel rsp-scf rsp-lay rsp-clu \
	  nrsp-lay
#
# programs
#
PROGRAMS = nr-scf.run nrsp-scf.run nrdlm-scf.run rsp-scf_ksave.run rsp-scf_save0.run \
rsp-scf_save.run rsp-lay_save.run rsp-clu.run nrsp-lay.run \
rsp-lay_save_mpi.run rsp-clu.run_mpi.run
#
# cross referrences
#
CROSSREF = nr-scf.crf nrsp-scf.crf rsp-scf.crf rsp-lay.crf rsp-clu.crf \
	   nrsp-lay.crf
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
$(MAINPATH)genlatt2d.o  $(MAINPATH)genlatt3d.o	\
$(MAINPATH)getg2d.o 	$(MAINPATH)getsfg.o 	 $(MAINPATH)gjinv.o \
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
$(MAINPATH)struc3d.o    $(MAINPATH)sublattpot.o \
$(MAINPATH)tcpaio.o 	$(MAINPATH)tau2d.o	 $(MAINPATH)vecrot.o \
$(MAINPATH)vecsrt2d.o 	$(MAINPATH)vecsrt3d.o \
$(MAINPATH)vector2d.o   $(MAINPATH)vector3d.o \
$(MAINPATH)ylag.o       $(MAINPATH)ylmc.o 	 $(MAINPATH)ytrafo.o \
$(MAINPATH)zero.o       $(MAINPATH)zmesh.o \
$(MAINPATH)clebsch.o    $(MAINPATH)gafill.o      $(MAINPATH)relmtrx.o
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
$(MAINPATH)genlatt2d.f 	$(MAINPATH)genlatt3d.f \
$(MAINPATH)getg2d.f     $(MAINPATH)getsfg.f      $(MAINPATH)gjinv.f \
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
$(MAINPATH)struc3d.f	 $(MAINPATH)sublattpot.f \
$(MAINPATH)tcpaio.f	$(MAINPATH)tau2d.f	 $(MAINPATH)vecrot.f \
$(MAINPATH)vecsrt2d.f 	$(MAINPATH)vecsrt3d.f \
$(MAINPATH)vector2d.f   $(MAINPATH)vector3d.f \
$(MAINPATH)ylag.f       $(MAINPATH)ylmc.f	 $(MAINPATH)ytrafo.f \
$(MAINPATH)zero.f       $(MAINPATH)zmesh.f \
$(MAINPATH)clebsch.f    $(MAINPATH)gafill.f      $(MAINPATH)relmtrx.f
#
MAINPRJ = \
$(MAINPATH)alloccheck.prj \
$(MAINPATH)alphamat.prj   $(MAINPATH)averagebpot.prj $(MAINPATH)blm.prj \
$(MAINPATH)bulkmad1.prj   $(MAINPATH)bulkmad2.prj \
$(MAINPATH)broyd.prj      $(MAINPATH)cerf.prj        $(MAINPATH)chcore.prj \
$(MAINPATH)coulint.prj    $(MAINPATH)csbf.prj        $(MAINPATH)delkam.prj \
$(MAINPATH)ibz.prj        $(MAINPATH)directg.prj     $(MAINPATH)dugo.prj \
$(MAINPATH)erf.prj        $(MAINPATH)error.prj       $(MAINPATH)exchange.prj \
$(MAINPATH)fac.prj        $(MAINPATH)findangles.prj  $(MAINPATH)fullbz.prj \
$(MAINPATH)gaunt.prj      $(MAINPATH)gauntcalc.prj \
$(MAINPATH)genlatt2d.prj  $(MAINPATH)genlatt3d.prj \
$(MAINPATH)getg2d.prj     $(MAINPATH)getsfg.prj      $(MAINPATH)gjinv.prj \
$(MAINPATH)gstore.prj     $(MAINPATH)ilkam.prj       $(MAINPATH)incgamma.prj \
$(MAINPATH)indekx.prj     $(MAINPATH)init.prj        $(MAINPATH)inouh.prj \
$(MAINPATH)inth.prj       $(MAINPATH)jobnam.prj      $(MAINPATH)kamstr.prj \
$(MAINPATH)kmesh2d.prj \
$(MAINPATH)latgen2d.prj   $(MAINPATH)latgen3d.prj \
$(MAINPATH)legpol.prj \
$(MAINPATH)legzero.prj    $(MAINPATH)ll.prj          $(MAINPATH)lmfill.prj \
$(MAINPATH)madcheck.prj \
$(MAINPATH)madelungb.prj  $(MAINPATH)madelungi.prj   $(MAINPATH)madelungv.prj \
$(MAINPATH)matops.prj     $(MAINPATH)newfl.prj       $(MAINPATH)packer.prj \
$(MAINPATH)phi.prj        $(MAINPATH)poisson.prj     $(MAINPATH)primkr2d.prj \
$(MAINPATH)rcore.prj      $(MAINPATH)resld.prj       $(MAINPATH)rgroup.prj \
$(MAINPATH)rotmat.prj     $(MAINPATH)rotmat2.prj \
$(MAINPATH)rotvec.prj     $(MAINPATH)rsimp.prj \
$(MAINPATH)scrinv.prj     $(MAINPATH)scrstr.prj      $(MAINPATH)setspk.prj \
$(MAINPATH)signum.prj     $(MAINPATH)speck.prj       $(MAINPATH)speck2.prj \
$(MAINPATH)spherh.prj     $(MAINPATH)sphrm.prj       $(MAINPATH)surfgf.prj \
$(MAINPATH)struc3d.prj    $(MAINPATH)sublattpot.prj \
$(MAINPATH)tcpaio.prj	  $(MAINPATH)tau2d.prj       $(MAINPATH)vecrot.prj \
$(MAINPATH)vecsrt2d.prj   $(MAINPATH)vecsrt3d.prj \
$(MAINPATH)vector2d.prj   $(MAINPATH)vector3d.prj \
$(MAINPATH)ylag.prj       $(MAINPATH)ylmc.prj	     $(MAINPATH)ytrafo.prj \
$(MAINPATH)zero.prj       $(MAINPATH)zmesh.prj \
$(MAINPATH)clebsch.prj    $(MAINPATH)gafill.prj      $(MAINPATH)relmtrx.prj
#
MAINSPG = \
$(MAINPATH)alloccheck.spg \
$(MAINPATH)alphamat.spg   $(MAINPATH)averagebpot.spg $(MAINPATH)blm.spg \
$(MAINPATH)bulkmad1.spg   $(MAINPATH)bulkmad2.spg \
$(MAINPATH)broyd.spg      $(MAINPATH)cerf.spg        $(MAINPATH)chcore.spg \
$(MAINPATH)coulint.spg    $(MAINPATH)csbf.spg        $(MAINPATH)delkam.spg \
$(MAINPATH)ibz.spg        $(MAINPATH)directg.spg     $(MAINPATH)dugo.spg \
$(MAINPATH)erf.spg        $(MAINPATH)error.spg       $(MAINPATH)exchange.spg \
$(MAINPATH)fac.spg        $(MAINPATH)findangles.spg  $(MAINPATH)fullbz.spg \
$(MAINPATH)gaunt.spg      $(MAINPATH)gauntcalc.spg \
$(MAINPATH)genlatt2d.spg  $(MAINPATH)genlatt3d.spg \
$(MAINPATH)getg2d.spg     $(MAINPATH)getsfg.spg      $(MAINPATH)gjinv.spg \
$(MAINPATH)gstore.spg     $(MAINPATH)ilkam.spg       $(MAINPATH)incgamma.spg \
$(MAINPATH)indekx.spg     $(MAINPATH)init.spg        $(MAINPATH)inouh.spg \
$(MAINPATH)inth.spg       $(MAINPATH)jobnam.spg      $(MAINPATH)kamstr.spg \
$(MAINPATH)kmesh2d.spg \
$(MAINPATH)latgen2d.spg   $(MAINPATH)latgen3d.spg \
$(MAINPATH)legpol.spg \
$(MAINPATH)legzero.spg    $(MAINPATH)ll.spg          $(MAINPATH)lmfill.spg \
$(MAINPATH)madcheck.spg \
$(MAINPATH)madelungb.spg  $(MAINPATH)madelungi.spg   $(MAINPATH)madelungv.spg \
$(MAINPATH)matops.spg     $(MAINPATH)newfl.spg       $(MAINPATH)packer.spg \
$(MAINPATH)phi.spg        $(MAINPATH)poisson.spg     $(MAINPATH)primkr2d.spg \
$(MAINPATH)rcore.spg      $(MAINPATH)resld.spg       $(MAINPATH)rgroup.spg \
$(MAINPATH)rotmat.spg     $(MAINPATH)rotmat2.spg \
$(MAINPATH)rotvec.spg     $(MAINPATH)rsimp.spg \
$(MAINPATH)scrinv.spg     $(MAINPATH)scrstr.spg      $(MAINPATH)setspk.spg \
$(MAINPATH)signum.spg     $(MAINPATH)speck.spg       $(MAINPATH)speck2.spg \
$(MAINPATH)spherh.spg     $(MAINPATH)sphrm.spg       $(MAINPATH)surfgf.spg \
$(MAINPATH)struc3d.spg    $(MAINPATH)sublattpot.spg \
$(MAINPATH)tcpaio.spg     $(MAINPATH)tau2d.spg      $(MAINPATH)vecrot.spg \
$(MAINPATH)vecsrt2d.spg   $(MAINPATH)vecsrt3d.spg \
$(MAINPATH)vector2d.spg   $(MAINPATH)vector3d.spg \
$(MAINPATH)ylag.spg       $(MAINPATH)ylmc.spg	     $(MAINPATH)ytrafo.spg \
$(MAINPATH)zero.spg       $(MAINPATH)zmesh.spg \
$(MAINPATH)clebsch.spg    $(MAINPATH)gafill.spg      $(MAINPATH)relmtrx.spg
#
MAINSMB = \
$(MAINPATH)alloccheck.smb \
$(MAINPATH)alphamat.smb   $(MAINPATH)averagebpot.smb $(MAINPATH)blm.smb \
$(MAINPATH)bulkmad1.smb   $(MAINPATH)bulkmad2.smb \
$(MAINPATH)broyd.smb      $(MAINPATH)cerf.smb        $(MAINPATH)chcore.smb \
$(MAINPATH)coulint.smb    $(MAINPATH)csbf.smb        $(MAINPATH)delkam.smb \
$(MAINPATH)ibz.smb        $(MAINPATH)directg.smb     $(MAINPATH)dugo.smb \
$(MAINPATH)erf.smb        $(MAINPATH)error.smb       $(MAINPATH)exchange.smb \
$(MAINPATH)fac.smb        $(MAINPATH)findangles.smb  $(MAINPATH)fullbz.smb \
$(MAINPATH)gaunt.smb      $(MAINPATH)gauntcalc.smb \
$(MAINPATH)genlatt2d.smb  $(MAINPATH)genlatt3d.smb \
$(MAINPATH)getg2d.smb     $(MAINPATH)getsfg.smb      $(MAINPATH)gjinv.smb \
$(MAINPATH)gstore.smb     $(MAINPATH)ilkam.smb       $(MAINPATH)incgamma.smb \
$(MAINPATH)indekx.smb     $(MAINPATH)init.smb        $(MAINPATH)inouh.smb \
$(MAINPATH)inth.smb       $(MAINPATH)jobnam.smb      $(MAINPATH)kamstr.smb \
$(MAINPATH)kmesh2d.smb \
$(MAINPATH)latgen2d.smb   $(MAINPATH)latgen3d.smb \
$(MAINPATH)legpol.smb \
$(MAINPATH)legzero.smb    $(MAINPATH)ll.smb          $(MAINPATH)lmfill.smb \
$(MAINPATH)madcheck.smb \
$(MAINPATH)madelungb.smb  $(MAINPATH)madelungi.smb   $(MAINPATH)madelungv.smb \
$(MAINPATH)matops.smb     $(MAINPATH)newfl.smb       $(MAINPATH)packer.smb \
$(MAINPATH)phi.smb        $(MAINPATH)poisson.smb     $(MAINPATH)primkr2d.smb \
$(MAINPATH)rcore.smb      $(MAINPATH)resld.smb       $(MAINPATH)rgroup.smb \
$(MAINPATH)rotmat.smb     $(MAINPATH)rotmat2.smb \
$(MAINPATH)rotvec.smb     $(MAINPATH)rsimp.smb \
$(MAINPATH)scrinv.smb     $(MAINPATH)scrstr.smb      $(MAINPATH)setspk.smb \
$(MAINPATH)signum.smb     $(MAINPATH)speck.smb       $(MAINPATH)speck2.smb \
$(MAINPATH)spherh.smb     $(MAINPATH)sphrm.smb       $(MAINPATH)surfgf.smb \
$(MAINPATH)struc3d.smb    $(MAINPATH)sublattpot.smb \
$(MAINPATH)tcpaio.smb     $(MAINPATH)tau2d.smb      $(MAINPATH)vecrot.smb \
$(MAINPATH)vecsrt2d.smb   $(MAINPATH)vecsrt3d.smb \
$(MAINPATH)vector2d.smb   $(MAINPATH)vector3d.smb \
$(MAINPATH)ylag.smb       $(MAINPATH)ylmc.smb	     $(MAINPATH)ytrafo.smb \
$(MAINPATH)zero.smb       $(MAINPATH)zmesh.smb \
$(MAINPATH)clebsch.smb    $(MAINPATH)gafill.smb      $(MAINPATH)relmtrx.smb
#
NONRELPATH = $(KPATH)nonrel/
NONRELOBJ = \
$(NONRELPATH)fullcpa.o   $(NONRELPATH)initia.o \
$(NONRELPATH)initzero.o  $(NONRELPATH)locquant.o  $(NONRELPATH)matrep.o \
$(NONRELPATH)readpot.o   $(NONRELPATH)phystau.o   $(NONRELPATH)schreq.o \
$(NONRELPATH)sphbas.o    $(NONRELPATH)tmat.o      $(NONRELPATH)tmatini.o \
$(NONRELPATH)wafu.o      $(NONRELPATH)moment.o    $(NONRELPATH)dens.o
#
NONRELSRC = \
$(NONRELPATH)fullcpa.f   $(NONRELPATH)initia.f \
$(NONRELPATH)initzero.f  $(NONRELPATH)locquant.f  $(NONRELPATH)matrep.f \
$(NONRELPATH)readpot.f   $(NONRELPATH)phystau.f   $(NONRELPATH)schreq.f \
$(NONRELPATH)sphbas.f    $(NONRELPATH)tmat.f      $(NONRELPATH)tmatini.f \
$(NONRELPATH)wafu.f      $(NONRELPATH)moment.f    $(NONRELPATH)dens.f
#
NONRELPRJ = \
$(NONRELPATH)fullcpa.prj   $(NONRELPATH)initia.prj \
$(NONRELPATH)initzero.prj  $(NONRELPATH)locquant.prj  $(NONRELPATH)matrep.prj \
$(NONRELPATH)readpot.prj   $(NONRELPATH)phystau.prj   $(NONRELPATH)schreq.prj \
$(NONRELPATH)sphbas.prj    $(NONRELPATH)tmat.prj      $(NONRELPATH)tmatini.prj \
$(NONRELPATH)wafu.prj      $(NONRELPATH)moment.prj    $(NONRELPATH)dens.prj
#
NONRELSPG = \
$(NONRELPATH)fullcpa.spg   $(NONRELPATH)initia.spg \
$(NONRELPATH)initzero.spg  $(NONRELPATH)locquant.spg  $(NONRELPATH)matrep.spg \
$(NONRELPATH)readpot.spg   $(NONRELPATH)phystau.spg   $(NONRELPATH)schreq.spg \
$(NONRELPATH)sphbas.spg    $(NONRELPATH)tmat.spg      $(NONRELPATH)tmatini.spg \
$(NONRELPATH)wafu.spg      $(NONRELPATH)moment.spg    $(NONRELPATH)dens.spg
#
NONRELSMB = \
$(NONRELPATH)fullcpa.smb   $(NONRELPATH)initia.smb \
$(NONRELPATH)initzero.smb  $(NONRELPATH)locquant.smb  $(NONRELPATH)matrep.smb \
$(NONRELPATH)readpot.smb   $(NONRELPATH)phystau.smb   $(NONRELPATH)schreq.smb \
$(NONRELPATH)sphbas.smb    $(NONRELPATH)tmat.smb      $(NONRELPATH)tmatini.smb \
$(NONRELPATH)wafu.smb      $(NONRELPATH)moment.smb    $(NONRELPATH)dens.smb
#
NRSCFPATH = $(KPATH)nr-scf/
NRSCFOBJ = \
$(NRSCFPATH)enprint.o \
$(NRSCFPATH)newpot.o   $(NRSCFPATH)pothandle.o  $(NRSCFPATH)printres.o \
$(NRSCFPATH)printscf.o $(NRSCFPATH)readgeom.o   $(NRSCFPATH)readini.o  \
$(NRSCFPATH)restart.o  $(NRSCFPATH)vgen.o
#
NRSCFSRC = \
$(NRSCFPATH)enprint.f \
$(NRSCFPATH)newpot.f   $(NRSCFPATH)pothandle.f  $(NRSCFPATH)printres.f \
$(NRSCFPATH)printscf.f $(NRSCFPATH)readgeom.f   $(NRSCFPATH)readini.f  \
$(NRSCFPATH)restart.f  $(NRSCFPATH)vgen.f
#
NRSCFPRJ = \
$(NRSCFPATH)enprint.prj \
$(NRSCFPATH)newpot.prj   $(NRSCFPATH)pothandle.prj  $(NRSCFPATH)printres.prj \
$(NRSCFPATH)printscf.prj $(NRSCFPATH)readgeom.prj   $(NRSCFPATH)readini.prj  \
$(NRSCFPATH)restart.prj  $(NRSCFPATH)vgen.prj
#
NRSCFSPG = \
$(NRSCFPATH)enprint.spg \
$(NRSCFPATH)newpot.spg   $(NRSCFPATH)pothandle.spg  $(NRSCFPATH)printres.spg \
$(NRSCFPATH)printscf.spg $(NRSCFPATH)readgeom.spg   $(NRSCFPATH)readini.spg  \
$(NRSCFPATH)restart.spg  $(NRSCFPATH)vgen.spg
#
NRSCFSMB = \
$(NRSCFPATH)enprint.smb \
$(NRSCFPATH)newpot.smb   $(NRSCFPATH)pothandle.smb  $(NRSCFPATH)printres.smb \
$(NRSCFPATH)printscf.smb $(NRSCFPATH)readgeom.smb   $(NRSCFPATH)readini.smb  \
$(NRSCFPATH)restart.smb  $(NRSCFPATH)vgen.smb
#
NRSPSCFPATH = $(KPATH)nrsp-scf/
NRSPSCFOBJ = \
$(NRSPSCFPATH)enprint.o \
$(NRSPSCFPATH)newpot.o   $(NRSPSCFPATH)pothandle.o  $(NRSPSCFPATH)printscf.o\
$(NRSPSCFPATH)printres.o $(NRSPSCFPATH)readini.o    $(NRSPSCFPATH)readgeom.o \
$(NRSPSCFPATH)restart.o  $(NRSPSCFPATH)vgen.o
#
NRSPSCFSRC = \
$(NRSPSCFPATH)enprint.f \
$(NRSPSCFPATH)newpot.f   $(NRSPSCFPATH)pothandle.f  $(NRSPSCFPATH)printscf.f\
$(NRSPSCFPATH)printres.f $(NRSPSCFPATH)readini.f    $(NRSPSCFPATH)readgeom.f \
$(NRSPSCFPATH)restart.f  $(NRSPSCFPATH)vgen.f
#
NRSPSCFPRJ = \
$(NRSPSCFPATH)enprint.prj \
$(NRSPSCFPATH)newpot.prj   $(NRSPSCFPATH)pothandle.prj  $(NRSPSCFPATH)printscf.prj\
$(NRSPSCFPATH)printres.prj $(NRSPSCFPATH)readini.prj    $(NRSPSCFPATH)readgeom.prj \
$(NRSPSCFPATH)restart.prj  $(NRSPSCFPATH)vgen.prj
#
NRSPSCFSPG = \
$(NRSPSCFPATH)enprint.spg \
$(NRSPSCFPATH)newpot.spg   $(NRSPSCFPATH)pothandle.spg  $(NRSPSCFPATH)printscf.spg\
$(NRSPSCFPATH)printres.spg $(NRSPSCFPATH)readini.spg    $(NRSPSCFPATH)readgeom.spg \
$(NRSPSCFPATH)restart.spg  $(NRSPSCFPATH)vgen.spg
#
NRSPSCFSMB = \
$(NRSPSCFPATH)enprint.smb \
$(NRSPSCFPATH)newpot.smb   $(NRSPSCFPATH)pothandle.smb  $(NRSPSCFPATH)printscf.smb\
$(NRSPSCFPATH)printres.smb $(NRSPSCFPATH)readini.smb    $(NRSPSCFPATH)readgeom.smb \
$(NRSPSCFPATH)restart.smb  $(NRSPSCFPATH)vgen.smb
#
NRDLMSCFPATH = $(KPATH)nrdlm-scf/
NRDLMSCFOBJ = \
$(NRDLMSCFPATH)enprint.o \
$(NRDLMSCFPATH)newpot.o   $(NRDLMSCFPATH)pothandle.o  $(NRDLMSCFPATH)printscf.o\
$(NRDLMSCFPATH)printres.o $(NRDLMSCFPATH)readini.o    $(NRDLMSCFPATH)readgeom.o \
$(NRDLMSCFPATH)restart.o  $(NRDLMSCFPATH)vgen.o \
$(NRDLMSCFPATH)cpacoord.o  $(NRDLMSCFPATH)fullcpa.o   $(NRDLMSCFPATH)initia.o \
$(NRDLMSCFPATH)initzero.o  $(NRDLMSCFPATH)locquant.o  $(NRDLMSCFPATH)matrep.o \
$(NRDLMSCFPATH)readpot.o   $(NRDLMSCFPATH)phystau.o   $(NRDLMSCFPATH)schreq.o \
$(NRDLMSCFPATH)sphbas.o    $(NRDLMSCFPATH)tmat.o      $(NRDLMSCFPATH)tmatini.o \
$(NRDLMSCFPATH)wafu.o      $(NRDLMSCFPATH)moment.o    $(NRDLMSCFPATH)dens.o
#
NRDLMSCFSRC = \
$(NRDLMSCFPATH)enprint.f \
$(NRDLMSCFPATH)newpot.f   $(NRDLMSCFPATH)pothandle.f  $(NRDLMSCFPATH)printscf.f\
$(NRDLMSCFPATH)printres.f $(NRDLMSCFPATH)readini.f    $(NRDLMSCFPATH)readgeom.f \
$(NRDLMSCFPATH)restart.f  $(NRDLMSCFPATH)vgen.f \
$(NRDLMSCFPATH)cpacoord.f  $(NRDLMSCFPATH)fullcpa.f   $(NRDLMSCFPATH)initia.f \
$(NRDLMSCFPATH)initzero.f  $(NRDLMSCFPATH)locquant.f  $(NRDLMSCFPATH)matrep.f \
$(NRDLMSCFPATH)readpot.f   $(NRDLMSCFPATH)phystau.f   $(NRDLMSCFPATH)schreq.f \
$(NRDLMSCFPATH)sphbas.f    $(NRDLMSCFPATH)tmat.f      $(NRDLMSCFPATH)tmatini.f \
$(NRDLMSCFPATH)wafu.f      $(NRDLMSCFPATH)moment.f    $(NRDLMSCFPATH)dens.f
#
RELPATH = $(KPATH)rel/
RELOBJ = \
$(RELPATH)bop.o         $(RELPATH)brmat.o \
$(RELPATH)enblloyd.o    $(RELPATH)denblloyd.o \
$(RELPATH)dirmag1-op.o  $(RELPATH)dirmag2-op.o \
$(RELPATH)gfill.o       $(RELPATH)initia.o      $(RELPATH)matlms.o \
$(RELPATH)magnet.o      $(RELPATH)matr.o        $(RELPATH)matr2.o \
$(RELPATH)matrep.o      $(RELPATH)matrot.o      $(RELPATH)matrot2.o \
$(RELPATH)moment.o      $(RELPATH)phystau.o \
$(RELPATH)readpot.o     $(RELPATH)replms.o 	$(RELPATH)replmsf.o \
$(RELPATH)sphbas.o      $(RELPATH)spzwafu.o 	$(RELPATH)matnr.o \
$(RELPATH)tinvmat.o     $(RELPATH)tmat.o        $(RELPATH)wafu.o
#
RELSRC = \
$(RELPATH)bop.f         $(RELPATH)brmat.f \
$(RELPATH)enblloyd.f    $(RELPATH)denblloyd.f \
$(RELPATH)dirmag1-op.f  $(RELPATH)dirmag2-op.f \
$(RELPATH)gfill.f       $(RELPATH)initia.f      $(RELPATH)matlms.f \
$(RELPATH)magnet.f      $(RELPATH)matr.f        $(RELPATH)matr2.f \
$(RELPATH)matrep.f      $(RELPATH)matrot.f      $(RELPATH)matrot2.f \
$(RELPATH)moment.f      $(RELPATH)phystau.f \
$(RELPATH)readpot.f     $(RELPATH)replms.f 	$(RELPATH)replmsf.f \
$(RELPATH)sphbas.f      $(RELPATH)spzwafu.f 	$(RELPATH)matnr.f \
$(RELPATH)tinvmat.f     $(RELPATH)tmat.f        $(RELPATH)wafu.f
#
RELLAYOBJ = \
$(RELPATH)fullcpa.o
#
RELLAYSRC = \
$(RELPATH)fullcpa.f
#
RELPRJ = \
$(RELPATH)bop.prj         $(RELPATH)brmat.prj \
$(RELPATH)enblloyd.prj    $(RELPATH)denblloyd.prj \
$(RELPATH)dirmag1-op.prj  $(RELPATH)dirmag2-op.prj  $(RELPATH)fullcpa.prj \
$(RELPATH)gfill.prj       $(RELPATH)initia.prj      $(RELPATH)matlms.prj \
$(RELPATH)magnet.prj      $(RELPATH)matr.prj        $(RELPATH)matr2.prj \
$(RELPATH)matrep.prj      $(RELPATH)matrot.prj      $(RELPATH)matrot2.prj \
$(RELPATH)moment.prj      $(RELPATH)phystau.prj \
$(RELPATH)readpot.prj     $(RELPATH)replms.prj 	    $(RELPATH)replmsf.prj \
$(RELPATH)sphbas.prj      $(RELPATH)spzwafu.prj     $(RELPATH)matnr.prj \
$(RELPATH)tinvmat.prj     $(RELPATH)tmat.prj        $(RELPATH)wafu.prj
#
RELSPG = \
$(RELPATH)bop.spg         $(RELPATH)brmat.spg \
$(RELPATH)enblloyd.spg    $(RELPATH)denblloyd.spg \
$(RELPATH)dirmag1-op.spg  $(RELPATH)dirmag2-op.spg  $(RELPATH)fullcpa.spg \
$(RELPATH)gfill.spg       $(RELPATH)initia.spg      $(RELPATH)matlms.spg \
$(RELPATH)magnet.spg      $(RELPATH)matr.spg        $(RELPATH)matr2.spg \
$(RELPATH)matrep.spg      $(RELPATH)matrot.spg      $(RELPATH)matrot2.spg \
$(RELPATH)moment.spg      $(RELPATH)phystau.spg \
$(RELPATH)readpot.spg     $(RELPATH)replms.spg      $(RELPATH)replmsf.spg \
$(RELPATH)sphbas.spg      $(RELPATH)spzwafu.spg     $(RELPATH)matnr.spg \
$(RELPATH)tinvmat.spg     $(RELPATH)tmat.spg        $(RELPATH)wafu.spg
#
RELSMB = \
$(RELPATH)bop.smb         $(RELPATH)brmat.smb \
$(RELPATH)enblloyd.smb    $(RELPATH)denblloyd.smb \
$(RELPATH)dirmag1-op.smb  $(RELPATH)dirmag2-op.smb  $(RELPATH)fullcpa.smb \
$(RELPATH)gfill.smb       $(RELPATH)initia.smb      $(RELPATH)matlms.smb \
$(RELPATH)magnet.smb      $(RELPATH)matr.smb        $(RELPATH)matr2.smb \
$(RELPATH)matrep.smb      $(RELPATH)matrot.smb      $(RELPATH)matrot2.smb \
$(RELPATH)moment.smb      $(RELPATH)phystau.smb \
$(RELPATH)readpot.smb     $(RELPATH)replms.smb      $(RELPATH)replmsf.smb \
$(RELPATH)sphbas.smb      $(RELPATH)spzwafu.smb     $(RELPATH)matnr.smb \
$(RELPATH)tinvmat.smb     $(RELPATH)tmat.smb        $(RELPATH)wafu.smb
#
RSPSCFPATH = $(KPATH)rsp-scf/
RSPSCF-ksave-OBJ = \
$(RSPSCFPATH)cpacoord_ksave.o \
$(RSPSCFPATH)dens.o        $(RSPSCFPATH)enprint.o \
$(RSPSCFPATH)initzero.o    $(RSPSCFPATH)localrot2.o   $(RSPSCFPATH)locquant.o \
$(RSPSCFPATH)newdir.o      $(RSPSCFPATH)newpot.o      $(RSPSCFPATH)pothandle.o \
$(RSPSCFPATH)printres.o    $(RSPSCFPATH)printscf.o    $(RSPSCFPATH)readgeom.o \
$(RSPSCFPATH)readini.o     $(RSPSCFPATH)restart.o     $(RSPSCFPATH)sorttmat.o \
$(RSPSCFPATH)tmatini.o     $(RSPSCFPATH)vgen.o
#
RSPSCF-save0-OBJ = \
$(RSPSCFPATH)cpacoord_save0.o \
$(RSPSCFPATH)dens.o        $(RSPSCFPATH)enprint.o \
$(RSPSCFPATH)initzero.o    $(RSPSCFPATH)localrot2.o   $(RSPSCFPATH)locquant.o \
$(RSPSCFPATH)newdir.o      $(RSPSCFPATH)newpot.o      $(RSPSCFPATH)pothandle.o \
$(RSPSCFPATH)printres.o    $(RSPSCFPATH)printscf.o    $(RSPSCFPATH)readgeom.o \
$(RSPSCFPATH)readini.o     $(RSPSCFPATH)restart.o     $(RSPSCFPATH)sorttmat.o \
$(RSPSCFPATH)tmatini.o     $(RSPSCFPATH)vgen.o
#
RSPSCF-save-OBJ = \
$(RSPSCFPATH)cpacoord_save.o \
$(RSPSCFPATH)dens.o        $(RSPSCFPATH)enprint.o \
$(RSPSCFPATH)initzero.o    $(RSPSCFPATH)localrot2.o   $(RSPSCFPATH)locquant.o \
$(RSPSCFPATH)newdir.o      $(RSPSCFPATH)newpot.o      $(RSPSCFPATH)pothandle.o \
$(RSPSCFPATH)printres.o    $(RSPSCFPATH)printscf.o    $(RSPSCFPATH)readgeom.o \
$(RSPSCFPATH)readini.o     $(RSPSCFPATH)restart.o     $(RSPSCFPATH)sorttmat.o \
$(RSPSCFPATH)tmatini.o     $(RSPSCFPATH)vgen.o
#
RSPSCFSRC = \
$(RSPSCFPATH)cpacoord_ksave.f	$(RSPSCFPATH)cpacoord_save0.f \
$(RSPSCFPATH)cpacoord_save.f \
$(RSPSCFPATH)dens.f	   $(RSPSCFPATH)enprint.f \
$(RSPSCFPATH)initzero.f    $(RSPSCFPATH)localrot2.f   $(RSPSCFPATH)locquant.f \
$(RSPSCFPATH)newdir.f      $(RSPSCFPATH)newpot.f      $(RSPSCFPATH)pothandle.f \
$(RSPSCFPATH)printres.f    $(RSPSCFPATH)printscf.f    $(RSPSCFPATH)readgeom.f \
$(RSPSCFPATH)readini.f     $(RSPSCFPATH)restart.f     $(RSPSCFPATH)sorttmat.f \
$(RSPSCFPATH)tmatini.f     $(RSPSCFPATH)vgen.f
#
RSPSCFPRJ = \
$(RSPSCFPATH)cpacoord.prj  $(RSPSCFPATH)dens.prj      $(RSPSCFPATH)enprint.prj \
$(RSPSCFPATH)initzero.prj  $(RSPSCFPATH)localrot2.prj $(RSPSCFPATH)locquant.prj \
$(RSPSCFPATH)newdir.prj    $(RSPSCFPATH)newpot.prj    $(RSPSCFPATH)pothandle.prj \
$(RSPSCFPATH)printres.prj  $(RSPSCFPATH)printscf.prj  $(RSPSCFPATH)readgeom.prj \
$(RSPSCFPATH)readini.prj   $(RSPSCFPATH)restart.prj   $(RSPSCFPATH)sorttmat.prj \
$(RSPSCFPATH)tmatini.prj   $(RSPSCFPATH)vgen.prj
#
RSPSCFSPG = \
$(RSPSCFPATH)cpacoord.spg  $(RSPSCFPATH)dens.spg      $(RSPSCFPATH)enprint.spg \
$(RSPSCFPATH)initzero.spg  $(RSPSCFPATH)localrot2.spg $(RSPSCFPATH)locquant.spg \
$(RSPSCFPATH)newdir.spg    $(RSPSCFPATH)newpot.spg    $(RSPSCFPATH)pothandle.spg \
$(RSPSCFPATH)printres.spg  $(RSPSCFPATH)printscf.spg  $(RSPSCFPATH)readgeom.spg \
$(RSPSCFPATH)readini.spg   $(RSPSCFPATH)restart.spg   $(RSPSCFPATH)sorttmat.spg \
$(RSPSCFPATH)tmatini.spg   $(RSPSCFPATH)vgen.spg
#
RSPSCFSMB = \
$(RSPSCFPATH)cpacoord.smb  $(RSPSCFPATH)dens.smb      $(RSPSCFPATH)enprint.smb \
$(RSPSCFPATH)initzero.smb  $(RSPSCFPATH)localrot2.smb $(RSPSCFPATH)locquant.smb \
$(RSPSCFPATH)newdir.smb    $(RSPSCFPATH)newpot.smb    $(RSPSCFPATH)pothandle.smb \
$(RSPSCFPATH)printres.smb  $(RSPSCFPATH)printscf.smb  $(RSPSCFPATH)readgeom.smb \
$(RSPSCFPATH)readini.smb   $(RSPSCFPATH)restart.smb   $(RSPSCFPATH)sorttmat.smb \
$(RSPSCFPATH)tmatini.smb   $(RSPSCFPATH)vgen.smb
#
RSPLAYPATH = $(KPATH)rsp-lay/
RSPLAY-ksave-OBJ = \
$(RSPLAYPATH)cpacoord_ksave.o \
$(RSPLAYPATH)dens.o \
$(RSPLAYPATH)enprint.o \
$(RSPLAYPATH)initzero.o \
$(RSPLAYPATH)localrot2.o \
$(RSPLAYPATH)locquant.o \
$(RSPLAYPATH)madcolay.o \
$(RSPLAYPATH)newdir.o \
$(RSPLAYPATH)newpot.o \
$(RSPLAYPATH)phystau2.o \
$(RSPLAYPATH)pothandle.o \
$(RSPLAYPATH)printres.o \
$(RSPLAYPATH)printscf.o \
$(RSPLAYPATH)readclu.o \
$(RSPLAYPATH)readgeom.o \
$(RSPLAYPATH)readini.o \
$(RSPLAYPATH)restart.o \
$(RSPLAYPATH)sorttmat.o \
$(RSPLAYPATH)tau2d_god.o \
$(RSPLAYPATH)tauhout.o \
$(RSPLAYPATH)thout.o \
$(RSPLAYPATH)tmatini.o \
$(RSPLAYPATH)unsct.o \
$(RSPLAYPATH)vgen.o \
$(RSPLAYPATH)vhout.o
#
RSPLAY-save0-OBJ = \
$(RSPLAYPATH)cpacoord_save0.o \
$(RSPLAYPATH)dens.o \
$(RSPLAYPATH)enprint.o \
$(RSPLAYPATH)initzero.o \
$(RSPLAYPATH)localrot2.o \
$(RSPLAYPATH)locquant.o \
$(RSPLAYPATH)madcolay.o \
$(RSPLAYPATH)newdir.o \
$(RSPLAYPATH)newpot.o \
$(RSPLAYPATH)phystau2.o \
$(RSPLAYPATH)pothandle.o \
$(RSPLAYPATH)printres.o \
$(RSPLAYPATH)printscf.o \
$(RSPLAYPATH)readclu.o \
$(RSPLAYPATH)readgeom.o \
$(RSPLAYPATH)readini.o \
$(RSPLAYPATH)restart.o \
$(RSPLAYPATH)sorttmat.o \
$(RSPLAYPATH)tau2d_god.o \
$(RSPLAYPATH)tauhout.o \
$(RSPLAYPATH)thout.o \
$(RSPLAYPATH)tmatini.o \
$(RSPLAYPATH)unsct.o \
$(RSPLAYPATH)vgen.o \
$(RSPLAYPATH)vhout.o
#
RSPLAY-save-OBJ = \
$(RSPLAYPATH)cpacoord_save.o \
$(RSPLAYPATH)dens.o \
$(RSPLAYPATH)enprint.o \
$(RSPLAYPATH)initzero.o \
$(RSPLAYPATH)localrot2.o \
$(RSPLAYPATH)locquant.o \
$(RSPLAYPATH)madcolay.o \
$(RSPLAYPATH)newdir.o \
$(RSPLAYPATH)newpot.o \
$(RSPLAYPATH)phystau2.o \
$(RSPLAYPATH)pothandle.o \
$(RSPLAYPATH)printres.o \
$(RSPLAYPATH)printscf.o \
$(RSPLAYPATH)readclu.o \
$(RSPLAYPATH)readgeom.o \
$(RSPLAYPATH)readini.o \
$(RSPLAYPATH)restart.o \
$(RSPLAYPATH)sorttmat.o \
$(RSPLAYPATH)tau2d_god.o \
$(RSPLAYPATH)tauhout.o \
$(RSPLAYPATH)thout.o \
$(RSPLAYPATH)tmatini.o \
$(RSPLAYPATH)unsct.o \
$(RSPLAYPATH)vgen.o \
$(RSPLAYPATH)vhout.o
#
RSPLAYSRC = \
$(RSPLAYPATH)cpacoord_ksave.f \
$(RSPLAYPATH)cpacoord_save0.f \
$(RSPLAYPATH)cpacoord_save.f \
$(RSPLAYPATH)dens.f \
$(RSPLAYPATH)enprint.f \
$(RSPLAYPATH)initzero.f \
$(RSPLAYPATH)localrot2.f \
$(RSPLAYPATH)locquant.f \
$(RSPLAYPATH)madcolay.f \
$(RSPLAYPATH)newdir.f \
$(RSPLAYPATH)newpot.f \
$(RSPLAYPATH)phystau2.f \
$(RSPLAYPATH)pothandle.f \
$(RSPLAYPATH)printres.f \
$(RSPLAYPATH)printscf.f \
$(RSPLAYPATH)readclu.f \
$(RSPLAYPATH)readgeom.f \
$(RSPLAYPATH)readini.f \
$(RSPLAYPATH)restart.f \
$(RSPLAYPATH)sorttmat.f \
$(RSPLAYPATH)tau2d_god.f \
$(RSPLAYPATH)tauhout.f \
$(RSPLAYPATH)thout.f \
$(RSPLAYPATH)tmatini.f \
$(RSPLAYPATH)unsct.f \
$(RSPLAYPATH)vgen.f \
$(RSPLAYPATH)vhout.f
#
RSPLAYPRJ = \
$(RSPLAYPATH)cpacoord_save.prj \
$(RSPLAYPATH)dens.prj     \
$(RSPLAYPATH)enprint.prj \
$(RSPLAYPATH)initzero.prj \
$(RSPLAYPATH)localrot2.prj \
$(RSPLAYPATH)locquant.prj \
$(RSPLAYPATH)madcolay.prj \
$(RSPLAYPATH)newdir.prj   \
$(RSPLAYPATH)newpot.prj   \
$(RSPLAYPATH)phystau2.prj \
$(RSPLAYPATH)pothandle.prj \
$(RSPLAYPATH)printres.prj \
$(RSPLAYPATH)printscf.prj \
$(RSPLAYPATH)readclu.prj \
$(RSPLAYPATH)readgeom.prj \
$(RSPLAYPATH)readini.prj  \
$(RSPLAYPATH)restart.prj  \
$(RSPLAYPATH)sorttmat.prj \
$(RSPLAYPATH)tau2d_god.prj \
$(RSPLAYPATH)tauhout.prj \
$(RSPLAYPATH)thout.prj \
$(RSPLAYPATH)tmatini.prj  \
$(RSPLAYPATH)unsct.prj \
$(RSPLAYPATH)vgen.prj \
$(RSPLAYPATH)vgen_bcast.prj \
$(RSPLAYPATH)vhout.prj
#
RSPLAYSPG = \
$(RSPLAYPATH)cpacoord_save.spg \
$(RSPLAYPATH)dens.spg     \
$(RSPLAYPATH)enprint.spg \
$(RSPLAYPATH)initzero.spg \
$(RSPLAYPATH)localrot2.spg \
$(RSPLAYPATH)locquant.spg \
$(RSPLAYPATH)madcolay.spg \
$(RSPLAYPATH)newdir.spg   \
$(RSPLAYPATH)newpot.spg   \
$(RSPLAYPATH)phystau2.spg \
$(RSPLAYPATH)pothandle.spg \
$(RSPLAYPATH)printres.spg \
$(RSPLAYPATH)printscf.spg \
$(RSPLAYPATH)readclu.spg \
$(RSPLAYPATH)readgeom.spg \
$(RSPLAYPATH)readini.spg  \
$(RSPLAYPATH)restart.spg  \
$(RSPLAYPATH)sorttmat.spg \
$(RSPLAYPATH)tau2d_god.spg \
$(RSPLAYPATH)tauhout.spg \
$(RSPLAYPATH)thout.spg \
$(RSPLAYPATH)tmatini.spg  \
$(RSPLAYPATH)unsct.spg \
$(RSPLAYPATH)vgen.spg \
$(RSPLAYPATH)vgen_bcast.spg \
$(RSPLAYPATH)vhout.spg
#
RSPLAYSMB = \
$(RSPLAYPATH)cpacoord_save.smb \
$(RSPLAYPATH)dens.smb     \
$(RSPLAYPATH)enprint.smb \
$(RSPLAYPATH)initzero.smb \
$(RSPLAYPATH)localrot2.smb \
$(RSPLAYPATH)locquant.smb \
$(RSPLAYPATH)madcolay.smb \
$(RSPLAYPATH)newdir.smb    \
$(RSPLAYPATH)newpot.smb  \
$(RSPLAYPATH)phystau2.smb \
$(RSPLAYPATH)pothandle.smb \
$(RSPLAYPATH)printres.smb \
$(RSPLAYPATH)printscf.smb \
$(RSPLAYPATH)readclu.smb \
$(RSPLAYPATH)readgeom.smb \
$(RSPLAYPATH)readini.smb  \
$(RSPLAYPATH)restart.smb  \
$(RSPLAYPATH)sorttmat.smb \
$(RSPLAYPATH)tau2d_god.smb \
$(RSPLAYPATH)tauhout.smb \
$(RSPLAYPATH)thout.smb \
$(RSPLAYPATH)tmatini.smb  \
$(RSPLAYPATH)unsct.smb \
$(RSPLAYPATH)vgen.smb \
$(RSPLAYPATH)vgen_bcast.smb \
$(RSPLAYPATH)vhout.smb
#
MPILAYOBJ = \
$(RSPLAYPATH)vgen_bcast.o
#
MPILAYSRC = \
$(RSPLAYPATH)vgen_bcast.f
#
RSPCLUPATH = $(KPATH)rsp-clu/
RSPCLU-OBJ = \
$(RSPCLUPATH)alphamat1.o    \
$(RSPCLUPATH)deltau.o       \
$(RSPCLUPATH)dens.o \
$(RSPCLUPATH)dtat.o \
$(RSPCLUPATH)ecoreimp.o \
$(RSPCLUPATH)enprint1.o \
$(RSPCLUPATH)gettaucl.o \
$(RSPCLUPATH)invdeltau.o \
$(RSPCLUPATH)initzero1.o \
$(RSPCLUPATH)localrot2imp.o \
$(RSPCLUPATH)locquant2.o \
$(RSPCLUPATH)madcoimp.o \
$(RSPCLUPATH)newdir1.o \
$(RSPCLUPATH)newpot1.o \
$(RSPCLUPATH)pothandleimp.o \
$(RSPCLUPATH)printres1.o \
$(RSPCLUPATH)printscf1.o \
$(RSPCLUPATH)readclu.o \
$(RSPCLUPATH)readgeom.o \
$(RSPCLUPATH)readini.o \
$(RSPCLUPATH)restart1.o \
$(RSPCLUPATH)sorttmat.o \
$(RSPCLUPATH)tauhin.o \
$(RSPCLUPATH)thin.o \
$(RSPCLUPATH)tmatiniimp.o \
$(RSPCLUPATH)vgen1.o \
$(RSPCLUPATH)vhin.o \
$(RSPCLUPATH)derivt.o \
$(RSPCLUPATH)gradm.o \
$(RSPCLUPATH)jij.o \
$(RSPCLUPATH)phase.o
#
RSPCLUSRC = \
$(RSPCLUPATH)alphamat1.f    \
$(RSPCLUPATH)deltau.f       \
$(RSPCLUPATH)dens.f \
$(RSPCLUPATH)dtat.f \
$(RSPCLUPATH)ecoreimp.f \
$(RSPCLUPATH)enprint1.f \
$(RSPCLUPATH)gettaucl.f \
$(RSPCLUPATH)invdeltau.f \
$(RSPCLUPATH)initzero1.f \
$(RSPCLUPATH)localrot2imp.f \
$(RSPCLUPATH)locquant2.f \
$(RSPCLUPATH)madcoimp.f \
$(RSPCLUPATH)newdir1.f \
$(RSPCLUPATH)newpot1.f \
$(RSPCLUPATH)pothandleimp.f \
$(RSPCLUPATH)printres1.f \
$(RSPCLUPATH)printscf1.f \
$(RSPCLUPATH)readclu.f \
$(RSPCLUPATH)readgeom.f \
$(RSPCLUPATH)readini.f \
$(RSPCLUPATH)restart1.f \
$(RSPCLUPATH)sorttmat.f \
$(RSPCLUPATH)tauhin.f \
$(RSPCLUPATH)thin.f \
$(RSPCLUPATH)tmatiniimp.f \
$(RSPCLUPATH)vgen1.f \
$(RSPCLUPATH)vhin.f \
$(RSPCLUPATH)derivt.f \
$(RSPCLUPATH)gradm.f \
$(RSPCLUPATH)jij.o \
$(RSPCLUPATH)phase.f
#
MPICLUOBJ = \
$(RSPCLUPATH)collect_data.o
#
MPICLUSRC = \
$(RSPCLUPATH)collect_data.f
#
RSPCLUPRJ = \
$(RSPCLUPATH)alphamat1.prj  \
$(RSPCLUPATH)deltau.prj     \
$(RSPCLUPATH)dens.prj     \
$(RSPCLUPATH)dtat.prj     \
$(RSPCLUPATH)ecoreimp.prj \
$(RSPCLUPATH)enprint1.prj \
$(RSPCLUPATH)gettaucl.prj \
$(RSPCLUPATH)invdeltau.prj \
$(RSPCLUPATH)initzero1.prj \
$(RSPCLUPATH)localrot2imp.prj \
$(RSPCLUPATH)locquant1.prj \
$(RSPCLUPATH)madcoimp.prj \
$(RSPCLUPATH)newdir1.prj   \
$(RSPCLUPATH)newpot1.prj   \
$(RSPCLUPATH)pothandleimp.prj \
$(RSPCLUPATH)printres1.prj \
$(RSPCLUPATH)printscf1.prj \
$(RSPCLUPATH)readclu.prj \
$(RSPCLUPATH)readgeom.prj \
$(RSPCLUPATH)readini.prj  \
$(RSPCLUPATH)restart1.prj  \
$(RSPCLUPATH)sorttmat.prj \
$(RSPCLUPATH)tauhin.prj \
$(RSPCLUPATH)thin.prj \
$(RSPCLUPATH)tmatiniimp.prj  \
$(RSPCLUPATH)vgen1.prj \
$(RSPCLUPATH)vhin.prj \
$(RSPCLUPATH)phase.prj \
$(RSPCLUPATH)collect_data.prj
#
RSPCLUSPG = \
$(RSPCLUPATH)alphamat1.spg  \
$(RSPCLUPATH)deltau.spg     \
$(RSPCLUPATH)dens.spg     \
$(RSPCLUPATH)dtat.spg \
$(RSPCLUPATH)ecoreimp.spg \
$(RSPCLUPATH)enprint1.spg \
$(RSPCLUPATH)gettaucl.spg \
$(RSPCLUPATH)invdeltau.spg \
$(RSPCLUPATH)initzeroq1.spg \
$(RSPCLUPATH)localrot2imp.spg \
$(RSPCLUPATH)locquant1.spg \
$(RSPCLUPATH)madcoimp.spg \
$(RSPCLUPATH)newdir1.spg   \
$(RSPCLUPATH)newpot1.spg   \
$(RSPCLUPATH)pothandleimp.spg \
$(RSPCLUPATH)printres1.spg \
$(RSPCLUPATH)printscf1.spg \
$(RSPCLUPATH)readclu.spg \
$(RSPCLUPATH)readgeom.spg \
$(RSPCLUPATH)readini.spg  \
$(RSPCLUPATH)restart1.spg  \
$(RSPCLUPATH)sorttmat.spg \
$(RSPCLUPATH)tauhin.spg \
$(RSPCLUPATH)thin.spg \
$(RSPCLUPATH)tmatiniimp.spg  \
$(RSPCLUPATH)vgen1.spg \
$(RSPCLUPATH)vhin.spg \
$(RSPCLUPATH)phase.spg \
$(RSPCLUPATH)collect_data.spg
#
RSPCLUSMB = \
$(RSPCLUPATH)alphamat1.smb  \
$(RSPCLUPATH)dens.smb     \
$(RSPCLUPATH)deltau.smb     \
$(RSPCLUPATH)dtat.smb \
$(RSPCLUPATH)ecoreimp.smb \
$(RSPCLUPATH)enprint1.smb \
$(RSPCLUPATH)gettaucl.smb \
$(RSPCLUPATH)invdeltau.smb \
$(RSPCLUPATH)initzero1.smb \
$(RSPCLUPATH)localrot2imp.smb \
$(RSPCLUPATH)locquant1.smb \
$(RSPCLUPATH)madcoimp.smb \
$(RSPCLUPATH)newdir1.smb    \
$(RSPCLUPATH)newpot1.smb  \
$(RSPCLUPATH)pothandleimp.smb \
$(RSPCLUPATH)printres1.smb \
$(RSPCLUPATH)printscf1.smb \
$(RSPCLUPATH)readclu.smb \
$(RSPCLUPATH)readgeom.smb \
$(RSPCLUPATH)readini.smb  \
$(RSPCLUPATH)restart1.smb  \
$(RSPCLUPATH)sorttmat.smb \
$(RSPCLUPATH)tauhin.smb \
$(RSPCLUPATH)thin.smb \
$(RSPCLUPATH)tmatiniimp.smb  \
$(RSPCLUPATH)vgen1.smb \
$(RSPCLUPATH)vhin.smb \
$(RSPCLUPATH)phase.smb \
$(RSPCLUPATH)collect_data.smb 
#
#--------------------------------------
RSPCLUDOS-OBJ = \
$(RSPCLUPATH)alphamat1.o    \
$(RSPCLUPATH)deltau.o       \
$(RSPCLUPATH)dens.o \
$(RSPCLUPATH)dtat.o \
$(RSPCLUPATH)ecoreimp.o \
$(RSPCLUPATH)enprint1.o \
$(RSPCLUPATH)gettaucl.o \
$(RSPCLUPATH)invdeltau.o \
$(RSPCLUPATH)initzero1.o \
$(RSPCLUPATH)localrot2imp.o \
$(RSPCLUPATH)locquant_dos.o \
$(RSPCLUPATH)madcoimp.o \
$(RSPCLUPATH)newdir1.o \
$(RSPCLUPATH)newpot1.o \
$(RSPCLUPATH)pothandleimp.o \
$(RSPCLUPATH)printres1.o \
$(RSPCLUPATH)printscf1.o \
$(RSPCLUPATH)readclu.o \
$(RSPCLUPATH)readgeom.o \
$(RSPCLUPATH)readini.o \
$(RSPCLUPATH)restart1.o \
$(RSPCLUPATH)sorttmat.o \
$(RSPCLUPATH)tauhin.o \
$(RSPCLUPATH)thin.o \
$(RSPCLUPATH)tmatiniimp.o \
$(RSPCLUPATH)vgen1.o \
$(RSPCLUPATH)vhin.o \
$(RSPCLUPATH)derivt.o \
$(RSPCLUPATH)gradm.o \
$(RSPCLUPATH)jij.o 
#
RSPCLUDOSSRC = \
$(RSPCLUPATH)alphamat1.f    \
$(RSPCLUPATH)deltau.f       \
$(RSPCLUPATH)dens.f \
$(RSPCLUPATH)dtat.f \
$(RSPCLUPATH)ecoreimp.f \
$(RSPCLUPATH)enprint1.f \
$(RSPCLUPATH)gettaucl.f \
$(RSPCLUPATH)invdeltau.f \
$(RSPCLUPATH)initzero1.f \
$(RSPCLUPATH)localrot2imp.f \
$(RSPCLUPATH)locquant_dos.f \
$(RSPCLUPATH)madcoimp.f \
$(RSPCLUPATH)newdir1.f \
$(RSPCLUPATH)newpot1.f \
$(RSPCLUPATH)pothandleimp.f \
$(RSPCLUPATH)printres1.f \
$(RSPCLUPATH)printscf1.f \
$(RSPCLUPATH)readclu.f \
$(RSPCLUPATH)readgeom.f \
$(RSPCLUPATH)readini.f \
$(RSPCLUPATH)restart1.f \
$(RSPCLUPATH)sorttmat.f \
$(RSPCLUPATH)tauhin.f \
$(RSPCLUPATH)thin.f \
$(RSPCLUPATH)tmatiniimp.f \
$(RSPCLUPATH)vgen1.f \
$(RSPCLUPATH)vhin.f \
$(RSPCLUPATH)derivt.f \
$(RSPCLUPATH)gradm.f \
$(RSPCLUPATH)jij.o
#
RSPCLUDOSPRJ = \
$(RSPCLUPATH)alphamat1.prj  \
$(RSPCLUPATH)deltau.prj     \
$(RSPCLUPATH)dens.prj     \
$(RSPCLUPATH)dtat.prj     \
$(RSPCLUPATH)ecoreimp.prj \
$(RSPCLUPATH)enprint1.prj \
$(RSPCLUPATH)gettaucl.prj \
$(RSPCLUPATH)invdeltau.prj \
$(RSPCLUPATH)initzero1.prj \
$(RSPCLUPATH)localrot2imp.prj \
$(RSPCLUPATH)locquant_dos.prj \
$(RSPCLUPATH)madcoimp.prj \
$(RSPCLUPATH)newdir1.prj   \
$(RSPCLUPATH)newpot1.prj   \
$(RSPCLUPATH)pothandleimp.prj \
$(RSPCLUPATH)printres1.prj \
$(RSPCLUPATH)printscf1.prj \
$(RSPCLUPATH)readclu.prj \
$(RSPCLUPATH)readgeom.prj \
$(RSPCLUPATH)readini.prj  \
$(RSPCLUPATH)restart1.prj  \
$(RSPCLUPATH)sorttmat.prj \
$(RSPCLUPATH)tauhin.prj \
$(RSPCLUPATH)thin.prj \
$(RSPCLUPATH)tmatiniimp.prj  \
$(RSPCLUPATH)vgen1.prj \
$(RSPCLUPATH)vhin.prj
#
RSPCLUDOSSPG = \
$(RSPCLUPATH)alphamat1.spg  \
$(RSPCLUPATH)deltau.spg     \
$(RSPCLUPATH)dens.spg     \
$(RSPCLUPATH)dtat.spg \
$(RSPCLUPATH)ecoreimp.spg \
$(RSPCLUPATH)enprint1.spg \
$(RSPCLUPATH)gettaucl.spg \
$(RSPCLUPATH)invdeltau.spg \
$(RSPCLUPATH)initzeroq1.spg \
$(RSPCLUPATH)localrot2imp.spg \
$(RSPCLUPATH)locquant_dos.spg \
$(RSPCLUPATH)madcoimp.spg \
$(RSPCLUPATH)newdir1.spg   \
$(RSPCLUPATH)newpot1.spg   \
$(RSPCLUPATH)pothandleimp.spg \
$(RSPCLUPATH)printres1.spg \
$(RSPCLUPATH)printscf1.spg \
$(RSPCLUPATH)readclu.spg \
$(RSPCLUPATH)readgeom.spg \
$(RSPCLUPATH)readini.spg  \
$(RSPCLUPATH)restart1.spg  \
$(RSPCLUPATH)sorttmat.spg \
$(RSPCLUPATH)tauhin.spg \
$(RSPCLUPATH)thin.spg \
$(RSPCLUPATH)tmatiniimp.spg  \
$(RSPCLUPATH)vgen1.spg \
$(RSPCLUPATH)vhin.spg
#
RSPCLUDOSSMB = \
$(RSPCLUPATH)alphamat1.smb  \
$(RSPCLUPATH)dens.smb     \
$(RSPCLUPATH)deltau.smb     \
$(RSPCLUPATH)dtat.smb \
$(RSPCLUPATH)ecoreimp.smb \
$(RSPCLUPATH)enprint1.smb \
$(RSPCLUPATH)gettaucl.smb \
$(RSPCLUPATH)invdeltau.smb \
$(RSPCLUPATH)initzero1.smb \
$(RSPCLUPATH)localrot2imp.smb \
$(RSPCLUPATH)locquant_dos.smb \
$(RSPCLUPATH)madcoimp.smb \
$(RSPCLUPATH)newdir1.smb    \
$(RSPCLUPATH)newpot1.smb  \
$(RSPCLUPATH)pothandleimp.smb \
$(RSPCLUPATH)printres1.smb \
$(RSPCLUPATH)printscf1.smb \
$(RSPCLUPATH)readclu.smb \
$(RSPCLUPATH)readgeom.smb \
$(RSPCLUPATH)readini.smb  \
$(RSPCLUPATH)restart1.smb  \
$(RSPCLUPATH)sorttmat.smb \
$(RSPCLUPATH)tauhin.smb \
$(RSPCLUPATH)thin.smb \
$(RSPCLUPATH)tmatiniimp.smb  \
$(RSPCLUPATH)vgen1.smb \
$(RSPCLUPATH)vhin.smb
#
#-----------------------------------------------
NRSPLAYPATH = $(KPATH)nrsp-lay/
NRSPLAY-OBJ = \
$(NRSPLAYPATH)cpacoord.o \
$(NRSPLAYPATH)enprint.o \
$(NRSPLAYPATH)madcolay.o \
$(NRSPLAYPATH)newpot.o \
$(NRSPLAYPATH)phystau2.o \
$(NRSPLAYPATH)pothandle.o \
$(NRSPLAYPATH)printres.o \
$(NRSPLAYPATH)printscf.o \
$(NRSPLAYPATH)readclu.o \
$(NRSPLAYPATH)readgeom.o \
$(NRSPLAYPATH)readini.o \
$(NRSPLAYPATH)restart.o \
$(NRSPLAYPATH)tau2d_god.o \
$(NRSPLAYPATH)tauhout.o \
$(NRSPLAYPATH)thout.o \
$(NRSPLAYPATH)unsct.o \
$(NRSPLAYPATH)vgen.o \
$(NRSPLAYPATH)vhout.o
#
NRSPLAYSRC = \
$(NRSPLAYPATH)cpacoord.f \
$(NRSPLAYPATH)enprint.f \
$(NRSPLAYPATH)madcolay.f \
$(NRSPLAYPATH)newpot.f \
$(NRSPLAYPATH)phystau2.f \
$(NRSPLAYPATH)pothandle.f \
$(NRSPLAYPATH)printres.f \
$(NRSPLAYPATH)printscf.f \
$(NRSPLAYPATH)readclu.f \
$(NRSPLAYPATH)readgeom.f \
$(NRSPLAYPATH)readini.f \
$(NRSPLAYPATH)restart.f \
$(NRSPLAYPATH)tau2d_god.f \
$(NRSPLAYPATH)tauhout.f \
$(NRSPLAYPATH)thout.f \
$(NRSPLAYPATH)unsct.f \
$(NRSPLAYPATH)vgen.f \
$(NRSPLAYPATH)vhout.f
#
NRSPLAYPRJ = \
$(NRSPLAYPATH)cpacoord.prj \
$(NRSPLAYPATH)enprint.prj \
$(NRSPLAYPATH)madcolay.prj \
$(NRSPLAYPATH)newpot.prj   \
$(NRSPLAYPATH)phystau2.prj \
$(NRSPLAYPATH)pothandle.prj \
$(NRSPLAYPATH)printres.prj \
$(NRSPLAYPATH)printscf.prj \
$(NRSPLAYPATH)readclu.prj \
$(NRSPLAYPATH)readgeom.prj \
$(NRSPLAYPATH)readini.prj  \
$(NRSPLAYPATH)restart.prj  \
$(NRSPLAYPATH)tau2d_god.prj \
$(NRSPLAYPATH)tauhout.prj \
$(NRSPLAYPATH)thout.prj \
$(NRSPLAYPATH)tmatini.prj  \
$(NRSPLAYPATH)unsct.prj \
$(NRSPLAYPATH)vgen.prj \
$(NRSPLAYPATH)vhout.prj
#
NRSPLAYSPG = \
$(NRSPLAYPATH)cpacoord.spg \
$(NRSPLAYPATH)enprint.spg \
$(NRSPLAYPATH)madcolay.spg \
$(NRSPLAYPATH)newpot.spg   \
$(NRSPLAYPATH)phystau2.spg \
$(NRSPLAYPATH)pothandle.spg \
$(NRSPLAYPATH)printres.spg \
$(NRSPLAYPATH)printscf.spg \
$(NRSPLAYPATH)readclu.spg \
$(NRSPLAYPATH)readgeom.spg \
$(NRSPLAYPATH)readini.spg  \
$(NRSPLAYPATH)restart.spg  \
$(NRSPLAYPATH)tau2d_god.spg \
$(NRSPLAYPATH)tauhout.spg \
$(NRSPLAYPATH)thout.spg \
$(NRSPLAYPATH)tmatini.spg  \
$(NRSPLAYPATH)unsct.spg \
$(NRSPLAYPATH)vgen.spg \
$(NRSPLAYPATH)vhout.spg
#
NRSPLAYSMB = \
$(NRSPLAYPATH)cpacoord.smb \
$(NRSPLAYPATH)dens.smb     \
$(NRSPLAYPATH)enprint.smb \
$(NRSPLAYPATH)madcolay.smb \
$(NRSPLAYPATH)newpot.smb  \
$(NRSPLAYPATH)phystau2.smb \
$(NRSPLAYPATH)pothandle.smb \
$(NRSPLAYPATH)printres.smb \
$(NRSPLAYPATH)printscf.smb \
$(NRSPLAYPATH)readclu.smb \
$(NRSPLAYPATH)readgeom.smb \
$(NRSPLAYPATH)readini.smb  \
$(NRSPLAYPATH)restart.smb  \
$(NRSPLAYPATH)tau2d_god.smb \
$(NRSPLAYPATH)tauhout.smb \
$(NRSPLAYPATH)thout.smb \
$(NRSPLAYPATH)tmatini.smb  \
$(NRSPLAYPATH)unsct.smb \
$(NRSPLAYPATH)vgen.smb \
$(NRSPLAYPATH)vhout.smb
#
#
     skkr-obj   = $(MAINOBJ) $(NONRELOBJ) $(NRSCFOBJ) $(NRSPSCFOBJ) $(NRDLMSCFOBJ) \
	$(RELOBJ) $(RELLAYOBJ) $(RSPSCF-ksave-OBJ) $(RSPSCF-save0-OBJ) $(RSPSCF-save-OBJ) \
        $(RSPLAY-save-OBJ) $(RSPCLU-OBJ) $(MPILAYYOBJ) $(MPICLUOBJ)
     nr-scf-obj = $(MAINOBJ) $(NONRELOBJ) $(NRSCFOBJ)
     nrsp-scf-obj = $(MAINOBJ) $(NONRELOBJ) $(NRSPSCFOBJ)
     nrdlm-scf-obj = $(MAINOBJ) $(NRDLMSCFOBJ)
     nrsp-lay-obj = $(MAINOBJ) $(NONRELOBJ) $(NRSPLAY-OBJ)
     rsp-scf_ksave-obj = $(MAINOBJ) $(RELOBJ) $(RSPSCF-ksave-OBJ)
     rsp-scf_save0-obj = $(MAINOBJ) $(RELOBJ) $(RSPSCF-save0-OBJ)
     rsp-scf_save-obj = $(MAINOBJ) $(RELOBJ) $(RSPSCF-save-OBJ)
     rsp-lay_save-obj = $(MAINOBJ) $(RELOBJ) $(RELLAYOBJ) $(RSPLAY-save-OBJ)
     rsp-lay_save_mpi-obj = $(MAINOBJ) $(RELOBJ) $(RELLAYOBJ) $(RSPLAY-save-OBJ) $(MPILAYOBJ)
     rsp-clu-obj = $(MAINOBJ) $(RELOBJ) $(RSPCLU-OBJ)
     rsp-clu_mpi-obj = $(MAINOBJ) $(RELOBJ) $(RSPCLU-OBJ) $(MPICLUOBJ)
     rsp-clu-dos-obj = $(MAINOBJ) $(RELOBJ) $(RSPCLUDOS-OBJ)
#
     skkr-src   = $(MAINSRC) $(NONRELSRC) $(NRSCFSRC) $(NRSPSCFSRC) $(NRDLMSCFSRC) \
	$(RELSRC) $(RSPSCFSRC)
     nr-scf-src = $(MAINSRC) $(NONRELSRC) $(NRSCFSRC)
     nrsp-scf-src = $(MAINSRC) $(NONRELSRC) $(NRSPSCFSRC)
     nrdlm-scf-src = $(MAINSRC) $(NRDLMSCFSRC)
     nrsp-lay-src = $(MAINSRC) $(NONRELSRC) $(NRSPLAYSRC)
     rsp-scf-src = $(MAINSRC) $(RELSRC) $(RSPSCFSRC)
     rsp-lay-src = $(MAINSRC) $(RELSRC) $(RSPLAYSRC)
     rsp-clu-src = $(MAINSRC) $(RELSRC) $(RSPCLUSRC)
     rsp-clu-dos-src = $(MAINSRC) $(RELSRC) $(RSPCLUDOSSRC)
#
     skkr-prj   = $(MAINPRJ) $(NONRELPRJ) $(NRSCFPRJ) $(NRSPSCFPRJ) \
	$(RELPRJ) $(RSPSCFPRJ)
     nr-scf-prj = $(MAINPRJ) $(NONRELPRJ) $(NRSCFPRJ)
     nrsp-scf-prj = $(MAINPRJ) $(NONRELPRJ) $(NRSPSCFPRJ)
     nrsp-lay-prj = $(MAINPRJ) $(NONRELPRJ) $(NRSPLAYPRJ)
     rsp-scf-prj = $(MAINPRJ) $(RELPRJ) $(RSPSCFPRJ)
     rsp-lay-prj = $(MAINPRJ) $(RELPRJ) $(RSPLAYPRJ)
     rsp-clu-prj = $(MAINPRJ) $(RELPRJ) $(RSPCLUPRJ)
     rsp-clu-dos-prj = $(MAINPRJ) $(RELPRJ) $(RSPCLUDOSPRJ)
#
     skkr-smb   = $(MAINSMB) $(NONRELSMB) $(NRSCFSMB) $(NRSPSCFSMB) \
	$(RELSMB) $(RSPSCFSMB)
     nr-scf-smb = $(MAINSMB) $(NONRELSMB) $(NRSCFSMB)
     nrsp-scf-smb = $(MAINSMB) $(NONRELSMB) $(NRSPSCFSMB)
     nrsp-lay-smb = $(MAINSMB) $(NONRELSMB) $(NRSPLAYSMB)
     rsp-scf-smb = $(MAINSMB) $(RELSMB) $(RSPSCFSMB)
     rsp-lay-smb = $(MAINSMB) $(RELSMB) $(RSPLAYSMB)
     rsp-clu-smb = $(MAINSMB) $(RELSMB) $(RSPCLUSMB)
     rsp-clu-dos-smb = $(MAINSMB) $(RELSMB) $(RSPCLUDOSSMB)
#
     skkr-spg   = $(MAINSPG) $(NONRELSPG) $(NRSCFSPG) $(NRSPSCFSPG) \
	$(RELSPG) $(RSPSCFSPG)
     nr-scf-spg = $(MAINSPG) $(NONRELSPG) $(NRSCFSPG)
     nrsp-scf-spg = $(MAINSPG) $(NONRELSPG) $(NRSPSCFSPG)
     nrsp-lay-spg = $(MAINSPG) $(NONRELSPG) $(NRSPLAYSPG)
     rsp-scf-spg = $(MAINSPG) $(RELSPG) $(RSPSCFSPG)
     rsp-lay-spg = $(MAINSPG) $(RELSPG) $(RSPLAYSPG)
     rsp-clu-spg = $(MAINSPG) $(RELSPG) $(RSPCLUSPG)
     rsp-clu-dos-spg = $(MAINSPG) $(RELSPG) $(RSPCLUDOSSPG)
#
# dependencies
#
# .SUFFIXES: .o .f
.SUFFIXES: .crf .run .o .prj .f
#
# some extra intuitive stuff for serial/MPI compile
rsp-lay_save_mpi.run rsp-clu_mpi.run: FC = $(FCPAR)
#
rsp-lay_save_mpi.run rsp-clu_mpi.run: FFLAGS += -DMPIP
#
.f.o:
	$(FC) -c $(FFLAGS) $<
	mv $(*F).o $(*D)
# $<   : source files of an out-of-date module
# $(*F): "file part" 
# $(*D): "directory part" of target
#
.f.prj:
	@cp param.h ../
	$(FU) $(FUFLAGS) $(FUPRO) -out=$* $<
	@rm -f ../param.h 
# $*   : file name (without the suffix) of the source file
# $@   : full target name of the current target
#
# activities: executables
#
all: $(PROGRAMS)
#
nr-scf.run: $(NRSCFPATH)skkr.o $(nr-scf-obj)
	@echo " "
	@echo "Linking nr-scf.run ..."
	@echo "-----------------------"
	@echo " "
	$(FC) -o $@ $(NRSCFPATH)skkr.o $(nr-scf-obj) $(LFLAGS)
	@echo " "
	@echo "done"
#
nrsp-scf.run: $(NRSPSCFPATH)skkr.o $(nrsp-scf-obj)
	@echo " "
	@echo "Linking nrsp-scf.run ..."
	@echo "-----------------------"
	@echo " "
	$(FC) -o $@ $(NRSPSCFPATH)skkr.o $(nrsp-scf-obj) $(LFLAGS)
	@echo " "
	@echo "done"
#
nrdlm-scf.run: $(NRDLMSCFPATH)skkr.o $(nrdlm-scf-obj)
	@echo " "
	@echo "Linking nrdlm-scf.run ..."
	@echo "-----------------------"
	@echo " "
	$(FC) -o $@ $(NRDLMSCFPATH)skkr.o $(nrdlm-scf-obj) $(LFLAGS)
	@echo " "
	@echo "done"
#
rsp-scf_ksave.run: $(RSPSCFPATH)skkr.o $(rsp-scf_ksave-obj)
	@echo " "
	@echo "Linking rsp-scf_ksave.run ..."
	@echo "-----------------------------"
	@echo " "
	$(FC) -o $@ $(RSPSCFPATH)skkr.o $(rsp-scf_ksave-obj) $(LFLAGS)
	@echo " "
	@echo "done"
#
rsp-scf_save0.run: $(RSPSCFPATH)skkr.o $(rsp-scf_save0-obj)
	@echo " "
	@echo "Linking rsp-scf_save0.run ..."
	@echo "-----------------------------"
	@echo " "
	$(FC) -o $@ $(RSPSCFPATH)skkr.o $(rsp-scf_save0-obj) $(LFLAGS)
	@echo " "
	@echo "done"
#
rsp-scf_save.run: $(RSPSCFPATH)skkr.o $(rsp-scf_save-obj)
	@echo " "
	@echo "Linking rsp-scf_save.run ..."
	@echo "----------------------------"
	@echo " "
	$(FC) -o $@ $(RSPSCFPATH)skkr.o $(rsp-scf_save-obj) $(LFLAGS)
	@echo " "
	@echo "done"
#
rsp-lay_save_mpi.run: $(RSPLAYPATH)skkr.o $(rsp-lay_save_mpi-obj)
	@echo " "
	@echo "Linking rsp-lay_save_mpi.run ..."
	@echo "----------------------------"
	@echo " "
	$(FC) -o $@ $(RSPLAYPATH)skkr.o $(rsp-lay_save_mpi-obj) $(LFLAGS)
	@echo " "
	@echo "done"
#
rsp-lay_save.run: $(RSPLAYPATH)skkr.o $(rsp-lay_save-obj)
	@echo " "
	@echo "Linking rsp-lay_save.run ..."
	@echo "----------------------------"
	@echo " "
	$(FC) -o $@ $(RSPLAYPATH)skkr.o $(rsp-lay_save-obj) $(LFLAGS)
	@echo " "
	@echo "done"
#
nrsp-lay.run: $(NRSPLAYPATH)skkr.o $(nrsp-lay-obj)
	@echo " "
	@echo "Linking nrsp-lay.run ..."
	@echo "----------------------------"
	@echo " "
	$(FC) -o $@ $(NRSPLAYPATH)skkr.o $(nrsp-lay-obj) $(LFLAGS)
	@echo " "
	@echo "done"
#
rsp-clu.run: $(RSPCLUPATH)skkr.o $(rsp-clu-obj)
	@echo " "
	@echo "Linking rsp-clu.run ..."
	@echo "----------------------------"
	@echo " "
	$(FC) -o $@ $(RSPCLUPATH)skkr.o $(rsp-clu-obj) $(LFLAGS)
	@echo " "
	@echo "done"
#
rsp-clu_mpi.run: $(RSPCLUPATH)skkr.o $(rsp-clu_mpi-obj)
	@echo " "
	@echo "Linking rsp-clu.run ..."
	@echo "----------------------------"
	@echo " "
	$(FC) -o $@ $(RSPCLUPATH)skkr.o $(rsp-clu_mpi-obj) $(LFLAGS)
	@echo " "
	@echo "done"
#
rsp-clu-dos.run: $(RSPCLUPATH)skkr_dos.o $(rsp-clu-dos-obj)
	@echo " "
	@echo "Linking rsp-clu.run ..."
	@echo "----------------------------"
	@echo " "
	$(FC) -o $@ $(RSPCLUPATH)skkr_dos.o $(rsp-clu-dos-obj) $(LFLAGS)
	@echo " "
	@echo "done"
#
# delete object files
#
del-all.o:
	@echo " "
	@echo "Deleting ALL object-files ..."
	@rm -f $(skkr-obj) $(rsp-clu_mpi-obj) $(NRSCFPATH)skkr.o $(NRSPSCFPATH)skkr.o \
	$(RSPSCFPATH)skkr.o $(RSPLAYPATH)skkr.o $(RSPCLUPATH)skkr.o $(NRSPLAYPATH)skkr.o \
	$(RSPCLUPATH)skkr.o $(RSPCLUDOS-OBJ) $(RSPCLUPATH)skkr_dos.o
	@echo "done"
	@echo " "
#
del-main.o:
	@echo " "
	@echo "Deleting object-files from main/ ..."
	@rm -f $(MAINOBJ)
	@echo "done"
	@echo " "
#
del-nonrel.o:
	@echo " "
	@echo "Deleting object-files from all nonrel. subdirectories"
	@rm -f $(NONRELOBJ) $(NRSCFOBJ) $(NRSPSCFOBJ) \
	$(NRSCFPATH)skkr.o $(NRSPSCFPATH)skkr.o
	@echo "done"
	@echo " "
#
del-rel.o:
	@echo " "
	@echo "Deleting object-files from all rel. subdirectories"
	@rm -f $(RELOBJ) $(RSPSCF-ksave-OBJ) $(RSPSCF-save0-OBJ) \
	$(RSPSCF-save-OBJ) $(RSPSCFPATH)skkr.o
	@echo "done"
	@echo " "
#
# cross references
#
allcrf: $(CROSSREF)
#
nr-scf.crf: $(NRSCFPATH)skkr.prj $(nr-scf-prj)
	@echo " "
	@echo "Making nr-scf.crf ..."
	@echo "----------------------"
	@echo " "
	$(FU) $(FUFLAGS) -out=$@ $(NRSCFPATH)skkr.prj $(nr-scf-prj)
	@echo " "
	@echo "done"
#
nrsp-scf.crf: $(NRSPSCFPATH)skkr.prj $(nrsp-scf-prj)
	@echo " "
	@echo "Making nrsp-scf.crf ..."
	@echo "----------------------"
	@echo " "
	$(FU) $(FUFLAGS) -out=$@ $(NRSPSCFPATH)skkr.prj $(nrsp-scf-prj)
	@echo " "
	@echo "done"
#
rsp-scf.crf: $(RSPSCFPATH)skkr.prj $(rsp-scf-prj)
	@echo " "
	@echo "Making rsp-scf.crf ..."
	@echo "----------------------"
	@echo " "
	$(FU) $(FUFLAGS) -out=$@ $(RSPSCFPATH)skkr.prj $(rsp-scf-prj)
	@echo " "
	@echo "done"
#
#
rsp-lay.crf: $(RSPLAYPATH)skkr.prj $(rsp-lay-prj)
	@echo " "
	@echo "Making rsp-lay.crf ..."
	@echo "----------------------"
	@echo " "
	$(FU) $(FUFLAGS) -out=$@ $(RSPLAYPATH)skkr.prj $(rsp-lay-prj)
	@echo " "
	@echo "done"
#
nrsp-lay.crf: $(NRSPLAYPATH)skkr.prj $(nrsp-lay-prj)
	@echo " "
	@echo "Making nrsp-lay.crf ..."
	@echo "----------------------"
	@echo " "
	$(FU) $(FUFLAGS) -out=$@ $(NRSPLAYPATH)skkr.prj $(nrsp-lay-prj)
	@echo " "
	@echo "done"
#
rsp-clu.crf: $(RSPCLUPATH)skkr.prj $(rsp-clu-prj)
	@echo " "
	@echo "Making rsp-clu.crf ..."
	@echo "----------------------"
	@echo " "
	$(FU) $(FUFLAGS) -out=$@ $(RSPCLUPATH)skkr.prj $(rsp-clu-prj)
	@echo " "
	@echo "done"
#
# delete project files and listings
#
del-all.prj:
	@echo " "
	@echo "Deleting ALL project-files ..."
	@rm -f $(skkr-prj) $(NRSCFPATH)skkr.prj $(NRSPSCFPATH)skkr.prj \
	$(RSPSCFPATH)skkr.prj $(RSPLAYPATH)skkr.prj \
        $(RSPCLUPATH)skkr.prj $(NRSPLAYPATH)skkr.prj
	@echo "done"
	@echo " "
#
del-all.lis:
	@echo " "
	@echo "Deleting ALL listing-files ..."
	@rm -f $(MAINPATH)*.lis $(NONRELPATH)*.lis $(NRSCFPATH)*.lis \
	$(NRSPSCFPATH)*.lis $(RELPATH)*.lis $(RSPSCFPATH)*.lis
	@echo "done"
	@echo " "
#
# source tool
#
spag-main:    
	$(FS) fig=$(PLFCONF) log=$@ 215=$(MAINPATH) 230=$(MAINPATH) \
	232=$(MAINPATH) $(MAINSRC)
#
spag-nonrel:    
	$(FS) fig=$(PLFCONF) log=$@ 215=$(NONRELPATH) 230=$(NONRELPATH) \
	232=$(NONRELPATH) $(NONRELSRC)
#
spag-nr-scf:    
	$(FS) fig=$(PLFCONF) log=$@ 215=$(NRSCFPATH) 230=$(NRSCFPATH) \
	232=$(NRSCFPATH) $(NRSCFPATH)skkr.f $(NRSCFSRC)
#
spag-nrsp-scf:    
	$(FS) fig=$(PLFCONF) log=$@ 215=$(NRSPSCFPATH) 230=$(NRSPSCFPATH) \
	232=$(NRSPSCFPATH) $(NRSPSCFPATH)skkr.f $(NRSPSCFSRC)
#
spag-rel:    
	$(FS) fig=$(PLFCONF) log=$@ 215=$(RELPATH) 230=$(RELPATH) \
	232=$(RELPATH) $(RELSRC)
#
spag-rsp-scf:    
	$(FS) fig=$(PLFCONF) log=$@ 215=$(RSPSCFPATH) 230=$(RSPSCFPATH) \
	232=$(RSPSCFPATH) $(RSPSCFPATH)skkr.f $(RSPSCFSRC)
#
# 215  : search list for INCLUDE files
# 230  : directory for output code file
# 232  : directory for symbol tables file
#
nr-scf.sum: 
	$(FX) fig=$(PLFCONF) sum=$@ to=gxchk-nr-scf \
	$(NRSCFPATH)skkr.smb $(nr-scf-smb)
#
nrsp-scf.sum: 
	$(FX) fig=$(PLFCONF) sum=$@ to=gxchk-nrsp-scf \
	$(NRSPSCFPATH)skkr.smb $(nrsp-scf-smb)
#
rsp-scf.sum: 
	$(FX) fig=$(PLFCONF) sum=$@ to=gxchk-rsp-scf \
	$(RSPSCFPATH)skkr.smb $(rsp-scf-smb)
#
# delete codes, symbol tables and log-files created by spag
#
del-all.spag:
	@echo " "
	@echo "Deleting EVERTHING created by SPAG ..."
	@rm -f spag-main spag-nonrel spag-nr-scf spag-nrsp-scf \
	spag-rsp-scf 
	@rm -f $(skkr-smb) $(NRSCFPATH)skkr.smb $(NRSPSCFPATH)skkr.smb \
	$(RSPSCFPATH)skkr.smb
	@rm -f gxchk-nr-scf gxchk-nrsp-scf gxchk-rsp-scf
	@rm -f $(skkr-spg) $(NRSCFPATH)skkr.spg $(NRSPSCFPATH)skkr.spg \
	$(RSPSCFPATH)skkr.spg
	@rm -f nr-scf.sum nrsp-scf.sum rsp-scf.sum
	@echo "done"
	@echo " "
#
del-smb.spag:
	@echo " "
	@echo "Cleaning SPAG-files ..."
	@rm -f spag-main spag-nonrel spag-nr-scf spag-nrsp-scf \
	spag-rsp-scf
	@rm -f $(skkr-smb) $(NRSCFPATH)skkr.smb $(NRSPSCFPATH)skkr.smb \
	$(RSPSCFPATH)skkr.smb
	@rm -f gxchk-nr-scf gxchk-nrsp-scf gxchk-rsp-scf
	@echo "done"
	@echo " "
#
# utils
#
skkr.tar:
	@echo " "
	@echo "Making tar-file: skkr.tar ..."
	@echo "-----------------------------"
	tar cvf $@ ./makefile ./plusFORT.config ./param.h ./param.h-save \
	./compile-all ./compile-rsp_save \
	./intro_skkr.ps ./nrsp-io.tar ./rsp-io.tar \
	$(skkr-src) \
	$(NRSCFPATH)skkr.f \
	$(NRSPSCFPATH)skkr.f \
	$(NRDLMSCFPATH)skkr.f \
	$(RSPSCFPATH)skkr.f
	@echo " "
	@echo "done"
#
skkr-spag.tar:
	@echo " "
	@echo "Making tar-file: skkr-spag.tar ..."
	@echo "----------------------------------"
	tar cvf $@ ./makefile ./param.h ./plusFORT.config \
        ./*.sum $(skkr-spg) $(NRSCFPATH)skkr.spg $(NRSPSCFPATH)skkr.spg \
	$(RSPSCFPATH)skkr.spg
	@echo " "
	@echo "done"
#
changespg:
	@echo " "
	@echo "Change spg-extension to f ..."
	@spg2f
	@echo "done"
	@echo " "
