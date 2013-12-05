##
## Copyright 2003-2008 ANSYS, Inc. 
## All Rights Reserved  
##

#----------------------------------------------------------------------#
# makefile for user defined functions.
#----------------------------------------------------------------------#


#----------------------------------------------------------------------#
# User modifiable section.
#----------------------------------------------------------------------#
SOURCES=fsi.c fsi_udf.c
FLUENT_INC=/home/software/ansys/ansys_121_64bit/v121/fluent

# Precompiled User Object files (for example .o files from .f sources)
USER_OBJECTS=../libprecice.a ../libpython2.7.so

#----------------------------------------------------------------------#
# Build targets (do not modify below this line).
#----------------------------------------------------------------------#
CC= gcc
SHELL= /bin/sh
RELEASE= 12.1.2
VERSION= `basename "\`pwd\`"`
ARCH= `expr "\`pwd\`" : '.*/\(.*\)/[23].*'`
ARCHC=$(ARCH)
TARGET= libudf.so
UDFDATA= udf_names.c
OBJECTS= $(UDFDATA:.c=.o) $(SOURCES:.c=.o)

INCLUDES= -I$(FLUENT_INC)/fluent$(RELEASE)/$(ARCH)/$(VERSION) \
          -I$(FLUENT_INC)/fluent$(RELEASE)/src \
          -I$(FLUENT_INC)/fluent$(RELEASE)/cortex/src \
          -I$(FLUENT_INC)/fluent$(RELEASE)/client/src \
          -I$(FLUENT_INC)/fluent$(RELEASE)/tgrid/src \
          -I$(FLUENT_INC)/fluent$(RELEASE)/multiport/src \
          -I$(FLUENT_INC)/fluent$(RELEASE)/multiport/mpi_wrapper/include \
          -I$(FLUENT_INC)/fluent$(RELEASE)/multiport/mpi_wrapper/src\
	  -I.

CFLAGS_ALPHA=		-call_shared -fast -w0 -std1
CFLAGS_HPUX10P8=	+z -O +w1 -Aa -D_HPUX_SOURCE 
CFLAGS_HPUX11=		+z -O +w1 -Aa -D_HPUX_SOURCE 
CFLAGS_HPUX11_64=	+z -O +w1 -Aa -D_HPUX_SOURCE +DD64 +DA2.0W
CFLAGS_HPUX11_ia64=	-O +DD64 +DSblended -Ae
CFLAGS_HIUX=
CFLAGS_IRIX65=		-KPIC -xansi -fullwarn -O -n32 -mips3
CFLAGS_IRIX65_MIPS4=	-KPIC -xansi -fullwarn -O -n32 -mips4
CFLAGS_IRIX65_MIPS4_64=	-KPIC -xansi -fullwarn -O3 -64 -mips4 -OPT:IEEE_arithmetic=3:roundoff=3:Olimit=8000 -LNO:fusion=2:fission=0:prefetch=2 -signed
CFLAGS_AIX43=		-O -qalign=natural
CFLAGS_POWER3=		-O3  -qlanglvl=ansi -qarch=pwr3  -qfloat=hssngl -qansialias -qnolm -qalign=natural
CFLAGS_POWER3_64=	-O -q64 -qnolm -qalign=natural
CFLAGS_ULTRA=		-KPIC -fast -v 
CFLAGS_ULTRA_64=	-KPIC -fast -v -xarch=v9
CFLAGS_FUJITSU_PP=	-Kfast_GP=2,prefetch,V9
CFLAGS_LNX86=		-fpic -shared -ansi -Wall -m32 -O	# options for CC= gcc
CFLAGS_LNAMD64=         -fPIC -shared -ansi -Wall -O -DPTR_RESTRICT=	# options for CC= gcc
CFLAGS_LNIA64=		-fpic -shared -ansi -Wall -O -DPTR_RESTRICT=    # gcc options
CFLAGS_GCC=		-fPIC -shared -ansi -Wall -O -DPTR_RESTRICT=
CFLAGS_GCC64=		-fpic -shared -ansi -Wall -O -m64 -DPTR_RESTRICT= 

LDFLAGS_ALPHA=		-shared -lm
LDFLAGS_HPUX10P8=	-b -lm
LDFLAGS_HPUX11=		-b -lm
LDFLAGS_HPUX11_64=	-b -lm
LDFLAGS_HPUX11_ia64=	-b -lm
LDFLAGS_HIUX=		-b -lm
LDFLAGS_IRIX65=		-shared -n32 -mips3 -lm
LDFLAGS_IRIX65_MIPS4=	-shared -n32 -mips4 -lm
LDFLAGS_IRIX65_MIPS4_64= -shared -64 -mips4 -lm
LDFLAGS_AIX43=
LDFLAGS_POWER3=
LDFLAGS_POWER3_64=	-b64 -lmass -lm
LDFLAGS_ULTRA=		-G -lm
LDFLAGS_ULTRA_64=	-G -lm
LDFLAGS_FUJITSU_PP=	-G -lm
LDFLAGS_LNX86=		-shared -lm -m elf_i386
LDFLAGS_LNAMD64=	-shared -lm
LDFLAGS_LNIA64=		-shared -lm
LDFLAGS_GCC=		-shared -lm 
LDFLAGS_GCC64=		-shared -lm

.c.o:
	@echo "# Compiling $@ because of $?"
	@rm -f $*.o
	$(CC) $(CFLAGS) $(INCLUDES) -c $*.c

default: $(UDFDATA)
	@$(MAKE) $(ARCHC)	

$(TARGET): makefile $(UDFDATA) $(OBJECTS)
	@echo "# Linking $@ because of $?"
	$(LD) $(LDFLAGS) $(OBJECTS) $(USER_OBJECTS) -o $@

$(UDFDATA): makefile $(SOURCES)
	@echo '# Generating $@ because of $?'
	@echo '/* This file generated automatically. */' > $@
	@echo '/*          Do not modify.            */' >> $@
	@echo '#include "udf.h"' >> $@
	@echo '#include "prop.h"' >> $@
	@echo '#include "dpm.h"' >> $@
	@sed -n 's/^ *\(DEFINE_[_A-Z]* *([^{]*)\)[{ ]*/extern \1;/p' $(SOURCES) >> $@
	@echo 'UDF_Data udf_data[] = {' >> $@
	@sed -n 's/^ *DEFINE_\([_A-Z]*\) *( *\([_a-zA-Z0-9]*\)[^{]*)[{ ]*/\{"\2", (void (*)())\2, UDF_TYPE_\1\},/p' $(SOURCES) >> $@
	@echo '};' >> $@
	@echo 'int n_udf_data = sizeof(udf_data)/sizeof(UDF_Data);' >> $@
	@echo '#include "version.h"' >> $@
	@echo 'void UDF_Inquire_Release(int *major, int *minor, int *revision)' >> $@
	@echo '{' >> $@
	@echo '  *major = RampantReleaseMajor;' >> $@
	@echo '  *minor = RampantReleaseMinor;' >> $@
	@echo '  *revision = RampantReleaseRevision;' >> $@
	@echo '}' >> $@

alpha:
	$(MAKE) $(TARGET) "CFLAGS=$(CFLAGS_ALPHA)" "LD=cc" "LDFLAGS=$(LDFLAGS_ALPHA)"

hpux10p8:
	$(MAKE) $(TARGET) "CFLAGS=$(CFLAGS_HPUX10P8)" "LDFLAGS=$(LDFLAGS_HPUX10P8)"

hpux11:
	$(MAKE) $(TARGET) "CFLAGS=$(CFLAGS_HPUX11)" "LDFLAGS=$(LDFLAGS_HPUX11)"

hpux11_64:
	$(MAKE) $(TARGET) "CFLAGS=$(CFLAGS_HPUX11_64)" "LDFLAGS=$(LDFLAGS_HPUX11_64)"

hpux11_ia64:
	$(MAKE) $(TARGET) "CFLAGS=$(CFLAGS_HPUX11_ia64)" "LDFLAGS=$(LDFLAGS_HPUX11_ia64)"

hiux::
	$(MAKE) $(TARGET) "CFLAGS=$(CFLAGS_HIUX)" "LDFLAGS=$(LDFLAGS_HIUX)"

irix65:
	$(MAKE) $(TARGET) "CFLAGS=$(CFLAGS_IRIX65)" "LDFLAGS=$(LDFLAGS_IRIX65)"

irix65_mips4:
	$(MAKE) $(TARGET) "CFLAGS=$(CFLAGS_IRIX65_MIPS4)" "LDFLAGS=$(LDFLAGS_IRIX65_MIPS4)"

irix65_mips4_64:
	$(MAKE) $(TARGET) "CFLAGS=$(CFLAGS_IRIX65_MIPS4_64)" "LDFLAGS=$(LDFLAGS_IRIX65_MIPS4_64)"

aix43:
	$(MAKE) $(TARGET) "CFLAGS=$(CFLAGS_AIX43)" "LDFLAGS=$(LDFLAGS_AIX43)" "LD=makeshared"

aix51 power3:
	$(MAKE) $(TARGET) "CFLAGS=$(CFLAGS_POWER3)" "LDFLAGS=$(LDFLAGS_POWER3)" "LD=makeshared"

aix51_64 power3_64:
	$(MAKE) $(TARGET) "CFLAGS=$(CFLAGS_POWER3_64)" "LDFLAGS=$(LDFLAGS_POWER3_64)" "LD=makeshared"

ultra:
	$(MAKE) $(TARGET) "CFLAGS=$(CFLAGS_ULTRA)" "LDFLAGS=$(LDFLAGS_ULTRA)"

ultra_64:
	$(MAKE) $(TARGET) "CFLAGS=$(CFLAGS_ULTRA_64)" "LDFLAGS=$(LDFLAGS_ULTRA_64)"

fujitsu_pp:
	$(MAKE) $(TARGET) "CC=/opt/FSUNf90/bin/fcc"  "CFLAGS=$(CFLAGS_FUJITSU_PP)" "LDFLAGS=$(LDFLAGS_FUJITSU_PP)"

lnx86 lnx86_rh8:
	$(MAKE) $(TARGET) "CFLAGS=$(CFLAGS_LNX86)" "LDFLAGS=$(LDFLAGS_LNX86)"

lnamd64:
	$(MAKE) $(TARGET) LD=mpicxx.mpich2 "CFLAGS=$(CFLAGS_LNAMD64)" "LDFLAGS=$(LDFLAGS_LNAMD64)"

lnia64:
	$(MAKE) $(TARGET) "CFLAGS=$(CFLAGS_LNIA64)" "LDFLAGS=$(LDFLAGS_LNIA64)" "LD=$(CC)"

gcc:
	if [ ! -f  sunmath.h ]; then  \
          ln  -s  /usr/include/math.h  sunmath.h ; \
        fi; \
	$(MAKE) $(TARGET) CC=gcc LD=gcc "CFLAGS=$(CFLAGS_GCC)" "LDFLAGS=$(LDFLAGS_GCC)"
	rm sunmath.h

gcc64:
	if [ ! -f  sunmath.h ]; then  \
          ln  -s  /usr/include/math.h  sunmath.h ; \
        fi; \
	$(MAKE) $(TARGET) CC=gcc LD=gcc "CFLAGS=$(CFLAGS_GCC64)" "LDFLAGS=$(LDFLAGS_GCC64)"
	rm sunmath.h

clean:
	@rm -f $(TARGET) $(UDFDATA) $(OBJECTS)
