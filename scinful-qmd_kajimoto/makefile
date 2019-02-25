F77     = ifort
OPT     = -fast -noautomatic -132
LIBES   =
INCLUDES = 
TERGET  = ./scinful-qmd

OBJS = \
       banker41.o     etol5.o        main51.o       s2j.o    \
       choos.o        gamma.o        mdp.o          scin.o   \
       coll00.o       ground.o       mfield00.o     sdml.o   \
       conv.o         input6.o       phitsdump4.o   sigma4.o \
       decay4.o       interact-n.o   qmddflt.o      trans2.o \
       dedx2.o        interact2.o    qmdinit.o      utl00.o  \
       domega.o       jqmd00.o       qmdout1.o      kalbach.o\
       effn2.o        kinema4.o      qmdout2.o               \
       energetic.o    legendre.o     random.o

.f.o	:
	$(F77) $(OPT) $(INCLUDES) -c $<

all	:$(OBJS)
	$(F77) -o $(TERGET) $(OPT) $(OBJS) $(LIBES)

clean	:
	rm -f $(OBJS)
