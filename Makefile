SYSTEM := $(shell uname -s)
ROOTI = -I${shell root-config --incdir}
ROOTL = ${shell root-config --glibs} -lSpectrum
DEBUGFLAG = -g3
GCC = g++ -fPIC
all: reader MCA_data.so calibrate apply_calibration analyse compare time_dep testDB

%.o: %.cxx %.h
	${GCC} ${DEBUGFLAG} $^ -c ${ROOTI} 

reader: reader.cxx MCA_dataDict.o MCA_data.o
	${GCC} ${DEBUGFLAG} $^ -o $@ ${ROOTI} ${ROOTL}

calibrate: calibrate.cxx MCA_dataDict.o MCA_data.o
	${GCC} ${DEBUGFLAG} $^ -o $@ ${ROOTI} ${ROOTL}

analyse: analyse.cxx MCA_dataDict.o MCA_data.o
	${GCC} ${DEBUGFLAG} $^ -o $@ ${ROOTI} ${ROOTL}

compare: compare.cxx MCA_dataDict.o MCA_data.o
	${GCC} ${DEBUGFLAG} $^ -o $@ ${ROOTI} ${ROOTL}

time_dep: time_dep.cxx MCA_dataDict.o MCA_data.o
	${GCC} ${DEBUGFLAG} $^ -o $@ ${ROOTI} ${ROOTL}

apply_calibration: apply_calibration.cxx MCA_dataDict.o MCA_data.o
	${GCC} ${DEBUGFLAG} $^ -o $@ ${ROOTI} ${ROOTL}

testDB: testDB.cxx MCA_dataDict.o MCA_data.o
	${GCC} ${DEBUGFLAG} $^ -o $@ ${ROOTI} ${ROOTL}

MCA_dataDict.cxx MCA_dataDict.h:	MCA_data.h MCA_dataLinkDef.h
	rootcint -f $@ -c $(ROOTI) -p $^

cnf_reader: cnf_reader.C
	${GCC} ${DEBUGFLAG} $^ -o $@ ${ROOTI} ${ROOTL}

MCA_data.so: MCA_dataDict.o MCA_data.o
ifeq ($(SYSTEM),Darwin)
	$(GCC) -dynamiclib -single_module -O $^ $(FFTSL) -o $@ ${ROOTL}
#	$(GCC) -bundle -flat_namespace -undefined suppress -o $@ -O $^
else    
	$(GCC) -shared -O $^ $(FFTSL) -o $@
endif

clean:
	rm -f *Dict* *.o *.so reader apply_calibration analyse compare time_dep
