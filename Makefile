GXX = g++
SLD = 
DLD = -lpthread
OBJ = pssw/ssw_cpp.cpp pssw/ssw.c sRNATarPredictor.cpp

qTar : ${OBJ}
	${GXX} -o qTar ${OBJ} ${DLD} ${SLD}

