
CFLAGS	         =
FFLAGS	         =
CPPFLAGS         =
FPPFLAGS         =
LOCDIR           = 
EXAMPLESC        = main.cpp steady.cpp obs.cpp
EXAMPLESF        = 
MANSEC           = KSP
CLEANFILES       = rhs.vtk solution.vtk
NP               = 1
OBJ=$(EXAMPLESC:.cpp=.o)

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

EXE = SteadyState
${EXE}: ${OBJ}  chkopts
	-${CLINKER} -o ${EXE} ${OBJ}  ${PETSC_KSP_LIB}
#	${RM} ${OBJ}

#EXE = steady
#${EXE}: ${EXE}.o  chkopts
#	-${CLINKER} -o ${EXE} ${EXE}.o  ${PETSC_KSP_LIB}
#	${RM} ${EXE}.o

touch:
	touch *.cpp *.h
