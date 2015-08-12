
CFLAGS	         =
FFLAGS	         =
CPPFLAGS         =
FPPFLAGS         =
LOCDIR           = 
EXAMPLESC        = steady.cpp ex23.c ex11.c
EXAMPLESF        = 
MANSEC           = KSP
CLEANFILES       = rhs.vtk solution.vtk
NP               = 1

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

steady: steady.o  chkopts
	-${CLINKER} -o steady steady.o  ${PETSC_KSP_LIB}
	${RM} steady.o

ex23: ex23.o  chkopts
	-${CLINKER} -o ex23 ex23.o  ${PETSC_KSP_LIB}
	${RM} ex23.o

ex11: ex11.o  chkopts
	-${CLINKER} -o ex11 ex11.o  ${PETSC_KSP_LIB}
	${RM} ex11.o
