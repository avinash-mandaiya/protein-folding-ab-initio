PROGRAM = NIR4
FILES.c = NIR4.c motifproj.c 
FILES.h = pfs.h 
FILES.o = ${FILES.c:.c=.o}
	
CC      = gcc
WFLAGS  = -Wall
SFLAGS  = -std=c11
OFLAGS  = -O
MFLAGS  = -fsanitize=address
CFLAGS  = $(MFLAGS) ${WFLAGS} ${SFLAGS} ${OFLAGS} 
LDFLAGS = -lm

all: ${PROGRAM} 

${PROGRAM}: ${FILES.o}
	${CC} -o $@ ${CFLAGS} ${FILES.o} ${LDFLAGS}

NIR4.o: ${FILES.h}
motifproj.o:    ${FILES.h}

ProID        = 1LB0
NumSol       = 10
EGS          = -200.0
CHstep       = 10
id           = t24
beta         = 0.5
maxiter      = 100000
iterstride   = 200
stoperr      = 0.0005
epsilon      = 0.0001
seed         = 2
RAD          = 15.0
HB1          = 7
HB2          = 5
epsRR        = 1.
epsRS        = 1.
epsSS        = 1.
ESCALE       = 1.0
BRAD         = 0.3

run: ${PROGRAM}
	$(addprefix ./,${PROGRAM}) ${ProID} ${NumSol} ${EGS} ${CHstep} ${id} ${beta} ${maxiter} ${iterstride} ${stoperr} ${epsilon} ${seed} ${RAD} ${HB1} ${HB2} ${epsRR} ${epsRS} ${epsSS} ${ESCALE} ${BRAD}

.PHONY: run 
