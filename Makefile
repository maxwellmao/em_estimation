include Makefile.config

EXEC = em_estimation

all: $(EXEC)

em_estimation: em_estimation.h em_estimation.cpp vec_mat_operation.h common.h
	$(CC) $(CXXFLAGS) $(CXXOPENMP) -o em_estimation em_estimation.cpp $(LDFLAGS) $(LIBS)

