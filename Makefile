# -*- Makefile -*-

arch = Linux_MPI
setup_file = setup/Make.$(arch)

include $(setup_file)


HPCG_DEPS = obj/CG.o \
	    obj/CG_ref.o \
	    obj/ComputeResidual.o \
	    obj/ExchangeHalo.o \
	    obj/GenerateGeometry.o \
	    obj/GenerateProblem.o \
	    obj/GenerateProblem_ref.o \
	    obj/MixedBaseCounter.o \
	    obj/SetupHalo.o \
	    obj/SetupHalo_ref.o \
	    obj/ComputeDotProduct.o \
	    obj/ComputeDotProduct_ref.o \
	    obj/mytimer.o \
	    obj/ComputeSPMV.o \
	    obj/ComputeSPMV_ref.o \
	    obj/ComputeSYMGS.o \
	    obj/ComputeSYMGS_ref.o \
	    obj/ComputeWAXPBY.o \
	    obj/ComputeWAXPBY_ref.o \
	    obj/ComputeMG_ref.o \
	    obj/ComputeMG.o \
	    obj/ComputeProlongation_ref.o \
	    obj/ComputeRestriction_ref.o \
	    obj/GenerateCoarseProblem.o \
	    obj/init.o \
	    obj/finalize.o

# These header files are included in many source files, so we recompile every file if one or more of these header is modified.
PRIMARY_HEADERS = src/Geometry.hpp src/SparseMatrix.hpp src/Vector.hpp src/CGData.hpp \
                  src/MGData.hpp src/hpcg.hpp

all: bin/xhpcg

bin/xhpcg: obj/main.o $(HPCG_DEPS)
	$(LINKER) $(LINKFLAGS) obj/main.o $(HPCG_DEPS) $(HPCG_LIBS) -o bin/xhpcg

clean:
	rm -f obj/*.o bin/xhpcg

.PHONY: all clean

obj/main.o: src/main.cpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -Isrc $< -o $@

obj/CG.o: src/CG.cpp src/CG.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -Isrc $< -o $@

obj/CG_ref.o: src/CG_ref.cpp src/CG_ref.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -Isrc $< -o $@

obj/ComputeResidual.o: src/ComputeResidual.cpp src/ComputeResidual.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -Isrc $< -o $@

obj/ExchangeHalo.o: src/ExchangeHalo.cpp src/ExchangeHalo.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -Isrc $< -o $@

obj/GenerateGeometry.o: src/GenerateGeometry.cpp src/GenerateGeometry.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -Isrc $< -o $@

obj/GenerateProblem.o: src/GenerateProblem.cpp src/GenerateProblem.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -Isrc $< -o $@

obj/GenerateProblem_ref.o: src/GenerateProblem_ref.cpp src/GenerateProblem_ref.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -Isrc $< -o $@

obj/MixedBaseCounter.o: src/MixedBaseCounter.cpp src/MixedBaseCounter.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -Isrc $< -o $@

obj/SetupHalo.o: src/SetupHalo.cpp src/SetupHalo.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -Isrc $< -o $@

obj/SetupHalo_ref.o: src/SetupHalo_ref.cpp src/SetupHalo_ref.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -Isrc $< -o $@

obj/ComputeDotProduct.o: src/ComputeDotProduct.cpp src/ComputeDotProduct.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -Isrc $< -o $@

obj/ComputeDotProduct_ref.o: src/ComputeDotProduct_ref.cpp src/ComputeDotProduct_ref.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -Isrc $< -o $@

obj/finalize.o: src/finalize.cpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -Isrc $< -o $@

obj/init.o: src/init.cpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -Isrc $< -o $@

obj/mytimer.o: src/mytimer.cpp src/mytimer.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -Isrc $< -o $@

obj/ComputeSPMV.o: src/ComputeSPMV.cpp src/ComputeSPMV.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -Isrc $< -o $@

obj/ComputeSPMV_ref.o: src/ComputeSPMV_ref.cpp src/ComputeSPMV_ref.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -Isrc $< -o $@

obj/ComputeSYMGS.o: src/ComputeSYMGS.cpp src/ComputeSYMGS.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -Isrc $< -o $@

obj/ComputeSYMGS_ref.o: src/ComputeSYMGS_ref.cpp src/ComputeSYMGS_ref.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -Isrc $< -o $@

obj/ComputeWAXPBY.o: src/ComputeWAXPBY.cpp src/ComputeWAXPBY.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -Isrc $< -o $@

obj/ComputeWAXPBY_ref.o: src/ComputeWAXPBY_ref.cpp src/ComputeWAXPBY_ref.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -Isrc $< -o $@

obj/ComputeMG_ref.o: src/ComputeMG_ref.cpp src/ComputeMG_ref.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -Isrc $< -o $@

obj/ComputeMG.o: src/ComputeMG.cpp src/ComputeMG.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -Isrc $< -o $@

obj/ComputeProlongation_ref.o: src/ComputeProlongation_ref.cpp src/ComputeProlongation_ref.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -Isrc $< -o $@

obj/ComputeRestriction_ref.o: src/ComputeRestriction_ref.cpp src/ComputeRestriction_ref.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -Isrc $< -o $@

obj/GenerateCoarseProblem.o: src/GenerateCoarseProblem.cpp src/GenerateCoarseProblem.hpp $(PRIMARY_HEADERS)
	$(CXX) -c $(CXXFLAGS) -Isrc $< -o $@


