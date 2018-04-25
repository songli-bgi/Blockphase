object=blockphase
Ofile=main.o Hap_Init.o GenoHap_Init.o Haplotype_graph.o Init_Type.o HMM_Algo.o Thread_Pool.o Run_Wizzard.o Arg_Wizzard.o anyarg.o Random_Wizzard.o Output_Wizzard.o sq_Pearson_cc.o freader.o Hap_Conjunct.o Gibbsgeno_Store.o EM_hfs.o Sample_Orders.o

flag= -O2
#MACRO= -DDEBUG
MACRO= 

GCC=/usr/bin/g++

all:$(Ofile)
	$(GCC) $(flag) -o $(object) $(Ofile) $(MACRO) -lpthread -lz
main.o:main.cpp options.h
	$(GCC) $(flag) -c main.cpp $(MACRO)
Hap_Init.o:Hap_Init.cpp
	$(GCC) $(flag) -c Hap_Init.cpp $(MACRO)
GenoHap_Init.o:GenoHap_Init.cpp
	$(GCC) $(flag) -c GenoHap_Init.cpp $(MACRO)
Haplotype_graph.o:Haplotype_graph.cpp
	$(GCC) $(flag) -c Haplotype_graph.cpp $(MACRO)
Init_Type.o:Init_Type.cpp
	$(GCC) $(flag) -c Init_Type.cpp $(MACRO)
HMM_Algo.o:HMM_Algo.cpp
	$(GCC) $(flag) -c HMM_Algo.cpp $(MACRO)
Thread_Pool.o:Thread_Pool.cpp
	$(GCC) $(flag) -c Thread_Pool.cpp $(MACRO)
Run_Wizzard.o:Run_Wizzard.cpp
	$(GCC) $(flag) -c Run_Wizzard.cpp $(MACRO)
Arg_Wizzard.o:Arg_Wizzard.cpp options.h
	$(GCC) $(flag) -c Arg_Wizzard.cpp $(MACRO)
anyarg.o:anyarg.cpp
	$(GCC) $(flag) -c anyarg.cpp $(MACRO)
Random_Wizzard.o:Random_Wizzard.cpp
	$(GCC) $(flag) -c Random_Wizzard.cpp $(MACRO)
Output_Wizzard.o:Output_Wizzard.cpp
	$(GCC) $(flag) -c Output_Wizzard.cpp $(MACRO)
sq_Pearson_cc.o:sq_Pearson_cc.cpp
	$(GCC) $(flag) -c sq_Pearson_cc.cpp $(MACRO)
freader.o:freader.cpp
	$(GCC) $(flag) -c freader.cpp $(MACRO)
Hap_Conjunct.o:Hap_Conjunct.cpp
	$(GCC) $(flag) -c Hap_Conjunct.cpp $(MACRO)
Gibbsgeno_Store.o:Gibbsgeno_Store.cpp
	$(GCC) $(flag) -c Gibbsgeno_Store.cpp $(MACRO)
EM_hfs.o:EM_hfs.cpp
	$(GCC) $(flag) -c EM_hfs.cpp $(MACRO)
Sample_Orders.o:Sample_Orders.cpp
	$(GCC) $(flag) -c Sample_Orders.cpp $(MACRO)
clean:
	-rm -f $(object) $(Ofile)
