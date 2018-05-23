# Blockphase
Blockphase is a software for accurate genotype calling and haplotype phasing. Our newly designed approach constructed two steps: first, we built a haplotype graph(HG) by randomly dividing the haplotypes into blocks with disjoint contiguous markers, treat the same distinct haplotypes as one which is the key feature of reducing redundancies states. Then, we used backward sampling to update one individual by its genotype likelihood. 

# Building Blockphase
Blockphase requires the Boost C++ library, please add them to your enviroment path:

export PATH=/path/to/boost/include:$PATH

export LD_LIBRARY_PATH=/path/to/boost/lib:$LD_LIBRARY_PATH

git clone https://github.com/songli-bgi/Blockphase.git

Get in to the "Blcokphase" directory. 

cd Blockphase

If you had beed set up the boost library path, you can type "make" in your command line

make

# Usage
You can get the help information by:

./blockphase -h/--help 

examples:

./blockphase -i $infile -o $outfile 

./blockphase -i $infile -o $outfile -b 5 -c 15 -s 20
