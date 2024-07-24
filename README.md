## Standard RAxML version 8.2.12

---

**When using RAxML please cite the following paper:**

_A. Stamatakis: "RAxML Version 8: A tool for Phylogenetic Analysis and Post-Analysis of Large Phylogenies"._
_Bioinformatics (2014) 30 (9): 1312-1313._

---

### Quick start under Linux:

standard, SSE3 or AVX version?

In general you should try to compile the SSE3 version that makes use of capabilities on relativiely recent processors (most Intel or AMD chips not older than 3-4 years).The SSE3 version will run about 40% faster than the non-SSE3 version.

If you have a more recently bought processor (within the last 1-2 years), please also try to compile the AVX version which once again uses the capabilities of modern processors better and can be 10-30% faster than the SSE3 version. AVX will work on the Intel i7 (sandy-bridge) series processors as well as on the very recently released AMD Bulldozer systems.

Compiling should work out of the box with all reasonably recent versions of the GNU gcc and Intel icc compilers. If you want to use icc replace gcc by icc in the Makefiles. Please direct all your RAxML questions to our google group (and only after having used the search function!): https://groups.google.com/forum/?hl=de&fromgroups#!forum/raxml

#### Sequential version:

type:

`make -f Makefile.gcc`

`rm *.o`

or

`make -f Makefile.SSE3.gcc`

`rm *.o`

or

`make -f Makefile.AVX.gcc`

`rm *.o`

#### Pthreads version:

type:

`make -f Makefile.PTHREADS.gcc`

`rm *.o`

or

`make -f Makefile.SSE3.PTHREADS.gcc`

`rm *.o`

or

`make -f Makefile.AVX.PTHREADS.gcc`

`rm *.o`

Coarse-grain MPI version:

type:

`make -f Makefile.MPI.gcc`

`rm *.o`

or

`make -f Makefile.SSE3.MPI.gcc`

`rm *.o`

or

`make -f Makefile.AVX.MPI.gcc`

`rm *.o`

#### Hybrid MPI/Pthreads version:

Before using this version, please read this paper here:

http://sco.h-its.org/exelixis/pubs/Exelixis-RRDR-2010-3.pdf

and look at these slides:

http://sco.h-its.org/exelixis/resource/doc/Phylo100225.pdf

type:

`make -f Makefile.HYBRID.gcc"`

`rm *.o`

or

`make -f Makefile.SSE3.HYBRID.gcc`

`rm *.o`

or

`make -f Makefile.AVX.HYBRID.gcc`

`rm *.o`
