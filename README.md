**MicroSCOPE** 

_MicroBooNE Searches for Curious Occurrences, Particles, and Exotics_

**Goal:**

To provide a unified analysis framework for BSM analyses, flexible to the needs of the analyser. Through an abstract base "module" class, which performs generic and modifiable methods on the data, and a "module manager" which handles the execution of these methods, running any number of modules sequentially can be performed through lightweight top-level executable scripts and configuration files.

**Requirements:**

ROOT 6.20 or higher, c++ 17 or higher (may work with older versions, but I haven't tested this)

**Building Instructions:**

Building of makefiles is done using cmake.

1) In the top directory `mkdir build`
2) `cd build`
3) `cmake ..`
4) `make`
