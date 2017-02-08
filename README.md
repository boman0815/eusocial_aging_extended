# eusocial_aging_extended
 version for simulation runs

## build script

```
#!/bin/bash
# Run this script in the folder with the source files.
# Perhaps you first need to:
#
#    chmod +x build
#
# This allows execution of the script.
# Then running it:
#
#    ./build
#
# This will produce an executable called 'my_exe'


# Load Intel C++ compiler
module load intel

# Compile with Intel compiler
icc -o my_exe -std=c++11 colony.cpp corand.cpp simulation.cpp individual.cpp utils.cpp linalg.cpp random.cpp










# Load GCC 5.1.0
# module load GCC/5.1.0
# Build the executable using GCC
# g++ -o my_exe -std=c++11 colony.cpp corand.cpp simulation.cpp individual.cpp utils.cpp linalg.cpp random.cpp
```
