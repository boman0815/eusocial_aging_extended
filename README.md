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
```

## `run`

Just call this one:

```
#!/bin/bash
# Run this script in the folder with the executable
# Perhaps you first need to:
#
#    chmod +x run
#
# This allows execution of this script.
# Then running it:
#
#    ./run
#
# This will put 20 (or so) jobs request in the queue

for i in `seq 1 20`;
do
  sbatch run_one
done
```

It will call `run_one`.

## `run_one`

```
#!/bin/bash
# Call this file using
#
#  sbatch run_one
#
# Do not do this as such, prefer
# using the 'run' script:
#
# ./run
#
# This 'run' script will call this one for you :-)
#

#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1 
#SBATCH --mem=10GB
#SBATCH --job-name=run_one
#SBATCH --output=run_one-%j.log
./my_exe
```

## Troubleshooting

### Making scripts execute

To make `my_file` an executable scuipt: 
```
chmod +x my_file
```

Where ``my_file` can be `run`, `build` or whatever.

Useful oneliner:

```
chmod +x run; chmod +x run_one; chmod +x build
```


### Windows line endings

File cannot be read, due to Windows line endings.

```
dos2unix my_file
```

Where ``my_file` can be `run`, `build` or whatever.

Useful oneliner:

```
dos2unix run; dos2unix run_one; dos2unix build
```

### Renaming text files

Executable scripts with a `.txt` extension would be considered misleading on GNU/Linux.
Here a one-liner to rename all these.

```
mv build.txt build; mv run.txt run; mv run_one.txt run_one
```
### preliminaries for the cluster

error: 

./my_exe: /usr/lib64/libstdc++.so.6: version GLIBCXX_3.4.19' not found (required by ./my_exe).

./my_exe: /usr/lib64/libstdc++.so.6: version CXXABI_1.3.8' not found (required by ./my_exe).

load module before running the executable.

```
module load GCC/5.1.0
```
