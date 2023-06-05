## General Directions for how to run Fortran functions in Python

### 1. Install f2py

### 2. Use Makefile or f2py commands to compile Fortran functions into shared object file
Note that this will need to include all object files and modules, which can clutter up your working directory.
I recommend putting all object/mod files into a separate directory.

Use commands: 
```
python3 -m numpy.f2py -h get_va.pyf -m get_va getva.f90
```
This creates the signature file: "get_va.pyf", which allows us to see the parameters
and types of each variable.
Then:
```
python3 -m numpy.f2py -c -m get_va get_va.pyf getva.f90
```

This creates the shared object file, which we can then import into Python.

This works for simple Fortran functions, but for more complicated ones, I 
would recommend using a Makefile. See the Makefile in this directory for an example.
### 3. Import the shared object library into Python

Use command:
```
from psipy.fortran import get_va
```
Use the name from the below command as the function to import
```
python3 -m numpy.f2py -c -m get_va get_va.pyf getva.f90
```
### 4. Call the Fortran functions from Python

```
get_va.getva(br, bt, bp, rho, nr, nt, np)
```
"get_va" is the shared library object you import, while "getva" is the 
name of the actual function or subroutine in the Fortran code.