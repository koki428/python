# R2D2_py
Python code for analyzing  the result of R2D2

# Upgrade resolution

if you want to use upgrade resolution procedure

- Oakforest-PACS

```
export FFLAGS=
f2py --fcompiler=ifort -m regrid -c --opt='-fPIC' regrid.f90
```

- Ubuntu 18.04 server

```
f2py --fcompiler=gfortran -m regrid -c --opt='-O3' regrid.f90
```
