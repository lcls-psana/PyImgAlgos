# pyimgalgos
A set of classes and methds to deal with LCLS data.

## Documentation
- Sphinx generated documentation: https://lcls-psana.github.io/pyimgalgos/
<!--- - GitHub Pages: https://github.com/lcls-psana/pyimgalgos/wiki --->

## Install in conda release directory

### Create conda release
See for detail Psana Developer Documentation [3] 
```
cd <my-conda-repo>
source conda_setup
```

### Clone package
**on pslogin:**
```
git clone https://github.com/lcls-psana/pyimgalgos.git
# or 
condarel --addpkg --name pyimgalgos --tag HEAD```
```

### Build 
```
scons
```
then run application(s)

## References
- [1] https://lcls-psana.github.io/pyimgalgos/
- [2] https://github.com/lcls-psana/pyimgalgos/wiki
- [3] https://confluence.slac.stanford.edu/display/PSDMInternal/Psana+Developer+Documentation

