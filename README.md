# dissociation_heom
compiling: 
f2py -m sparsity -c sparsity.f90 --fcompiler=intelem
f2py -m hierarchy_bose -c hierarchy_bose.f90 --fcompiler=intelem
pgf90 -o propagation vals.cuf propagation.cuf

run:
python main.py
./propagation
