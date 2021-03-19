# provie -U for Hubbard U calculations. Good for any TX2 compounds The first 2n-2 arguments are poscars and outcars of corresponding doped TX2 compounds, the last 2 arguments are bulk TX2 poscar and outcar. poscar always before outcar. nfile=n (n-1 no of doped structures and 1 bulk)

python ../src/formation_energy.py doped_TMD/POSCAR doped_TMD/OUTCAR bulk_TMD/POSCAR bulk_TMD/OUTCAR --mu 0 --mat V --nfile 2 -typ TX2 -U 0 --font 16 --legend V

