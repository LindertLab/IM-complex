gmx_gpu grompp -f minim.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 3
gmx_gpu mdrun -v -deffnm em -nb gpu
gmx_gpu grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 3
gmx_gpu mdrun -deffnm nvt -nb gpu
gmx_gpu grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 3
gmx_gpu mdrun -deffnm npt -nb gpu
gmx_gpu grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_10.tpr -maxwarn 3
gmx_gpu mdrun -deffnm md_10 -nb gpu