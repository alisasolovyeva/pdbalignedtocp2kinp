"""
The purpose of this script is to generate cp2k input files, where
different ligands are placed into a binding site of a fixed protein.

The pdb codes are aligned using Biopython (https://biopython.org/).
The input files correspond to truncated pdb files, where a protein
environment is the same. In addition, cp2k input templates and 
submission template are required.

The following output will be created:
    Directories named after a fixed protein with subdirectories named
    after pdb codes from which ligands are taken.
    Each "ligand" subdirectory contains two subdirectories with input
    files for cp2k DFT calculation for ligand and PM6 calculation of 
    the complex.

The script is written by Alisa Solovyeva in 2018.

"""

import os, sys
import Bio.PDB
import numpy

mypath=os.getcwd()

aligned_to = ["6frf","6fqo","5nlk"]
pdb_code   = ["6frf","5nlk"]
ligands    = ['92E', 'UL4', 'YE5', 'XZ8', '5QN', '8Q6', 'QPR', '0BC', 
        'E3T', 'UT0', 'PKU', 'E3B', '4I8', 'F31', '5QR', '5J5', 
        '7MX', '5XS', 'E2T','E3H']
parts      = ["ligand","complex"]

for item in aligned_to:
    pdb_filename = '{0}{1}'.format(item,".pdb")
    ref_model = Bio.PDB.PDBParser().get_structure(item, pdb_filename)
    ref_model = ref_model[0]
    for align in pdb_code:
        if item!=align:
            coords = []
            alt_filename = '{0}{1}'.format(align,".pdb")
            pdb_out_filename = '{0}_{1}_{2}{3}'.format(align,'aligned_to',
                    item,'.pdb')
            alt_model = Bio.PDB.PDBParser().get_structure(align, alt_filename)
            alt_model=alt_model[0]
            ref_atoms = []
            alt_atoms = []
            altligat  = []
            for (ref_chain, alt_chain) in zip(ref_model, alt_model) :   
                for ref_res in ref_chain:
                    for alt_res in alt_chain:
                        if (ref_res.resname == alt_res.resname) and (ref_res.id == alt_res.id):
                            for (at, altat) in zip(ref_res, alt_res):
                                coords.append('{0}  {1} {2} {3} {4}'.format(
                                    (list(at.get_id()))[0], (at.get_coord())[0],
                                    (at.get_coord())[1],(at.get_coord())[2],'\n'))
                                if (list(ref_res.id))[0] != 'H':
                                    ref_atoms.append(at)
                                    alt_atoms.append(altat)
                for alt_res in alt_chain:
                    for lig in ligands:
                        if alt_res.resname == lig:
                            for ligat in alt_res:
                                altligat.append(ligat)
           #-----------------align these paired atom lists------------- 
            super_imposer = Bio.PDB.Superimposer()
            super_imposer.set_atoms(ref_atoms, alt_atoms)
            super_imposer.apply(alt_model.get_atoms())
            io=Bio.PDB.PDBIO()
            io.set_structure(alt_model)
            io.save(pdb_out_filename)

            super_imposer.set_atoms(ref_atoms, alt_atoms)
            super_imposer.apply(altligat)
            ligand_coords=[]
            for atom in altligat:
                ligcoord=atom.get_coord()
                ligatname=atom.get_id()
                ligand_coords.append('{0}  {1} {2} {3}'.format((list(ligatname))[0],
                    ligcoord[0],ligcoord[1],ligcoord[2]))
                coords.append('{0}  {1} {2} {3} {4}'.format((list(ligatname))[0],
                    ligcoord[0],ligcoord[1],ligcoord[2],'\n'))                
            aname=[]
            x=[]
            y=[]
            z=[]
           #-----------------move to center of 55 A cell--------------- 
            for line in coords:
                aname.append(line.split()[0])
                x.append(float(line.split()[1]))
                y.append(float(line.split()[2]))
                z.append(float(line.split()[3]))
            min_x = float(min(x))
            min_y = float(min(y))
            min_z = float(min(z))
            max_x = float(max(x))
            max_y = float(max(y))
            max_z = float(max(z))
            mol_d = ((max_x-min_x)+(max_y-min_y)+(max_z-min_z))/3.0
            new_x=[]
            new_y=[]
            new_z=[]
            abstand=(float(55.0)-mol_d)/2.0
            for j,k,h in zip(x,y,z):
                new_x.append((float(j)+(abstand-min_x)))
                new_y.append((float(k)+(abstand-min_y)))
                new_z.append((float(h)+(abstand-min_z)))
            coord_in_center = []
            for el in range(len(aname)):
                coord_in_center.append('{0}  {1}  {2}  {3}'.format(aname[el], 
                    new_x[el], new_y[el], new_z[el]))
            #----------------create directory--------------------------
            cp2kpath=[]
            for n in parts:
                cp2kcomplex_path = "{0}/{1}_{2}/{3}/{4}/".format(mypath,'aligned_to',item,align,n)
                cp2kpath.append("{0}/{1}_{2}/{3}/{4}/".format(mypath,'aligned_to',item,align,n))
                directory = os.path.dirname(cp2kcomplex_path)
                if not os.path.exists(directory):
                    try:
                        os.makedirs(directory)
                    except OSError:
                        print ("Creation of the directory %s failed" % directory)
                    else:
                        print ("Successfully created the directory %s" % directory)
            #---------------write cp2k complex input------------------
            template=open("pm6.inp","r")
            inp=open('{0}/{1}{2}'.format(cp2kpath[1],align,'.inp'),'w')
            for line in template:
                if line.split()[0] =='PROJECT':
                    inp.write('{0}  {1} {2}'.format('  PROJECT',align,'\n'))
                elif line.split()[0] =='&COORD':
                    inp.write('{0}'.format(line))
                    for l in coord_in_center:
                        inp.write('{0} {1}'.format(l,'\n'))
                elif line.split()[0] =='CHARGE' and (align=='4tqn' or align=='5ep7' or align=='5mpk' or align=='5mpn'):
                    inp.write('{0}'.format('    CHARGE -2\n'))
                elif line.split()[0] =='CHARGE' and (align=='5eng'):
                    inp.write('{0}'.format('    CHARGE 1\n'))
                elif line.split()[0] =='CHARGE' and (align=='6frf' or align=='6fqo' or align=='6fqu' or align=='5nlk' or align=='5h85' or align=='5mqg' or align=='5owk' or align=='5mqk'):
                    inp.write('{0}'.format('    CHARGE -1\n'))
                else:
                    inp.write('{0}'.format(line))
            template.close()
            inp.close()
            runinp=open("run","r")
            cmxrun=open('{0}/{1}'.format(cp2kpath[1],'run'),'w')
            for line in runinp:
                if len(line.split( ))>1 and line.split( )[1]=='--job-name=':
                    cmxrun.write('{0}{1}{2}{3}'.format('#SBATCH --job-name=',align,'complex','pm6 \n'))
                elif len(line.split( ))>1 and line.split( )[1]=='--nodes=18':
                    cmxrun.write('{0}'.format('#SBATCH --nodes=4 \n'))
                elif len(line.split( ))>1 and line.split()[0]=='srun':
                    cmxrun.write('{0} {1}{2} {3} {4}{5}'.format('srun -n $SLURM_NTASKS --ntasks-per-node=$SLURM_NTASKS_PER_NODE -c $SLURM_CPUS_PER_TASK cp2k.psmp',align,'.inp','>',align,'.out \n'))
                else:
                    cmxrun.write('{0}'.format(line))
            cmxrun.close()        
            runinp.close()
            #---------------write cp2k ligand input------------------
            if not os.path.exists('{0}/{1}/'.format(cp2kpath[0],'dft_sp')):
                try:
                    os.makedirs('{0}/{1}/'.format(cp2kpath[0],'dft_sp'))
                except OSError:
                    print ("Creation of the directory failed")
                else:
                    print ("Successfully created dir")
            #---------------write cp2k complex input------------------
            dftinp=open("dft.inp","r")
            liginp=open('{0}/{1}/{2}{3}'.format(cp2kpath[0],'dft_sp',align,'.inp'),'w')
            for line in dftinp:
                if line.split()[0] =='PROJECT':
                    liginp.write('{0}  {1} {2}'.format('  PROJECT',align,'\n'))
                elif line.split()[0] =='&COORD':
                    liginp.write('{0}'.format(line))
                    for coord in ligand_coords:
                        liginp.write('{0}  {1} {2} {3} {4}'.format(coord.split()[0],
                            (float(coord.split()[1])+(abstand-min_x)),
                            (float(coord.split()[2])+(abstand-min_y)),
                            (float(coord.split()[3])+(abstand-min_z)),'\n'))
                elif line.split()[0] =='CHARGE' and (align=='4tqn' or align=='5ep7' or align=='5mpk' or align=='5mpn'):
                    liginp.write('{0}'.format('    CHARGE -1\n'))
                elif line.split()[0] =='CHARGE' and (align=='5eng'):
                    liginp.write('{0}'.format('    CHARGE 2\n'))
                elif line.split()[0] =='CHARGE' and (align=='6frf' or align=='6fqo' or align=='6fqu' or align=='5nlk' or align=='5h85' or align=='5mqg' or align=='5owk' or align=='5mqk'):
                    liginp.write('{0}'.format('    CHARGE 0\n'))
                else:
                    liginp.write('{0}'.format(line))

            liginp.close()
            dftinp.close()
            runinp=open("run","r")
            ligrun=open('{0}/{1}/{2}'.format(cp2kpath[0],'dft_sp','run'),'w')
            for line in runinp:
                if len(line.split( ))>1 and line.split( )[1]=='--job-name=':
                    ligrun.write('{0}{1}{2}{3}'.format('#SBATCH --job-name=',align,'ligand','dft \n'))
                elif len(line.split( ))>1 and line.split()[0]=='srun':
                    ligrun.write('{0} {1}{2} {3} {4}{5}'.format('srun -n $SLURM_NTASKS --ntasks-per-node=$SLURM_NTASKS_PER_NODE -c $SLURM_CPUS_PER_TASK cp2k.psmp',align,'.inp','>',align,'.out \n'))
                else:
                    ligrun.write('{0}'.format(line))
            runinp.close()
            ligrun.close()
