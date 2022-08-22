import numpy as np

#Inputs
K = 100 #Base stiffnes
div = 20 #Number of subdivisions
r_ = (0.9,1.1) #distance for create bonds
potfile = "potfile_pert.dat"
chimapfile = "chimap_pert.dat"

with open(chimapfile,'w') as chif:
    chif.write('#chimap generated with input_gen.py\n\nnum_mat 1\n\n')

    chif.write('PROPERTIES: chi material type\n')
    for i in range(div+1):
        chif.write(f'{i/div} homo {i+1}\n')

with open(potfile,'w') as potf:
    potf.write('#potfile generated with input_gen.py\n\npair_style zero 1.0\npair_coeff * *\n\n')
    
    #create groups
    for i in range(div+1):
        potf.write(f'group {i+1} type {i+1}    #chi equal {i/div}\n')
    potf.write(f'group {div+2} type {div+2}    #non-opt\n')

    #define coefficents and create bonds
    potf.write('\nbond_style harmonic\n')
    bond_count = 0
    for i in range(div+1):
        for j in range(i,div+1):
            potf.write(f'bond_coeff {bond_count+1} {max(0.001,K*(i/div)**2*(j/div)**2):.4g} 1\n')
            bond_count += 1
    #non-opt bonds
    for i in range(div+1):
            potf.write(f'bond_coeff {bond_count+1} {max(0.001,K*(i/div)**2*(j/div)**2):.4g} 1\n')
            bond_count += 1
    potf.write(f'bond_coeff {bond_count+1} {K:.4g} 1\n')
    bond_count = 0
    for i in range(div+1):
        for j in range(i,div+1):
            potf.write(f'create_bonds many {i+1} {j+1} {bond_count+1} {r_[0]} {r_[1]}\n')
            bond_count += 1
    #non-opt bonds
    for i in range(div+1):
            potf.write(f'create_bonds many {i+1} {div+2} {bond_count+1} {r_[0]} {r_[1]}\n')
            bond_count += 1
    potf.write(f'create_bonds many {div+2} {div+2} {bond_count+1} {r_[0]} {r_[1]}\n')
    bond_count += 1
    potf.write('\n')

    #delete groups
    for i in range(div+1):
        potf.write(f'group {i+1} delete\n')

    #print stats
    potf.write(f'\n#{div+2} particle types\n#{bond_count} bonds created')
