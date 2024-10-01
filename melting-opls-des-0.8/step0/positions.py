
"""

Autora: Gabriela Barreto Correa

"""


""" PACOTES """
import numpy as np

""" INPUTS """
N1 = int(500*1) #numero total de atomos de Cl
N2 = int(500*3) #numero total de atomos de CA
N3 = int(500*1) #numero total de atomos de CS
N4 = int(500*1) #numero total de atomos de CW
N5 = int(500*1) #numero total de atomos de N
N6 = int(500*1) #numero total de atomos de O
Nt =  int(500 + 500*21) #numero total de atomos
nframes = 1

""" INICIANDO MATRIZES E VETORES """
rv1 = np.zeros((nframes, N1, 3), dtype=np.float_)
rv2 = np.zeros((nframes, N2, 3), dtype=np.float_)
rv3 = np.zeros((nframes, N3, 3), dtype=np.float_)
rv4 = np.zeros((nframes, N4, 3), dtype=np.float_)
rv5 = np.zeros((nframes, N5, 3), dtype=np.float_)
rv6 = np.zeros((nframes, N6, 3), dtype=np.float_)

""" LEITURA DO CENTRO DE MASSA DAS MOLECULAS """
f = open("chcl.pdb", "r")

for frame in range(1):
    i1 = 0
    i2 = 0
    i3 = 0
    i4 = 0
    i5 = 0
    i6 = 0
    for n in range(1):
        f.readline()
    for n in range(Nt):
        a1, a2, ida, a4, a5, a6, rx, ry, rz, b1, b2 = f.readline().split() 
        if (int(ida) == 10): 
            rv1[frame, i1, :] = float(rx), float(ry), float(rz)
            i1 += 1
        if (int(ida) == 1): 
            rv2[frame, i2, :] = float(rx), float(ry), float(rz)
            i2 += 1
        elif (int(ida) == 2):
            rv3[frame, i3, :] = float(rx), float(ry), float(rz)
            i3 += 1   
        elif (int(ida) == 3):
            rv4[frame, i4, :] = float(rx), float(ry), float(rz)
            i4 += 1       
        elif (int(ida) == 7):
            rv5[frame, i5, :] = float(rx), float(ry), float(rz)
            i5 += 1      
        elif (int(ida) == 8):
            rv6[frame, i6, :] = float(rx), float(ry), float(rz)
            i6 += 1      
f.close()

f = open('positions_Cl.txt', 'w')
for n in range(N1):
    f.write('create_atoms \t %u \t single \t %f \t %f \t %f \t remap yes \t units box\n' % (11, rv1[0,n,0], rv1[0,n,1], rv1[0,n,2]))
f.close()

f = open('positions_CA.txt', 'w')
for n in range(N2):
    f.write('create_atoms \t %u \t single \t %f \t %f \t %f \t remap yes \t units box\n' % (12, rv2[0,n,0], rv2[0,n,1], rv2[0,n,2]))
f.close()

f = open('positions_CS.txt', 'w')
for n in range(N3):
    f.write('create_atoms \t %u \t single \t %f \t %f \t %f \t remap yes \t units box\n' % (13, rv3[0,n,0], rv3[0,n,1], rv3[0,n,2]))
f.close()

f = open('positions_CW.txt', 'w')
for n in range(N4):
    f.write('create_atoms \t %u \t single \t %f \t %f \t %f \t remap yes \t units box\n' % (14, rv4[0,n,0], rv4[0,n,1], rv4[0,n,2]))
f.close()

f = open('positions_N.txt', 'w')
for n in range(N5):
    f.write('create_atoms \t %u \t single \t %f \t %f \t %f \t remap yes \t units box\n' % (15, rv5[0,n,0], rv5[0,n,1], rv5[0,n,2]))
f.close()

f = open('positions_O.txt', 'w')
for n in range(N6):
    f.write('create_atoms \t %u \t single \t %f \t %f \t %f \t remap yes \t units box\n' % (16, rv6[0,n,0], rv6[0,n,1], rv6[0,n,2]))
f.close()





