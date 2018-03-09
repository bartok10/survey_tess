import numpy as np
############# function to select adjacent voids #############
def overlap(void):
    #v_0 = int(void[0])
    if len(void) > 3.:
        j = 1
        cont = 3
        voids = []
        while int(void[j]) != 0:
            v_j = np.zeros(int(void[j]))
            for i in range(len(v_j)):
                voids.append(int(void[cont+i]))
            cont = cont + len(v_j) + 2
            j = j + len(v_j) + 2
    elif len(void) == 3.:
        voids = -1
    return voids
################################################################
#
gal = np.loadtxt('../test_st/gal_pos.txt')
vol = np.loadtxt('../test_st/cat_zobov.ascii.vol')
zones = np.loadtxt('../test_st/cat_zobov.zone')
v_mem = []
v_vols = []
for i in range(int(max(zones))+1):
    v_mem.append(gal[zones==i])
    v_vols.append(vol[zones==i])
f = open('/home/jarmijo/test_st/cat_zobov.void','r')
af  = f.read().splitlines()
f.close()
apf = []
for i in range(len(af)):
    apf.append(af[i].split())
stvoids = []
stvol = []
vovp = []
for i in range(1,len(apf)):
    vovp.append(overlap(apf[i]))
for i in range(len(vovp)):
    if vovp[i] != -1 :
        tmpv=[];tmpvl=[]
        for j in range(len(vovp[i])):
            tmpv.append(v_mem[vovp[i][j]])
            tmpvl.append(v_vols[vovp[i][j]])
        tmpa = np.concatenate(tmpv,axis=0)
        tmpb = np.concatenate(tmpvl,axis=0)
        stvoids.append(np.concatenate([v_mem[i],tmpa]))
        stvol.append(np.concatenate([v_vols[i],tmpb]))
    elif vovp[i] == -1 :
        stvoids.append(v_mem[i])
        stvol.append(v_vols[i])
        
