import numpy as np
import os
import sys
import subprocess
import shutil

def zapis(NI,NJ,R_in,eps,max_iter,L,delta_r,h1,N_sec,P0,order_acc):
	with open('Input.txt', 'w', encoding='utf-8') as f:
		f.write(str(NI)+'	!NI' + '\n')
		f.write(str(NJ)+'	!NJ'+'\n')
		f.write(str(R_in)+'	!R_in'+'\n')
		f.write(str(eps)+'	!eps'+'\n')
		f.write(str(max_iter)+'	!max_iter'+'\n')
		f.write(str(L)+'	!L'+'\n')
		f.write(str(delta_r)+'	!delta_r'+'\n')
		f.write(str(h1)+'	!h1'+'\n')
		f.write(str(N_sec)+'	!N_sec'+'\n')
		f.write(str(P0)+'	!P0'+'\n')
		f.write(str(order_acc)+'	!order_acc'+'\n')


def viv(process):
    for line in process.stdout:
        print(line, end='', flush=True)
    process.wait()

def zapysk(NI,NJ,R_in,eps,max_iter,L,delta_r,h1,N_sec,P0,order_acc):
    zapis(NI,NJ,R_in,eps,max_iter,L,delta_r,h1,N_sec,P0,order_acc)
    process = subprocess.Popen(['./q'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1)
    viv(process)


#--------------базовый вариант----------------
def baza(NI=41,NJ=41,R_in=0.1,eps=1e-8,max_iter=50000,L=10,delta_r=0.1,h1=2,N_sec=4,P0=0,order_accs=[1,2],outdir='base'):
    os.makedirs(outdir, exist_ok=True)
    for order_acc in order_accs:
        os.makedirs(os.path.join(outdir, str(order_acc)), exist_ok=True)
    with open('force.plt', 'w', encoding='utf-8'):
        pass
    subprocess.run(['gfortran', 'main.f90','function.f90', 'solver.f90', '-o', './q'], check=True)
    for order_acc in order_accs:
        zapysk(NI,NJ,R_in,eps,max_iter,L,delta_r,h1,N_sec,P0,order_acc)
        shutil.copy('Input.txt', os.path.join(outdir,str(order_acc), 'Input.txt'))
        for name in ('residuals.plt', 'fields.plt'):
            shutil.move(os.path.join('output',name), os.path.join(outdir,str(order_acc), name))
        shutil.move('force.plt', os.path.join(outdir,'force.plt'))
#-------------------------------------------------

#-------------------параметрическое исследование-----------
def parametr(NI=41,NJ=41,R_in=0.1,eps=1e-8,max_iter=50000,L=10,delta_r=0.1,h1=2,N_sec=4,P0=0,order_acc=2,outdir='parametr',
    delta_rs=[0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3],N_secs=[1,2,3,4,5,6,7,8,9,10]):
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(os.path.join(outdir, 'delta_r'), exist_ok=True)
    os.makedirs(os.path.join(outdir, 'N_sec'), exist_ok=True)
    with open('force.plt', 'w', encoding='utf-8') as f:
        pass
    subprocess.run(['gfortran', 'main.f90','function.f90', 'solver.f90', '-o', './q'], check=True)
    for N_s in N_secs:
        zapysk(NI,NJ,R_in,eps,max_iter,L,delta_r,h1,N_s,P0,order_acc)
    shutil.move('force.plt', os.path.join(outdir,'N_sec','force.plt'))
    for delta in delta_rs:
        zapysk(NI,NJ,R_in,eps,max_iter,L,delta,h1,N_sec,P0,order_acc)
    shutil.move('force.plt', os.path.join(outdir,'delta_r','force.plt'))



#---------------------------------------------------

#---------------аналитика--------------------
def analitika(NI=41,NJ=41,R_in=0.1,eps=1e-8,max_iter=50000,L=1,delta_r=0.1,h1=2,N_sec=100,P0=0,order_acc=2,outdir='analitic'):
    os.makedirs(outdir, exist_ok=True)
    subprocess.run(['gfortran', 'main.f90','function.f90', 'solver.f90', '-o', './q'], check=True)
    zapysk(NI,NJ,R_in,eps,max_iter,L,delta_r,h1,N_sec,P0,order_acc)
    for name in ('residuals.plt', 'fields.plt'):
        shutil.move(os.path.join('output',name), os.path.join(outdir, name))
#-------------------------------------------------

#---------------порядок точности------------------
def poradok(NI=41,NJ=41,R_in=0.1,eps=1e-8,max_iter=50000,L=10,delta_r=0.1,h1=2,N_sec=4,P0=0,order_acc=2,
    outdir='poradok', outname='porad.plt', order_accs=[1,2], ni_fixed=21,nj_fixed=21, ni_list=[11,21,41],nj_list=[11,21,41]):
    os.makedirs(os.path.join(outdir, 'res'), exist_ok=True)
    os.makedirs(os.path.join(outdir, 'plt'), exist_ok=True)
    with open('force.plt', 'w', encoding='utf-8') as f:
        pass

    subprocess.run(['gfortran', 'main.f90','function.f90', 'solver.f90', '-o', './q'], check=True)

    for order_acc in order_accs:
        with open('force.plt', 'a', encoding='utf-8') as f:
            f.write('order_acc = ' + str(order_acc) + '\n')
        with open('force.plt', 'a', encoding='utf-8') as f:
            f.write('NI = ' + str(ni_fixed) + '\n')
        for NJ in nj_list:
            zapysk(ni_fixed,NJ,R_in,eps,max_iter,L,delta_r,h1,N_sec,P0,order_acc)

        with open('force.plt', 'a', encoding='utf-8') as f:
            f.write('nj = ' + str(nj_fixed) + '\n')
        for NI in ni_list:
            zapysk(NI,nj_fixed,R_in,eps,max_iter,L,delta_r,h1,N_sec,P0,order_acc)

    shutil.copy('force.plt', os.path.join(outdir, outname))
#------------------------------------------------------------
