from function.function import *

os.makedirs('input', exist_ok=True)
os.makedirs('output', exist_ok=True)

#baza(NI=41,NJ=41,R_in=0.1,eps=1e-8,max_iter=50000,L=10,delta_r=0.3,h1=2)

#analitika(NI=41,NJ=41,R_in=0.1,eps=1e-8,max_iter=50000,L=1,delta_r=0.3,h1=2,N_sec=100,P0=0,order_acc=2)

poradok(nj_fixed=41, ni_list=[11,21,41,81,161,321], ni_fixed=41,nj_list=[11,21,41,81,161,321],delta_r=0.3)

#parametr(NI=81,NJ=81,R_in=0.1,eps=1e-8,max_iter=500000,L=10,delta_r=0.3,h1=2,N_sec=4,order_acc=2,
#    delta_rs=[0,0.025,0.05,0.075,0.1,0.15,0.2,0.25,0.3],N_secs=[1,2,3,4,5,6,7,8,9,10])
