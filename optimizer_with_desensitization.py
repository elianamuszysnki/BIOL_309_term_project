#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 22 14:54:36 2020

@author: elianamuszynski
"""

import matplotlib.pyplot as plt
from scipy.optimize import minimize

#Instructions:

# (1) To obtain 2P simulation using Table 1 Parameters: line 312 = False 
# (2) To obtain 2P simulation w/ Nedler Mead Optimization (w/o desensitzation): line 312 = True
# (3) To obtain 2P simulation w/ Nedler Mead Optimization (w/ desensitzation) :line 312 = True and uncomment line 260


#parameters to esimtate 

# unaltered betwen WT / KO 
tauCa_init = 30.1 #milliseconds
kCa_init = 1.19
ko_init =0.0009  #/ms
tin_init = 3.0  #ms
Ph_init =0.17

kenter_init = 1
trecovery_init = 200

# WT
WT_kmax_init = 0.026#/ms 
WT_kr_init= 4.05
WT_nr_init= 1 
WT_Kf_init = 2.4 
WT_nf_init = 1.15

# KO
KO_kmax_init = 0.00825 #/ms  
KO_kr_init= 0.65
KO_nr_init= 3.92 
KO_Kf_init = 19.4
KO_nf_init = 1.5

WT_array_PPR_5= []
WT_array_PPR_10= []
WT_array_PPR_20= []
WT_array_PPR_50= []

KO_array_PPR_5= []
KO_array_PPR_10= []
KO_array_PPR_20= []
KO_array_PPR_50= []

#WT experimental data

WT_5hz = [(1,1.00),
            (2,1.35),
            (3,1.08),
            (4,1.36),
            (5,1.13),
            (6,1.20),
            (7,1.24),
            (8,1.30),
            (9,1.25),
            (10,1.21)]

WT_10hz = [(1,1.00),
            (2,1.40),
            (3,1.31),
            (4,1.26),
            (5,1.19),
            (6,1.25),
            (7,1.23),
            (8,1.30),
            (9,1.34),
            (10,1.23)]

WT_20hz = [(1,1.00),
            (2,1.96),
            (3,1.83),
            (4,1.61),
            (5,1.49),
            (6,1.50),
            (7,1.60),
            (8,1.52),
            (9,1.51),
            (10,1.58)]

WT_50hz = [(1,1.00),
            (2,2.60),
            (3,2.78),
            (4,1.96),
            (5,1.56),
            (6,1.45),
            (7,1.44),
            (8,1.47),
            (9,1.50),
            (10,1.40)]

#KO experimental data

KO_5hz =  [(1,1.00),
            (2,0.69),
            (3,0.60),
            (4,0.65),
            (5,0.67),
            (6,0.65),
            (7,0.66),
            (8,0.68),
            (9,0.73),
            (10,0.68)]



KO_10hz =  [(1,1.00),
            (2,0.73),
            (3,0.68),
            (4,0.63),
            (5,0.64),
            (6,0.66),
            (7,0.67),
            (8,0.64),
            (9,0.68),
            (10,0.66)]

KO_20hz =  [(1,1.00),
            (2,0.67),
            (3,0.61),
            (4,0.58),
            (5,0.57),
            (6,0.58),
            (7,0.63),
            (8,0.62),
            (9,0.54),
            (10,0.59)]

KO_50hz =  [(1,1.00),
            (2,0.82),
            (3,0.58),
            (4,0.55),
            (5,0.52),
            (6,0.56),
            (7,0.58),
            (8,0.55),
            (9,0.56),
            (10,0.51)]

def calculate_XYZ(X,Y,Z,Pr,Krec,t,tAP,tin):
   
    #dX
    dX = Krec * Z
    if t % tAP == 0 :
        dX += -Pr*X

    #dY
    dY = - Y / tin
    if t % tAP == 0 :
        dY += Pr*X
   
   
    #dZ
    dZ = Y / tin - Krec * Z
   
    X += dX
    Y += dY
    Z += dZ
   
    return (X,Y,Z)

def calculate_R(R, t, tAP, input_strength, kenter, trecovery):
    
    D = 1 - R
    dR = D / trecovery 
    if t % tAP == 0:
        dR -= kenter * R * input_strength
        
    R += dR
        
    return R

def run_simulation_reutrn_error(is_KO, x_vector, experimental_data, frequency):
    
    
    
    tauCa = x_vector[0] #milliseconds
    kCa = x_vector[1]
    ko = x_vector[2]#/ms
    tin = x_vector[3]
    Ph = x_vector[4]
    
    kenter = x_vector[5]
    trecovery = x_vector[6]
    
    if is_KO: #is true

        kmax = x_vector[7]
        kr = x_vector[8]
        nr = x_vector[9]
        Kf = x_vector[10]
        nf = x_vector[11]
        
    else:
        
        kmax = x_vector[12] 
        kr = x_vector[13]
        nr = x_vector[14]
        Kf = x_vector[15]
        nf = x_vector[16]
    
    Ca = 0.000001 #first impulse has already occured at t=0

   
    tAP= 1000 / frequency
    t_sim = int(tAP * 10) #show course of time over 10 pulses, match scale
    
    Poh=0.55 # high release probability --> constant
    Pol=0.03 # low release probability --> constant
    
    
    X1=Ph # Initial value for high sites
    Y1=0
    Z1=0
    
    X2=1 - Ph # Initial value for low sites
    Y2=0
    Z2=0
    
    R=1 #ready receptors
    
    function_array_PPR=[]
   

    #run simulation   
   
    for t in range (0,t_sim) :
        #Pr
        Pr1 = (1-Poh) / (1+ (Kf/Ca)**nf) + Poh
        Pr2 = (1-Pol) / (1+ (Kf/Ca)**nf) + Pol
       
        #Krec
        Krec = ko + (kmax - ko) / (1 + (kr / Ca) ** nr)
       
        #dCa
        dCa = -Ca/tauCa*(1/(1+kCa/Ca)) # more Ca larger the slope down (decay rate)
        if t % tAP == 0 :
            dCa +=1
            
        input_strength = (Pr1 * X1 + Pr2 * X2) #weighted average of probability
        
        if t % tAP == 0:
            PPR = input_strength*R / (Poh * Ph + Pol * (1 - Ph))   # strength of activation * proportion of active receptors left
            function_array_PPR.append(PPR)
       
       
        X1,Y1,Z1 = calculate_XYZ(X1, Y1, Z1, Pr1, Krec, t, tAP, tin)
        X2,Y2,Z2 = calculate_XYZ(X2, Y2, Z2, Pr2, Krec, t ,tAP, tin)
        
        #COMMENT THIS OUT FOR NO-DESENSITIZATION
        # R = calculate_R(R, t, tAP, input_strength, kenter, trecovery) 
        
        Ca += dCa
        
    #calculate sum of squared erros
        
    error_sum = 0
    
    for pulse in range (0,10):    
        
        difference_squared = (function_array_PPR[pulse] - experimental_data[pulse][1]) ** 2
        error_sum += difference_squared


    # print (error_sum)    
    return error_sum, function_array_PPR
 
def cumulative_error(x_vector):
    global WT_array_PPR_5
    global WT_array_PPR_10
    global WT_array_PPR_20
    global WT_array_PPR_50
    
    global KO_array_PPR_5
    global KO_array_PPR_10
    global KO_array_PPR_20
    global KO_array_PPR_50

    WT_error_sum_5, WT_array_PPR_5 = run_simulation_reutrn_error(False, x_vector, WT_5hz, 5)
    WT_error_sum_10, WT_array_PPR_10 = run_simulation_reutrn_error(False, x_vector, WT_10hz, 10)
    WT_error_sum_20, WT_array_PPR_20 = run_simulation_reutrn_error(False, x_vector, WT_20hz, 20)
    WT_error_sum_50, WT_array_PPR_50 = run_simulation_reutrn_error(False, x_vector, WT_50hz, 50)
    
    KO_error_sum_5, KO_array_PPR_5 = run_simulation_reutrn_error(True, x_vector, KO_5hz, 5)
    KO_error_sum_10, KO_array_PPR_10 = run_simulation_reutrn_error(True, x_vector, KO_10hz, 10)
    KO_error_sum_20, KO_array_PPR_20 = run_simulation_reutrn_error(True, x_vector, KO_20hz, 20)
    KO_error_sum_50, KO_array_PPR_50 = run_simulation_reutrn_error(True, x_vector, KO_50hz, 50)
    
    cumulative_error_sum =  WT_error_sum_5 + WT_error_sum_10 + WT_error_sum_20  + WT_error_sum_50 \
                            + KO_error_sum_5 +KO_error_sum_10 + KO_error_sum_20 + KO_error_sum_50
    
    print (cumulative_error_sum)
    return cumulative_error_sum



x_vector_init = [tauCa_init, kCa_init, ko_init, tin_init, Ph_init, kenter_init, trecovery_init, #index 0-6
                 KO_kmax_init , KO_kr_init, KO_nr_init, KO_Kf_init, KO_nf_init, #index 7-11
                 WT_kmax_init , WT_kr_init, WT_nr_init, WT_Kf_init, WT_nf_init,   #index 12-16       
                ]

#True to run optimizer, False to use Table 1 parameters (MAKE SURE YOU COMMENT OUT DESENSITIZATION FIRST)
if True: 

    solution = minimize(cumulative_error, 
                        x_vector_init,
                        method = 'Nelder-Mead', 
                        bounds = [(0,None), (0,None), (0,None),(0,None),(0,None),(0,None), (0,None), (0,None),(0,None),(0,None),
                                  (0,None), (0,None)],
                        tol = 0.0001)
    
    print(solution.x)
    print ("COMMON")
    print ("tauCa :" + str(solution.x[0]))
    print ("kCa :" + str(solution.x[1]))
    print ("ko :" + str(solution.x[2]))
    print ("tin :" + str(solution.x[3]))
    print ("Ph :" + str(solution.x[4]))
    print ("k enter :" + str(solution.x[5]))
    print ("t rec:" + str(solution.x[6]))
    
    print ("KO")
    print ("kmax :" + str(solution.x[7]))
    print ("kr :" + str(solution.x[8]))
    print ("nr :" + str(solution.x[9]))
    print ("kf :" + str(solution.x[10]))
    print ("nf :" + str(solution.x[11]))
    
    print("WT")
    print ("kmax :" + str(solution.x[12]))
    print ("kr :" + str(solution.x[13]))
    print ("nr :" + str(solution.x[14]))
    print ("kf :" + str(solution.x[15]))
    print ("nf :" + str(solution.x[16]))


else: 
     cumulative_error(x_vector_init)   #unitless, just looking to make small, relative measurement
     


def frequency_plot (axs, WT_exp, KO_exp, WT_calc, KO_calc, frequency, ymin,ymax):

    axs.plot([item[1] for item in WT_exp],  'bx' , label = " WT experimental")
    axs.plot([item[1] for item in KO_exp], 'bx' , label = " KO experimental")
    axs.plot(WT_calc, 'r--' )
    axs.plot(KO_calc, 'r--')
    plt.ylim((ymin,ymax))
    axs.set_title (frequency) 

fig, axs = plt.subplots(2, 2)
       
frequency_plot(axs[0, 0], WT_5hz, KO_5hz, WT_array_PPR_5, KO_array_PPR_5, "5 Hz", 0.5, 1.5)
frequency_plot(axs[0, 1], WT_10hz, KO_10hz, WT_array_PPR_10, KO_array_PPR_10, "10 Hz", 0.5, 2)
frequency_plot(axs[1, 0], WT_20hz, KO_20hz, WT_array_PPR_20, KO_array_PPR_20, "20 Hz", 0.5, 2.5)
frequency_plot(axs[1, 1], WT_50hz, KO_50hz, WT_array_PPR_50, KO_array_PPR_50, "50 Hz", 0.5, 4)


# plt.legend()
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0) 
plt.savefig('response ratio.png', dpi=300)  
plt.show()


