#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  6 20:44:08 2020

@author: glisantplasa
"""

#CODE USED TO GET VALUES FROM THE SIMPLIFIED MODEL OF X IN Syt7 Mathematical Model:
#CODE FROM ELIANA M. & GLISANT P.

import matplotlib.pyplot as plt


def get_other_x(Fr):
    Ca = 0.01 #first impulse has already occured at t=0
    tauCa = 30.1 #milliseconds

    kCa= 1.19

    
    tAP= 1000/Fr
    t_sim = int(tAP*10) #show course of time over 10 pulses, match scale
    array_Ca =[]

    Poh=0.55 # high release probability 
    Pol=0.03 # low release probability 
    Ph=0.17 #fraction high release sites
    Pl=0.83 #fraction low release sites

    Po = Poh*Ph + Pol*Pl #weighted average of release probability naught
    
    #For WT Syt 7
    Kf = 2.4 #WT Syt 7
    nf = 1.15
    array_Pr =[]

    ko=0.0009 #/ms
    kmax = 0.026 #/ms  (WT syt 7)
    kr= 4.05 # WT
    nr= 1 #WT
    array_Krec = []
    
    """For Syt7 KO
    Kf = 19.4 
    nf = 1.5
    array_Pr =[]

    ko=0.0009 #/ms
    kmax = 0.00825 #/ms 
    kr= 0.65 # WT
    nr= 3.92 #WT
    array_Krec = []
    """
    X=1.0 # for PPR response, eq 10
    Y=0
    Z=0

    other_X=[]
    
    tin = 3.0 #ms

    PPR = 1
    array_PPR_easy= []

    for t in range (0,t_sim) :
        
        #Pr
        if t % tAP == 0 :
            Pr = (1-Po) / (1+ (Kf/Ca)**nf) + Po
            array_Pr.append(Pr)
        
        
        #Krec
        Krec = ko + (kmax -ko)/(1+(kr/Ca)**nr)
        
        
        array_Ca.append(Ca)
        
        array_Krec.append(Krec)
        
        
        #dCa
        dCa = -Ca/tauCa*(1/(1+kCa/Ca)) 
        if t % tAP == 0 :
            dCa +=1
        
        if t % tAP == 0:
            PPR = (Pr*X)/Po 
            
            array_PPR_easy.append(PPR) 
        
        
        if t%tAP==0:
            X=(X*(1-Pr))
           
            other_X.append(X)
        else:
            dX= ((1-X))*Krec/(((tin*(0.3473))))
            X+=(dX)
            other_X.append(X)
        
        Ca += dCa
       
        
   
    
    
    #plt.plot(other_X, label = "My X")
    plt.plot(array_PPR_easy, label = "PPR")
    plt.legend()
    plt.show()

   
    return array_PPR_easy
    #return other_X
