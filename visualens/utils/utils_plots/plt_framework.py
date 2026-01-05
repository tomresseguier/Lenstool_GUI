#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 17:09:43 2023

@author: Tom
"""


from matplotlib import pyplot as plt

def plt_framework(full_tick_framework=False, ticks='in', image=False, width=None, figsize=10, drawscaler=1., tickscaler=0.8, minor_ticks=True) :
    if width=='full' :
        figsize = 10
    if width=='half' :
        figsize = 6
    if width=='quarter' :
        figsize = 4
    fontsize = 9
    
    figsize = figsize/drawscaler
    fontsize = fontsize/drawscaler
    
    font = {'size':fontsize, 'family':'serif'} #font={'size':18, 'family':'serif'}
    plt.rc('font', **font)
    plt.rcParams['ytick.labelsize'] = 'large'
    plt.rcParams['xtick.labelsize'] = 'large'
    plt.rcParams['axes.labelsize'] = fontsize * 4/3 #24
    
    plt.rcParams['figure.figsize'] = (figsize,figsize) #(16,16) or (12,9)
    plt.rcParams['figure.dpi'] = 950 / figsize
    #plt.rcParams['figure.figsize']=(12,12)
    #plt.rcParams['figure.figsize']=(20,20)
    #plt.rcParams['figure.dpi']=30
    if minor_ticks :
        plt.rcParams['xtick.minor.visible']='True'
        plt.rcParams['ytick.minor.visible']='True'
    else :
        plt.rcParams['xtick.minor.visible']='False'
        plt.rcParams['ytick.minor.visible']='False'
    
    if full_tick_framework :
        plt.rcParams['xtick.top']=True
        plt.rcParams['ytick.right']=True
        
        if ticks=='in' :
            plt.rcParams['xtick.direction']='in'
            plt.rcParams['ytick.direction']='in'
        if ticks=='out' :
            plt.rcParams['xtick.direction']='out'
            plt.rcParams['ytick.direction']='out'
        
        plt.rcParams['xtick.major.top']=True
        plt.rcParams['xtick.minor.top']=True
        plt.rcParams['ytick.major.right']=True
        plt.rcParams['ytick.minor.right']=True
        
        plt.rcParams['xtick.major.size']=10 *tickscaler/drawscaler
        plt.rcParams['xtick.minor.size']=6 *tickscaler/drawscaler
        plt.rcParams['ytick.major.size']=10 *tickscaler/drawscaler
        plt.rcParams['ytick.minor.size']=6 *tickscaler/drawscaler
        
        plt.rcParams['xtick.major.width']=1.75 *tickscaler/drawscaler
        plt.rcParams['xtick.minor.width']=1.2 *tickscaler/drawscaler
        plt.rcParams['ytick.major.width']=1.75 *tickscaler/drawscaler
        plt.rcParams['ytick.minor.width']=1.2 *tickscaler/drawscaler
    
    plt.rcParams['savefig.bbox'] = 'tight'
    plt.rcParams['savefig.format'] = 'pdf'
    plt.rcParams['figure.autolayout'] = True
    plt.rcParams['savefig.pad_inches'] = 0
    
    if image :
        plt.rcParams['savefig.dpi']=150
    else :
        plt.rcParams['savefig.dpi']=300
    
    #plt.rcParams['text.usetex']=True
