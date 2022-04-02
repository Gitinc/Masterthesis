import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import os
import glob
import pandas as pd


def plotSolution():
    soll = np.array([])
    ist = np.array([])
    x = np.array([])
    y = np.array([])
    z = np.array([])
    completeName = os.path.join("/home/noah.jaitner/Local/Masterarbeit/bin/Release/", "narma.dat")
    with open(completeName) as my_file:
        for line in my_file:
            rows = line.split(';')
            rows[0] = float(rows[0].replace('\n',''))
            #rows[1] = float(rows[1].replace('\n',''))
            #rows[2] = float(rows[2].replace('\n',''))
            x = np.append(x,rows[0])
            #y = np.append(y,rows[1])
            #z = np.append(z,rows[2])
            #ist = np.append(ist,rows[1])
    #plt.plot(soll[:200])
    #fig = plt.figure()
    #ax = fig.add_subplot(111,projection='3d')
    plt.plot(x[:1000])
    #plt.plot(y)
    #plt.plot(ist[:200])
    plt.show()

def plotNRMSEwithSTD():
    NRMSEMackey1 = np.array([])
    STDMackey1 = np.array([])
    NRMSELorenz1 = np.array([])
    STDLorenz1 = np.array([])
    NRMSELorenz5 = np.array([])
    STDLorenz5 = np.array([])
    NRMSEMackey5 = np.array([])
    STDMackey5 = np.array([])
    NRMSENarma10 = np.array([])
    STDNarma10 = np.array([])
    dataLorenz1 = os.path.join("/home/noah.jaitner/Local/Masterarbeit/bin/Release/", "lorenz_MeanWithSTD_M1_N1_Nv30_g0.001_J02_Pre1_Delay10_Delay20_K0.05_Rtik0.000005.txt")
    dataMackey1 = os.path.join("/home/noah.jaitner/Local/Masterarbeit/bin/Release/","mackey_MeanWithSTD_M2_N2_Nv30_g1_J00_Pre5_Delay10_Delay20_K0.05_Rtik0.000005.txt")
    dataLorenz5 = os.path.join("/home/noah.jaitner/Local/Masterarbeit/bin/Release/","lorenz_MeanWithSTD_M2_N2_Nv30_g0.001_J02_Pre1_Delay10_Delay20_K0.05_Rtik0.000005.txt")
    dataMackey5 = os.path.join("/home/noah.jaitner/Local/Masterarbeit/bin/Release/","mackey_MeanWithSTD_M1_N1_Nv30_g1_J00_Pre5_Delay10_Delay20_K0.05_Rtik0.000005.txt")
    dataNarma10 = os.path.join("/home/noah.jaitner/Local/Masterarbeit/bin/Release/","narma_MeanWithSTD_M2_N2_Nv30_g1.8_J00.4_Pre0_Delay10_Delay210_K0.02_Rtik0.000005.txt")
    for line in open(dataMackey1, 'r'):
        values = [float(s) for s in line.split()]
        NRMSEMackey1 = np.append(NRMSEMackey1,values[0])
        STDMackey1 = np.append(STDMackey1,values[1])
    for line in open(dataLorenz1, 'r'):
        values = [float(s) for s in line.split()]
        NRMSELorenz1 = np.append(NRMSELorenz1,values[0])
        STDLorenz1 = np.append(STDLorenz1,values[1])
    for line in open(dataMackey5, 'r'):
        values = [float(s) for s in line.split()]
        NRMSEMackey5 = np.append(NRMSEMackey5,values[0])
        STDMackey5 = np.append(STDMackey5,values[1])
    for line in open(dataLorenz5, 'r'):
        values = [float(s) for s in line.split()]
        NRMSELorenz5 = np.append(NRMSELorenz5,values[0])
        STDLorenz5 = np.append(STDLorenz5,values[1])
    for line in open(dataNarma10, 'r'):
        values = [float(s) for s in line.split()]
        NRMSENarma10 = np.append(NRMSENarma10,values[0])
        STDNarma10 = np.append(STDNarma10,values[1])
    K = np.linspace(0,0.2,100)
    plt.errorbar(K,NRMSEMackey1,STDMackey1, label = "2 Memory Cells")
    #plt.errorbar(K,NRMSELorenz1,STDLorenz1)
    plt.errorbar(K,NRMSEMackey5,STDMackey5, label="1 Memory Cell")
    #plt.errorbar(K,NRMSELorenz5,STDLorenz5)
    #plt.errorbar(K,NRMSENarma10,STDNarma10)
    plt.title("Mackey 5 Step Prediction")
    plt.legend()
    plt.xlabel("D")
    plt.ylabel("NRMSE")
    plt.savefig("nrmseMackey2MemoryCells.png")
    plt.show()

def plotLorenz():
    lorenzFiles = glob.glob('/home/noah.jaitner/Local/Masterarbeit/bin/Release/lorenz/lorenz.dat_MeanWithSTD_M2_N2_Nv30_g0.012_J00.85_Pre1_D*.txt')
    NRMSE = np.empty(shape=(len(lorenzFiles),80))
    STD = np.empty(shape=(len(lorenzFiles),80))
    gValue = [None]*len(lorenzFiles)
    NRMSE_1 = np.array([])
    STD_1 = np.array([])
    for i in range(len(lorenzFiles)):
        completeName = lorenzFiles[i]
        with open(completeName) as my_file:
            gValue[i] = completeName.split("_")[8]
            for line in my_file:
                rows = line.split(' ')
                rows[1] = float(rows[1])
                rows[0] = float(rows[0])
                NRMSE_1 = np.append(NRMSE_1[:79],rows[0])
                STD_1 = np.append(STD_1[:79],rows[1])
            NRMSE[i] = NRMSE_1
            STD[i] = STD_1
            NRMSE_1 = np.array([])
            STD_1 = np.array([])
    K = np.linspace(0.01,1.0,80)
    for i in range(len(lorenzFiles)):
        plt.errorbar(K,NRMSE[i],STD[i]/2, label= gValue[i])
    plt.legend()
    plt.xlabel("K")
    plt.ylabel("NRMSE")
    plt.savefig("lorenz.png")
    plt.show()

def plotNarma():
    lorenzFiles = glob.glob('/home/noah.jaitner/Local/Masterarbeit/bin/Release/narma/narma.dat_MeanWithSTD_M2_N2_Nv*_g1.3_J00_Pre1_D0.txt')
    NRMSE = np.empty(shape=(len(lorenzFiles),80))
    STD = np.empty(shape=(len(lorenzFiles),80))
    gValue = [None]*len(lorenzFiles)
    NRMSE_1 = np.array([])
    STD_1 = np.array([])
    for i in range(len(lorenzFiles)):
        completeName = lorenzFiles[i]
        with open(completeName) as my_file:
            gValue[i] = completeName.split("_")[8]
            for line in my_file:
                rows = line.split(' ')
                rows[1] = float(rows[1])
                rows[0] = float(rows[0])
                NRMSE_1 = np.append(NRMSE_1,rows[0])
                STD_1 = np.append(STD_1,rows[1])
            NRMSE[i] = NRMSE_1
            STD[i] = STD_1
            NRMSE_1 = np.array([])
            STD_1 = np.array([])
    K = np.linspace(0.01,1.0,100)
    for i in range(len(lorenzFiles)):
        plt.errorbar(K,NRMSE[i],STD[i]/2, label= gValue[i])
    plt.legend()
    plt.xlabel("K")
    plt.ylabel("NRMSE")
    plt.savefig("narma.png")
    plt.show()

def plotMackey():
    lorenzFiles = glob.glob('/home/noah.jaitner/Local/Masterarbeit/bin/Release/mackey/mackey.dat_MeanWithSTD_M2_N2_Nv30_g1_J00_Pre1_D*.txt')
    NRMSE = np.empty(shape=(len(lorenzFiles),80))
    STD = np.empty(shape=(len(lorenzFiles),80))
    gValue = [None]*len(lorenzFiles)
    NRMSE_1 = np.array([])
    STD_1 = np.array([])
    for i in range(len(lorenzFiles)):
        completeName = lorenzFiles[i]
        with open(completeName) as my_file:
            gValue[i] = completeName.split("_")[8]
            for line in my_file:
                rows = line.split(' ')
                rows[1] = float(rows[1])
                rows[0] = float(rows[0])
                NRMSE_1 = np.append(NRMSE_1[:79],rows[0])
                STD_1 = np.append(STD_1[:79],rows[1])
            NRMSE[i] = NRMSE_1
            STD[i] = STD_1
            NRMSE_1 = np.array([])
            STD_1 = np.array([])
    K = np.linspace(0.01,.8,80)
    for i in range(len(lorenzFiles)):
        plt.errorbar(K,NRMSE[i],STD[i]/2, label= gValue[i])
    plt.legend()
    plt.xlabel("K")
    plt.ylabel("NRMSE")
    plt.savefig("mackey.png")
    plt.show()

def colorPlotDelayNv(prediction, system = "m"):
    if system == "m":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/mackey/mackey_MeanWithSTD_M2_N2_Nv'
        title = "Mackey " + prediction + " Step Prediction"
    if system == "l":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/lorenz/lorenz_MeanWithSTD_M2_N2_Nv'
        title = "Lorenz " + prediction + " Step Prediction"
    if system == "n":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/narma/narma_MeanWithSTD_M2_N2_Nv'
        title = "Narma10"
    j = 0
    NRMSE_all = np.empty(shape=(50,21))
    for i in range(4,54):
        NRMSE = np.array([])
        if system == "m" and prediction == '1':
            complete_FilePath = FilePath +  str(i) + '_g1_J00_Pre' + prediction + '_Delay10_Delay20_K0.05_Rtik0.000005.txt'
        if system == "l" and prediction == '1':
            complete_FilePath = FilePath +  str(i) + '_g0.001_J02_Pre' + prediction + '_Delay10_Delay20_K0.05_Rtik0.000005.txt'
        if system == "m" and prediction == '5':
            complete_FilePath = FilePath +  str(i) + '_g1_J00_Pre' + prediction + '_Delay10_Delay20_K0.05_Rtik0.000005.txt'
        if system == "l" and prediction == '5':
            complete_FilePath = FilePath +  str(i) + '_g0.001_J02_Pre' + prediction + '_Delay10_Delay20_K0.05_Rtik0.000005.txt'
        if system == "n":
            complete_FilePath = FilePath +  str(i) + '_g1.8_J00.4_Pre0_Delay10_Delay20_K0.02_Rtik0.000005.txt'
        try:
            for line in open(complete_FilePath, 'r'):
                values = [float(s) for s in line.split()]
                NRMSE = np.append(NRMSE,values[0])
        except FileNotFoundError:
            NRMSE = np.append(NRMSE,0.95*np.ones(21))
        NRMSE[NRMSE > 1] = 1
        NRMSE_all[j] = NRMSE[:21]
        j = j+1
    Nv = np.linspace(4,54,50)
    Delay = np.linspace(0,20,21)
    fig, ax = plt.subplots()
    c = ax.pcolormesh(Delay,Nv,NRMSE_all)
    ax.set_title(title)
    ax.xaxis.set_ticks(np.arange(0, 22, 2))
    ax.yaxis.set_ticks(np.arange(4, 55, 4))
    plt.xlabel("Delay1")
    plt.ylabel("Nv")
    fig.colorbar(c, ax=ax)
    plt.savefig(system + prediction + "DelayNv.png")
    plt.show()

def colorPlotDelayNoise(prediction, system = "m"):
    if system == "m":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/mackey/mackey_MeanWithSTD_M2_N2_Nv'
    if system == "l":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/lorenz/lorenz_MeanWithSTD_M2_N2_Nv'
    if system == "n":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/narma/narma_MeanWithSTD_M2_N2_Nv'
    j = 0
    NRMSE_all = np.empty(shape=(50,40))
    for i in range(4,54):
        NRMSE = np.array([])
        if system == "m" and prediction == '1':
            complete_FilePath = FilePath +  str(i) + '_g1_J00_Pre' + prediction + '_Delay13_Delay21_K0.05_Rtik0.000005.txt'
            title = "Mackey 1 Step Prediction"
        if system == "l" and prediction == '1':
            complete_FilePath = FilePath +  str(i) + '_g0.001_J02_Pre' + prediction + '_Delay11_Delay23_K0.05_Rtik0.000005.txt'
            title = "Lorenz X 1 Step Prediction"
        if system == "m" and prediction == '5':
            complete_FilePath = FilePath +  str(i) + '_g1_J00_Pre' + prediction + '_Delay16_Delay214_K0.05_Rtik0.000005.txt'
            title = "Mackey 5 Step Prediction"
        if system == "l" and prediction == '5':
            complete_FilePath = FilePath +  str(i) + '_g0.001_J02_Pre' + prediction + '_Delay11_Delay23_K0.05_Rtik0.000005.txt'
            title = "Lorenz 5 Step Prediction"
        if system == "n":
            complete_FilePath = FilePath +  str(i) + '_g1.8_J00.4_Pre0_Delay11_Delay210_K0.02_Rtik0.000005.txt'
            title = "Narma10"
        try:
            for line in open(complete_FilePath, 'r'):
                values = [float(s) for s in line.split()]
                NRMSE = np.append(NRMSE,values[0])
        except FileNotFoundError:
            NRMSE = np.append(NRMSE,np.ones(40))
        NRMSE[NRMSE > 1] = 1
        NRMSE_all[j] = NRMSE[:40]
        j = j+1
    Nv = np.linspace(4,54,50)
    D = np.linspace(0,0.08,40)
    fig, ax = plt.subplots()
    c = ax.pcolormesh(D, Nv, NRMSE_all)
    ax.set_title(title)
    fig.colorbar(c, ax=ax)
    plt.savefig(system + ".png")
    plt.show()

def colorPlotDelay(prediction, storung, system = "m"):
    if system == "m":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/mackey/mackey_MeanWithSTD_M2_N2_Nv30_g1_J00_Pre'
        title = "Mackey " + prediction + " Step Prediction"
    if system == "l":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/lorenz/lorenz_MeanWithSTD_M2_N2_Nv30_g0.001_J02_Pre'
        title = "Lorenz " + prediction + " Step Prediction"
    if system == "n":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/narma/narma_MeanWithSTD_M2_N2_Nv30_g1.8_J00.4_Pre'
        title = "Narma10"
    j = 0
    NRMSE_all = np.empty(shape=(20,20))
    for i in range(20):
        NRMSE = np.array([])
        if system == "m":
            complete_FilePath = FilePath + prediction + '_D' + storung + '_Delay2'+ str(i) + '_K0.05_Rtik0.000005.txt'
        if system == "l":
            complete_FilePath = FilePath + prediction + '_D' + storung + '_Delay2'+ str(i) + '_K0.05_Rtik0.000005.txt'
        if system == "n":
            complete_FilePath = FilePath + '0' + '_D' + storung + '_Delay2'+ str(i) + '_K0.02_Rtik0.000005.txt'
        for line in open(complete_FilePath, 'r'):
            values = [float(s) for s in line.split()]
            NRMSE = np.append(NRMSE,values[0])
        NRMSE[NRMSE > 1] = 1
        NRMSE_all[j] = NRMSE[:20]
        j = j+1
    test_df = pd.DataFrame(NRMSE_all)
    fig, ax = plt.subplots()
    c = ax.pcolormesh(NRMSE_all)
    plt.xlabel("Delay1")
    plt.ylabel("Delay2")
    ax.set_title(title)

    ax.xaxis.set_ticks(np.arange(0, 21, 2))
    ax.yaxis.set_ticks(np.arange(0, 21, 2))

    fig.colorbar(c, ax=ax)
    plt.savefig(system + prediction + storung + "Delay.png")
    print(np.argwhere(NRMSE_all == np.min(NRMSE_all)))
    plt.show()

def colorPlotKgMackey(prediction):
    FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/mackey/mackey_MeanWithSTD_M2_N2_Nv30_g'
    NRMSE_all = np.empty(shape=(20,100))
    j = 0
    for i in range(1,21):
        NRMSE = np.array([])
        complete_FilePath = ''.join([FilePath,"{:.1f}".format(i/20),"_J00_Pre",prediction,'_Delay10_Delay20_K000_Rtik0.000005.txt'])
        try:
            for line in open(complete_FilePath, 'r'):
                values = [float(s) for s in line.split()]
                NRMSE = np.append(NRMSE,values[0])
        except FileNotFoundError:
            NRMSE = np.append(NRMSE,np.ones(40))
        NRMSE[NRMSE > 1] = 1
        NRMSE[NRMSE < 0] = -1
        NRMSE_all[j] = NRMSE[:1]
        j = j+1
    K = np.linspace(0,1,100)
    g = np.linspace(0.1,2.0,20)
    fig, ax = plt.subplots()
    print(np.where(NRMSE_all == np.min(NRMSE_all)))
    c = ax.pcolormesh(K,g,NRMSE_all)
    plt.xlabel("K")
    plt.ylabel("g")
    fig.colorbar(c, ax=ax)
    plt.show()

def colorPlotKgLorenz(prediction):
    FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/lorenz/lorenz_MeanWithSTD_M2_N2_Nv30_g'
    NRMSE_all = np.empty(shape=(24,100))
    j = 0
    for i in range(1,24):
        NRMSE = np.array([])
        complete_FilePath = ''.join([FilePath,"{:.3f}".format(i/200),"_J01_Pre",prediction,'_Delay10_Delay20_K000_Rtik0.000005.txt'])
        try:
            for line in open(complete_FilePath, 'r'):
                values = [float(s) for s in line.split()]
                NRMSE = np.append(NRMSE,values[0])
        except FileNotFoundError:
            NRMSE = np.append(NRMSE,np.ones(40))
        NRMSE[NRMSE > 1] = 1
        NRMSE[NRMSE < 0] = -1
        NRMSE_all[j] = NRMSE[:1]
        j = j+1
    K = np.linspace(0,1,100)
    g = np.linspace(0.02,0.132,66)
    fig, ax = plt.subplots()
    c = ax.pcolormesh(K,g,NRMSE_all)
    plt.xlabel("K")
    plt.ylabel("g")
    fig.colorbar(c, ax=ax)
    plt.show()

def colorPlotKgNarma(prediction):
    FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/narma/narma_MeanWithSTD_M2_N2_Nv30_g'
    NRMSE_all = np.empty(shape=(20,100))
    j = 0
    for i in range(1,21):
        NRMSE = np.array([])
        complete_FilePath = ''.join([FilePath,"{:.1f}".format(i/20),"_J00.4_Pre",prediction,'_Delay10_Delay210_K000_Rtik0.000005.txt'])
        try:
            for line in open(complete_FilePath, 'r'):
                values = [float(s) for s in line.split()]
                NRMSE = np.append(NRMSE,values[0])
        except FileNotFoundError:
            NRMSE = np.append(NRMSE,np.ones(40))
        NRMSE[NRMSE > 1] = 1
        NRMSE[NRMSE < 0] = -1
        NRMSE_all[j] = NRMSE[:1]
        j = j+1
    K = np.linspace(0,1,100)
    g = np.linspace(0.1,2.0,20)
    fig, ax = plt.subplots()
    c = ax.pcolormesh(K,g,NRMSE_all)
    plt.xlabel("K")
    plt.ylabel("g")
    fig.colorbar(c, ax=ax)
    plt.show()

def ParameterOptimization(prediction, system = "m"):
    gPoints = 0
    if system == "m":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/Parameter/mackey/mackey_MeanWithSTD_M4_N4_Nv30_g'
        title = "Mackey " + prediction + " Step Prediction"
        gPoints = 40
    if system == "l":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/Parameter/lorenz/lorenz_MeanWithSTD_M4_N4_Nv30_g'
        title = "Lorenz " + prediction + " Step Prediction"
        gPoints = 12
    if system == "n":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/Parameter/narma/narma_MeanWithSTD_M4_N4_Nv30_g'
        title = "Narma10"
        gPoints = 40
    j = 0
    NRMSE_all = np.empty(shape=(gPoints,100))
    line_lowest = np.empty(shape=(100))
    for i in range(gPoints):
        NRMSE = np.array([])
        if system == "m":
            complete_FilePath = FilePath + "{:.1f}".format(0.1 + i/10) + '_J00_Pre' + prediction + '_Delay10_Delay20_K0.02_Rtik0.000005.txt'
            g = np.linspace(0.1,4.0,gPoints)
            I_min = 0.4
        if system == "l":
            complete_FilePath = FilePath + "{:.3f}".format(0.001 + i/1000) + '_J01_Pre' + prediction + '_Delay10_Delay20_K0.02_Rtik0.000005.txt'
            g = np.linspace(0.001,0.013,gPoints)
            I_min = -15
        if system == "n":
            complete_FilePath = FilePath + "{:.1f}".format(0.1 + i/10) + '_J00.4_Pre' + prediction + '_Delay10_Delay210_K0.02_Rtik0.000005.txt'
            g = np.linspace(0.1,4.0,gPoints)
            I_min = 0.1
        for line in open(complete_FilePath, 'r'):
            values = [float(s) for s in line.split()]
            NRMSE = np.append(NRMSE,values[0])
        NRMSE[NRMSE > 1] = 1
        NRMSE_all[j] = NRMSE
        j = j+1
    J_0 = np.linspace(0,2.5,100)
    fig, ax = plt.subplots()
    c = ax.pcolormesh(J_0,g,NRMSE_all)
    ax.set_title(title)
    plt.xlabel("J_0")
    plt.ylabel("g")
    fig.colorbar(c, ax=ax)
    plt.show()

def ParameterOptimizationOverK(N, prediction, system = "m"):
    gPoints = 0
    if system == "m":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/Parameter/mackey/mackey_MeanWithSTD_M' + N +'_N' + N + '_Nv30_g'
        title = "Mackey " + prediction + " Step Prediction"
        gPoints = 40
        g = np.linspace(0.1,4.0,gPoints)
        if prediction == '5':
            J = "0.075"
        else:
            J = "0.275"
    if system == "l":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/Parameter/lorenz/lorenz_MeanWithSTD_M' + N +'_N' + N + '_Nv30_g'
        title = "Lorenz " + prediction + " Step Prediction"
        gPoints = 12
        g = np.linspace(0.001,0.013,gPoints)
    if system == "n":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/Parameter/narma/narma_MeanWithSTD_M' + N +'_N' + N + '_Nv30_g'
        title = "Narma10"
        gPoints = 40
        g = np.linspace(0.1,4.0,gPoints)
    j = 0
    NRMSE_all = np.empty(shape=(gPoints,100))
    for i in range(gPoints):
        NRMSE = np.array([])
        if system == "m":
            complete_FilePath = FilePath + "{:.1f}".format(0.1 + i/10) + '_J0' + J +'_Pre' + prediction + '_Delay10_Delay20_K000_Rtik0.000005.txt'
        if system == "l":
            complete_FilePath = FilePath + "{:.3f}".format(0.001 + i/1000) + '_J00.1_Pre' + prediction + '_Delay10_Delay20_K000_Rtik0.000005.txt'
        if system == "n":
            complete_FilePath = FilePath + "{:.1f}".format(0.1 + i/10) + '_J00.1_Pre' + prediction + '_Delay10_Delay210_K000_Rtik0.000005.txt'
        for line in open(complete_FilePath, 'r'):
            values = [float(s) for s in line.split()]
            NRMSE = np.append(NRMSE,values[0])
        NRMSE[NRMSE > 1] = 1
        NRMSE_all[j] = NRMSE
        j = j+1
    fig, ax = plt.subplots()
    K = np.linspace(0,1,100)
    c = ax.pcolormesh(K,g,NRMSE_all)
    ax.set_title(title)
    plt.xlabel("K")
    plt.ylabel("g")
    fig.colorbar(c, ax=ax)
    plt.show()

def ParameterOptimizationOverK2d(N, prediction, g, system = "m"):
    gPoints = 0
    if system == "m":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/Parameter/mackey/mackey_MeanWithSTD_M' + N +'_N' + N + '_Nv30_g'
        title = "Mackey " + prediction + " Step Prediction"
        gPoints = 40
        if prediction == '5':
            J = "0.075"
        else:
            J = "0.275"
    if system == "l":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/Parameter/lorenz/lorenz_MeanWithSTD_M' + N +'_N' + N + '_Nv30_g'
        title = "Lorenz " + prediction + " Step Prediction"
        gPoints = 12
    if system == "n":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/Parameter/narma/narma_MeanWithSTD_M' + N +'_N' + N + '_Nv30_g'
        title = "Narma10"
        gPoints = 40
    j = 0
    NRMSE = np.array([])
    STD = np.array([])
    if system == "m":
        complete_FilePath = FilePath + g + '_J0' + J +'_Pre' + prediction + '_Delay10_Delay20_K000_Rtik0.000005.txt'
    if system == "l":
        complete_FilePath = FilePath + g + '_J00.1_Pre' + prediction + '_Delay10_Delay20_K000_Rtik0.000005.txt'
    if system == "n":
        complete_FilePath = FilePath + g + '_J00.1_Pre' + prediction + '_Delay10_Delay210_K000_Rtik0.000005.txt'
    for line in open(complete_FilePath, 'r'):
        values = [float(s) for s in line.split()]
        NRMSE = np.append(NRMSE,values[0])
        STD = np.append(STD,values[1])
    NRMSE[NRMSE > 1] = 1
    x = np.linspace(0,1,100)
    plt.errorbar(x,NRMSE,STD/2)
    plt.title(title)
    plt.xlabel("K")
    plt.ylabel("NRMSE")
    plt.show()

def Noise2d(prediction, dnla, system = "m"):
    J = '0.075'
    if prediction == '1':
        J = '0.275'
    if system == "m":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/other/mackey/mackey_MeanWithSTD_M31_N31_Nv30_g3.5_J0'+ J +'_Pre'
        title = "Mackey " + prediction + " Step Prediction"
    if system == "l":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/other/lorenz/lorenz_MeanWithSTD_M31_N31_Nv30_g0.01_J00.1_Pre'
        title = "Lorenz " + prediction + " Step Prediction"
    if system == "n":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/other/narma/narma_MeanWithSTD_M4_N4_Nv30_g0.5_J00.1_Pre'
        title = "Narma10"
    j = 0
    NRMSE = np.array([])
    STD = np.array([])
    if system == "m":
        complete_FilePath = FilePath + prediction + '_Delay10_Delay20_K0_Rtik0.000005_DNLA'+ dnla +'.txt'
    if system == "l":
        complete_FilePath = FilePath + prediction + '_Delay10_Delay20_K0_Rtik0.000005_DNLA'+ dnla +'.txt'
    if system == "n":
        complete_FilePath = FilePath + '0_Delay10_Delay210_K0_Rtik0.000005_DNLA'+ dnla +'.txt'
    for line in open(complete_FilePath, 'r'):
        values = [float(s) for s in line.split()]
        NRMSE = np.append(NRMSE,values[0])
        STD = np.append(STD,values[1])
    NRMSE[NRMSE > 1] = 1
    x = np.linspace(0,0.2,100)
    plt.errorbar(x,NRMSE,STD/2)
    plt.title(title)
    plt.xlabel("K")
    plt.ylabel("NRMSE")
    plt.show()

def NoiseNLA(prediction, system = "m"):
    if system == "m":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/other/mackey/mackey_MeanWithSTD_M4_N4_Nv30_g3.5_J00.075_Pre'
        title = "Mackey " + prediction + " Step Prediction"
    if system == "l":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/other/lorenz/lorenz_MeanWithSTD_M4_N4_Nv30_g0.01_J00.1_Pre'
        title = "Lorenz " + prediction + " Step Prediction"
    if system == "n":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/other/narma/narma_MeanWithSTD_M4_N4_Nv30_g0.5_J00.1_Pre'
        title = "Narma10"
    j = 0
    NRMSE_all = np.empty(shape=(100,100))
    for i in range(100):
        NRMSE = np.array([])
        if system == "m":
            complete_FilePath = FilePath + prediction + '_Delay10_Delay20_K0.05_Rtik0.000005_DNLA'+ "{:.3f}".format(i*0.002)+'.txt'
        if system == "l":
            complete_FilePath = FilePath + prediction + '_Delay10_Delay20_K0.06_Rtik0.000005_DNLA'+ "{:.3f}".format(i*0.002)+'.txt'
        if system == "n":
            complete_FilePath = FilePath + '0_Delay10_Delay210_K0.8_Rtik0.000005_DNLA'+ "{:.3f}".format(i*0.002)+'.txt'
        for line in open(complete_FilePath, 'r'):
            values = [float(s) for s in line.split()]
            NRMSE = np.append(NRMSE,values[0])
        NRMSE[NRMSE > 1] = 1
        NRMSE_all[j] = NRMSE
        j = j+1
    fig, ax = plt.subplots()
    c = ax.pcolormesh(NRMSE_all)
    plt.xlabel("D")
    plt.ylabel("DNLA")
    ax.set_title(title)
    fig.colorbar(c, ax=ax)
    plt.show()

def N4NvDelay(prediction, system = "m"):
    if system == "m":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/DelayNv4N/mackey/mackey_MeanWithSTD_M4_N4_Nv'
        title = "Mackey " + prediction + " Step Prediction"
    if system == "l":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/DelayNv4N/lorenz/lorenz_MeanWithSTD_M4_N4_Nv'
        title = "Lorenz " + prediction + " Step Prediction"
    if system == "n":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/DelayNv4N/narma/narma_MeanWithSTD_M4_N4_Nv'
        title = "Narma10"
    j = 0
    NRMSE_all = np.empty(shape=(39,20))
    for i in range(2,41):
        NRMSE = np.array([])
        if system == "m":
            complete_FilePath = FilePath + str(i) + '_g3.5_J00.275_Pre' + prediction + '_Delay10_Delay20_K0_Rtik0.000005_DNLA0.txt'
        if system == "l":
            complete_FilePath = FilePath + str(i) + '_g0.01_J00.1_Pre' + prediction + '_Delay10_Delay20_K0_Rtik0.000005_DNLA0.txt'
        if system == "n":
            complete_FilePath = FilePath + str(i) + '_g0.5_J00.1_Pre0_Delay10_Delay20_K0_Rtik0.000005_DNLA0.txt'
        try:
            for line in open(complete_FilePath, 'r'):
                values = [float(s) for s in line.split()]
                NRMSE = np.append(NRMSE,values[0])
        except FileNotFoundError:
            NRMSE = np.append(NRMSE,np.ones(20))
        NRMSE[NRMSE > 1] = 1
        NRMSE_all[j] = NRMSE[:20]
        j = j+1
    fig, ax = plt.subplots()
    Nv = np.linspace(2,41,39)
    Delay = np.linspace(0,21,20)
    c = ax.pcolormesh(Delay, Nv, NRMSE_all)
    plt.xlabel("Nv")
    plt.ylabel("Delay")
    ax.set_title(title)
    fig.colorbar(c, ax=ax)
    plt.show()

if __name__ == "__main__":
    #plotSolution()
    #plotStoerung()
    #colorPlotKgNarma("0")
    #colorPlotKgMackey("1")
    #colorPlotKgMackey("5")
    #colorPlotKgLorenz("1")
    #colorPlotKgLorenz("5")
    #plotNRMSEwithSTD()
    #plotLorenz()
    #plotNarma()
    #plotMackey()
    #colorPlotDelayNv("1","m")
    #colorPlotDelayNv("0","n")
    #colorPlotDelayNv("5","m")
    #colorPlotDelayNv("1","l")
    #colorPlotDelayNv("5","l")
    #colorPlotDelayNoise("1","m")
    #colorPlotDelayNoise("1","l")
    #colorPlotDelayNoise("0","n")
    #colorPlotDelayNoise("5","m")
    #colorPlotDelayNoise("5","l")
    #colorPlotDelay("1","0.1","m")
    #colorPlotDelay("1","0.1","l")
    #colorPlotDelay("0","0.1","n")
    #colorPlotDelay("5","0.1","m")
    #colorPlotDelay("5","0.1","l")
    #StatematrixCalc()
    #ParameterOptimization("1","m")
    #ParameterOptimization("5","m")
    #ParameterOptimization("1","l")
    #ParameterOptimization("5","l")
    #ParameterOptimization("0","n")
    #ParameterOptimizationOverK("4","1","m")
    #ParameterOptimizationOverK("4","5","m")
    #ParameterOptimizationOverK("4","1","l")
    #ParameterOptimizationOverK("4","5","l")
    #ParameterOptimizationOverK("4","0","n")
    #ParameterOptimizationOverK2d("4","5","0.5","m")
    #ParameterOptimizationOverK2d("4","5","0.012","l")
    #ParameterOptimizationOverK2d("4","1","0.5","m")
    #ParameterOptimizationOverK2d("4","1","0.012","l")
    #ParameterOptimizationOverK2d("4","0","1.0","n")
    #NoiseNLA('5','m')
    #NoiseNLA('1','l')
    #NoiseNLA('5','l')
    #NoiseNLA('0','n')
    Noise2d('5','0','m')
    Noise2d('5','0','l')
    Noise2d('1','0','m')
    Noise2d('1','0','l')
    Noise2d('0','0','n')
    #N4NvDelay('1','m')
    #N4NvDelay('1','l')
    #N4NvDelay('5','l')
    #N4NvDelay('0','n')

