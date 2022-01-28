import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import os
import glob


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
    NRMSEMackey = np.array([])
    STDMackey = np.array([])
    NRMSELorenz = np.array([])
    STDLorenz = np.array([])
    dataLorenz = os.path.join("/home/noah.jaitner/Local/Masterarbeit/bin/Release/", "lorenz.dat_MeanWithSTD_M2_N2_Nv40_g0.012_J00.80_Pre1_D0.0.txt")
    dataMackey = os.path.join("/home/noah.jaitner/Local/Masterarbeit/bin/Release/","lorenz.dat_MeanWithSTD_M2_N2_Nv40_g0.012_J00.80_Pre1_D0.0.txt")
    for line in open(dataMackey, 'r'):
        values = [float(s) for s in line.split()]
        NRMSEMackey = np.append(NRMSEMackey,values[0])
        STDMackey = np.append(STDMackey,values[1])
    for line in open(dataLorenz, 'r'):
        values = [float(s) for s in line.split()]
        NRMSELorenz = np.append(NRMSELorenz,values[0])
        STDLorenz = np.append(STDLorenz,values[1])
    K = np.linspace(0.01,1.0,100)
    plt.errorbar(K[:20],NRMSEMackey[:20],STDMackey[:20]/2,label = "Mackey")
    plt.errorbar(K[:20],NRMSELorenz[:20],STDLorenz[:20]/2, label = "Lorenz")
    #plt.errorbar(K,NRMSE_training,STD_training, label = "Training")
    plt.legend()
    plt.xlabel("K")
    plt.ylabel("NRMSE")
    plt.savefig("5step.png")
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
    
def colorPlotDelayNv(system = "m"):
    if system == "m":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/mackey/mackey_MeanWithSTD_M2_N2_Nv'
    if system == "l":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/lorenz/lorenz_MeanWithSTD_M2_N2_Nv'
    if system == "n":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/narma/narma_MeanWithSTD_M2_N2_Nv'
    j = 0
    NRMSE_all = np.empty(shape=(44,20))
    for i in range(6,50):
        NRMSE = np.array([])
        if system == "m":
            complete_FilePath = FilePath +  str(i) + '_g1_J00_Pre1_D0.05_Delay21_K0.05_Rtik0.000005.txt'
        if system == "l":
            complete_FilePath = FilePath +  str(i) + '_K0.05_Rtik0.000005.txt'
        if system == "n":
            complete_FilePath = FilePath + str(i) + '_g1.8_J00.4_Pre0_D0.05_Delay210_K0.02_Rtik0.000005.txt'
        for line in open(complete_FilePath, 'r'):
            values = [float(s) for s in line.split()]
            NRMSE = np.append(NRMSE,values[0])
        NRMSE[NRMSE > 1] = 1
        NRMSE_all[j] = NRMSE[:20]
        j = j+1
    Nv = np.linspace(2,44,44)
    D = np.linspace(0,20,20)
    fig, ax = plt.subplots()
    c = ax.pcolormesh(D, Nv, NRMSE_all)
    ax.set_title(system)
    fig.colorbar(c, ax=ax)
    plt.savefig(system + ".png")
    plt.show()
    
    
def ParameterOptimization(prediction, system = "m"):
    if system == "m":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/mackey/mackey.dat_MeanWithSTD_M2_N2_Nv20_g-0.01_J0'
    if system == "l":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/lorenz/lorenz.dat_MeanWithSTD_M2_N2_Nv20_g-0.001_J0'
    if system == "n":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/narma/narma.dat_MeanWithSTD_M2_N2_Nv20_g-1_J0'
    j = 0
    NRMSE_all = np.empty(shape=(20,100))
    for i in range(0,20):
        NRMSE = np.array([])
        if system == "m":
            complete_FilePath = FilePath + str(i*0.1) + '_Pre'+ prediction + '_Delay0_K0.05.txt'
        if system == "l":
            complete_FilePath = FilePath + str(i*0.1) + '_Pre'+ prediction + '_Delay0_K0.05.txt'
        if system == "n":
            complete_FilePath = FilePath + str(i*0.1) + '_Pre'+ prediction + '_Delay0_K0.05.txt'
        try:
            for line in open(complete_FilePath, 'r'):
                values = [float(s) for s in line.split()]
                NRMSE = np.append(NRMSE,values[0])
        except FileNotFoundError:
            NRMSE = np.append(NRMSE,0.95*np.ones(100))
        NRMSE[NRMSE > 1] = 1
        NRMSE_all[j] = NRMSE
        j = j+1
    J0 = np.linspace(0,2,20)
    g = np.linspace(0,1,100)
    fig, ax = plt.subplots()
    c = ax.pcolormesh(g, J0, NRMSE_all)
    ax.set_title(system + prediction)
    fig.colorbar(c, ax=ax)
    plt.show()
  
def colorPlotDelay(prediction, storung, system = "m"):
    if system == "m":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/mackey/mackey_MeanWithSTD_M2_N2_Nv30_g1_J00_Pre'
    if system == "l":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/lorenz/lorenz_MeanWithSTD_M2_N2_Nv30_g0.001_J02_Pre'
    if system == "n":
        FilePath = '/home/noah.jaitner/Local/Masterarbeit/bin/Release/narma/narma_MeanWithSTD_M2_N2_Nv30_g1.8_J00.4_Pre' 
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
    Nv = np.linspace(0,20,20)
    D = np.linspace(0,20,20)
    fig, ax = plt.subplots()
    c = ax.pcolormesh(D, Nv, NRMSE_all)
    ax.set_title(system + prediction)
    fig.colorbar(c, ax=ax)
    plt.savefig(system + prediction + ".png")
    plt.show()
    
def StatematrixCalc():
    complete_FilePath_1 = "/home/noah.jaitner/Local/Masterarbeit/bin/Release/Statematrix.txt"
    complete_FilePath = "/home/noah.jaitner/Local/Masterarbeit/bin/Release/Statematrix_0.txt"
    with open(complete_FilePath) as textFile:
        D_0 = [line.split() for line in textFile]
    with open(complete_FilePath_1) as textFile:
        D_1 = [line.split() for line in textFile]    
        
    D_0 = np.array(D_0).astype(np.float)
    D_1 = np.array(D_1).astype(np.float)
    
    Differenz = abs(D_0-D_1)/D_0

    Differenz = np.nan_to_num(Differenz)
    
    summe = 0
    
    for i in range(len(Differenz)):
        for j in range(len(Differenz[0])):
            summe = summe + Differenz[i][j]
    
    Mittelwert = summe/(len(Differenz)*len(Differenz[0]))
    
    print(Mittelwert)
    
if __name__ == "__main__":
    #plotSolution()
    #plotStoerung()
    #plotNRMSEwithSTD()
    #plotLorenz()
    #plotNarma()
    #plotMackey()
    #colorPlotDelayNv('m')
    #colorPlotDelayNv('n')
    colorPlotDelay("1","0.5","m")
    colorPlotDelay("1","0.5","l")
    colorPlotDelay("0","0.5","n")
    colorPlotDelay("5","0.5","m")
    colorPlotDelay("5","0.5","l")
    #StatematrixCalc()
    #ParameterOptimization("1","m")
    #ParameterOptimization("1","l")
    #ParameterOptimization("5","m")
    #ParameterOptimization("5","l")
    #ParameterOptimization("1","n")
    #ParameterOptimization("5","n")
