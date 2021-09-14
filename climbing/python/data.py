#!/usr/bin/env python
import numpy as np
from matplotlib import pyplot as plt
import numpy as np

def RopeData():
    data = np.genfromtxt("../data/RopeData1.txt")
    return (data.T)

if __name__ == '__main__':
    fig=plt.figure(figsize=(7,21))
    data = RopeData()
    print(data)

    plt.plot(data[0],data[1],'r')
    plt.ylabel('F$_{y}$ [N]')
    plt.xlabel('Time [s]')
    plt.show()
