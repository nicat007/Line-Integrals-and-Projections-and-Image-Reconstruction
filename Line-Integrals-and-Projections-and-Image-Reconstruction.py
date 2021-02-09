'''
Created on Nov 28, 2019

@author: Nicat
'''

import numpy as np
import scipy.io as sio
#import sympy as sp
from skimage.transform import radon, rescale
from pylab import *
from matplotlib import pyplot as plt


def lineIntegral(table, M, imgMtrx):    ## line integral for forward projection
    length = len(table)
    p = 0
    for l in range(length-1):
        dist = math.sqrt((table[l+1][0]-table[l][0])**2 + (table[l+1][1]-table[l][1])**2)
        midx = (table[l][0]+table[l+1][0])/2
        midy = (table[l][1]+table[l+1][1])/2
        rowdata = int(M/2 - math.floor(midy)-1)
        columndata = int(M/2 + math.ceil(midx)-1)

        p = p+imgMtrx[rowdata][columndata] * dist  # calculated projection

    return p



if __name__ == '__main__':
    img = sio.loadmat('/home/nicat/Desktop/Nicat/Biomedcal Area/EE 415/filt-back-proj-master/square.mat')   ## load a matrix
    img = img['square']
    #img = [[1, 2],[3, 4]]
    #img = [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 15, 16]]
    # img = [[1, 2], [3, 4], [5, 6],[7,8],[9,10],[11,12],[13,14]]
    # img = rescale(img, scale=0.4, mode='reflect', multichannel=False)



    ## getting the necessary constant and arrays

    size = np.shape(img)

    M = size[0] # size of the image
    n_beams = 100  # number of beams
    step_size = np.pi / 180
    diag = M * np.sqrt(2)  # diagonal size of .mat matrix
    interval = diag / (n_beams - 1)

    t_range = np.arange(-diag / 2, diag / 2 + interval, interval)

    theta_range = np.arange(0, np.pi , step_size)

    x_range = np.arange(-M / 2, M / 2 + 1, 1)
    y_range = np.arange(-M / 2, M / 2 + 1, 1)
    y_table = [[], []]
    x_table = [[], []]

    proj = np.zeros([len(t_range), len(theta_range)])
    rows = np.zeros([len(t_range), len(theta_range)])
    columns = np.zeros([len(t_range), len(theta_range)])
    distance = np.zeros([len(x_range), 1])

    ### for nested for loops to calculate x-y table for forward projection:
    theta_j = 0
    for theta in theta_range:
        t_i = 0
        for t in t_range:
            if (np.sin(theta)==1):
                y_table = [[x, t] for x in x_range if abs(t) <= M / 2]  # calc. & saving in interrsctn points wrt x val.
            elif (np.sin(theta)!=1 and np.sin(theta) > 1.23e-15):  # /0 gives infinity which is not required
                y_table = [[x, float((t - x * np.cos(theta)) / np.sin(theta))] for x in x_range if abs((t - x * np.cos(theta)) / np.sin(theta)) <= M / 2]  # calculating and saving in intersection points accoording x values

            if (np.cos(theta)==1):
                x_table = [[t, y] for y in y_range if abs(t) <= M / 2]  # calc. & saving in intersec. points wrt x val.

            elif (np.cos(theta) != 1 and abs(np.cos(theta))>1.23e-15):  # /0 gives infinity which is not required
                x_table = [[float((t - y * np.sin(theta)) / np.cos(theta)), y] for y in y_range if abs((t - y * np.sin(theta)) / np.cos(theta)) <= M / 2]  # calculating and saving in intersection points according y values


            ###  concatenating tables
            xy_table = [[], []]
            try:
                xy_table = np.concatenate((x_table, y_table))
            except:
                if x_table == [[],[]] and y_table != [[],[]]:
                    xy_table = y_table
                elif y_table == [[],[]] and x_table != [[],[]]:
                    xy_table = x_table

            if(theta_j == 21):
                s=3

            ## removing the repeated rows
            if np.count_nonzero(xy_table) != 0:

                xy_table = np.unique(xy_table, axis=0)

                proj[t_i][theta_j] = lineIntegral(xy_table, M, img)  # calculating the projection for each theta and t values

            t_i = t_i + 1

        theta_j = theta_j + 1

    plt.title("Forward projection values")
    plt.imshow(proj)

    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(8,4.5))


    ax1.set_title('Original image')
    ax1.imshow(img, cmap=plt.cm.Greys_r)







    #                                                                     ## backprojection

    # filtering:

    # bp_filtered = [[]]
    # fft_p = np.transpose(np.fft.fft2(np.transpose(proj)))
    #
    # for tt in np.arange(len(theta_range)-1):
    #     proj[:][tt] = np.multiply((sio.signal.blackmanharris(n_beams, sym=True)), fft_p[:][tt])


    ## getting the reconstructed
    newImage=np.zeros([M,M])


    theta_j = 0
    for theta in theta_range:
        t_i = 0
        for t in t_range:
            if (np.sin(theta)==1):  # /0 gives infinity which is not required
                y_table = [[x, t] for x in x_range if abs(t) <= M / 2]  # calc. & saving in interrsctn points wrt x val.
            elif (np.sin(theta)!=1 and np.sin(theta) > 1.23e-15):  # /0 gives infinity which is not required
                y_table = [[x, float((t - x * np.cos(theta)) / np.sin(theta))] for x in x_range if abs((t - x * np.cos(theta)) / np.sin(theta)) <= M / 2]  # calculating and saving in intersection points accoording x values

            if (np.cos(theta)==1):  # /0 gives infinity which is not required
                x_table = [[t, y] for y in y_range if abs(t) <= M / 2]  # calc. & saving in intersec. points wrt x val.

            elif (np.cos(theta) != 1 and abs(np.cos(theta))> 1.23e-15):  # /0 gives infinity which is not required
                x_table = [[float((t - y * np.sin(theta)) / np.cos(theta)), y] for y in y_range if abs((t - y * np.sin(theta)) / np.cos(theta)) <= M / 2]  # calculating and saving in intersection points according y values

            # if len(x_table) != 0:
            #     x_table = x_table[~np.all(x_table == 0, axis=1)]
            # if len(y_table) != 0:
            #
            #     y_table = y_table[~np.all(y_table == 0, axis=1)]

            xy_table = [[], []]
            try:
                xy_table = np.concatenate((x_table, y_table))
            except:
                if np.count_nonzero(x_table) != 0:              # and len(xy_table) != 0:
                # if len(x_table) != 0 and np.count_nonzero(x_table) != 0:  # and len(xy_table) != 0:
                    xy_table = x_table
                elif np.count_nonzero(y_table) != 0:             # and len(xy_table) != 0:
                # elif len(y_table) != 0 and np.count_nonzero(y_table) != 0:  # and len(xy_table) != 0:
                    xy_table = y_table


            if np.count_nonzero(xy_table) != 0:
            # if len(xy_table) != 0 and np.count_nonzero(xy_table) != 0:

                xy_table = np.unique(xy_table, axis=0)

                ##backimg[][] = lineIntegral2(xy_table, M,proj, t_i, theta_j )



            length = len(xy_table)
            for l in range(length - 1):
                dist = math.sqrt((xy_table[l + 1][0] - xy_table[l][0]) ** 2 + (xy_table[l + 1][1] - xy_table[l][1]) ** 2)
                midx = (xy_table[l][0] + xy_table[l + 1][0]) / 2
                midy = (xy_table[l][1] + xy_table[l + 1][1]) / 2
                rowdata = int(M / 2 - math.floor(midy) - 1)
                columndata = int(M / 2 + math.ceil(midx) - 1)


                newImage[rowdata][columndata] = newImage[rowdata][columndata] + proj[t_i][theta_j] * dist  # calculated projection

            t_i = t_i + 1

        theta_j = theta_j + 1

    ax2.set_title("Reconstructed image")
    ax2.imshow(newImage, cmap=plt.cm.Greys_r)
    fig.tight_layout()
    plt.show()
