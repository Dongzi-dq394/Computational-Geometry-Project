import numpy as np
import matplotlib.pyplot as plt

import sys


# This function generates data points in rectangle shape
# -------------------------------------------------------------- #
# Parameters:
# number: Total number of data points we need.
# x_left/x_right: The left/right boundary of these data points.
# y_bottom/y_up: The bottom/up boundary of these data points.
# 
# Return(s):
# data_list: Array with [x_i, y_i] inside
# -------------------------------------------------------------- #
# Among these points: no three points colinear and no same x-coord
def GenerateDataRect(number, x_left, x_right, y_bottom, y_up):
    if number<=0 or (x_left>=x_right) or (y_bottom>=y_up):
        sys.exit("Invalid parameters for generating data points.")
    
    data_list = []
    x_coord = set()

    while number>0:
        number -= 1

        # map the [0, 1] to [x_left, x_right] and [y_bottom, y_up]
        x, y = np.random.rand(), np.random.rand()
        x = round(x_left + x * (x_right - x_left), 2)
        y = round(y_bottom + y * (y_up - y_bottom), 2)

        flag = True
        # To avoid two points with same x-coord
        if x in x_coord:
            flag = False
        
        # To avoid colinear case
        if flag:
            for i in range(len(data_list)-1):
                for j in range(i+1, len(data_list)):
                    x1, y1 = data_list[i]
                    x2, y2 = data_list[j]
                    s1, s2 = (y1-y2)/(x1-x2), (y-y1)/(x-x1)
                    if s1==s2:
                        flag = False
                        break
                if not flag:
                    break
        
        if not flag:
            number += 1
        else:
            x_coord.add(x)
            data_list.append([x, y])
    
    return data_list


# This function plots the data points itself.
# -------------------------------------------------------------- #
# Parameters:
# P: The 2D data points in general position.
# -------------------------------------------------------------- #
def PlotData(P):
    plt.figure(figsize=(12, 5))
    for x, y in P:
        plt.scatter(x, y, c='black', s=25)
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.show()

def PlotGroupData(sub_P):
    plt.figure(figsize=(12, 5))
    colors = ['red', 'green', 'orange', 'blue', 'indigo']
    ind = 0
    for pp in sub_P:
        for x, y in pp:
            plt.scatter(x, y, c=colors[ind%len(colors)], s=25)
        ind += 1
    plt.show()

# This function plots the data points and its Convex Hull.
# -------------------------------------------------------------- #
# Parameters:
# P: The 2D data points in general position.
# V: Array list for those points lying on the Convex Hull.
# -------------------------------------------------------------- #
def PlotDataWithHull(P, V):
    plt.figure(figsize=(12, 5))
    for x, y in P:
        plt.scatter(x, y, c='black', s=25)
    # V.append(V[0])
    for i in range(len(V)-1):
        x1, y1 = V[i]
        x2, y2 = V[i+1]
        plt.plot([x1, x2], [y1, y2], '--', c='red')
    plt.plot([V[-1][0], V[0][0]], [V[-1][1], V[0][1]], '--', c='red')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.show()


def PlotSubHull(P, sub_hull):
    plt.figure(figsize=(12, 5))
    for x, y in P:
        plt.scatter(x, y, c='black', s=25)
    colors = ['red', 'green', 'orange', 'blue', 'indigo']
    ind = 0
    for V in sub_hull:
        # V.append(V[0])
        for i in range(len(V)-1):
            x1, y1 = V[i]
            x2, y2 = V[i+1]
            plt.plot([x1, x2], [y1, y2], c=colors[ind%len(colors)])
        plt.plot([V[-1][0], V[0][0]], [V[-1][1], V[0][1]], c=colors[ind%len(colors)])
        ind += 1
    plt.show()

def PlotTest(P, sub_hull, bighull, q_set):
    plt.figure(figsize=(12, 5))
    for x, y in P:
        plt.scatter(x, y, c='black', s=25)
    colors = ['red', 'green', 'orange', 'blue', 'indigo']
    ind = 0
    for V in sub_hull:
        # V.append(V[0])
        for i in range(len(V)-1):
            x1, y1 = V[i]
            x2, y2 = V[i+1]
            plt.plot([x1, x2], [y1, y2], c=colors[ind%len(colors)])
        plt.plot([V[-1][0], V[0][0]], [V[-1][1], V[0][1]], c=colors[ind%len(colors)])
        ind += 1
    # for (x, y), _, _ in q_set:
    #     plt.scatter(x, y, c='blue', s=40, marker='^')
    # for x, y in bighull:
    #     plt.scatter(x, y, c='black', s=30)
    for i in range(1, len(bighull)-1):
        plt.plot([bighull[i][0], bighull[i+1][0]], [bighull[i][1], bighull[i+1][1]], '--', c='black')
    xs, ys = bighull[-1]
    # for (x, y), _, _ in q_set:
    #     if x==xs and y==ys:
    #         plt.plot([bighull[-2][0], x], [bighull[-2][1], y], '--', c='pink')
    #     else:
    #         plt.plot([bighull[-2][0], x], [bighull[-2][1], y], '--', c='black')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.show()