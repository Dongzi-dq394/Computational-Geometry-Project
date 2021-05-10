import numpy as np
import matplotlib.pyplot as plt

import sys
import math
from tqdm import tqdm
import time
from util import *


# These are used to make decision for the binary search on wrapping
CCW = 1 # Counter-Clock-Wise orientation
CW = -1 # Clock-Wise orientation
COL = 0 # colinear for three points (excluded in our code)


# This function computes the orientation for three given points.
# -------------------------------------------------------------- #
# Parameters:
# p1/p2/p3: The 2D data points in general position.
#
# Return(s):
# -1: p1->p2->p3 is CW; 1: p1->p2->p3 is CCW; 0: Colinear
# -------------------------------------------------------------- #
def orientation(p1, p2, p3):
    x1, y1 = p1[0], p1[1]
    x2, y2 = p2[0], p2[1]
    x3, y3 = p3[0], p3[1]
    # We compute this val by comparing two slopes.
    val = (y2 - y1) * (x3 - x2) - (x2 - x1) * (y3 - y2)
    # val = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
    if val==0:
        return COL # Colinear
    return CW if val>0 else CCW


# This function computes the tangent point in Hull drawn from p,
# which is useful on merging step for the final big Convex Hull.
# We use Binary Search and orientation for implementation.
# -------------------------------------------------------------- #
# Parameters:
# p: The 2D data point in general position.
# Hull: Array of 2D points representing the hull (hull points).
#
# Return(s):
# The index of the point to which the tangent is drawn from p.
# -------------------------------------------------------------- #
def BinarySearchTangent(Hull, p):
    # print(Hull)
    l, r = 0, len(Hull)
    l_before = orientation(p, Hull[0], Hull[r-1])
    l_after = orientation(p, Hull[0], Hull[(l+1)%len(Hull)])
    while l<r:
        m = (l+r)//2
        m_before = orientation(p, Hull[m], Hull[(m-1)%len(Hull)])
        m_after = orientation(p, Hull[m], Hull[(m+1)%len(Hull)])
        m_side = orientation(p, Hull[l], Hull[m])
        if m_before!=CW and m_after!=CW:
            return m
        elif (((m_side == CCW) and (l_after == CW or l_before == l_after)) or (m_side == CW and m_before == CW)):
            r = m
        else:
            l = m+1
        # l_before = -m_after
        l_before = orientation(p, Hull[l], Hull[r-1])
        l_after = orientation(p, Hull[l], Hull[(l+1)%len(Hull)])
    return l


# This function computes the Convex Hull in O(nlogn) time.
# -------------------------------------------------------------- #
# Parameters:
# P: The 2D data points in general position.
#
# Return(s):
# vertices: List of points lying on the Convex Hull in ccw order
#           from the point with smallest x-coord.
# -------------------------------------------------------------- #
def GrahamScan(P):
    vertices = []
    arr = [item+[i] for i, item in enumerate(P)]
    arr.sort(key=lambda x: x[0])
        
    # Maintain a stack to know which vertex should be popped out
    # We start from the leftmost point building the lower hull
    stack = [arr[0]]
    for i in range(1, len(arr)):
        while len(stack)>=2:
            x1, y1, _ = stack[-2]
            x2, y2, _ = stack[-1]
            nx, ny, _ = arr[i]
            ow, nw = (y1-y2)/(x1-x2), (ny-y2)/(nx-x2)
            if ow>nw:
                stack.pop()
            else:
                break
        stack.append(arr[i])
    # we add index of all vertex on the convex hull orderly to the vertices
    for x, y, ind in stack:
        vertices.append([x, y])
        
    # We do this again but from the rightmost side to build up the upper hull
    stack = [arr[-1]]
    for i in range(len(arr)-2, -1, -1):
        while len(stack)>=2:
            x1, y1, _ = stack[-2]
            x2, y2, _ = stack[-1]
            nx, ny, _ = arr[i]
            ow, nw = (y1-y2)/(x1-x2), (ny-y2)/(nx-x2)
            if ow>nw:
                stack.pop()
            else:
                break
        stack.append(arr[i])
    # The first one and the last one have already been added into self.vertices
    for i in range(1, len(stack)-1):
        vertices.append([stack[i][0], stack[i][1]])
    
    return vertices


# The main function for the Chan's algorithm.
# -------------------------------------------------------------- #
# Parameters:
# P: A set of n points in general position.
# m: The largest possible size of each subsets after partitioning.
# H: The iteration number finding the vertices on the Convex Hull.
#    We guess this number from 1 trickly.
# -------------------------------------------------------------- #
def Hull2D(P, m, H, test):
    if test:
        PlotData(P)

    l = len(P)

    # Partition the P and put them into sub_P
    sub_P = []
    start = 0
    while l>m:
        sub_P.append(P[start:start+m])
        l -= m
        start += m
    # If the remaining size is smaller than 3, we merge it into the last one
    # otherwise, we regard it as a single one sub set
    if l>0:
        if l>=3:
            sub_P.append(P[start:])
        else:
            sub_P[-1] += P[start:]

    if test:
        PlotGroupData(sub_P)
    
    # Apply Graham's Scan Algorithm on each sub set and put results in sub_hull
    sub_hull = []
    for partition in sub_P:
        sub_hull.append(GrahamScan(partition))
    
    if test:
        PlotSubHull(P, sub_hull)

    # We find the leftmost point among all points in P.
    # In the same time, we record the hull index of this point
    # and the point index in this hull.
    lp, temp = None, float('Inf')
    hullind = posind = None
    for i, hh in enumerate(sub_hull):
        for j, (x, y) in enumerate(hh):
            if x<temp:
                temp = x
                lp = [x, y]
                hullind = i
                posind = j

    # This is the big Convex Hull for the entire P
    # We first put two points inside, where the first one is useful for 
    # the later calculation; and the second one is the first point in the hull.
    bighull = [[lp[0], lp[1]+1], [lp[0], lp[1]]]
    for _ in range(H):
        # Select previous two points for calculation
        p_k, p_k_1 = bighull[-1], bighull[-2]

        vector = (p_k_1[0]-p_k[0], p_k_1[1]-p_k[1])
        dis = math.sqrt(vector[0]**2 + vector[1]**2)
        
        # We compute the tangent drawn from p_k in each sub hull
        # from sub_hull and put them into q_set.
        q_set = []
        for i, hh in enumerate(sub_hull):
            # If this subhull contains the previous best point,
            # then we just choose the next point as the tangent.
            if hullind == i:
                q_set.append((hh[(posind+1)%len(hh)], i, (posind+1)%len(hh)))
            # Otherwise, we use Binary Search method for calculation.
            else:
                ind = BinarySearchTangent(hh, p_k)
                q_set.append((hh[ind], i, ind))
        
        # We compare each point in q_set to select the best one by
        # computing the cos value among each three points.
        cos_val = 1
        p_next = None
        # print(p_k_1, p_k, q_set)
        for q, i, ind in q_set:
            if not (q[0]==p_k[0] and q[1]==p_k[1]):
                vector1 = (q[0]-p_k[0], q[1]-p_k[1])
                cur_cos = (vector1[0]*vector[0]+vector1[1]*vector[1])/(dis*math.sqrt(vector1[0]**2+vector1[1]**2))
                # We update the hullind and posind.
                if cur_cos<cos_val:
                    cos_val = cur_cos
                    p_next = q
                    hullind = i
                    posind = ind
        
        # If p_next is same to the bighull[1], which is lp,
        # then we know that we have already selected all the points
        # we need in this big Convex Hull. So we return it.
        if p_next[0]==bighull[1][0] and p_next[1]==bighull[1][1]:
            # bighull.append(p_next)
            return bighull[1:]
        
        # Otherwise, we add this new point into bighull and continue.
        bighull.append(p_next)
        
        if test:
            PlotTest(P, sub_hull, bighull, q_set)
    return []


def Chan2D(P, test):
    n = len(P)
    t = 1
    while True:
        m = H = min(1<<(1<<t), n)
        L = Hull2D(P, m, H, test)
        if L:
            return L
        t += 1


# This function the results by Chan's Algorithm with Graham's Scan.
# -------------------------------------------------------------- #
# Parameters:
# V1: The hull point array generated using Graham's Scan algo.
# V2: The hull point array generated using Chan's algo.
# 
# Return(s): True if two point sets are totally same, False otherwise
# -------------------------------------------------------------- #
def check(V1, V2):
    if len(V1)!=len(V2):
        return False
    for (x1, y1), (x2, y2) in zip(V1, V2):
        if x1!=x2 or y1!=y2:
            return False
    return True


# Main function starts from here:
if __name__=='__main__':
    testmode = True if int(sys.argv[1])==1 else False
    logic = int(sys.argv[2])
    
    if logic not in [1,2,3]:
        sys.exit("Wrong params input.")

    # Confirm the correctness of Chan's algorithm
    if logic==1:
        T = testcases = 100
        res = 0
        while T>0:
            T -= 1
            P = GenerateDataRect(150, 0, 100, 10, 40)
            # PlotData(P)
            V_confirm = GrahamScan(P)
            # PlotDataWithHull(P, V)
            V = Chan2D(P, testmode)
            # PlotDataWithHull(P, V)
            if check(V_confirm, V):
                res += 1
        print("Among %d test cases, there are %d wrong example." % (testcases, testcases-res))

    # illustrate the whole process of Chan's algorithm
    elif logic==2:
        number = int(sys.argv[3])
        limit = int(sys.argv[4])
        while True:
            P = GenerateDataRect(number, 0, 100, 10, 40)
            if len(GrahamScan(P))<=limit:
                break
        V = Chan2D(P, testmode)
        PlotDataWithHull(P, V)
        
    # Compute the total time used of two different methods
    else:
        T = 100
        t_graham = t_chan = 0
        for number in range(10, 151, 10):
            print("The number of points: %d" % number)
            for _ in tqdm(range(T)):
                P = GenerateDataRect(number, 0, 100, 10, 40)
                t_start_graham = time.time()
                GrahamScan(P)
                t_end_graham = time.time()
                t_start_chan = time.time()
                Chan2D(P, testmode)
                t_end_chan = time.time()
                t_graham += t_end_graham - t_start_graham
                t_chan += t_end_chan - t_start_chan
            print("The total time for Graham's Scan Algorithm is: %.6fs." % t_graham)
            print("The total time for Chan's Algorithm is: %.6fs." % t_chan)