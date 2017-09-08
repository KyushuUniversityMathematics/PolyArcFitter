#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Poly-arc fitter


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.linalg import solve
import argparse
import sys

# コマンドライン引数解析
parser = argparse.ArgumentParser(description='poly arc fitting of ordered points')
parser.add_argument('dataset', help='Path to data file')
parser.add_argument('--deviation_tolerance', '-d', type=float, default=0.1,
                        help='maximum allowed deviation from the original points')
parser.add_argument('--corner_angle', '-c', type=float, default=20.0,
                        help='minimum angle in degrees to be classified as a corner')
parser.add_argument('--arc_limit', '-a', type=float, default=300,
                        help='points separated more than this value will be connected by a line segment')
parser.add_argument('--visualise', '-v', action='store_true',
                        help='visualise the result using matplotlib')
parser.add_argument('--show_circle', action='store_true',
                        help='draw full circles rather than arcs')
args = parser.parse_args()

cos_theta = np.cos(np.pi*args.corner_angle/180)

sys.stderr.write('corner angle: {}, deviation tolerance: {}, arc length limit: {}, dataset: {}\n'.format(args.corner_angle,args.deviation_tolerance,args.arc_limit,args.dataset))


# 角(角度が大きい点と隣の点と離れている点)を見つけ、そのインデックスのリストを返す
def find_corners(polyline):
    #角度の大きい点
    p =polyline[1:,:]-polyline[:-1,:]
    norm_p=np.linalg.norm(p,axis=1)
    cosine =np.sum(p[:-1]*p[1:],axis=1)/(norm_p[:-1]*norm_p[1:])
    (angled , )=np.where(cosine < cos_theta)
    angled = angled +1
    #離れている点
    (distant, )=np.where(norm_p > args.arc_limit)
    distant_end = distant + 1
    #点をまとめてソート
    corners = [0,len(polyline)-1]
    corners.extend(angled.tolist()+distant.tolist()+distant_end.tolist())
    corners = list(set(corners))
    corners.sort()
    sys.stderr.write('Number of corners: {}\n'.format(len(corners)))
    return corners

#polylineを一つの線分で近似できるか判定し、できる場合は開始地点、終了地点、最大エラーを返す
def segment_fit(polyline):
    vs = polyline[0]
    ve = polyline[-1]
    d = ve - vs
    d = d / np.linalg.norm(d)
    max_error = np.max(np.abs(np.cross(d, polyline - vs)))
    if  max_error < args.deviation_tolerance:
        return 'line', vs, ve, max_error
    else:
        return None

#polylineを一つの円弧で近似できるか判定し、
#近似できる場合はその中心,半径,開始地点、終了地点、開始角、終了角、最大エラーを返す
def arc_fit(polyline):
    o, R = approximate_arc(polyline)
    max_error = np.max(np.abs(np.sqrt(np.sum((polyline - o)**2,axis=1))-R))
    if max_error>args.deviation_tolerance:
        return None
    # starting and ending angles in degrees
    ds, de, dm = polyline[0] - o, polyline[-1] - o, polyline[len(polyline)//2] - o
    th1, th2 = getArgsSet(ds,de,dm)
    return 'arc', o, R, polyline[0], polyline[-1], th1, th2, max_error

#polylineを近似する円弧の中心、半径を返す
def approximate_arc(polyline):
    n = len(polyline)
    xs, ys = polyline[0]
    xe, ye = polyline[-1]
    x_sum,y_sum= np.sum(polyline,axis=0)
    xx_sum,yy_sum =np.sum(polyline**2,axis=0)
    xy_sum = np.sum(np.prod(polyline,axis=1))
    xxx_sum,yyy_sum =np.sum(polyline**3,axis=0)
    xxy_sum=np.sum((polyline**2)[:,0]*polyline[:,1],axis=0)
    xyy_sum=np.sum(polyline[:,0]*(polyline**2)[:,1],axis=0)
    
    # only upper triangular part is necessary
    A = np.array([[2*xx_sum, 2*xy_sum, 2*x_sum, xs, xe],
                  [2*xy_sum, 2*yy_sum, 2*y_sum, ys, ye],
                  [2*x_sum, 2*y_sum, 2*n, 1, 1],
                  [xs, ys, 1, 0, 0],
                  [xe, ye, 1, 0, 0]])
    
    bb = np.array([-2 * (xxx_sum + xyy_sum),
                  -2 * (xxy_sum + yyy_sum),
                  -2 * (xx_sum + yy_sum),
                  -(xs ** 2 + ys ** 2),
                  -(xe ** 2 + ye ** 2)])
    
    # 対称行列なので scipy の方が早いけれど
    #sol = solve(A, bb, overwrite_a=True, overwrite_b=True, assume_a = 'sym')
    # scipy 使いたくなければ、上の代わりに下を
    sol = np.linalg.solve(A, bb)
    
    a, b, c, _, _ = sol
    
    return np.array([-a/2, -b/2]), np.sqrt(a**2/4 + b**2/4 - c)

# polyline を、エラーが deviation_tolerance に収まるように、適宜分割しつつ円弧・線分補間 
def segments_arcs_fit(polyline, seg_arc_list):
    while True:
        end, seg_arc = binary_search(polyline)
        seg_arc_list.append(seg_arc)
        if end >= len(polyline)-1:
            break
        polyline = polyline[end:]

# 分割点を二分探索
def binary_search(polyline):
    high = len(polyline)
    low = 2
    cur = high
    seg_arc = 'line', polyline[0], polyline[1], 0.0
    while True:
        try_seg_arc = segment_fit(polyline[:cur])
        if try_seg_arc:  # try line segment first
            low = cur
            seg_arc = try_seg_arc
        else:
            try_seg_arc = arc_fit(polyline[:cur])
            if try_seg_arc:
                low = cur
                seg_arc = try_seg_arc
            else:
                high = cur
        cur = (low + high + 1) // 2
        if high-low <= 1:
            #print("OK",try_seg_arc[0], cur, low, high)
            return low-1, seg_arc

# datum を poly-arc で近似して、arc のリストを返す
def polyarc_fit(datum):
    corners = find_corners(datum)
    seg_arc_list = []
    for i in range(len(corners)-1):
        #print("corner: ",corners[i],corners[i+1])
        segments_arcs_fit(datum[corners[i]:corners[i+1]+1], seg_arc_list)    
    return seg_arc_list

# 表示
def draw(datum,seg_arc_list):
    # draw original points
    plt.plot(datum[:, 0], datum[:, 1], 'x', alpha=0.3, c='b')
    #plt.plot(datum[0,0],datum[0,1], '>', alpha=0.5, c='r')
    #plt.plot(datum[-1,0],datum[-1,1], '<', alpha=0.5, c='r')
    # draw the simplified poly curve    
    ax = plt.axes()
    max_error = 0
    num_lines = 0
    num_arcs = 0
    for seg_arc in seg_arc_list:
        if seg_arc[0] == 'line':
            num_lines += 1
            _, vs, ve, part_error = seg_arc
            max_error = max(part_error,max_error)
            # draw a line segment
            plt.plot([vs[0],ve[0]], [vs[1],ve[1]], alpha=0.5, c='g')
        elif seg_arc[0] == 'arc':
            num_arcs += 1
            _, o, R, vs, ve, th1, th2, part_error = seg_arc
            max_error = max(part_error,max_error)
            # draw circle and arc
            circ = patches.Circle(xy=o, radius=R, ec='y', fill=False, alpha=0.5)
            if th2>th1:
                arc = patches.Arc(xy=o, width=2*R, height=2*R, ec='r', theta1=th1*180/np.pi, theta2=th2*180/np.pi)
            else:
                arc = patches.Arc(xy=o, width=2*R, height=2*R, ec='r', theta1=th2*180/np.pi, theta2=th1*180/np.pi)
            ax.add_patch(arc)
            if args.show_circle:
                ax.add_patch(circ)
        # draw the start-end points
        plt.plot(vs[0], vs[1], 'v', alpha=0.5, c='k')
        plt.plot(ve[0], ve[1], '^', alpha=0.5, c='k')
        
    sys.stderr.write('Maximum deviation: {} \n # of segments and arcs: {}, # of lines {}, # of arcs {}\n'.format(max_error,len(seg_arc_list),num_lines,num_arcs))
    # show plot window
    plt.axes().set_aspect('equal', 'datalim')    
    plt.show()

## 角度の関係
#(x,y):対象とするベクトルが (vx,vy):基準とするベクトルの右手にあるか
def isPointRightHand(x, y, vx, vy):
    return(vy*x - vx*y > 0)

#弧の向きが時計回りか判定
#vs:始点 ve:終点 vm:中継点
def isCircleClockwise(ds,de,dm):
    if isPointRightHand(de[0], de[1], ds[0], ds[1]):
        return (isPointRightHand(dm[0], dm[1], ds[0], ds[1]) and isPointRightHand(de[0], de[1], dm[0], dm[1]))
    else:
        return (isPointRightHand(dm[0], dm[1], ds[0], ds[1]) or isPointRightHand(de[0], de[1], dm[0], dm[1]))

#弧の向きを加味した始点と終点のラジアン角を返す
#vs:始点 ve:終点 vm:中継点
def getArgsSet(ds,de,dm):
    vs_arg = np.arctan2(ds[1], ds[0])
    ve_arg = np.arctan2(de[1], de[0])
    vm_arg = np.arctan2(dm[1], dm[0])
    
    if isCircleClockwise(ds,de,dm):
        if isPointRightHand(dm[0], dm[1], ds[0], ds[1]):
            if vs_arg < 0 and ve_arg > 0:
                vs_arg += 2*np.pi
        else:
            if vs_arg < 0:
                vs_arg += 2*np.pi
            elif ve_arg > 0:
                ve_arg -= 2*np.pi
    else:
        if isPointRightHand(de[0], de[1], ds[0], ds[1]):
            if ve_arg < 0:
                ve_arg += 2*np.pi
            elif vs_arg > 0:
                vs_arg -= 2*np.pi
        else:
            if vs_arg > 0 and ve_arg < 0:
                ve_arg += 2*np.pi
    return (vs_arg, ve_arg)


# フォーマット済みテキストを出力
def formated_text(datum,seg_arc_list):
    max_error = 0
    is_new_line_segment = True
    for seg_arc in seg_arc_list:
        if seg_arc[0] == 'line':
            _, vs, ve, part_error = seg_arc
            max_error = max(part_error,max_error)
            is_new_line_segment = False
            if is_new_line_segment:
                print('{}\t{}'.format(vs[0],vs[1]))
            print('{}\t{}'.format(ve[0],ve[1]))
        elif seg_arc[0] == 'arc':
            _, o, R, vs, ve, th1, th2, part_error = seg_arc
            max_error = max(part_error,max_error)
            is_new_line_segment = True
            print('{}\t{}\t{}\t{}\t{}'.format(R,o[0],o[1],th1,th2))
    sys.stderr.write('Maximum deviation: {}\n'.format(max_error))

##
if __name__ == '__main__':
    datum = np.loadtxt(args.dataset, delimiter='\t')
    seg_arc_list = polyarc_fit(datum)

    if args.visualise:
        draw(datum,seg_arc_list)
    else:
        formated_text(datum,seg_arc_list)
