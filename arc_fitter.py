#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Poly-arc fitter
# by Akira Hirakawa, Chihiro Matsufuji, Yosuke Onitsuka, Shizuo Kaji
# This work is partially supported by the Kyushu university IMI short term research programme
# "Evaluation technique for a three-dimensional geometric modeling and software development"
# held in Sep. 2017


from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
#from scipy.linalg import solve
import argparse
import sys

# throughout the code "polyline" means an ordered sequence of 2D points

# least square fit the polyline by an arc
def least_sq_fit(polyline):
    X = np.hstack(( polyline, np.ones((len(polyline), 1)) ))
    Y = np.vstack((X[0], X[-1]))    
    p = np.sum(polyline**2, axis=1)
    q = np.array([p[0], p[-1]])    
    A = np.vstack(( np.hstack((np.dot(X.T, X), Y.T)),
                    np.hstack((Y, np.zeros((2,2)) )) ))    
    bb = -np.hstack((np.dot(X.T, p), q))    
    # with scipy: only upper triangular part is necessary
    # a, b, c, _, _ = solve(A, bb, overwrite_a=True, overwrite_b=True, assume_a = 'sym')
    # without scipy
    a, b, c, _, _ = np.linalg.solve(A, bb)
    o = -np.array([a, b])/2    # centre
    R = np.sqrt(a**2/4 + b**2/4 - c)   # radius
    max_error = np.max(np.abs(np.sqrt(np.sum((polyline - o)**2,axis=1))-R))

    if R > args.max_radius or max_error > args.deviation_tolerance:
        return None
    # compute starting and ending angles in degrees
    ds, de, dm = polyline[0] - o, polyline[-1] - o, polyline[len(polyline)//2] - o
    # make three points counter-clockwise
    if np.cross(dm-ds, de-dm) <= 0:
        ds, de = de, ds
    th1, th2 = getArgsForCCW(ds,de)
    return 'arc', o, R, th1, th2, max_error

# middle centre fit the polyline by an arc; counter-clockwise or clockwise
def middle_centre_fit(polyline):
    return under_arc_fit_with_max(polyline, polyline[0], polyline[-1])\
            or under_arc_fit_with_max(polyline, polyline[-1], polyline[0])


# dictionary of fitting functions
arc_fit_function = {
    'least_sq': least_sq_fit,
    'middle': middle_centre_fit,
}

# command line arguments
parser = argparse.ArgumentParser(description='poly arc fitting of ordered points')
parser.add_argument('dataset', help='Path to data file')
parser.add_argument('--arc_limit', '-a', type=float, default=300,
                        help='points separated more than this value will be connected by a line segment')
parser.add_argument('--deviation_tolerance', '-d', type=float, default=0.1,
                        help='maximum allowed deviation (in Eucledian distance) from the original points')
parser.add_argument('--corner_points', '-cp', default='',
                        help='Path to a csv file containing indices of corner points')
parser.add_argument('--corner_angle', '-ca', type=float, default=20.0,
                        help='angle in degrees greater than this will be classified as a corner')
parser.add_argument('--max_radius', '-r', type=float, default=float('inf'),
                        help='Maximum allowed radius of circular arcs.')
parser.add_argument('--visualise', '-v', action='store_true',
                        help='visualise the result using matplotlib')
parser.add_argument('--show_circle', action='store_true',
                        help='draw full circles rather than arcs')
parser.add_argument('--fit', '-f', choices=arc_fit_function.keys(), default='least_sq',
                        help='fitting method')
parser.add_argument('--max_ratio', '-mr', type=float, default=10,
                        help='points with the distance ratio to adjacent points greater than this value will be classified as a corner.')
args = parser.parse_args()

# cosine of the corner angle
cos_theta = np.cos(np.pi*args.corner_angle/180)

# print the settings
sys.stderr.write('method: {}, corner angle: {}, deviation tolerance: {}, arc length limit: {}, max radius: {},  extreme limit: {}, dataset: {}\n'.format(args.fit,args.corner_angle,args.deviation_tolerance,args.arc_limit,args.max_radius,args.max_ratio,args.dataset))


# detect corner points, which the final polyarc must pass
def find_corners(polyline):
    # find large angles
    p = polyline[1:,:]-polyline[:-1,:] # tangent vector
    norm_p=np.linalg.norm(p,axis=1)
    cosine = np.sum(p[:-1]*p[1:],axis=1)/(norm_p[:-1]*norm_p[1:])
    (angled, )=np.where(cosine < cos_theta)
    angled = angled +1    # i-th point is between p[i-1] and p[i]
    # find distant points
    (distant, )=np.where(norm_p > args.arc_limit)
    distant_end = distant + 1
    # find three consecutive points with a big distant ratio
    f1 = norm_p[1:]/norm_p[:-1]
    f2 = norm_p[:-1]/norm_p[1:]
    ratio = np.max(np.vstack((f1,f2)),axis=0)
    (extreme,) = np.where(ratio > args.max_ratio)
    extreme = extreme + 1
    # first and last points are corners
    corners = [0,len(polyline)-1]
    # user specified corner points
    if args.corner_points:
        corners.extend(np.loadtxt(args.corner_points, delimiter=",", dtype=np.int32))
    # amalgamate all corner points
    corners.extend(angled.tolist()+distant.tolist()+distant_end.tolist()+extreme.tolist())
    corners = list(set(corners)) # remove dup
    corners.sort()
    sys.stderr.write('Number of corners: {}\n'.format(len(corners)))
    return corners

# returns the line connecting the first and the last points on the polyline
# if the maximum Eucledian distance between the points on the polyline and the line is below the tolerance
# TODO: compute the distance with the line "segment" rather than the line
def segment_fit(polyline):
    vs, ve = polyline[0], polyline[-1]
    max_error = np.max(np.abs(np.cross((ve-vs)/np.linalg.norm(ve-vs), polyline - vs)))
    return ('line', vs, ve, max_error) if max_error < args.deviation_tolerance else None

# sub-routine for under_arc_fit_with_max
# determine the allowed interval for t
# more precisely, solution to abs( sqrt(x**2 + (y-t)**2) - sqrt(d**2 + t**2) ) = e
def get_boundary(x, y, d, e):
    a = e*np.sqrt((d**2 - e**2 + x**2 + y**2)**2 - (2*d*x)**2)
    b = y*(d**2 + e**2 - x**2 - y**2)
    c = 2*(e-y)*(e+y)
    return (a+b)/c, (-a+b)/c

# check if polyline can be approximated by an arc connecting vs and ve counter-clockwisely    
def under_arc_fit_with_max(polyline, vs, ve, eps=1e-10):
    # set coordinates so that vs -> (-d,0), ve ->(d,0)
    d = np.linalg.norm(ve-vs)/2
    e = args.deviation_tolerance
    lo = (vs+ve)/2
    # axis vectors
    x_axis = (ve-vs)/np.linalg.norm(ve-vs)
    y_axis = np.array([-x_axis[1], x_axis[0]])
    # fail if radius is too big
    if d > args.max_radius:
        return None
    
    # range of the centre location (0,t) to meet the max_radius requirement
    limit = np.sqrt(args.max_radius**2 - d**2)
    low, high = -limit, limit
    
    # consider only those points that are farther from end points than tolerance (*)
    poly_t = polyline[(np.linalg.norm(polyline-vs, axis=1) > e) & (np.linalg.norm(polyline-ve, axis=1) > e)]
    poly_t -= lo    
    # coordinate transform for points
    x = np.dot(poly_t, x_axis)
    y = np.dot(poly_t, y_axis)
    
    # there are two arcs connecting the end points; one goes above x-axix, the other under it.
    # here we consider only the latter.
    if np.any(y >= e-eps):
        return None
    
    # find indices of points which lie between x=-d and x=d (inner) and outside the range (outer)
    inner_ind = np.abs(x) < d
    outer_ind = np.abs(x) >= d

    # find indices of points which is close to x-axis (close)
    close_ind = np.abs(y) <= e - eps

    # for inner-close points, the range of t is [t1,inf]
    close_inner = close_ind & inner_ind
    if np.any(close_inner):
        t1, _ = get_boundary(x[close_inner], y[close_inner], d, e)  # always t1>t2
        low = max(low, np.max(t1))
    
    # for outer-close points, the range of t is [-inf,t2]
    close_outer = close_ind & outer_ind
    if np.any(close_outer):
        _, t2 = get_boundary(x[close_outer], y[close_outer], d, e)  # always t1>t2
        high = min(high, np.min(t2))
    
    # note that due to (*), we only have above two cases for close points


    under_fit = np.abs(y+e) < eps
    
    # when a point's neighbour touches x-axis, the equation degenerates to a linear one. Handle this case.
    under_inner = under_fit & inner_ind    
    if np.any(under_inner):
        x_masked = x[under_inner]
        low = max(low, np.max(((d**2 - x_masked**2)**2 - (2*d*e)**2)/(4*e*(d**2 - x_masked**2))))

    under_outer = under_fit & outer_ind
    if np.any(under_outer):
        x_masked = x[under_outer]
        high = min(high, np.min(((d**2 - x_masked**2)**2 - (2*d*e)**2)/(4*e*(d**2 - x_masked**2))))
    
    # the "usual case" when points are below x-axis
    far_ind = (y <= -e-eps)
    if np.any(far_ind):
        t1, t2 = get_boundary(x[far_ind], y[far_ind], d, e)
        low = max(low, np.max(t1))
        high = min(high, np.min(t2))  
    
    # fail if the range is empty
    if low > high:
        return None
    
    t = (low+high)/2
    o = lo + t*y_axis
    R = np.linalg.norm(o-vs)
    
    max_error = np.max(np.abs(np.sqrt(np.sum((polyline - o)**2,axis=1))-R))
    if max_error>args.deviation_tolerance:
        return None

    th1, th2 = getArgsForCCW(vs-o, ve-o)
    
    return 'arc', o, R, th1, th2, max_error

# angles for the counter-clockwise arc connecting ds de: note th1<=th2
def getArgsForCCW(ds, de):
    th1 = np.arctan2(ds[1], ds[0])
    th2 = np.arctan2(de[1], de[0])
    th2 += np.ceil((th1 - th2)/(2*np.pi))*2*np.pi
    if th2 >= 2*np.pi:
        th1 -= 2*np.pi
        th2 -= 2*np.pi
    return th1, th2 

# main loop to divide and approximater polyline
def segments_arcs_fit(polyline):
    seg_arc_list = []
    while True:
        end, seg_arc = binary_search(polyline)
        seg_arc_list.append(seg_arc)
        if end >= len(polyline)-1:
            return seg_arc_list
        polyline = polyline[end:]

# find knots by binary search
def binary_search(polyline):
    high = len(polyline)
    low = 2
    cur = high
    seg_arc = 'line', polyline[0], polyline[1], 0.0
    while True:
        try_seg_arc = np.any(polyline[0] != polyline[-1]) and \
            (segment_fit(polyline[:cur]) or arc_fit_function[args.fit](polyline[:cur]))
        if try_seg_arc:
            low = cur
            seg_arc = try_seg_arc
        else:
            high = cur
        cur = (low + high + 1) // 2
        if high-low <= 1:
            return low-1, seg_arc

# draw the result
def draw(datum,seg_arc_list):
    # draw original points
    plt.plot(datum[:, 0], datum[:, 1], 'x', alpha=0.3, c='b')
    #plt.plot(datum[0,0],datum[0,1], '>', alpha=0.5, c='r')
    #plt.plot(datum[-1,0],datum[-1,1], '<', alpha=0.5, c='r')
    # draw the polyarc   
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
            _, o, R, th1, th2, part_error = seg_arc
            max_error = max(part_error,max_error)
            # draw an arc
            arc = patches.Arc(xy=o, width=2*R, height=2*R, ec='r', theta1=th1*180/np.pi, theta2=th2*180/np.pi)
            ax.add_patch(arc)
            if args.show_circle:
                circ = patches.Circle(xy=o, radius=R, ec='y', fill=False, alpha=0.5)
                ax.add_patch(circ)
            vs = o + np.array([R*np.cos(th1),R*np.sin(th1)])
            ve = o + np.array([R*np.cos(th2),R*np.sin(th2)])
        # draw the start-end points
        plt.plot(vs[0], vs[1], 'v', alpha=0.5, c='k')
        plt.plot(ve[0], ve[1], '^', alpha=0.5, c='k')
        
    sys.stderr.write('Maximum deviation: {} \n # of segments and arcs: {}, # of lines {}, # of arcs {}\n'.format(max_error,len(seg_arc_list),num_lines,num_arcs))
    # show plot window
    plt.axes().set_aspect('equal', 'datalim')    
    plt.show()


# format output
def formated_text(datum,seg_arc_list):
    max_error = 0
    is_new_line_segment = True
    for seg_arc in seg_arc_list:
        if seg_arc[0] == 'line':
            _, vs, ve, part_error = seg_arc
            if is_new_line_segment:
                print('{}\t{}'.format(vs[0],vs[1]))
            is_new_line_segment = False
            print('{}\t{}'.format(ve[0],ve[1]))
        elif seg_arc[0] == 'arc':
            _, o, R, th1, th2, part_error = seg_arc
            is_new_line_segment = True
            print('{}\t{}\t{}\t{}\t{}'.format(R,o[0],o[1],th1,th2))
        max_error = max(part_error,max_error)
    sys.stderr.write('Maximum deviation: {}\n'.format(max_error))

## start here
if __name__ == '__main__':
    datum = np.loadtxt(args.dataset, delimiter='\t')
    corners = find_corners(datum)
    seg_arc_list = [item for i in range(len(corners)-1) for item in segments_arcs_fit(datum[corners[i]:corners[i+1]+1])]
    if args.visualise:
        draw(datum,seg_arc_list)
    else:
        formated_text(datum,seg_arc_list)