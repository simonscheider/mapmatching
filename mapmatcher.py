"""
-------------------------------------------------------------------------------
# Name:        mapMatcher
# Purpose:      This python script allows map matching (matching of track points to a network)
#               in arcpy using a Hidden Markov model with
#               probabilities parameterized based on spatial + network distances.
#               Follows the ideas in Newson, Krumm (2009):
#               "Hidden markov Map Matching through noise and sparseness"
#
#               Example usage under '__main__'
#
# Author:      Simon Scheider
#
# Created:     17/02/2017
# Copyright:   (c) simon 2017
# Licence:     <your licence>

The code is written in Python 2.7 and depends on:

* arcpy (ships with ArcGIS and its own Python 2.7)
* networkx (# python pip install networkx (https://networkx.github.io))
    (note: requires installing GDAL first, which can be obtained as a wheel from
    http://www.lfd.uci.edu/~gohlke/pythonlibs/ and then installed with pip locally:
    python pip install GDAL-2.1.3-cp27-cp27m-win32.whl
    )

#-------------------------------------------------------------------------------
"""

__author__      = "Simon Scheider"
__copyright__   = ""

from math import exp, sqrt
import os
import arcpy
arcpy.env.overwriteOutput = True
import networkx as nx


def mapMatch(track, segments, decayconstantNet = 30, decayConstantEu = 10, maxDist = 50):
    """
    The main method. Based on the Viterbi algorithm for Hidden Markov models,
    see https://en.wikipedia.org/wiki/Viterbi_algorithm.
    It gets trackpoints and segments, and returns the most probable segment path (a list of segments) for the list of points.
    Inputs:
        @param track = a shape file (filename) representing a track, can also be unprojected (WGS84)
        @param segments = a shape file of network segments, should be projected (in meter) to compute Euclidean distances properly (e.g. GCS Amersfoord)
        @param decayconstantNet (optional) = the network distance (in meter) after which the match probability falls under 0.34 (exponential decay). (note this is the inverse of lambda)
        @param decayConstantEu (optional) = the Euclidean distance (in meter) after which the match probability falls under 0.34 (exponential decay). (note this is the inverse of lambda)
        @param maxDist (optional) = the Euclidean distance threshold (in meter) for taking into account segments candidates.

    note: depending on the type of movement, optional parameters need to be fine tuned to get optimal results.
    """

    #this array stores, for each point in a track, probability distributions over segments, together with the (most probable) predecessor segment taking into account a network distance
    V = [{}]

    #get track points, build network graph (graph, endpoints, lengths) and get segment info from arcpy
    points = getTrackPoints(track, segments)
    r = getSegmentInfo(segments)
    endpoints = r[0]
    lengths = r[1]
    graph = getNetworkGraph(segments,lengths)

    #init first point
    sc = getSegmentCandidates(points[0], segments, decayConstantEu, maxDist)
    for s in sc:
        V[0][s] = {"prob": sc[s], "prev": None, "path": None}
    # Run Viterbi when t > 0
    for t in range(1, len(points)):
        V.append({})
        #Store previous segment candidates
        lastsc = sc
        #Get segment candidates and their a-priori probabilities (based on Euclidean distance for current point t)
        sc = getSegmentCandidates(points[t], segments, decayConstantEu, maxDist)
        for s in sc:
            max_tr_prob = 0
            prev_ss = None
            path = None
            for prev_s in lastsc:
                #determine the most probable transition probability from previous candidates to s and get the corresponding network path
                n = getNetworkTransP(prev_s, s, graph, endpoints, decayconstantNet)
                np = n[0] #This is the network transition probability
                tr_prob = V[t-1][prev_s]["prob"]*np
                #this selects the most probable predecessor candidate and the path to it
                if tr_prob > max_tr_prob:
                    max_tr_prob = tr_prob
                    prev_ss = prev_s
                    path = n[1]
            #The final probability of a candidate is the product of a-priori and network transitional probability
            max_prob =  sc[s] * max_tr_prob
            V[t][s] = {"prob": max_prob, "prev": prev_ss, "path": path}

    #print V

    #opt is the result: a list of (matched) segments [s1, s2, s3,...] in the exact order of the point track: [p1, p2, p3,...]
    opt = []

    # get the highest probability at the end of the track
    max_prob = max(value["prob"] for value in V[-1].values())
    previous = None

    # Get most probable ending state and its backtrack
    for st, data in V[-1].items():
        if data["prob"] == max_prob:
            opt.append(st)
            previous = st
            break

    # Follow the backtrack till the first observation to fish out most probable states and corresponding paths
    for t in range(len(V) - 2, -1, -1):
        #Get the path between last and most probable previous segment and add it to the resulting path
        path = V[t + 1][previous]["path"]
        opt[0:0] =(path if path !=None else [])
        #Insert the previous segment
        opt.insert(0, V[t + 1][previous]["prev"])
        previous = V[t + 1][previous]["prev"]
    pointstr= [str(g.firstPoint.X)+' '+str(g.firstPoint.Y) for g in points]
    optstr= [str(i) for i in opt]
    print 'The path for points ['+' '.join(pointstr)+'] is: '
    print '[' + ' '.join(optstr) + '] with highest probability of %s' % max_prob

    return opt

def exportPath(opt, trackname):
    """
    This exports the list of segments into a shapefile, a subset of the loaded segment file, including all attributes
    """
    qr =  '"OBJECTID" IN ' +str(tuple(opt))
    outname = os.path.splitext(trackname)[0]+'_path'
    arcpy.SelectLayerByAttribute_management('segments_lyr',"NEW_SELECTION", qr)
    try:
        if arcpy.Exists(outname):
            arcpy.Delete_management(outname)
        arcpy.FeatureClassToFeatureClass_conversion('segments_lyr', arcpy.env.workspace, outname)
    except Exception:
        e = sys.exc_info()[1]
        print(e.args[0])

        # If using this code within a script tool, AddError can be used to return messages
        #   back to a script tool.  If not, AddError will have no effect.
        arcpy.AddError(e.args[0])


def getPDProbability(dist, decayconstant = 10):
    """
    The probability that given a certain distance between points and segments, the point is on the segment
    This needs to be parameterized
    Turn difference into a probability with exponential decay function
    """
    p = 1 if dist == 0 else round(1/exp(dist/decayconstant),4)
    return p

def getSegmentCandidates(point, segments, decayConstantEu, maxdist=50):
    """
    Returns closest segment candidates with a-priori probabilities.
    Based on maximal spatial distance of segments from point.
    """
    p = point.firstPoint #get the coordinates of the point geometry
    print "Neighbors of point "+str(p.X) +' '+ str(p.Y)+" : "
    #Select all segments within max distance
    arcpy.Delete_management('segments_lyr')
    arcpy.MakeFeatureLayer_management(segments, 'segments_lyr')
    arcpy.SelectLayerByLocation_management ("segments_lyr", "WITHIN_A_DISTANCE", point, maxdist)
    candidates = {}
    #Go through these, compute distances, probabilities and store them as candidates
    cursor = arcpy.da.SearchCursor('segments_lyr', ["OBJECTID", "SHAPE@"])
    for row in cursor:
        feat = row[1]
        #compute the spatial distance
        dist = point.distanceTo(row[1])
        #compute the corresponding probability
        candidates[row[0]] = getPDProbability(dist, decayConstantEu)
    del row
    del cursor
    print str(candidates)
    return candidates


def getNDProbability(dist,decayconstant = 30):
    """
    The probability that given a certain network distance between segments, one is the successor of the other in a track
    This needs to be parameterized
    Turn difference into a probability  with exponential decay function
    """
    p = 1 if dist == 0 else  round(1/exp(dist/decayconstant),2)
    return p

def getNetworkTransP(s1, s2, graph, endpoints, decayconstantNet):
    """
    Returns transition probability of going from segment s1 to s2, based on network distance of segments, as well as corresponding path
    """
    subpath = None
    if s1 == s2:
        dist = 0
    else:
        #Obtain edges (tuples of endpoints) for segment identifiers
        s1_edge = endpoints[s1]
        s2_edge = endpoints[s2]

        #This determines segment endpoints of the two segment that are closest to each other
        minpair = [0,0,100000]
        for i in range(0,2):
            for j in range(0,2):
                d = round(pointdistance(s1_edge[i],s2_edge[j]),2)
                if d<minpair[2]:
                    minpair = [i,j,d]
        s1_point = s1_edge[minpair[0]]
        s2_point = s2_edge[minpair[1]]

        if s1_point == s2_point:
            dist = 0
        else:
            if (not graph.has_node(s1_point)) or (not graph.has_node(s2_point)):
                print "node not in segment graph!"
            try:
                #This computes a shortes path (using segment length) on a graph where segment endpoints are nodes and segments are (undirected) edges
                dist = nx.shortest_path_length(graph, s1_point, s2_point, weight='length')
                path = nx.shortest_path(graph, s1_point, s2_point, weight='length')
                #get path edges
                path_edges = zip(path,path[1:])
                #print "edges: "+str(path_edges)
                subpath = []
                # get object ids for path edges
                for e in path_edges:
                    oid = graph.edge[e[0]][e[1]]["OBJECTID"]
                    subpath.append(oid)
                #print "oid path:"+str(subpath)
            except nx.NetworkXNoPath:
                print 'no path available, assume a large distance'
                dist = 100
    #print "network distance between "+str(s1) + ' and '+ str(s2) + ' = '+str(dist)
    return (getNDProbability(dist,decayconstantNet),subpath)

def pointdistance(p1, p2):
    dist = sqrt((p1[0]-p2[0])**2 +(p1[1]-p2[1])**2)
    return dist

def getTrackPoints(track, segments):
    """
    Turns track shapefile into a list of point geometries, reprojecting to the planar RS of the network file
    """
    trackpoints = []
    for row in arcpy.da.SearchCursor(track, ["SHAPE@"]):
        #make sure track points are reprojected to network reference system (should be planar)
        geom = row[0].projectAs(arcpy.Describe(segments).spatialReference)
        trackpoints.append(row[0])
    return trackpoints

def getNetworkGraph(segments,segmentlengths):
    """
    Builds a networkx graph from the network file, inluding segment length taken from arcpy.
    It selects the largest connected component of the network (to prevent errors from routing between unconnected parts)
    """
    g = nx.read_shp(segments)
    #This selects the largest connected component of the graph
    sg = list(nx.connected_component_subgraphs(g.to_undirected()))[0]
    print "graph size (excluding unconnected parts): "+str(len(g))
    # Get the length for each road segment and append it as an attribute to the edges in the graph.
    for n0, n1 in sg.edges_iter():
        oid = sg[n0][n1]["OBJECTID"]
        sg.edge[n0][n1]['length'] = segmentlengths[oid]
    return sg

def getSegmentInfo(segments):
    """
    Builds a dictionary for looking up endpoints of network segments (needed only because networkx graph identifies edges by nodes)
    """
    cursor = arcpy.da.SearchCursor(segments, ["OBJECTID", "SHAPE@"])
    endpoints = {}
    segmentlengths = {}
    for row in cursor:
          endpoints[row[0]]=((row[1].firstPoint.X,row[1].firstPoint.Y), (row[1].lastPoint.X, row[1].lastPoint.Y))
          segmentlengths[row[0]]= row[1].length
    del row
    del cursor
    print "Number of segments: "+ str(len(endpoints))
    return (endpoints,segmentlengths)


if __name__ == '__main__':

    #Test using the shipped data example
    arcpy.env.workspace = 'C:/Users/simon/Documents/GitHub/mapmatching'
    opt = mapMatch('testTrack.shp', 'testSegments.shp')
    #outputs testTrack_path.shp
    exportPath(opt, 'testTrack.shp')
