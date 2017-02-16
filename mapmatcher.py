#-------------------------------------------------------------------------------
# Name:        mapMatcher
# Purpose:      This python script allows map matching in arcpy using a Hidden Markov model with
#               probabilities parameterized based on spatial + network distances. Follows the ideas in Newson, Krumm (2009):
#               Hidden markov Map Matching through noise and sparseness
#
# Author:      Simon Scheider
#
# Created:     25/10/2016
# Copyright:   (c) simon 2016
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from math import exp
import os
from math import sqrt
import arcpy
arcpy.env.workspace = 'C:/Users/simon/Documents/GitHub/mapmatching'
arcpy.env.overwriteOutput = True
import networkx as nx
import numpy as np





#A simple map matching algorithm in Python, based on the idea of the Viterbi Algorithm/Hidden Markov Model, but allowing on-the fly distances in a GIS

def getPDProbability(dist, decayconstant = 10):
    #The probability that given a certain distance between points and segments, the point is on the segment
    #This needs to be parameterized
    #test: turn index difference into a probability (exponential function)
    p = 1 if dist == 0 else round(1/exp(dist/decayconstant),4)
    return p

def getSegmentCandidates(point, segments, maxdist=50):
    #Returns closest segment candidates with probabilities. Based on maximal spatial distance of segments from point, returning a probability for each segment
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
        candidates[row[0]] = getPDProbability(dist)
    del row
    del cursor
    print str(candidates)
    return candidates
    #test: calculates distances for points based on list index

def getNDProbability(dist,decayconstant = 30):
    #The probability that given a certain network distance between segments, one is the successor of the other in a track
    #This needs to be parameterized
    #test: turn index difference into a probability (exponential function)
    p = 1 if dist == 0 else  round(1/exp(dist/decayconstant),2)
    return p

def getNetworkTransP(s1, s2, graph, endpoints):
    #Returns transition probability of going from segment s1 to s2, based on network distance of segments
    if s1 == s2:
        dist = 0
    else:
        s1_edge = endpoints[s1]
        s2_edge = endpoints[s2]

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
                dist = nx.shortest_path_length(graph, s1_point, s2_point, weight='length')
            except nx.NetworkXNoPath:
                print 'no path available, assume a large distance'
                dist = 9
    print "network distance between "+str(s1) + ' and '+ str(s2) + ' = '+str(dist)
    return getNDProbability(dist)

def pointdistance(p1, p2):
    dist = sqrt((p1[0]-p2[0])**2 +(p1[1]-p2[1])**2)
    return dist

def getTrackPoints(track):
    trackpoints = []
    for row in arcpy.da.SearchCursor(track, ["SHAPE@"]):
        #make sure tracks are reprojected to network reference system (should be planar)
        geom = row[0].projectAs(arcpy.Describe('testSegments.shp').spatialReference)
        trackpoints.append(row[0])
    return trackpoints

def getNetworkGraph(segments,segmentlengths):
    g = nx.read_shp(segments)
    #This selects the largest connected component of the graph
    sg = list(nx.connected_component_subgraphs(g.to_undirected()))[0]
    print "road network size: "+str(len(g))
    # Get the length for each road segment.
    for n0, n1 in sg.edges_iter():
        oid = sg[n0][n1]["OBJECTID"]
        sg.edge[n0][n1]['length'] = segmentlengths[oid]
    return sg

def getSegmentInfo(segments):
    cursor = arcpy.da.SearchCursor(segments, ["OBJECTID", "SHAPE@"])
    endpoints = {}
    segmentlengths = {}
    for row in cursor:
          endpoints[row[0]]=((row[1].firstPoint.X,row[1].firstPoint.Y), (row[1].lastPoint.X, row[1].lastPoint.Y))
          segmentlengths[row[0]]= row[1].length
    del row
    del cursor
    print "Number of segments: "+ str(len(segmentendpoints))
    return (endpoints,segmentlengths)




#This is the actual method. It gets two lists of point and segment identifiers, and returns the most probable segment path for the list of points

def mapMatch(points, segments):
    #Based on the Viterbi algorithm, see https://en.wikipedia.org/wiki/Viterbi_algorithm
    #this array stores, for each point in a track, probability distributions over segments, together with the (most probable) predecessor segment
    V = [{}]
    #initiate network graph and get segment info about segments from arcpy
    endpoints = getSegmentInfo(segments)[0]
    lengths = getSegmentInfo(segments)[1]
    graph = getNetworkGraph(segments,lengths)

    #init first point
    sc = getSegmentCandidates(points[0], segments)
    for s in sc:
        V[0][s] = {"prob": sc[s], "prev": None}
    # Run Viterbi when t > 0
    for t in range(1, len(points)):
        V.append({})
        #Store previous segment candidates
        lastsc = sc
        #Get segment candidates for current point t
        sc = getSegmentCandidates(points[t], segments)
        #print str(sc)
        for s in sc:
            max_tr_prob = 0
            prev_ss = None
            for prev_s in lastsc:
                 #determine the most probable transition probability from previous candidates to s
                tr_prob = V[t-1][prev_s]["prob"]*getNetworkTransP(prev_s, s, graph, endpoints)
                if tr_prob > max_tr_prob:
                    max_tr_prob = tr_prob
                    prev_ss = prev_s
            max_prob = max_tr_prob * sc[s]
            V[t][s] = {"prob": max_prob, "prev": prev_ss}

    print V

    #opt is the result: a list of (matched) segments [s1, s2, s3,...] in the exact order of the point track: [p1, p2, p3,...]
    opt = []

    # get the highest probability at the end of the track
    max_prob = max(value["prob"] for value in V[-1].values())
    previous = None

    # Get most probable state and its backtrack
    for st, data in V[-1].items():
        if data["prob"] == max_prob:
            opt.append(st)
            previous = st
            break

    # Follow the backtrack till the first observation
    for t in range(len(V) - 2, -1, -1):
        #print V[t + 1][previous]["prev"]
        opt.insert(0, V[t + 1][previous]["prev"])
        previous = V[t + 1][previous]["prev"]
    pointstr= [str(g.firstPoint.X)+' '+str(g.firstPoint.Y) for g in points]
    optstr= [str(i) for i in opt]
    print 'The path for points ['+' '.join(pointstr)+'] is: [' + ' '.join(optstr) + '] with highest probability of %s' % max_prob
    return opt

def exportPath(opt, trackname):
    qr =  '"OBJECTID" IN ' +str(tuple(opt))
    print qr
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




points = getTrackPoints('testTrack.shp')
opt = mapMatch(points, 'testSegments.shp')
exportPath(opt, 'testTrack.shp')
