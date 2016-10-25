#-------------------------------------------------------------------------------
# Name:        mapMatching
# Purpose:      This python script allows map matching in arcpy using a Hidden Markov model with
#               probabilities parameterized based on spatial + network distances. Follows the ideas in Newson, Krumm (2009):
#               Hidden markov Map Matching through noise and sparseness
#
# Author:      Simon Groen, Simon Scheider
#
# Created:     25/10/2016
# Copyright:   (c) simon 2016
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from math import exp



#A simple map matching algorithm in Python, based on the idea of the Viterbi Algorithm/Hidden Markov Model, but allowing on-the fly distances in a GIS

def getPDProbability(dist):
    #The probability that given a certain distance between points and segments, the point is on the segment
    #This needs to be parameterized
    #test: turn index difference into a probability (exponential function)
    p = 1 if dist == 0 else round(1/exp(dist),4)
    return p

def getSegmentCandidates(point, segments, maxdist=10):
    #Returns closest segment candidates with probabilities. Based on maximal spatial distance of segments from point, returning a probability for each segment

    #test: calculates distances for points based on list index
    print point
    points = ['p1', 'p2', 'p3', 'p4', 'p5', 'p6']
    #segments = ['s1', 's2', 's3', 's4', 's5', 's6']
    candidates = {}
    for j in range(len(points)):
        if points[j] == point:
            #print j
            for i in range(len(segments)):
                #print i
                diff = abs(j-i)
                candidates[segments[i]] = getPDProbability(diff)
            break
    print candidates
    return candidates

def getNDProbability(dist):
    #The probability that given a certain network distance between segments, one is the successor of the other in a track
    #This needs to be parameterized
    #test: turn index difference into a probability (exponential function)
    p = 1 if dist == 0 else  round(1/exp(dist),2)
    return p

def getNetworkTransP(s1, s2, segments):
    #Returns transition probability of going from segment s1 to s2, based on network distance of segments

    #test: calculates distances for segments based on list index
    points = ['p1', 'p2', 'p3', 'p4', 'p5', 'p6']
    #segments = ['s1', 's2', 's3', 's4', 's5', 's6']
    diff = 0
    for j in range(len(segments)):
        if segments[j] == s1:
            for i in range(len(segments)):
                if segments[i] == s2:
                    diff = abs(j-i)

    return getNDProbability(diff)



#This is the actual method. It gets two lists of point and segment identifiers, and returns the most probable segment path for the list of points

def mapMatch(points, segments):
    #this array stores, for each point, probability distributions over segments, together with the (most probable) predecessor segment
    V = [{}]
    #init first point
    sc = getSegmentCandidates(points[0], segments)
    #print sc
    for s in sc:
        V[0][s] = {"prob": sc[s], "prev": None}

    # Run Viterbi when t > 0
    for t in range(1, len(points)):
        V.append({})
        #Store previous segment candidates
        lastsc = sc
        #Get segment candidates for current point t
        sc = getSegmentCandidates(points[t], segments)
        for s in sc:
            max_tr_prob = 0
            prev_ss =None
            for prev_s in lastsc:
                 #determine the most probable transition probability from previous candidates to s
                tr_prob = V[t-1][prev_s]["prob"]*getNetworkTransP(prev_s, s, segments)
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

    print 'The path for points ['+' '.join(points)+'] is: [' + ' '.join(opt) + '] with highest probability of %s' % max_prob
    return opt



#test

points = ['p1', 'p2', 'p3', 'p4', 'p5', 'p6']
segments = ['s1', 's2', 's3', 's4', 's5', 's6']
#segments = list(reversed(segments))

mapMatch(points, segments)




























def viterbi(obs, states, start_p, trans_p, emit_p):
#taken from https://en.wikipedia.org/wiki/Viterbi_algorithm
    V = [{}]

    for st in states:
        V[0][st] = {"prob": start_p[st] * emit_p[st][obs[0]], "prev": None}

    # Run Viterbi when t > 0
    for t in range(1, len(obs)):
        V.append({})
        for st in states:
            max_tr_prob = max(V[t-1][prev_st]["prob"]*trans_p[prev_st][st] for prev_st in states)
            for prev_st in states:
                if V[t-1][prev_st]["prob"] * trans_p[prev_st][st] == max_tr_prob:
                    max_prob = max_tr_prob * emit_p[st][obs[t]]
                    V[t][st] = {"prob": max_prob, "prev": prev_st}
                    break

    for line in dptable(V):
        print line

    opt = []

    # The highest probability

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
        opt.insert(0, V[t + 1][previous]["prev"])
        previous = V[t + 1][previous]["prev"]


    print 'The steps of states are ' + ' '.join(opt) + ' with highest probability of %s' % max_prob


def dptable(V):
    # Print a table of steps from dictionary
    yield " ".join(("%12d" % i) for i in range(len(V)))
    for state in V[0]:
        yield "%.7s: " % state + " ".join("%.7s" % ("%f" % v[state]["prob"]) for v in V)

#simple example:
states = ('Healthy', 'Fever')
observations = ('normal', 'cold', 'dizzy')
start_probability = {'Healthy': 0.6, 'Fever': 0.4}
transition_probability = {
   'Healthy' : {'Healthy': 0.7, 'Fever': 0.3},
   'Fever' : {'Healthy': 0.4, 'Fever': 0.6}
   }

emission_probability = {
   'Healthy' : {'normal': 0.5, 'cold': 0.4, 'dizzy': 0.1},
   'Fever' : {'normal': 0.1, 'cold': 0.3, 'dizzy': 0.6}
   }

#viterbi(observations,states,start_probability,transition_probability,emission_probability)
