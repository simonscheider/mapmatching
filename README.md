# Map Matcher

This python script (mapmatcher.py) allows map matching (matching of track points to a network)
in arcpy using a Hidden Markov model with
probabilities parameterized based on spatial + network distances.
Follows the ideas in Newson, Krumm (2009):
"Hidden markov Map Matching through noise and sparseness"

Author:      Simon Scheider

Created:     17/02/2017
   

## The code is written in Python 2.7 and depends on:

* arcpy (ships with ArcGIS and its own Python 2.7)
* [networkx](https://networkx.github.io) 

`python pip install networkx`

* note: requires installing GDAL first, which can be obtained as a windows wheel from [here](http://www.lfd.uci.edu/~gohlke/pythonlibs/)
     and then installed with pip locally:

`python pip install GDAL-2.1.3-cp27-cp27m-win32.whl`
    

## Usage
Example:

`arcpy.env.workspace = 'C:/Users/simon/Documents/GitHub/mapmatching'`

`opt = mapMatch('testTrack.shp', 'testSegments.shp')`

`#outputs testTrack_path.shp`

`exportPath(opt, 'testTrack.shp')`

Results are shown here:

<img src="https://github.com/simonscheider/mapmatching/blob/master/example.PNG" width="500" />


The main method is _mapMatch_. Based on the Viterbi algorithm for Hidden Markov models,
see https://en.wikipedia.org/wiki/Viterbi_algorithm, it gets trackpoints and segments, and returns the most probable segment path (a list of segments) for the list of points.

### Method _mapMatch_:
* @param **track** = a shape file (filename) representing a track, can be unprojected (WGS84)
        
* @param **segments** = a shape file of network segments, should be _projected_ (in meter) to compute Euclidean distances properly (e.g. GCS Amersfoord)
        
* @param _decayconstantNet_ (optional) = the network distance (in meter) after which the match probability falls under 0.34 ([exponential decay](https://en.wikipedia.org/wiki/Exponential_decay)). Default is 30 meters.
        
* @param _decayConstantEu_ (optional) = the Euclidean distance (in meter) after which the match probability falls under 0.34 (exponential decay). Default is 10 meters.
        
* @param _maxDist_ (optional) = the Euclidean distance threshold (in meter) for taking into account segments candidates. Default is 50 meters.

* result = delivers back a path (a list of segment ids)

#### Note: 
Depending on the type of movement, optional parameters need to be fine tuned to get optimal results.

### Method _exportPath_ :
exports the path into a shape file






