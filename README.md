# Map Matcher

This python script allows map matching (matching of tracking points to a network)
in arcpy using a Hidden Markov model with
probabilities parameterized based on spatial + network distances.
Follows the ideas in Newson, Krumm (2009):
"Hidden markov Map Matching through noise and sparseness"

Author:      Simon Scheider

Created:     16/03/2017
   

## Installation

To install as a toolbox in ArcGIS, see [mapMaptch.pyt](#mapmatchpyt-arcgis-python-toolbox)

The code is written in Python 2.7 and depends on:

* arcpy (ships with ArcGIS and its own Python 2.7)
* [networkx](https://networkx.github.io) 

`python pip install networkx`

* note: requires installing GDAL first, which can be obtained as a windows wheel from [here](http://www.lfd.uci.edu/~gohlke/pythonlibs/)
     and then installed with pip locally:

`python pip install GDAL-2.1.3-cp27-cp27m-win32.whl`

To install the mapmatcher Python module, simply download and execute this windows executable:
- [mapmatching/dist/mapmatcher-1.0.win32.exe](https://github.com/simonscheider/mapmatching/blob/master/dist/mapmatcher-1.0.win32.exe).    

## Usage

Example:

`from mapmatcher import mapmatcher`

`arcpy.env.workspace = 'C:/Users/simon/Documents/GitHub/mapmatching'`

`opt = mapmatcher.mapMatch('testTrack.shp', 'testSegments.shp')`

`#outputs testTrack_pth.shp`

`mapmatcher.exportPath(opt, 'testTrack.shp')`

The last method saves a new shape file named _testTrack_pth.shp_ in the current arcpy workspace, containing a sequence of segments to which the track was mapped.

Results are shown here:

<img src="https://github.com/simonscheider/mapmatching/blob/master/example.PNG" width="500" />


The main method is _mapMatch_. Based on the Viterbi algorithm for Hidden Markov models,
see https://en.wikipedia.org/wiki/Viterbi_algorithm, it gets trackpoints and segments, and returns the most probable segment path (a list of segments) for the list of points.

### Method _mapMatch_:
* @param **track** = a shape file (filename) with point geometries representing a track, can be unprojected (WGS84). The order of points in this file should reflect the temporal order.
        
* @param **segments** = a shape file of network segments, should be _projected_ (in meter) to compute Euclidean distances properly (e.g. GCS Amersfoord). _Note_: To compute network distances, the script turns this network into a graph using [networkx](https://networkx.github.io), based on coincidence of segment end points. It is therefore important that logically connected segments are also geometrically connected (no geometrical errors). Other than this, the script does have other requirements for the network.
        
* @param _decayconstantNet_ (optional) = the network distance (in meter) after which the match probability falls under 0.34 ([exponential decay](https://en.wikipedia.org/wiki/Exponential_decay)). Default is 30 meters. This distance parameter depends on the intervals between successing points in the track.
        
* @param _decayConstantEu_ (optional) = the Euclidean distance (in meter) after which the match probability falls under 0.34 (exponential decay). Default is 10 meters. This distance parameter depends on the measurement accuracy of tracking points.
        
* @param _maxDist_ (optional) = the Euclidean distance threshold (in meter) for taking into account segments candidates. Default is 50 meters. Depends also on measurement accuracy of track points.

* result = delivers back a path (a list of segment ids).

#### Note: 
Depending on the type of movement, optional parameters need to be fine tuned to get optimal results. For example, when tracking frequency is very slow, then track points are far apart, and then _decayconstantNet_ needs to be increased accordingly.

### Method _exportPath_ :
Exports the path into a  shape file named _segments_pth.shp_ inside the current ArcGIS workspace.


# mapMatch.pyt (ArcGIS Python toolbox)

To use the Python method as an ArcGIS toolbox, you need to do the following:

1. In your ArcGIS Python version (e.g. Folder `C:\Python27\ArcGIS10.3\Lib\site-packages`), install required modules for GDAL and networkx in a cmd window:

- if you have not installed it yet, install pip (http://pip.readthedocs.io/en/latest/installing/). 
- Download a suitable GDAL wheel from [here](http://www.lfd.uci.edu/~gohlke/pythonlibs/). Then execute:
- `python pip install GDAL-2.1.3-cp27-cp27m-win32.whl`
- `python pip install networkx`

2. Install mapmatcher Python module by downloading and executing the windows executable [mapmatching/dist/mapmatcher-1.0.win32.exe](https://github.com/simonscheider/mapmatching/blob/master/dist/mapmatcher-1.0.win32.exe). Make sure you select exactly the Python installation that ships with your ArcGIS as a target folder.

3. Download the ArcGIS Python toolbox [mapMatch.pyt](https://github.com/simonscheider/mapmatching/blob/master/mapMatch.pyt), together with meta data files [mapMatch.mapMatch.pyt.xml](https://github.com/simonscheider/mapmatching/blob/master/mapMatch.mapMatch.pyt.xml) and [mapMatch.pyt.xml](https://github.com/simonscheider/mapmatching/blob/master/mapMatch.pyt.xml) and drop it anywhere on your computer.

4. Now you can open the toolbox by clicking on it inside an ArcGIS Catalog Window:
<img src="https://github.com/simonscheider/mapmatching/blob/master/mapMatch.PNG" width="500" />

The tool saves a new shape file named _NameofInputTrack_pth.shp_ inside the current ArcGIS workspace that contains the path of segments to which the track was mapped. When executing, make sure the network is as small as possible to speed up.
