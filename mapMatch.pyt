
import arcpy
from mapmatcher import mapmatcher

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "mapMatch toolbox"
        self.alias = "mapMatcher"

        # List of tool classes associated with this toolbox
        self.tools = [mapMatch]


class mapMatch(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "mapMatch"
        self.description = ""
        '''This tool allows map matching (matching of track points to a network)
        in arcpy using a Hidden Markov model with
               probabilities parameterized based on spatial + network distances.
               Follows the ideas in Newson, Krumm (2009):
               "Hidden markov Map Matching through noise and sparseness"
               '''
        self.canRunInBackground = False

    def getParameterInfo(self):
        #Define parameter definitions

        # Input Features parameter
        in_track = arcpy.Parameter(
            displayName="Input track (shp)",
            name="in_track",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        #in_track.filter.type = 'FeatureClass'
        #in_track.filter.list = ["Point"]

        # Input Features parameter
        in_roads = arcpy.Parameter(
            displayName="Input road network (shp)",
            name="in_roads",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")

        #in_roads.filter.type = 'FeatureClass'
       # in_roads.filter.list = ["Polyline"]

        # parameter 1
        decayconstantNet = arcpy.Parameter(
            displayName="Decay constant on network (in meter)",
            name="decayconstantNet",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        decayconstantNet.value = "30"

        # parameter 2
        decayConstantEu = arcpy.Parameter(
            displayName="Decay constant for Eucl. distances (in meter)",
            name="decayconstantEu",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        decayConstantEu.value = "10"

        # parameter 2
        maxDist = arcpy.Parameter(
            displayName="maximum Eucl. distance for selecting segment candidates (in meter)",
            name="maxDist",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input")

        maxDist.value = "50"

##        # Derived Output Features parameter
##        out_features = arcpy.Parameter(
##            displayName="Output Features",
##            name="out_features",
##            datatype="GPFeatureLayer",
##            parameterType="Derived",
##            direction="Output")

        #out_features.parameterDependencies = [in_features.name]
        #out_features.schema.clone = True

        parameters = [in_track, in_roads, decayconstantNet, decayConstantEu, maxDist]

        return parameters

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        in_track  = parameters[0].valueAsText
        in_roads   = parameters[1].valueAsText
        decayconstantNet = parameters[2].valueAsText
        decayConstantEu = parameters[3].valueAsText
        maxDist = parameters[4].valueAsText

        opt = mapmatcher.mapMatch(in_track, in_roads,  decayconstantNet, decayConstantEu,maxDist)
        mapmatcher.exportPath(opt, in_track)



