import openrouteservice
import xml.etree.ElementTree as ET
import folium
import geopy
import sys

MPATH = "DistanceMatrix.txt"
BUILDING_DATA = "BuildingData.xml"
DIM = 5
RITCOORDS = [43.08395011861679, -77.67580747604372]

def main():
    print(sys.argv)
    client = openrouteservice.Client(key = "5b3ce3597851110001cf6248bdc117dcf33b439eb7cac98d1e1329b6") # TODO add API key
    global CAMPUS # Global because of scope issues. It's never edited once created anyway.
    CAMPUS = buildRIT(BUILDING_DATA)
    show_best_route(sys.argv[1], sys.argv[2], client)

class BuildingGroup:
    """
    "Group" of building classes. Contains buildings in a list in the subbuildings field.
    """
    def __init__(self, name, subbuildings, entrances = []):
        self.name = name
        self.subbuildings = subbuildings
        self.entrances = entrances

    def __string__(self):
        return self.name

    def distance_entrance(self, buildinggroup, client):
        """
        Returns the CLOSEST distance from this building to another, and the corresponding
        entrances. We do this in one function to avoid wasting computation. To do this we'll
        use the distance matrix API as it uses less calls.
        """
        # Figure out if we already have the distance stored
        # Get the indices for the group
        start = self.getIndex()
        end = buildinggroup.getIndex()

        val = getMatrixVal(start, end, MPATH)
        if val:
            # TODO make sure matrix file protocol order is consistent
            source = self.entrances[val[1]]
            destination = buildinggroup.entrances[val[2]]
            return val[0], (source, destination)
        
        else:
            # This code runs if the entry doesn't exist
        
            # Build lambda function for getting list of coords
            getCoords = lambda e : e.coords
            destcoords = list(map(getCoords, buildinggroup.entrances))
            sourcecoords = list(map(getCoords, self.entrances))
            coords = sourcecoords + destcoords

            # Making indices for the source and destinations in the list of coords we pass
            # to the API.
            sourcesidx = list(range(len(sourcecoords)))
            destinationsidx = list(range(len(sourcecoords), len(sourcecoords) + len(destcoords)))

            # Distance Matrix API call
            input("Making API call... press enter to continue")
            matrixString = openrouteservice.distance_matrix.distance_matrix(client, coords, profile = "foot-walking", destinations = destinationsidx, sources = sourcesidx, metrics = ["distance"])

            # Processing the API return

            # Distance matrix
            distances = matrixString["distances"]

            # Finding the best two entrances

            sourcedist = []
        
            for source in distances:
                mindest = min(source)
                sourcedist.append((mindest, source.index(mindest)))

            sourceidx = sourcedist.index(min(sourcedist))
            destidx = sourcedist[sourceidx][1]
            
            source = self.entrances[sourceidx]
            destination = buildinggroup.entrances[destidx]
            finaldistance = distances[sourceidx][destidx]

            # Write this new value to the file
            writeDistMatrix(start, end, (finaldistance, sourceidx, destidx), MPATH)

            # Returning the results

            return finaldistance, (source, destination)

    def getBestSequence(self, buildinggroup, client):
        """
        Use Dijkstra's Algorithm to get the best sequence of buildings to get to the destination. We'll return
        this as a list.
        """
        distinfo = lambda e : [e, 9999999, None]
        
        nodes = list(map(distinfo, CAMPUS.copy())) # Fresh list for the algorithm. Represents unvisited nodes.
        visited = []
        bestlist = []
        nodes[self.getIndex()][1] = 0 # Set starting vertex to zero
        current = nodes[self.getIndex()]
        start = current
        current[1] = 0
        prev = None
        done = False
        while not done:
            # Consider all unvisited neighbors
            minnode = [0,99999999999, None] # For finding minimum
            for node in nodes:
                if node == current:
                    continue
                # Check if node was visited
                nodedistance = self.distance_entrance(node[0], client)[0]
                newdistance = current[1] + nodedistance
                if newdistance < node[1]:
                    node[1] = newdistance
                # Get min value if needed
                if newdistance < minnode[1]:
                    minnode = node
            nodes.remove(current)
            current[2] = prev
            prev = current
            visited.append(current)
            if current[0] == buildinggroup:
                done = True
            else:
                current = minnode

        # Backtrack and get best route
        path = []
        while current:
            #print(current)
            path.insert(0, current[0])
            current = current[2]

        return path        

    def getIndex(self):
        return CAMPUS.index(self)

class Building:
    def __init__(self, name, ID, coords):
        self.name = name
        self.ID = ID
        self.coords = coords
    def __string__(self):
        return self.name
    def distance(self, building, client):
        """
        Returns the CLOSEST distance from this building to another. To do this we'll
        use the distance matrix API as it uses less calls.
        """
        pass
    def findParent(self, name, ID):
        """
        Figure out which buildinggroup object this building belongs to.
        """
        # Loop through each buildinggroup
        parent = None
        for group in CAMPUS:
            for building in group:
                if building == self:
                    parent = group
                    break

        return parent

class Entrance:
    def __init__(self, name, coords = (0, 0)):
        self.name = name #TODO Write a routine for entrance names (e.g. south)
        self.coords = coords
        
    def __string__(self):
        return self.name
        

def buildRIT(file: str):
    """
    Returns a list with all the building groups and their subclasses
    """
    tree = ET.parse(file)
    root = tree.getroot()
    buildinggrouptags = root.findall("BUILDINGGROUP") # All the buildinggroup tag roots in xml
    buildinggroups = [] # Starting empty list. We'll append to this and return it at the end
    for grouproot in buildinggrouptags:
        buildinggroups.append(buildGroup(grouproot))

    return buildinggroups

def buildGroup(grouproot):
    """
    Returns a buildingGroup object. Should contain sub-buildings and a list of entrances.
    """

    groupname = grouproot.attrib["name"]
    buildingtags = grouproot.findall(".//BUILDING") # All the building tag roots in xml
    subbuildings = [] # Starts empty, at the end we'll use it to construct the building group.
    entrancetags = grouproot.findall(".//ENTRANCE") # All the entrances in xml
    entrances = [] # Starts empty, we'll append to it to create the entrance list

    for building in buildingtags:
        # Create the building objects and add them to the list
        name = building.find("NAME").text
        ID = building.find("ID").text
        coords = tuple(map(float, building.find("COORDS").text.split(",")))
        newb = Building(name, ID, coords)
        subbuildings.append(newb)


    for idx, entrance in enumerate(entrancetags):
        # Create the entrance objects and add them to the list
        name = groupname + " entrance " + str(idx)
        coords = tuple(map(float, entrance.text.split(",")))
        ent = Entrance(name, coords)
        entrances.append(ent)


    # Throw this stuff into the group initialilzer

    grp = BuildingGroup(groupname, subbuildings, entrances)

    return grp

# TODO consider putting the matrix protocol stuff in another file

def getMatrix(file):
    matrix_file = open(file)
    matrix = matrix_file.read()
    matrix_file.close()
    return eval(matrix) # This is kind of unsafe... I don't have a better way yet though

def storeMatrix(file, matrix):
    matrix_file = open(file, "w") # This should clear the file before writing
    matrix_file.write(str(matrix))
    matrix_file.close()

def writeDistMatrix(src, dest, val, file):
    """
    Write to the distance matrix cache
    """
    # Load the matrix
    matrix = getMatrix(file)
    # Write the value
    matrix[src][dest] = val
    # Store the matrix
    storeMatrix(file, matrix)

def getMatrixVal(src, dest, file):
    """
    Check if a value exists in the distance matrix cache at the given location
    """
    # Load the matrix
    matrix = getMatrix(file)
    # Return the value
    return matrix[src][dest]

def clearDistMatrix(file):
    """
    Empty the distance matrix cache. DIM is a global variable for the size of the matrix 
    """

    # Make the empty matrix
    matrix = [[False] * DIM for x in range(DIM)]
    # Write to file
    storeMatrix(file, matrix)

def getGeoJSON(grouplist, client):
    """
    Returns a GPX file based on the list of buildinggroups provided
    """

    coords = []
    jsonlist = []
    for idx in range(0, len(grouplist) - 1):
        group = grouplist[idx]
        nextgroup = grouplist[idx + 1]
        _, coordspair = group.distance_entrance(nextgroup, client)
        coords1 = list(coordspair[0].coords)
        coords2 = list(coordspair[1].coords)
        coords = [coords1, coords2]
        gpxstring = str(openrouteservice.directions.directions(client, coords, profile = "foot-walking", format = "geojson")).replace("'",'"')
        jsonlist.append(gpxstring)
        #coords.extend([coords1, coords2])

    return jsonlist

def display_routes(center, jsonlist, startcoords, endcoords):

    m = folium.Map(location = center, zoom_start = 16.48)

    folium.CircleMarker(location = startcoords, radius = 2, weight = 5).add_to(m)
    folium.CircleMarker(location = endcoords, radius = 2, weight = 5).add_to(m)

    for track in jsonlist:
        folium.GeoJson(track).add_to(m)
    
    m.show_in_browser()

def show_best_route(building1str, building2str, client):
    # TODO add markers for each building using address lookup
    group1, group2 = None, None
    
    for group in CAMPUS:
        for building in group.subbuildings:
            if building1str.lower() in building.name.lower():
                point1 = building.coords[::-1]
                group1 = group
            if building2str.lower() in building.name.lower():
                point2 = building.coords[::-1]
                group2 = group

    seq = group1.getBestSequence(group2, client)
    jsonlist = getGeoJSON(seq, client)
    display_routes(RITCOORDS, jsonlist, point1, point2)

def getBuildingCoords(ID):
    locator = geopy.geocoders.Nominatim(user_agent = "myGeocoder")
    location = locator.geocode(str(ID) + " Lomb Memorial Drive")
    return (location.latitude, location.longitude)

if __name__ == "__main__":
    main()
