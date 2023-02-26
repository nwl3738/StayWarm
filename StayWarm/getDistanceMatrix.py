import openrouteservice
import xml.etree.ElementTree as ET
import pprint

def main():
    client = openrouteservice.Client(key = "")
    print(getDistance("gym", "sau", "BuildingData.xml", client))

def getMatrices(file, client):
    """
    Responsible for reading the xml file data and getting a list of coords.
    Then calls openrouteservice to get the coords into a text file, which can
    then be parsed by the main program.
    """
    tree = ET.parse(file)
    root = tree.getroot()
    coordstags = tree.findall(".//ENTRANCE")
    #coords = [list(map(float, tag.text.split(","))) for tag in coordstags] # Finally, a good excuse to use map()
    print(coords)

    #matrixout = open("DistanceMatrix.txt", "w")
    #matrixout.write(matrixString)
    #matrixout.close()

def getDistance(buildinggroup1: str, buildinggroup2: str, file: str, client):
    # Gets a distance matrix between two building groups
    tree = ET.parse(file)
    root = tree.getroot()
    groups = root.findall("BUILDINGGROUP")

    # Find the two building groups we're looking for
    entrances = []
    sizes = []
    for element in groups:
        if element.attrib["name"] == buildinggroup1 or element.attrib["name"] == buildinggroup2:
            new_elements = element.findall(".//ENTRANCE")
            print("ELEMENTS: ", new_elements)
            entrances += new_elements
            sizes.append(len(new_elements))
            print(sizes)
            

    print(entrances)
    coords = [list(map(float, tag.text.split(","))) for tag in entrances] # Finally, a good excuse to use map()
    print(coords)
    sources = list(range(sizes[0]))
    destinations = list(range(sizes[0], sizes[0] + sizes[1]))
    print(sources, destinations)

    #matrixString = openrouteservice.distance_matrix.distance_matrix(client, coords, profile = "foot-walking", destinations = destinations, sources = sources, metrics = ["distance"])

    distances = matrixString["distances"]

    # Find the lowest distance pair

    sourcedist = []
    
    for source in distances:
        mindest = min(source)
        sourcedist.append((mindest, source.index(mindest)))

    sourceidx = sourcedist.index(min(sourcedist))
    destidx = sourcedist[sourceidx][1]
    
    bestpairidx = (sourceidx, destidx)
    #print(bestpairidx)


    #pp = pprint.PrettyPrinter()
    #pp.pprint(matrixString)

    return bestpairidx, distances[bestpairidx[0]][bestpairidx[1]]
    
if __name__ == "__main__":
    main()
