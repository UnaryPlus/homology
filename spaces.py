from homology import SimplicialSpace

def sphere(n):
    faces = []
    for k in range(n + 2):
        face = list(range(k)) + list(range(k + 1, n + 2))
        faces.append(face)
    return SimplicialSpace(faces)

def torus():
    faces = ["ABF", "BFG", "BCG", "ACG", "DFG", "DEG", "AEG", "AEF", "ABD", "BDE", "BCE", "CEF", "CDF", "ACD"]
    return SimplicialSpace(faces)

def klein_bottle():
    faces = ["AEF", "ABF", "BFG", "BCG", "CDG", "ACD", "DEF", "DFH", "FGH", "EGH", "DEG", "ABD", "BDH", "BCH", "ACH", "AEH"]
    return SimplicialSpace(faces)