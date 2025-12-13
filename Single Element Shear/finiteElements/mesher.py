import gmsh
import sys

fname = "mesh.msh"

elemOrder = 1
if sys.argv[1] == "hex8":
    elemOrder = 1
if sys.argv[1] == "hex20":
    elemOrder = 2

gmsh.initialize()
gmsh.model.add(fname)

nx = 4
ny = 25
nz = 2

# nx = 1
# ny = 1
# nz = 1

Lx = 1
Ly = 10
Lz = 0.2

# Lx = 1
# Ly = 1
# Lz = 1

# Add points for main surface extrude
p1 = gmsh.model.geo.addPoint(0, 0, 0)
p2 = gmsh.model.geo.addPoint(Lx, 0, 0)
p3 = gmsh.model.geo.addPoint(Lx, Ly, 0)
p4 = gmsh.model.geo.addPoint(0, Ly, 0)

# Lines between points
l1 = gmsh.model.geo.addLine(p1, p2)
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p4)
l4 = gmsh.model.geo.addLine(p4, p1)

# Curve Loop for extrusion/surface
cl1 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
s1 = gmsh.model.geo.addPlaneSurface([cl1])

# Transfinite Curve - Sets discretization
# Opposing sides need to have the same number of nodes
# N = (n + 1) Points on curve for n Elements
gmsh.model.geo.mesh.setTransfiniteCurve(l1, nx+1)
gmsh.model.geo.mesh.setTransfiniteCurve(l3, nx+1)
gmsh.model.geo.mesh.setTransfiniteCurve(l2, ny+1)
gmsh.model.geo.mesh.setTransfiniteCurve(l4, ny+1)

gmsh.model.geo.mesh.setTransfiniteSurface(s1, "Alternate", [l1, l2, l3, l4])
gmsh.model.geo.mesh.setRecombine(2, s1)  # Recombine quads to hexas

# Extrude Surface defined by Curve Loop 1,
# Extrude for Lz, with nz subdivisions,
# Recombine quads to hexas
ov = gmsh.model.geo.extrude([(2, cl1)], 0, 0, Lz, [nz], recombine=True)

# Named Selections (Surfs)
mesh = gmsh.model.geo.addPhysicalGroup(ov[1][0], [ov[1][1]])  # Mesh
back = gmsh.model.geo.addPhysicalGroup(2, [s1])  # Back
front = gmsh.model.geo.addPhysicalGroup(ov[0][0], [ov[0][1]])  # Front
bot = gmsh.model.geo.addPhysicalGroup(ov[2][0], [ov[2][1]])  # Bottom
right = gmsh.model.geo.addPhysicalGroup(ov[3][0], [ov[3][1]])  # Right
top = gmsh.model.geo.addPhysicalGroup(ov[4][0], [ov[4][1]])  # Top
left = gmsh.model.geo.addPhysicalGroup(ov[5][0], [ov[5][1]])  # Left

gmsh.model.setPhysicalName(2, top, "fixed_top")
gmsh.model.setPhysicalName(2, bot, "fixed_bot")

# CAD Elements sync
gmsh.model.geo.synchronize()

gmsh.option.setNumber("Mesh.Smoothing", 20)
gmsh.model.mesh.generate(dim=3)  # 3D Mesh

# setOrder(1) - Generate Linear Elements
# setOrder(2) - Generate Serendipity Elements
gmsh.option.setNumber('Mesh.SecondOrderIncomplete', 1)  # Important!
gmsh.model.mesh.setOrder(elemOrder)

gmsh.option.setNumber("Geometry.PointNumbers", 1)
gmsh.option.setColor("Geometry.Color.Points", 255, 165, 0)
gmsh.option.setColor("General.Color.Text", 255, 255, 255)
gmsh.option.setColor("Mesh.Color.Points", 255, 0, 0)

gmsh.write(fname)
gmsh.finalize()

gmsh.initialize()
gmsh.open(fname)

entities = gmsh.model.getEntities()
entities.reverse()

file = open("nodeCoords.txt", "w")
file.close()
file = open("hexaCon.txt", "w")
file.close()
file = open("physicalGroups.txt", "w")
file.write("{:s},{:d}".format("inner", 0) + "\n")
file.close()

for e in entities:
    eType = gmsh.model.getType(e[0], e[1])
    if eType == "Discrete volume":
        # print("Entity " + str(e) + " of type " + eType)
        nodeTagsElem, nodeCoords, nodeParams = gmsh.model.mesh.getNodes(
            e[0], e[1])
        elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(
            e[0], e[1])
        tags, coords = gmsh.model.mesh.getNodesForPhysicalGroup(e[0], e[1])
        for t in elemTypes:
            name, dim, order, num, par, _ = gmsh.model.mesh.getElementProperties(
                t)
            # print(" - Element type: " + name + ", order " + str(order))

        x = coords[0::3]
        y = coords[1::3]
        z = coords[2::3]

        pTags = gmsh.model.getPhysicalGroups()
        for pGroup in range(len(pTags)):
            name = gmsh.model.getPhysicalName(
                pTags[pGroup][0], pTags[pGroup][1])
            if name:
                nodeTags, coords = gmsh.model.mesh.getNodesForPhysicalGroup(
                    pTags[pGroup][0], pTags[pGroup][1])
                # print("Physical Group:", name, "| Contains",
                #   len(nodeTags), "nodes, with tags:")
                # print(nodeTags)
                with open("physicalGroups.txt", "a") as file:
                    file.write("{:s},{:d}".format(name, pGroup) + "\n")

        with open("nodeCoords.txt", "a") as file:
            file.write("{:d}\n".format(len(x)))
            for ti, xi, yi, zi in zip(tags, x, y, z):
                pGroup_ind = 0
                for pGroup in range(len(pTags)):
                    name = gmsh.model.getPhysicalName(
                        pTags[pGroup][0], pTags[pGroup][1])
                    if name:
                        nodeTags, coords = gmsh.model.mesh.getNodesForPhysicalGroup(
                            pTags[pGroup][0], pTags[pGroup][1])
                        if ti in nodeTags:
                            pGroup_ind = pGroup
                        # print("Physical Group:",name,"| Contains",len(nodeTags),"nodes, with tags:")
                        # print(nodeTags)
                file.write("{:d},{:.6f},{:.6f},{:.6f},{:d}\n".format(
                    ti, xi, yi, zi, pGroup_ind))

        numElem = sum(len(i) for i in elemTags)
        elemNodeTags = elemNodeTags[0].reshape(numElem, num)  # type: ignore

        with open("hexaCon.txt", "a") as file:
            file.write("{:d},{:d},{:d},{:d},{:d}\n".format(numElem, num, nx, ny, nz))  # type: ignore
            for i in range(numElem):
                file.write(','.join("{:d}".format(k)
                           for k in elemNodeTags[i]) + "\n")

gmsh.finalize()
