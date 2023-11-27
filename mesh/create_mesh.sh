#!/bin/bash
gmsh -2 mesh.geo
ElmerGrid 14 2 mesh.msh -autoclean -order 3 2 1
mv mesh/mesh.nodes mesh/mesh.nodes.orig
awk -f meshdeform.awk mesh/mesh.nodes.orig > mesh/mesh.nodes
ElmerGrid 2 5 mesh 

