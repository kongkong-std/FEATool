# FEATool

A finite element analysis tool developed by *C programming language*, currently a *serial version*. The program contains 8 main parts:<br>
1. **Mesh process**
2. **Element type**
3. **Weak form**
4. **Element sftiffness matrix and load vector**
5. **Global stiffness matrix and load vector**
6. **Boundary conditions process**
7. **Linear system assembling**
8. **Solver**

## Mesh process
The mesh data is used as an input file. This project uses unstructured meshes to store mesh information, including mesh node information (which contains node numbers and node coordinates) and mesh element information (which includes element numbers and the nodes contained within each element, thus describing the topological relationships between the mesh nodes).

>2D mesh, triangular element
>> 3 nodes in each element

>3D mesh, tetrahedral element
>> 4 nodes in each element

>node data
>> 3D coordinates, (x, y, z)