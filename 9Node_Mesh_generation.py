import numpy as np
import matplotlib.pyplot as plt
import math as math
import csv
from copy import deepcopy

# load the nodes and element data of 8-noded element.

node_file = 'nodes.dat'
element_file = 'elements.dat'
nodes = np.loadtxt(node_file,delimiter=",")
elements = np.loadtxt(element_file,delimiter=",")

# find the shape of nodes and elements matrix.

nodes_shape = np.shape(nodes)
elements_shape = np.shape(elements)
a1 = nodes_shape[0]
a2 = nodes_shape[1]
b1 = elements_shape[0]
b2 = elements_shape[1]

# initialize the shape function array and assign the values to the variables.

shape_function = np.zeros(8)
xi = 0
omega = 0


# Assign the values to the shape function array.

shape_function[0] = (1/4)*(1-xi)*(1-omega)*(-xi-omega-1)
shape_function[1] = (1/4)*(1+xi)*(1-omega)*(+xi-omega-1)
shape_function[2] = (1/4)*(1+xi)*(1+omega)*(+xi+omega-1)
shape_function[3] = (1/4)*(1-xi)*(1+omega)*(-xi+omega-1)
shape_function[4] = (1/2)*(1-xi**2)*(1-omega)
shape_function[5] = (1/2)*(1+xi)*(1-omega**2)
shape_function[6] = (1/2)*(1-xi**2)*(1+omega)
shape_function[7] = (1/2)*(1-xi)*(1-omega**2)

# initilize the new nodes and elements matrix.

New_nodes = np.zeros((a1+b1,a2))
New_elements = np.zeros((b1,b2+1))

# assign the same nodes matrix values to the new node matrix.

for i in range(a1):
    for j in range(a2):
        New_nodes[i,j] = nodes[i,j]

# assign the same element matrix values to the new element matrix.

for i in range(b1):
    for j in range(b2):
        New_elements[i,j] = elements[i,j]      

# loop to assign the addition nodes to the new nodes matrix.
# loop to add the nodes number to the new element matrix.

for i in range(a1,a1+b1,1):
    New_nodes[i,:] = 0
    New_nodes[i,0] = i+1
    element_number = i - a1 + 1

# loop to find coordinate of the new ninth node of the element. 

    for j in range(1,8+1,1):
        node_number = elements[element_number-1,j]
        node_number = np.int(node_number)
        New_nodes[i,1] = New_nodes[i,1] + shape_function[j-1]*nodes[node_number-1,1]        
        New_nodes[i,2] = New_nodes[i,2] + shape_function[j-1]*nodes[node_number-1,2]
    New_elements[element_number-1,9] = i+1
    
# convert new element matrix into integer matrix from folat matrix.

New_elements = New_elements.astype(int)  
with open('New_nodes.dat', 'w') as file:
    for i in range(New_nodes.shape[0]):

# write the matrix column vice to convert the first column of new node matrix.
# into the integer since it contains the node number.

        file.write(str(int(New_nodes[i, 0])) + ',' + str(New_nodes[i, 1]) + ',' + str(New_nodes[i, 2]) + '\n')
        
with open('New_elements.dat', 'w', newline='') as file1:
    writer = csv.writer(file1)
    writer.writerows(New_elements)    
