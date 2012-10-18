#! /usr/bin/env python

import numpy as np

class Mesh(object):
    """
    Representation of a finite element mesh. Essentially a struct.

    Required Arguments and Properties
    ---------------------------------
    * elements : integer numpy array of the global node numbers of each
      element.
    * nodes : double precision numpy array of global node coordinates.
    * node_collections : a dictionary correlating names with node numbers.

    Properties
    ----------
    * interior_nodes : an array of the node numbers correlated with the
      'interior' node collection.
    """
    def __init__(self, elements, nodes, node_collections):
        self.elements = elements
        self.nodes = nodes
        self.node_collections = node_collections

        # for convenience in enforcing homogeneous boundary conditions label
        # the interior nodes.
        self.interior_nodes = np.fromiter(
            node_collections['interior'].__iter__(),
            count = len(node_collections['interior']),
            dtype=np.int)

    def estimate_nnz(self):
        """
        Estimate the number of nonzero entries present in some IJV-format
        sparse matrix constructed from inner products on this collection of
        elements.
        """
        return self.elements.shape[1]**2 * self.elements.shape[0]

class ArgyrisMesh(object):
    """
    Class to build an Argyris mesh from a quadratic mesh. Can handle a mesh
    with multiple boundary conditions.

    The algorithm is as follows:

    1. Treat the current midpoint nodes as the normal derivative basis
    functions.
    2. Extract the corner nodes as a separate array. Associate each corner node
    with five new nodes stacked at the same location.
    3. Update nodal coordinates and fix the element order.

    Required Arguments
    ------------------
    * mesh : a parsed .mesh file (has elements, nodes, and edge_collections)

    Properties
    ----------
    * elements : a numpy array listing the node numbers of every element.
    * edges_by_midpoint : a dictionary associating each element with a certain
      edge (indexed by the normal derivative basis function number)
    * node_collections : a list of ArgyrisNodeCollection objects.
    * nodes : a numpy array of node coordinates.

    Methods
    -------
    * save_files : save the mesh in a format compatible to the existing QGE
    code.
    """
    def __init__(self, mesh):
        self.elements = np.zeros((mesh.elements.shape[0], 21), dtype=np.int)
        self.elements[:,0:6] = mesh.elements

        # solve a lot of orientation problems later by ensuring that the corner
        # nodes are in sorted order.
        for element in self.elements:
            self._fix_element_order(element[0:6])

        # stack new nodes.
        self.stacked_nodes = \
            {node_number : np.array(range(mesh.elements.max() + 1 + 5*count,
                                          mesh.elements.max() + 1 + 5*count + 5))
             for count, node_number in enumerate(np.unique(self.elements[:,0:3]))}

        for element in self.elements:
            element[6:11]  = self.stacked_nodes[element[0]]
            element[11:16] = self.stacked_nodes[element[1]]
            element[16:21] = self.stacked_nodes[element[2]]

        # update the edges by elements.
        self.edges_by_midpoint = \
            {midpoint : (element_number + 1, k) for k in range(1,4) for
             element_number, midpoint in enumerate(self.elements[:,2+k])}

        self._fix_argyris_node_order()

        # set coordinates for the new nodes.
        self.nodes = np.zeros((self.elements.max(), 2))
        self.nodes[0:len(mesh.nodes),:] = mesh.nodes
        for stacked_node, new_nodes in self.stacked_nodes.items():
            self.nodes[new_nodes - 1] = mesh.nodes[stacked_node - 1]

        # Construct the edge collections.
        self._build_node_collections(mesh)

    def save_files(self):
        """
        Save the following data for compatibility with the QGE code:

            nodes.txt    : all nodal coordinates
            elements.txt : the element array for Argyris

        and for each collection of nodes with key NAME:

            NAME_edge_elements.txt : all edge tuples (end, end, midpoint)
            NAME_all.txt : all numbers of nodes in the collection.
        """
        np.savetxt('nodes.txt', self.nodes)
        np.savetxt('elements.txt', self.elements, fmt="%d")

        for collection in self.node_collections:
            collection.save_files()

    def _fix_element_order(self, element):
        """
        Ensure that the corners of the input quadratic element are in increasing
        order. For example, convert
        1 3 2 4 6 5
        to
        1 2 3 5 6 4
        """
        if element[0] > element[1]:
            element[0], element[1] = element[1], element[0]
            element[4], element[5] = element[5], element[4]
        if element[1] > element[2]:
            element[2], element[1] = element[1], element[2]
            element[3], element[5] = element[5], element[3]
        if element[0] > element[2]:
            element[2], element[0] = element[0], element[2]
            element[3], element[4] = element[4], element[3]
        if element[0] > element[1]:
            element[0], element[1] = element[1], element[0]
            element[4], element[5] = element[5], element[4]

    def _build_node_collections(self, mesh):
        """
        Handle the edges by building a list of ArgyrisNodeCollection objects.
        This is done by extracting the information regarding corner nodes and
        midpoints from the original edge data and saving the interior nodes as
        everything that was not a boundary node.
        """
        self.node_collections = []
        interior_function_values = set(mesh.elements[:, 0:3].flatten())
        interior_normal_derivatives = set(mesh.elements[:, 3:6].flatten())

        for border_name, collection in mesh.edge_collections.items():
            # save left points of edges.
            function_values = {x[0] for x in collection}

            normal_derivatives = {x[2] for x in collection}
            edges = {tuple(x[0:3]) for x in collection}

            self.node_collections.append(ArgyrisNodeCollection(function_values,
                normal_derivatives, edges, self, name = border_name))

            interior_function_values.difference_update(function_values)
            interior_normal_derivatives.difference_update(normal_derivatives)

        self.node_collections.append(ArgyrisNodeCollection(
            interior_function_values, interior_normal_derivatives, [], self,
            name = 'interior'))

    def _fix_argyris_node_order(self):
        """
        Fix the node orderings from the constructed format

            [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21]

        to the usual Argyris format of

            [1 2 3 7 8 12 13 17 18 9 10 11 14 15 16 19 20 21 4 6 5].
        """

        normal_derivatives1 = self.elements[:,3].copy()
        normal_derivatives2 = self.elements[:,5].copy()
        normal_derivatives3 = self.elements[:,4].copy()

        first_nodes  = self.elements[:,6:11].copy()
        second_nodes = self.elements[:,11:16].copy()
        third_nodes  = self.elements[:,16:21].copy()

        self.elements[:,18]    = normal_derivatives1
        self.elements[:,19]    = normal_derivatives2
        self.elements[:,20]    = normal_derivatives3

        self.elements[:,3:5]   = first_nodes[:,0:2]
        self.elements[:,9:12]  = first_nodes[:,2:5]

        self.elements[:,5:7]   = second_nodes[:,0:2]
        self.elements[:,12:15] = second_nodes[:,2:5]

        self.elements[:,7:9]   = third_nodes[:,0:2]
        self.elements[:,15:18] = third_nodes[:,2:5]

class ArgyrisNodeCollection(object):
    """
    Contains information about a group of nodes in an Argyris Mesh and any
    relevant edge data.

    Required Arguments
    ------------------
    * function_values : set of basis function numbers that approximate function
      values on the Argyris mesh.
    * normal_derivatives : set of the node numbers corresponding to normal
      derivative basis functions.
    * edges : set of tuples corresponding to (endpoint, endpoint, midpoint)
    * mesh : the relevant Argyris mesh.

    Optional Arguments
    ------------------
    * name : prefix on the output files. Defaults to 'inner'.
    """
    def __init__(self, function_values, normal_derivatives,
                 edges, mesh, name = 'inner'):
        self.function_values = function_values
        self.normal_derivatives = normal_derivatives
        self.name = name

        self.stacked_nodes = {node : mesh.stacked_nodes[node] for node in
                              self.function_values}
        self._edge_elements = [mesh.edges_by_midpoint[edge[-1]] + edge
                               for edge in edges]

    def save_files(self):
        """
        Save the data to text files; place all node numbers in the collection
        in one file and all information on edge elements in another.
        """
        if self._edge_elements: # don't save if there are no edge elements.
            np.savetxt(self.name + '_edge_elements.txt',
                       np.asarray(self._edge_elements, dtype=np.int), "%d")

        # Use list comprehensions because they do the same thing in python2.7
        # and python3.*; *.values() became an iterator in python3000.
        np.savetxt(self.name + '_all.txt',
                   np.unique(np.hstack([x for x in self.stacked_nodes.values()] +
                                       [x for x in self.stacked_nodes.keys()] +
                                       [x for x in self.normal_derivatives])),
                   "%d")

    def __str__(self):
        """For interactive debugging use."""
        return ("Node collection name: " + self.name + "\n" +
        "function values:\n" + str(self.function_values) + "\n" +
        "normal derivatives:\n" + str(self.normal_derivatives) + "\n" +
        "edge elements:\n" + str(self._edge_elements))
