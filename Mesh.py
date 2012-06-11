#! /usr/bin/env python

import numpy as np

class Mesh(object):
    """
    Representation of a finite element mesh. contains three pieces:

    elements: an array of integers, where each row contains node numbers of
    some element.

    nodes: array of global coordinates of nodes.

    node_collections: a dictionary of names ('inner', 'ocean', etc) paired
    to sets of node numbers.
    """
    def __init__(self, elements, nodes, node_collections):
        self.elements = elements
        self.nodes = nodes
        self.node_collections = node_collections

class ArgyrisMesh(object):
    """
    Class to build an Argyris mesh from a quadratic mesh. Can handle a mesh
    with multiple boundary conditions.

    The algorithm is as follows:

    1. Treat the current midpoint nodes as the normal derivative basis
    functions.
    2. For each corner of each element, see if nodes have been stacked on
    previously. If not, use the class-scope variable _new_index to add five
    new nodes at the current corner. Update the appropriate node container.

    Required Arguements
    -------------------
    * node_collections : list of the ArgyrisNodeCollection objects formed
      from the quadratic mesh.

    * original_elements : integer numpy array of the global node numbers of
      each element in the quadratic mesh.

    * original_nodes : 2xN numpy array of node coordinates on the quadratic
      mesh.

    Properties
    ----------
    * elements : a numpy array listing the node numbers of every element.

    * node_collections : a list of ArgyrisNodeCollection objects.

    * nodes : a numpy array of node coordinates.

    Methods
    -------
    * save_QGE_files : save the mesh in a format compatible to the existing QGE
    code.
    """
    def __init__(self, node_collections, original_elements, original_nodes):
        self.node_collections = node_collections
        self.elements = np.zeros((original_elements.shape[0], 21), dtype=np.int)
        self.elements[:,0:6] = original_elements
        self._new_index = original_elements.max() + 1

        # solve a lot of orientation problems later by ensuring that the corner
        # nodes are in sorted order.
        for element in self.elements:
            self._fix_element_order(element[0:6])

        # examine each element and add on appropriate Argyris nodes.
        for element in self.elements:
            for corner_number in range(0,3):
                self._test_and_update(element, corner_number)

        # construct the nodal coordinates array.
        self.nodes = np.zeros((self.elements.max(),2))
        self.nodes[0:original_nodes.shape[0], :] = original_nodes
        for collection in self.node_collections:
            for pair in collection.stacked_nodes.iteritems():
                for new_node in pair[1]:
                    self.nodes[new_node - 1, :] = original_nodes[pair[0] - 1, :]

    def _test_and_update(self, element, corner_number):
        """
        Update a corner node to have 5 more nodes stacked on it. Check if this
        node has been previously edited; if so, load the new nodes and node
        numbers. Otherwise increment the global counter.
        """
        for collection in self.node_collections:
            if element[corner_number] in collection.function_values:
                first = 6 + corner_number*5
                last  = 11 + corner_number*5
                try:
                    # if we already added new nodes at this corner, use those.
                    element[first:last] = \
                        collection.stacked_nodes[element[corner_number]]
                except KeyError:
                    # otherwise save the coordinate and numbers.
                    collection.update(self._new_index,
                                           element[corner_number])
                    element[first:last] = \
                        collection.stacked_nodes[element[corner_number]]
                    self._new_index += 5
                break

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

    def save_QGE_files(self):
        """
        Save the following data for compatibility with the QGE code:

            nodes.txt    : all nodal coordinates
            elements.txt : the element array for Argyris
            unodes.txt   : nodes corresponding to function values

        and for each border in edge_collections with key NAME, as well as the
        interior nodes:

            NAME_dx.txt       : nodes approximating x-derivatives
            NAME_dy.txt       : nodes approximating y-derivatives
            NAME_normal.txt   : nodes approximating normal derivatives
            NAME_function.txt : nodes approximating function values
            NAME_all.txt      : all nodes in the collection.
        """
        # save node indicies containing function values.
        u_nodes = np.unique(self.elements[:,0:3].flatten())

        # fix the numbering on the argyris.elements to match QGE code.
        elements = self.elements.copy()
        normal_derivatives1 = elements[:,3].copy()
        normal_derivatives2 = elements[:,4].copy()
        normal_derivatives3 = elements[:,5].copy()

        first_nodes  = elements[:,6:11].copy()
        second_nodes = elements[:,11:16].copy()
        third_nodes  = elements[:,16:21].copy()

        elements[:,18]    = normal_derivatives1
        elements[:,19]    = normal_derivatives3
        elements[:,20]    = normal_derivatives2

        elements[:,3:5]   = first_nodes[:,0:2]
        elements[:,9:12]  = first_nodes[:,2:5]

        elements[:,5:7]   = second_nodes[:,0:2]
        elements[:,12:15] = second_nodes[:,2:5]

        elements[:,7:9]   = third_nodes[:,0:2]
        elements[:,15:18] = third_nodes[:,2:5]

        np.savetxt('nodes.txt', self.nodes)
        np.savetxt('elements.txt', elements, fmt="%d")
        np.savetxt('unodes.txt', u_nodes, fmt="%d")

        # save the information stored in the node collections as well.
        for collection in self.node_collections:
            collection.write_to_files()

class ArgyrisNodeCollection(object):
    """
    Contains information about a group of nodes in an Argyris Mesh.
    """
    def __init__(self, function_values, normal_derivatives, name = 'inner'):
        self.function_values = function_values
        self.normal_derivatives = normal_derivatives
        self.name = name

        self.dx = set()
        self.dy = set()
        self.stacked_nodes = dict()

    def update(self, index, original_node):
        if original_node not in self.stacked_nodes.keys():
            self.dx.add(index)
            self.dy.add(index + 1)
            self.stacked_nodes[original_node] = np.array(
                range(index, index + 5), dtype=np.int)
        else:
            raise KeyError("Attempted to update a node a second time.")

    def write_to_files(self):
        """
        Write the node numbers to several text files based on type
        (normal derivatives, dxs, dys, function values)
        """
        np.savetxt(self.name + '_dx.txt', np.array(list(self.dx)), fmt="%d")
        np.savetxt(self.name + '_dy.txt', np.array(list(self.dy)), fmt="%d")
        np.savetxt(self.name + '_normal.txt',
                   np.array(list(self.normal_derivatives)), fmt="%d")
        np.savetxt(self.name + '_function.txt',
                   np.array(list(self.function_values)), fmt="%d")
        np.savetxt(self.name + '_all.txt',
                   np.sort(np.hstack(self.stacked_nodes.values() +
                                     self.stacked_nodes.keys() +
                                     list(self.normal_derivatives))), "%d")

    def __str__(self):
        return ("Node collection name: " + self.name + "\n" +
        "function values:\n" + str(self.function_values) + "\n" +
        "normal derivatives:\n" + str(self.normal_derivatives) + "\n" +
        "x derivatives:\n" + str(self.dx) + "\n" +
        "y derivatives:\n" + str(self.dy) + "\n")
