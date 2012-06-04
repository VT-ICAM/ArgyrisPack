#! /usr/bin/env python
'''
File: Mesh.py
Author: Dave Wells
'''
import numpy as np

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

    * original_elements : integer array of the global node numbers of each
      element in the quadratic mesh.
    """
    def __init__(self, node_collections, original_elements):
        self.node_collections = node_collections
        self.elements = np.zeros((original_elements.shape[0], 21), dtype=np.int)
        self.elements[:,0:6] = original_elements

        self._new_index = original_elements.max() + 1

        # examine each element and add on appropriate Argyris nodes.
        for element in self.elements:
            for corner_number in range(0,3):
                self._test_and_update(element, corner_number)

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
        np.savetxt(self.name + 'dx.txt', np.array(list(self.dx)), fmt="%d")
        np.savetxt(self.name + 'dy.txt', np.array(list(self.dy)), fmt="%d")
        np.savetxt(self.name + 'normal.txt',
                   np.array(list(self.normal_derivatives)), fmt="%d")
        np.savetxt(self.name + 'function.txt',
                   np.array(list(self.function_values)), fmt="%d")
        np.savetxt(self.name + 'all.txt',
                   np.sort(np.hstack(self.stacked_nodes.values() +
                                     self.stacked_nodes.keys() +
                                     list(self.normal_derivatives))), "%d")
