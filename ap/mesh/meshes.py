#! /usr/bin/env python
import random
import numpy as np
from collections import namedtuple
import ap.mesh.meshtools as meshtools
import ap.mesh.parsers as parsers

def mesh_factory(*args, **kwargs):
    """
    Parse a finite element mesh representation and then convert it to a
    Mesh or ArgyrisMesh object.

    Required Arguments
    ------------------
    * mesh_files      : text files comprising the finite element mesh.

    Keyword Arguments
    -----------------
    * argyris        : boolean to specify if the mesh should have
                       additional nodes added to transform it in to an
                       Argyris mesh. Defaults to False.

    * order          : Mesh element order. The nodes will be renumbered
                       appropriately (1 for linears, 2 for quadratics).
                       Defaults to None. This is not implemented yet.

    * projection     : function that projects nodes. Defaults to None
                       (no projection)

    * borders        : a dictionary correlating names with GMSH 'Physical
      Line' attributes. For example,

          borders = {'open' : (1,2)}

      will correlate edges on Physical Lines 1 and 2 with the 'open' edge
      collection.

    * default_border : the default edge collection for any edges that
                       are not in a special_border collection. Defaults
                       to 'land'.
    """
    parsed_mesh = parsers.parser_factory(*args)
    if 'argyris' in kwargs:
        keywords = kwargs.copy()
        del keywords['argyris']
        return ArgyrisMesh(parsed_mesh, **keywords)
    elif 'Argyris' in kwargs:
        keywords = kwargs.copy()
        del keywords['Argyris']
        return ArgyrisMesh(parsed_mesh, **keywords)
    else:
        return Mesh(parsed_mesh, **kwargs)

class Mesh(object):
    """
    Representation of a finite element mesh.

    Required Arguments
    ------------------
    * parsed_mesh : Something that has the same interface as a MeshParser
                    (has fields elements, nodes, and edge_collections)

    Optional Arguments
    ------------------
    * borders        : A dictionary correlating names with a tuple of
      GMSH physical line numbers. for example:

          borders = {'no_flow' : (1,2,3), 'coast' : (4,5,6)}

    * default_border : the name corresponding to the default edge
      collection. Defaults to "land".

    * projection     : function for transforming nodes (say from 3D to
      2D); for example,

          projection = lambda x : x[0:2]

      will project the nodes down to the XY plane. Defaults to the
      identity function.

    Properties
    ----------
    * elements         : element connectivity matrix.

    * nodes            : coordinates of nodes.

    * edge_collections : a dictionary relating the border names to the
      edge tuples that fall along that border. If possible, the last
      number in the tuple is the geometrical item number that the edge
      falls upon from GMSH. Otherwise it is -1. For example,

        print(t.edge_collections)
        => {'land': set([(3, 4, 7, 3), (4, 1, 8, 4), (2, 3, 6, 2),
            (1, 2, 5, 1)])}

    * boundary_nodes   : Set containing the node numbers of nodes on the
                         boundary.

    * interior_nodes   : Set containing the node numbers of nodes in the
                         interior.

    Methods
    -------
    * get_nnz() : Calculate the number of nonzero entries in a typical
                  finite element matrix (e.g. stiffness matrix) based on
                  the total number of inner products. This will be
                  exactly the value of the length of one of the
                  tripplet-form vectors.

    * savetxt(prefix="") : Save the mesh as text files
      1. prefix + nodes.txt
      2. prefix + elements.txt
      3. prefix + interior_nodes.txt
      4. prefix + boundary_nodes.txt

      and, additionally, for each edge collection save
      prefix + name + _edges.txt.
    """
    def __init__(self, parsed_mesh, borders=dict(), default_border="land",
                 projection=lambda x : x):
        self.elements = parsed_mesh.elements

        self.nodes = meshtools.project_nodes(projection, parsed_mesh.elements,
                                             parsed_mesh.nodes)

        self.edge_collections = \
            meshtools.organize_edges(parsed_mesh.edges, borders=borders,
                                     default_border=default_border)

        if max(map(len, self.edge_collections.values())) == 0:
            self.edge_collections = \
                {default_border :
                 set(meshtools.extract_boundary_edges(self.elements))}

        if len(np.unique(self.elements)) != self.nodes.shape[0]:
            self._fix_unused_nodes()

        self.boundary_nodes = {}
        interior_nodes = set(range(1, len(self.nodes)+1))
        for name, edge_collection in self.edge_collections.items():
            self.boundary_nodes[name] = \
                np.fromiter(set(node for edge in edge_collection
                                for node in edge[0:-1]), int)
            interior_nodes -= set(self.boundary_nodes[name])

        self.interior_nodes = np.fromiter(interior_nodes, int)

    def get_nnz(self):
        """
        Estimate the number of nonzero entries present in some IJV-format
        sparse matrix constructed from inner products on this collection of
        elements.
        """
        return self.elements.shape[1]**2 * self.elements.shape[0]

    def savetxt(self, prefix=""):
        if prefix:
            prefix += "_"
        np.savetxt(prefix + "nodes.txt", self.nodes)
        np.savetxt(prefix + "elements.txt", self.elements, fmt="%d")
        np.savetxt(prefix + "interior_nodes.txt", self.interior_nodes, fmt="%d")
        np.savetxt(prefix + "boundary_nodes.txt", self.boundary_nodes, fmt="%d")

        for name, collection in self.edge_collections.items():
            np.savetxt(prefix + name + '_edges.txt', [t for t in collection],
                       fmt='%d')

    def _fix_unused_nodes(self):
        """
        GMSH has a bug where it saves non-mesh nodes (that is, nodes that
        are not members of any element) to some files. Get around that
        issue by deleting the extra nodes and renumbering accordingly.
        """
        number_of_mesh_nodes = len(np.unique(self.elements))
        old_to_new = dict(zip(np.unique(self.elements),
                              range(1, number_of_mesh_nodes + 1)))

        new_to_old = {new_node : old_node
                      for (old_node, new_node) in old_to_new.items()}
        new_elements = np.array([[old_to_new[node] for node in element]
                                  for element in self.elements])
        new_nodes = np.array([self.nodes[new_to_old[new_node_number] - 1]
                              for new_node_number in new_to_old.keys()])

        try:
            edge_size = {3 : 2, 6 : 3, 10 : 4, 15 : 5}[self.elements.shape[1]]
        except KeyError:
            raise ValueError("Unsupported mesh type")
        new_edge_collections = dict()
        for key, collection in self.edge_collections.items():
            new_edge_collections[key] = set()
            for edge in collection:
                # geometrical information available
                if len(edge) == edge_size + 1:
                    new_edge_collections[key].add(tuple([old_to_new[node]
                        for node in edge[0:-1]] + [edge[-1]]))
                # geometrical information not available
                elif len(edge) == edge_size:
                    new_edge_collections[key].add(tuple([old_to_new[node]
                        for node in edge]))
                else:
                    raise ValueError("Mismatch between size of mesh and" +
                                     " size of edges")

        self.edge_collections = new_edge_collections
        self.elements = new_elements
        self.nodes = new_nodes

ArgyrisEdge = namedtuple('ArgyrisEdge', ['element_number', 'edge_type', 'edge'])

class ArgyrisMesh(object):
    """
    Class to build an Argyris mesh from a parsed mesh. Can handle a mesh
    with multiple boundary conditions.

    The algorithm is as follows:

    1. Treat the current midpoint nodes as the normal derivative basis
    functions.
    2. Extract the corner nodes as a separate array. Associate each corner node
    with five new nodes stacked at the same location.
    3. Update nodal coordinates and fix the element order.

    Required Arguments
    ------------------
    * mesh : a parsed mesh (inherits from the MeshParser class)

    Properties
    ----------
    * elements : a numpy array listing the node numbers of every element.
    * edges_by_midpoint : a dictionary associating each element with a certain
      edge (indexed by the normal derivative basis function number)
    * node_collections : a list of ArgyrisNodeCollection objects.
    * nodes : a numpy array of node coordinates.

    Methods
    -------
    * savetxt: save the mesh in multiple text files.
    """
    def __init__(self, parsed_mesh, borders=dict(), default_border="land",
                 projection=lambda x : x):
        if parsed_mesh.elements.shape[1] != 6:
            raise NotImplementedError("Support for changing mesh order is not "
                                      + "implemented.")
            # parsed_mesh = meshtools.change_order(parsed_mesh, 2)

        lagrange_mesh = Mesh(parsed_mesh, borders=borders,
                             default_border=default_border, projection=projection)

        # if not projected, try to flatten as a last resort.
        if (lagrange_mesh.nodes.shape[1] == 3 and
            np.all(lagrange_mesh.nodes[:,2] == lagrange_mesh.nodes[0,2])):
            lagrange_mesh.nodes = lagrange_mesh.nodes[:,0:2]

        if lagrange_mesh.nodes.shape[1] != 2:
            raise ValueError("Requires a 2D mesh; try a different projection.")

        self.elements = np.zeros((parsed_mesh.elements.shape[0],21), dtype=np.int)
        self.elements[:,0:6] = lagrange_mesh.elements

        # solve a lot of orientation problems later by ensuring that the corner
        # nodes are in sorted order.
        for element in self.elements:
            self._sort_corners_increasing(element[0:6])

        # stack the extra basis function nodes on the corners.
        max_lagrange_mesh = parsed_mesh.elements.max() + 1;
        self.stacked_nodes = \
            {node_number : np.arange(max_lagrange_mesh + 5*count,
                                     max_lagrange_mesh + 5*count + 5)
             for count, node_number in enumerate(np.unique(self.elements[:,0:3]))}

        for element in self.elements:
            element[6:11]  = self.stacked_nodes[element[0]]
            element[11:16] = self.stacked_nodes[element[1]]
            element[16:21] = self.stacked_nodes[element[2]]

        self._fix_argyris_node_order()

        # update the edges by elements.
        self.edges_by_midpoint = dict()
        edge_type_to_nodes = {1 : (0,1,18), 2 : (0,2,19), 3 : (1,2,20)}
        for element_number, element in enumerate(self.elements):
            for edge_type in range(1,4):
                (i,j,k) = edge_type_to_nodes[edge_type]
                edge = ArgyrisEdge(element_number = element_number + 1,
                                   edge_type = edge_type,
                                   edge = (element[i], element[j], element[k]))
                if element[17 + edge_type] in self.edges_by_midpoint:
                    if (self.edges_by_midpoint[element[17 + edge_type]].edge
                        != edge.edge):
                        raise ValueError("Mesh is not consistent")
                else:
                    self.edges_by_midpoint[element[17 + edge_type]] = edge

        # set coordinates for the new nodes.
        self.nodes = np.zeros((self.elements.max(), 2))
        self.nodes[0:lagrange_mesh.nodes.shape[0],:] = lagrange_mesh.nodes
        for stacked_node, new_nodes in self.stacked_nodes.items():
            self.nodes[new_nodes - 1] = self.nodes[stacked_node - 1]

        # Construct the edge collections.
        self._build_node_collections(lagrange_mesh)

    def savetxt(self, prefix=""):
        """
        Save the following text files:

            nodes.txt    : nodal coordinates
            elements.txt : element connectivity matrix

        and for each collection of nodes with key NAME:

            NAME_edges.txt : all edge tuples (end, end, midpoint)
            NAME_all.txt : all numbers of nodes in the collection.
        """
        if prefix:
            prefix += "_"
        np.savetxt(prefix + 'nodes.txt', self.nodes)
        np.savetxt(prefix + 'elements.txt', self.elements, fmt="%d")

        for collection in self.node_collections:
            prefix = prefix[0:-1]
            collection.savetxt(prefix)

    def _sort_corners_increasing(self, element):
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

    def _build_node_collections(self, lagrange_mesh):
        """
        Handle the edges by building a list of ArgyrisNodeCollection
        objects. This is done by extracting the information regarding
        corner nodes and midpoints from the lagrange edge data and saving
        the interior nodes as everything that was not a boundary node.
        """
        self.node_collections = []
        interior_function_values = set(lagrange_mesh.elements[:, 0:3].flatten())
        interior_normal_derivatives = set(lagrange_mesh.elements[:, 3:6].flatten())

        for border_name, collection in lagrange_mesh.edge_collections.items():
            # save left points of edges.
            function_values = {x[0] for x in collection}
            normal_derivatives = {x[2] for x in collection}

            self.node_collections.append(ArgyrisNodeCollection(function_values,
                normal_derivatives, collection, self, name = border_name))

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

        self.edges = [mesh.edges_by_midpoint[edge[-2]] for edge in edges]

    def savetxt(self, prefix = ""):
        """
        Save the data to text files; place all node numbers in the collection
        in one file and all information on edges in another.
        """
        if prefix:
            prefix += "_"
        if self.edges: # don't save if there are no edges.
            edge_array = np.array([[edge.element_number, edge.edge_type,
                                    edge.edge[0], edge.edge[1], edge.edge[2]]
                for edge in self.edges])
            np.savetxt(prefix + self.name + "_edges.txt", edge_array, "%d")

        # Use list comprehensions because they do the same thing in python2.7
        # and python3.*; *.values() became an iterator in python3000.
        np.savetxt(prefix + self.name + "_all.txt",
                   np.unique(np.hstack([x for x in self.stacked_nodes.values()] +
                                       [x for x in self.stacked_nodes.keys()] +
                                       [x for x in self.normal_derivatives])),
                   "%d")

    def __str__(self):
        """For interactive debugging use."""
        return ("Node collection name: " + self.name + "\n" +
        "function values:\n" + str(self.function_values) + "\n" +
        "normal derivatives:\n" + str(self.normal_derivatives) + "\n" +
        "edges:\n" + str(self._edge_elements))
