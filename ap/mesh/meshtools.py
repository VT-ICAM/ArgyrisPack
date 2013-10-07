#! /usr/bin/env python
"""Various algorithms for meshes that don't fit anywhere else."""
import itertools as it
import numpy as np
import ap.mesh.parsers as parsers


def sorted_edges(element, argyris=True):
    """
    Return the edges (as sorted (start, end) tuples) of a given element.

    Required Arguments
    ------------------
    * element : array-like container of the nodes comprising a finite element.
                Assumed to be in GMSH or Argyris order.

    Optional Arguments
    ------------------
    * argyris : Set to True (default) to treat the input as an Argyris
                element.
    """
    if argyris:
        return [tuple(sorted((element[i], element[j])))
                for (i, j) in [(0, 1), (0, 2), (1, 2)]]
    else:
        raise NotImplementedError


def extract_boundary_edges(elements):
    """
    Some meshes do not specify edge information. Extract it from the
    connectivity matrix.

    Required Arguments
    ------------------
    * elements : element connectivity matrix. Assumes Lagrange elements
                 and GMSH-style ordering.

    Output
    ------
    The output is a list of edge tuples. For example:

        print(extract_boundary_edges(mesh.elements))

        => [(1, 5, 7, -1), (2, 10, 12, -1), (3, 15, 17, -1),
            (4, 20, 22, -1), (5, 6, 8, -1), (6, 2, 9, -1),
            (10, 11, 13, -1), (11, 3, 14, -1), (15, 16, 18, -1),
            (16, 4, 19, -1), (20, 21, 23, -1), (21, 1, 24, -1)]

        print(mesh.edges)

        => [(1, 5, 7, 1), (5, 6, 8, 1), (6, 2, 9, 1), (2, 10, 12, 2),
            (10, 11, 13, 2), (11, 3, 14, 2), (3, 15, 17, 3),
            (15, 16, 18, 3), (16, 4, 19, 3), (4, 20, 22, 4),
            (20, 21, 23, 4), (21, 1, 24, 4)]

    where the last index, rather than being an edge label (based on
    GMSHs line attribute), is -1.
    """
    # Keep track of the original edges. At the end return the non-sorted edges.
    original_order = dict()
    sorted_edges = set()

    side_nodes = int((len(elements[0])-3)/3)
    for element in elements:
        # Guarantee uniqueness of edges by sorting the nodes. At the end map
        # the sorted versions back to the original versions.
        local_edges = [(element[i], ) + (element[j], ) +
                       tuple(element[3+i*side_nodes:3+(i+1)*side_nodes])
                       + (-1, ) for (i, j) in [(0, 1), (1, 2), (2, 0)]]

        local_sorted_edges = [tuple(sorted(t)) for t in local_edges]
        original_order.update(zip(local_sorted_edges, local_edges))

        for edge in local_sorted_edges:
            if edge in sorted_edges:
                sorted_edges.remove(edge)
            else:
                sorted_edges.add(edge)

    return [original_order[t] for t in sorted_edges]


def project_nodes(projection, elements, original_nodes,
                  attempt_flatten=False):
    """
    Given a projection and components of a finite element mesh, project
    the nodes with the supplied function. For quadratics, recalculate
    the location of the midpoint nodes afterwards. Assumes GMSH ordering
    of nodes for quadratics.

    Required Arguments
    ------------------
    * elements : element connectivity matrix. Assumes Lagrange elements
                 and GMSH-style ordering.
    * original_nodes : nodal coordinates corresponding to elements.

    Optional Arguments
    ------------------
    * attempt_flatten : try to (instead of applying the projection) drop
                        the last dimension.

    Output
    ------
    A numpy array of the projected nodes.
    """
    if attempt_flatten:
        if np.all(original_nodes[:, -1] == original_nodes[0, -1]):
            nodes = np.ascontiguousarray(original_nodes[:, 0:-1])
            return nodes

    nodes = np.array([projection(node) for node in original_nodes])

    # Do nothing for linears: there are no midpoints to fix.
    if elements.shape[1] == 3:
        pass
    # fix quadratics.
    elif elements.shape[1] == 6:
        for i in range(2):
            for (j, k) in [(0, 1), (1, 2), (2, 0)]:
                nodes[elements[:, j + 3] - 1, i] = \
                    0.5*(nodes[elements[:, j] - 1, i] +
                         nodes[elements[:, k] - 1, i])
    return nodes


def change_order(mesh, order):
    """
    Change the order of the elements in a mesh.
    """
    # upgrade linears to quadratics.
    reference_edges = [(0, 1), (1, 2), (2, 0)]
    if mesh.elements.shape[1] == 3 and order == 2:
        edges = {tuple(sorted((element[i], element[j])))
                 for (i, j) in reference_edges for element in mesh.elements}
        new_elements     = np.ones((mesh.elements.shape[0], 6), dtype=np.int)
        edge_to_midpoint = dict(zip(edges,
                                    it.count(start=len(mesh.nodes) + 1)))
        for element_number, element in enumerate(mesh.elements):
            local_edges = [tuple(sorted((element[i], element[j])))
                           for (i, j) in reference_edges]
            for i, edge in enumerate(local_edges):
                new_elements[element_number, i + 3] = edge_to_midpoint[edge]

        new_elements[:, 0:3] = mesh.elements
        new_nodes = np.ones((new_elements.max(), mesh.nodes.shape[1]))*np.nan
        new_nodes[0:mesh.nodes.shape[0]] = mesh.nodes
        for (i, j) in reference_edges:
            new_nodes[np.unique(new_elements[:, 3 + i]) - 1] = 0.5*(
                mesh.nodes[new_elements[:, i] - 1, :]
                + mesh.nodes[new_elements[:, j] - 1, :])
        new_edges = []
        for edge in mesh.edges:
            midpoint = edge_to_midpoint[tuple(sorted(edge[0:2]))]
            new_edge = tuple([edge[0], edge[1], midpoint, edge[2]])
            new_edges.append(new_edge)
    else:
        raise NotImplementedError("Unsupported mesh order conversion")
    return parsers.ParseArrays(new_elements, new_nodes, edges=new_edges)


def organize_edges(edges, borders=None, default_border='land'):
    """
    Organize edges in to various collections specified by borders.

    Required Arguments
    ------------------
    * edges          : List of edges of the form

        (node_number, node_number, ... , edge_type)

                       The sorting process requires edge_type to be the last
                       entry in each tuple.

    Optional Arguments
    ------------------
    * borders        : a dictionary correlating border numbers with names. For
                       example,

          borders = {'ocean_in' : 1, 'ocean_out' : 2}

    * default_border : the border to place all other edges in. Defaults to
                      'land'.
    """
    if borders is None:
        borders = dict()
    if default_border in borders:
        raise ValueError("Specific border and default border share same name")

    edge_collections = dict()
    edge_collections[default_border] = set(edges)
    for border, labels in borders.items():
        edge_collections[border] = {t for t in edge_collections[default_border]
                                    if t[-1] in labels}
        edge_collections[default_border] -= edge_collections[border]
    return edge_collections
