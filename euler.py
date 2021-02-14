import graph
import sys

log = False

def euler_assembly(kmers):
    graph = kmer_graph(kmers)
    path = eulerian_path(graph)
    assembly = assembly_from_path(graph, path)
    return assembly

class SpectralAssemblyGraph(graph.OutEdgePoppingDirectedGraph, 
                            graph.VertexLabeledDirectedGraph):
    pass

def kmer_graph(kmers):
    # determine the (k-1)-mers, which will be the vertices
    k_minus_1_mers = set()
    for kmer in kmers:
        k_minus_1_mers.update([kmer[:-1], kmer[1:]])
    ordered_k_minus_1_mers = sorted(k_minus_1_mers)
    
    # create a graph with (k-1)-mers as vertices and k-mers as edges
    g = SpectralAssemblyGraph(len(ordered_k_minus_1_mers))
    # label the vertices
    for i, k_minus_1_mer in enumerate(ordered_k_minus_1_mers):
        g.set_vertex_label(i, k_minus_1_mer)
    # add the edges
    for kmer in kmers:
        g.add_edge(g.vertex_index(kmer[:-1]), g.vertex_index(kmer[1:]))
    
    return g

def eulerian_path(g):
    # find unbalanced vertices
    fake_edge = unbalanced_vertices(g)

    # create a fake edge if necessary to guarantee the existence of an Eulerian cycle
    if fake_edge:
        if log:
            print("Adding fake edge:", fake_edge, file=sys.stderr)
            if g.has_edge(*fake_edge):
                print("Fake edge is a duplicate of a real edge", file=sys.stderr)
        g.add_edge(*fake_edge)
    else:
        if log:
            print("Graph is already balanced", file=sys.stderr)

    # find an Eulerian cycle in the graph
    path = eulerian_cycle(g)

    # if we added a "fake" edge, remove it
    if fake_edge:
        for i in range(0, len(path) - 1):
            if path[i: i + 2] == fake_edge:
                path = path[i + 1:] + path[1: i + 1]
                break

    return path

def unbalanced_vertices(g):
    """Returns the indices of the pair of vertices that are unbalanced.
    
    The first vertex of the pair will have one less outgoing than incoming
    edge and the second vertex of the pair has one less incoming than
    outgoing edge. If the graph is already balancd, None is returned."""
    unbalanced_vertices = [i for i in range(g.num_vertices())
                           if g.outdegree(i) != g.indegree(i)]
    if len(unbalanced_vertices) == 0:
        return None
    elif len(unbalanced_vertices) == 2:
        degree_diffs = [g.outdegree(i) - g.indegree(i) for i in unbalanced_vertices]
        if degree_diffs == [-1, 1]:
            return unbalanced_vertices
        elif degree_diffs == [1, -1]:
            return list(reversed(unbalanced_vertices))
        else:
            raise Exception("Unbalanced vertices do not have compatible degree differences")
    else:
        raise Exception("More than two unbalanced vertices")    

def eulerian_cycle(g):
    # Find an initial cycle starting with the first vertex
    cycle = single_cycle(g, 0)

   # While there are still unused edges, find subcycles and merge into main cycle
    while g.num_edges() > 0:
        if log:
            print("Current cycle:", cycle, file=sys.stderr)
        
        # Find index of first vertex in cycle that has outgoing edges
        for i in range(len(cycle)):
            if g.outdegree(cycle[i]) > 0:
                break        
        if log:
            print("Finding a cycle from position", i, file=sys.stderr)
        next_cycle = single_cycle(g, cycle[i])
        cycle[i:i] = next_cycle[:-1]

    return cycle

def single_cycle(g, start_vertex):
    cycle = [start_vertex]
    while True:
        curr_vertex, next_vertex = g.pop_edge(cycle[-1])
        cycle.append(next_vertex)
        if next_vertex == start_vertex:
            return cycle

def assembly_from_path(g, path):
    if log:
        print(path, file=sys.stderr)
    vertex_last_bases = [g.vertex_label(i)[-1] for i in path]
    return g.vertex_label(path[0])[:-1] + ''.join(vertex_last_bases)