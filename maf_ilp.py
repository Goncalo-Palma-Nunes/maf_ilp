from docplex.mp.model import Model

class Tree:
    """A rooted phylogenetic binary tree
    """

    class Node:
        def __init__(self, value, parent=None, sibling=None, children=[]):
            if not isinstance(value, int):
                raise ValueError("Node value must be an integer.")
            if not (parent is None or isinstance(parent, Tree.Node)):
                raise ValueError("Parent must be a Tree.Node instance or None.")
            if not (sibling is None or isinstance(sibling, Tree.Node)):
                raise ValueError("Sibling must be a Tree.Node instance or None.")
            if not isinstance(children, list) or not all(isinstance(child, Tree.Node) for child in children) or len(children) > 2:
                raise ValueError("Children must be a list of Tree.Node instances with at most two children.")
            
            self.value = value
            self.parent = parent
            self.children = children
            self.sibling = sibling

        def val(self):
            return self.value

        def __str__(self):
            return f"Node({self.value})"

        def __eq__(self, other):
            return (self.value == other.value) 

        def __hash__(self):
            return hash(self.value)

        def is_sibling(self, other):
            if not isinstance(other, Tree.Node):
                raise ValueError("Other must be a Tree.Node instance.")
            return self.parent == other.parent and self.parent is not None and self != other

        def is_root(self):
            return self.parent is None

        def is_leaf(self):
            return len(self.children) == 0

        def _depth(self):
            d, node = 0, self
            while node.parent is not None:
                node = node.parent
                d += 1
            return d

        def lca(self, other):
            if not isinstance(other, Tree.Node):
                raise ValueError("Other must be a Tree.Node instance.")
            u, v = self, other
            du, dv = u._depth(), v._depth()
            edges = set()
            while du > dv:
                edges.add(Tree.Edge(u.parent, u))
                u = u.parent; du -= 1
            while dv > du:
                edges.add(Tree.Edge(v.parent, v))
                v = v.parent; dv -= 1
            while u is not v: # identity, NOT == (eq is value-based, lines 26-27)
                if u is None or v is None:
                    raise ValueError("Nodes are not in the same tree.")
                edges.add(Tree.Edge(u.parent, u))
                edges.add(Tree.Edge(v.parent, v))
                u, v = u.parent, v.parent
            return (u, edges)


        def compute_lca_path(self, u):
            path = self.lca(u)
            path_edges = path[1]
            return path_edges

        def compute_triple_lca_path(self, u, v):
            non_sibling = u if self.is_sibling(v) else v
            path = self.compute_lca_path(non_sibling)
            path.add(Tree.Edge(self.sibling.parent, self.sibling))
            return path


    class Edge:

        def __init__(self, node1, node2):
            if not isinstance(node1, Tree.Node) or not isinstance(node2, Tree.Node):
                raise ValueError("Nodes must be Tree.Node instances.")
            
            self.node1 = node1
            self.node2 = node2

        def to(self):
            return self.node2

        def from_(self):
            return self.node1

        def __hash__(self):
            return hash((self.node1.val(), self.node2.val()))

        def __str__(self):
            return f"({self.node1}, {self.node2})"

        def __eq__(self, other):
            if not isinstance(other, Tree.Edge):
                return NotImplemented
            return (self.node1 == other.node1 and self.node2 == other.node2)

        def __repr__(self):
            return f"({self.node1}, {self.node2})"

    def __init__(self, edges, nodes, root, leaves):
        self.edges = edges
        self.nodes = nodes
        self.root = root
        self.leaves = leaves

        self._wire_pointers()
        self.validateTree()

    def _wire_pointers(self):
        all_nodes = set(self.nodes) | self.get_nodes()
        for n in all_nodes:                      # reset stale state
            n.parent, n.children, n.sibling = None, [], None
        for edge in self.edges:
            child, par = edge.to(), edge.from_()
            child.parent = par
            par.children.append(child)
        for n in all_nodes:                      # link siblings (same parent only)
            if len(n.children) == 2:
                n.children[0].sibling = n.children[1]
                n.children[1].sibling = n.children[0]

    def get_edges(self):
        return self.edges

    def get_nodes(self):
        nodes = set()
        for edge in self.edges:
            nodes.update([edge.from_(), edge.to()])
        return nodes

    def get_edges_incident_to_node(self, node):
        if not isinstance(node, Tree.Node):
            raise ValueError("Node must be a Tree.Node instance.")

        return [edge for edge in self.edges if node == edge.to()]

    def get_edges_outgoing_from_node(self, node):
        if not isinstance(node, Tree.Node):
            raise ValueError("Node must be a Tree.Node instance.")

        return [edge for edge in self.edges if node == edge.from_()]

    def validateTree(self):
        if not all(isinstance(edge, Tree.Edge) for edge in self.edges):
            raise ValueError("All edges must be instances of Tree.Edge.")
        if not type(self.edges) in [list, tuple, set]:
            raise ValueError("Edges must be a list, tuple, or set of Tree.Edge instances.")
        if not isinstance(self.root, Tree.Node) or self.root is None or not self.root.is_root():
            raise ValueError("Root must be a Tree.Node instance and must be the root of the tree.")


        # Check if the root is in the nodes
        if self.root not in self.get_nodes():
            raise ValueError("Root must be one of the nodes in the tree.")

        # Check if all edges are valid
        for edge in self.edges:
            if edge.from_() not in self.get_nodes() or edge.to() not in self.get_nodes():
                raise ValueError(f"Edge {edge} contains nodes that are not in the tree.")

        # Check for cycles and connectivity
        visited = set()
        stack = [self.root]

        while stack:
            node = stack.pop()
            if node in visited:
                raise ValueError("The tree contains a cycle.")
            visited.add(node)
            for edge in self.get_edges_outgoing_from_node(node):
                stack.append(edge.to())

        if visited != self.get_nodes():
            raise ValueError("The tree is not connected.")

        # Check if the tree is binary (each node has at most two children)
        for node in self.get_nodes():
            if len(self.get_edges_outgoing_from_node(node)) > 2:
                raise ValueError(f"Node {node} has more than two children, so the tree is not binary.")

        # Check if the tree is rooted (the root has no incoming edges)
        if len(self.get_edges_incident_to_node(self.root)) > 0:
            raise ValueError("The root node has incoming edges, so the tree is not rooted.")

        # Check if all nodes except the root have exactly one incoming edge
        for node in self.get_nodes():
            if node != self.root and len(self.get_edges_incident_to_node(node)) != 1:
                raise ValueError(f"Node {node} does not have exactly one incoming edge, so the tree is not rooted.")

        for leaf_key in self.leaves:
            leaf = self.leaves[leaf_key]
            if not isinstance(self.leaves, dict):
                raise ValueError("Leaves must be a dictionary of Tree.Node instances.")
            if not isinstance(leaf, Tree.Node):
                raise ValueError("All leaves must be instances of Tree.Node.")
            if leaf not in self.get_nodes():
                raise ValueError(f"Leaf {leaf} is not in the tree.")
            if not leaf.is_leaf():
                raise ValueError(f"Leaf {leaf} is not a leaf node.")

        return True

    def __str__(self):
        return f"Tree with edges: {self.edges}"

#####################################################
#####            Auxiliary Functions            #####
#####################################################

def is_compatible(leaf_i, leaf_j, leaf_k, tree1, tree2):
    """
    Checks if the triple (leaf_i, leaf_j, leaf_k) is compatible between two trees.
    :param leaf_i: The first leaf label.
    :param leaf_j: The second leaf label.
    :param leaf_k: The third leaf label.
    :param tree1: The first tree.
    :param tree2: The second tree.
    :return: True if the sibling relationships are coherent with respect to the two triples in both trees, False otherwise.
    """
    leaf_i_t1 = tree1.leaves[leaf_i]
    leaf_j_t1 = tree1.leaves[leaf_j]
    leaf_k_t1 = tree1.leaves[leaf_k]

    leaf_i_t2 = tree2.leaves[leaf_i]
    leaf_j_t2 = tree2.leaves[leaf_j]
    leaf_k_t2 = tree2.leaves[leaf_k]
    return (leaf_i_t1.sibling == leaf_i_t2.sibling and leaf_i_t1.sibling in (leaf_j_t1, leaf_k_t1)) or \
           (leaf_j_t1.sibling == leaf_j_t2.sibling and leaf_j_t1.sibling in (leaf_i_t1, leaf_k_t1)) or \
           (leaf_k_t1.sibling == leaf_k_t2.sibling and leaf_k_t1.sibling in (leaf_i_t1, leaf_j_t1))


def compute_incompatible_triples(tree1, tree2):
    """
    Computes the incompatible triples between two trees.
    :param tree1: The first tree.
    :param tree2: The second tree.
    :return: A list of incompatible triples.
    """
    num_leaves = len(tree1.leaves)
    incompatible_triples = []
    for i in range(num_leaves):
        for j in range(i + 1, num_leaves):
            for k in range(j + 1, num_leaves):
                leaf_i = tree1.leaves[i]
                leaf_j = tree1.leaves[j]
                leaf_k = tree1.leaves[k]
                if not (leaf_i.parent == leaf_j.parent or
                        leaf_i.parent == leaf_k.parent or
                        leaf_j.parent == leaf_k.parent):
                    continue

                # Check if the triple (leaf_i, leaf_j, leaf_k) is incompatible between the two trees
                if not is_compatible(i, j, k, tree1, tree2):
                    incompatible_triples.append(leaf_i.compute_triple_lca_path(leaf_j, leaf_k))
                    
    return incompatible_triples

def compute_incompatible_paths(tree1, tree2):
    """
    Computes the incompatible paths between two trees.
    :param tree1: The first tree.
    :param tree2: The second tree.
    :return: A list of incompatible paths.
    """
    num_leaves = len(tree1.leaves)
    incompatible_paths = []
    for i in range(num_leaves):
        for j in range(i + 1, num_leaves):
            if i == j:
                continue
            leaf_i_t1 = tree1.leaves[i]
            leaf_j_t1 = tree1.leaves[j]
            leaf_j_t2 = tree2.leaves[j]
            leaf_i_t2 = tree2.leaves[i]

            for k in range(num_leaves):
                if k == i or k == j:
                    continue
                for l in range(num_leaves):
                    if l == i or l == j or l == k:
                        continue
                    leaf_k_t1 = tree1.leaves[k]
                    leaf_l_t1 = tree1.leaves[l]
                    leaf_k_t2 = tree2.leaves[k]
                    leaf_l_t2 = tree2.leaves[l]


                    # Check if the paths (leaf_i, leaf_j) and (leaf_k, leaf_l) are incompatible between the two trees
                    i_j_tree1_path = leaf_i_t1.compute_lca_path(leaf_j_t1)
                    k_l_tree1_path = leaf_k_t1.compute_lca_path(leaf_l_t1)
                    i_j_tree2_path = leaf_i_t2.compute_lca_path(leaf_j_t2)
                    k_l_tree2_path = leaf_k_t2.compute_lca_path(leaf_l_t2)

                    t1_intersection = i_j_tree1_path.intersection(k_l_tree1_path)
                    t2_intersection = i_j_tree2_path.intersection(k_l_tree2_path)
                    if t2_intersection and not t1_intersection:
                        incompatible_paths.append(i_j_tree1_path.union(k_l_tree1_path))

    return incompatible_paths


# Tree 1
NUM_NODES = 9
NUM_EDGES = 8
NUM_LEAVES = 5
nodes_t1 = [Tree.Node(i) for i in range(0, NUM_NODES)]
leaves_t1 = {i: nodes_t1[i] for i in range(0, NUM_LEAVES)}
edges_t1 = [
    Tree.Edge(nodes_t1[5], nodes_t1[0]), Tree.Edge(nodes_t1[5], nodes_t1[1]), 
    Tree.Edge(nodes_t1[6], nodes_t1[2]), Tree.Edge(nodes_t1[6], nodes_t1[5]),
    Tree.Edge(nodes_t1[7], nodes_t1[3]), Tree.Edge(nodes_t1[7], nodes_t1[6]),
    Tree.Edge(nodes_t1[8], nodes_t1[4]), Tree.Edge(nodes_t1[8], nodes_t1[7])
    ]
edge_to_index_t1 = {edge: i for i, edge in enumerate(edges_t1)}
root = nodes_t1[NUM_NODES - 1]
tree1 = Tree(edges=edges_t1, nodes=nodes_t1, root=root, leaves=leaves_t1)
nodes_t2 = [Tree.Node(i) for i in range(0, NUM_NODES)]
leaves_t2 = {i: nodes_t2[i] for i in range(0, NUM_LEAVES)}
root_t2 = nodes_t2[NUM_NODES - 1]
edges_t2 = [
    Tree.Edge(nodes_t2[5], nodes_t2[0]), Tree.Edge(nodes_t2[5], nodes_t2[3]), 
    Tree.Edge(nodes_t2[6], nodes_t2[5]), Tree.Edge(nodes_t2[6], nodes_t2[2]),
    Tree.Edge(nodes_t2[7], nodes_t2[6]), Tree.Edge(nodes_t2[7], nodes_t2[1]),
    Tree.Edge(nodes_t2[8], nodes_t2[4]), Tree.Edge(nodes_t2[8], nodes_t2[7])
    ]
tree2 = Tree(edges=edges_t2, nodes=nodes_t2, root=root_t2, leaves=leaves_t2)




if __name__ == '__main__':
    #####################################################
    #####           Maximum Agreement Forest        #####
    #####                   ILP Model               #####
    #####################################################

    NUM_EDGES = len(edges_t1)

    m = Model(name='example_model')

    #edges_cut = [m.binary_var(name=f'edges_cut_{i + 1}') for i in range(NUM_EDGES)]
    edges_cut = [m.binary_var(name=f'edges_cut_{(edges_t1[i]).from_().val()},{edges_t1[i].to().val()}') for i in range(NUM_EDGES)]
    incompatible_triples = compute_incompatible_triples(tree1=tree1, tree2=tree2)
    incompatible_paths = compute_incompatible_paths(tree1=tree1, tree2=tree2)

    # Constraint 1 (Incompatible Triples)
    for idx, incompatible_triple in enumerate(incompatible_triples, 1):
        m.add_constraint(sum(edges_cut[edge_to_index_t1[i]] for i in incompatible_triple) >= 1, ctname=f'incompatible_triple_{idx}')

    # Constraint 2 (Incompatible Paths)
    for idx, incompatible_path in enumerate(incompatible_paths, 1):
        m.add_constraint(sum(edges_cut[edge_to_index_t1[i]] for i in incompatible_path) >= 1, ctname=f'incompatible_path_{idx}')

    # A Maximum Agreement Forest (MAF) minimizes the number of edges cut (in other words,
    # it maximizes the number of components in which two trees agree)
    m.minimize(sum(edges_cut))

    m.print_information()
    sol = m.solve()
    if sol is None:
        raise RuntimeError("Model could not be solved (is a solver installed/available?)")
    m.print_solution()
