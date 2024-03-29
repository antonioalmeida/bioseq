class Tree:

    def __init__(self, val=None, dist = 0, left = None, right = None, species='mekie'):
        self.value = val
        self.distance = dist
        self.left = left
        self.right = right
        self.species = species

    def get_cluster(self):
        res = []
        if self.value:
            res.append(self.value)
        else:
            if self.left:
                res.extend(self.left.get_cluster())
            if self.right:
                res.extend(self.right.get_cluster())
        return res

    def print_tree(self):
        print("> Root")
        self._print(0, "    ")
    
    def _print(self, level, side):
        tabs = ""
        for i in range(level): tabs += "\t"
        if self.value:
            print(tabs, side, "> Leaf(value=" + str(self.value) + ")")
        else:
            print(tabs, side, "> Node(length={%.4f})" % self.distance)
            if self.left: 
                self.left._print(level+1, "    ")
            if self.right: 
                self.right._print(level+1, "    ")

    def size(self):
        ''' size of the tree: returns two values
        - number of internal nodes of the tree
        - number of leaves'''
        numleaves = 0
        numnodes = 0
        if self.value:
            numleaves = 1
        else: 
            if (self.left != None):
                resl = self.left.size()
            else: resl = (0,0)
            if (self.right != None):  
                resr = self.right.size() 
            else: resr = (0,0)
            numnodes += (resl[0] + resr[0] + 1)
            numleaves += (resl[1] + resr[1])
        return numnodes, numleaves

    def exists_leaf(self, leafnum):
        if self.value:
            return self.value == leafnum
        else:
            if self.left:
                res_l = self.left.exists_leaf(leafnum)
            if self.right:
                res_r = self.right.exists_leaf(leafnum)
            return res_l or res_r
    
    def common_ancestor(self, l1, l2):
        ''' Return simplest tree that contains l1, l2
        '''
        if self.left and self.right:
            if self.left.exists_leaf(l1) and self.left.exists_leaf(l2):
                return self.left.common_ancestor(l1, l2)
            if self.right.exists_leaf(l1) and self.right.exists_leaf(l2):
                return self.right.common_ancestor(l1, l2)

            if self.left.exists_leaf(l1) and not self.left.exists_leaf(l2):
                return self
            if self.left.exists_leaf(l2) and not self.left.exists_leaf(l1):
                return self

    def distance_leaves(self, l1, l2):
        ca = self.common_ancestor(l1,l2)
        return ca.distance*2

    def phylo_tree(self, childattr='children', nameattr='name', indent='', last='updown'):

        name = str(self.value) if self.value else int(self.distance)*'--'

        children = lambda node: [node.left, node.right] if node.left != None else []
        nb_children = lambda node: sum(nb_children(child) for child in children(node)) + 1
        size_branch = {child: nb_children(child) for child in children(self)}

        """ Creation of balanced lists for "up" branch and "down" branch. """
        up = sorted(children(self), key=lambda node: nb_children(node))
        down = []
        while up and sum(size_branch[node] for node in down) < sum(size_branch[node] for node in up):
            down.append(up.pop())

        """ Printing of "up" branch. """
        for child in up:     
            next_last = 'up' if up.index(child) is 0 else ''
            next_indent = '{0}{1}{2}'.format(indent, ' ' if 'up' in last else '│', ' ' * len(name))
            child.phylo_tree(childattr, nameattr, next_indent, next_last)

        """ Printing of current node. """
        if last == 'up': start_shape = '┌'
        elif last == 'down': start_shape = '└'
        elif last == 'updown': start_shape = ' '
        else: start_shape = '├'

        if up: end_shape = '┤'
        elif down: end_shape = '┐'
        else: end_shape = ''

        print('{0}{1}{2}{3}'.format(indent, start_shape, name, end_shape))

        """ Printing of "down" branch. """
        for child in down:
            next_last = 'down' if down.index(child) is len(down) - 1 else ''
            next_indent = '{0}{1}{2}'.format(indent, ' ' if 'down' in last else '│', ' ' * len(name))
            child.phylo_tree(childattr, nameattr, next_indent, next_last)

