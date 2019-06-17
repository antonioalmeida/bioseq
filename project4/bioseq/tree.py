class Tree:

    def __init__(self, val=0, dist = 0, left = None, right = None):
        self.value = val
        self.distance = dist
        self.left = left
        self.right = right

    def get_cluster(self):
        res = []
        if self.value != -1:
            res.append(self.value)
        else:
            if self.left:
                res.extend(self.left.get_cluster())
            if self.right:
                res.extend(self.right.get_cluster())
        return res

    def print_tree(self):
        self.print_tree_rec(0, "Root")
    
    def print_tree_rec (self, level, side):
        tabs = ""
        for i in range(level): tabs += "\t"
        if self.value != -1:
            print(tabs, side, " - value:", self.value)
        else:
            print(tabs, side, "- Dist.: ", self.distance)
            if (self.left != None): 
                self.left.print_tree_rec(level+1, "Left")
            if (self.right != None): 
                self.right.print_tree_rec(level+1, "Right")

    def size(self):
        ''' size of the tree: returns two values
        - number of internal nodes of the tree
        - number of leaves'''
        numleaves = 0
        numnodes = 0
        if self.value != -1:
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
        if self.value >= 0:
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

def test():              
    a = Tree(1)
    b = Tree(2)
    c = Tree(3)
    d = Tree(4)
    e = Tree(-1, 2.0, b, c)
    f = Tree(-1, 1.5, d, a)
    g = Tree(-1, 4.5, e, f)
    g.print_tree()
    print(g.get_cluster())

    # testing exercise 3
    print(g.size())
    print(g.exists_leaf(1))
    print(g.exists_leaf(5))
    g.common_ancestor(1,4).print_tree()
    g.common_ancestor(1,2).print_tree()
    print(g.distance_leaves(1,4))
    print(g.distance_leaves(1,2))

if __name__ == '__main__':
    test()
