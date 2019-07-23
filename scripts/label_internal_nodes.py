from ete3 import Tree


# Load a tree structure from a newick file.
tree = Tree('bd_1.nw')
print(tree)
#tree = Tree('(((A,B),C),D);')

edge = 0
for node in tree.traverse():
   if not node.is_leaf():
      node.name = "Node_%d" %edge
      edge += 1

print(tree.write(format=8))
print(tree)

# We can also write into a file
tree.write(format=2, outfile="new_tree.nw")

# (((A,B)NODE_2,C)NODE_1,D);
#
