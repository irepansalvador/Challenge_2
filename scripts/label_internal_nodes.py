from ete3 import Tree
import sys


# Load a tree structure from a newick file.
print(sys.argv[1]);
tree = Tree(sys.argv[1])
print(tree)

out = sys.argv[1] + "_intNodes.nw";
print(out);
#tree = Tree('(((A,B),C),D);')

edge = 0
for node in tree.traverse():
   if not node.is_leaf():
      node.name = "Node_%d" %edge
      edge += 1

print(tree.write(format=8))
print(tree)

# We can also write into a file
tree.write(format=2, outfile=out)

# (((A,B)NODE_2,C)NODE_1,D);
#
