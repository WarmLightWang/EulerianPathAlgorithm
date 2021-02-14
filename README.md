# The Eulerian path algorithm for spectral assembly

Euler_assemble that takes as input a set of k-mers and outputs a shortest superstring that has a k-mer spectrum identical to the set of input k-mers. Using the Eulerian path finding algorithm to solve this problem.

We will assume that:

1.we are assembling a single-stranded DNA sequence
2.there exists such a superstring
3.𝑘>1 

There exist valid k-mer spectra for which the corresponding k-mer graph has two unbalanced vertices, and the “fake edge” one would add to balance those vertices already exists in the graph. To make the Eulerian cycle algorithm work, one still needs to add an extra “fake edge” in this case, which results in two identical edges, making the graph a “multigraph”. Using a data structure that can accomodate this case.
