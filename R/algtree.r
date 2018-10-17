
# Algebraic manipulating of trees.
# v 0.5

library("treetop")
library("igraph")
library("ape")
library("treespace")
library("phytools")
library("adegenet")
library("adegraphics")
library("ggplot2")
library("blockcluster")



source("tree_algebra/R/NtoS.R")
source("tree_algebra/R/charOperations.R")
source("tree_algebra/R/unlabelledFuncs.R")

folder = "tree_algebra/R/Composition/T_";
flder = "tree_algebra/R/Transformations/";

# Need to load igraph and gmp for nthPhylo factory(), phylobase for phylo4()
# phylo object tree corresponding to n, n is the root

ta <- nthPhylo(8) 
tb <- nthPhylo(9)
tc <- nthPhylo(10)

plotlabels(ta)

# Get distance between two trees
distunlab(ta, tb)

# Prune a tree to the desired size, randomly uniformly dropping tips (reduction, 6 > 4)
t1 <- rtree(6)
t2 <- prunetosize(t1, 4)

# A vector of labels for the tree including tips and internal nodes
treelabels(rtree(10)) # the max is the root

find the adddition, multiplication, morphs!
play with the length of branch

# all paths:

library(ape)

num <- 5
set.seed(1)
t = rtree(num, FALSE)
plot(t, show.tip.label=FALSE)
nodelabels(); # shows internal node numbers
tiplabels();  # shows numerical tip values (name collision w/ above?!)

paths <- nodepath(t);

for (i in 1:length(paths))
{
	cat ("path: ", paths[[i]], " sum: ", sum(paths[[i]]), "\n")
}

#------------------------------
# Comparisons
#------------------------------

# rtree() takes the number of tpis as input
# the trees are not necessarly equal:

t1 <- rtree(6);
t2 <- rtree(6);
t3 <- rtree(6);
save_tree_plot(t1, "t1");
save_tree_plot(t2, "t2");
save_tree_plot(t3, "t3");
distunlab(t1, t2);
distunlab(t2, t3);
distunlab(t1, t3);

# However, nthPhylo() takes the root as input, same input give the same tree
ta <- nthPhylo(8)
tb <- nthPhylo(9)
tc <- nthPhylo(8)
distunlab(ta, tb)
distunlab(ta, tc) # = 0

# t1$tip.label # tips

save_tree_plot <- function(t, name, folder, title)
{
    fname = paste(folder, name,".pdf", sep = "");
    pdf(file = fname, width = 6, height = 8, family = "Helvetica");
    plot(t, main=title, show.tip.label=FALSE, font = 1, label.offset = 0.4);
    nodelabels(); # shows internal node numbers
    tiplabels(col="black");#, bg=mycol);
    edgelabels(round(t$edge.length, 3), bg="white", frame = "none", adj = c(0, -0.2));
    dev.off();
}

# Reproducing figure S1, done in FigureS5/
for (i in 1:5)
{
    save_tree_plot(nthPhylo(Z.to.N(i)), i);
}

#------------------------------------------------
# Checking convexity of the distance
# Doesnt always work!
#------------------------------------------------

testit <- function(a, b, c)
{
    ta <- nthPhylo(a);
    tb <- nthPhylo(b);
    tc <- nthPhylo(c);

    # Metric in [Colijn2017, end of page 5]
    d_ab <- abs(N.to.Qp(a)$p - N.to.Qp(b)$p);
    d_bc <- abs(N.to.Qp(b)$p - N.to.Qp(c)$p);
    d_ac <- abs(N.to.Qp(a)$p - N.to.Qp(c)$p);

    cat(d_ab + d_bc, "\n");
    cat(d_ac, "\n\n");

    save_tree_plot(ta, a)
    save_tree_plot(tb, b)
    save_tree_plot(tc, c)
}
#------------------------------------------------

# Figure S1
t <- nthPhylo(Z.to.N(-3) );
save_tree_plot(t, "-3")

# N.to.Qp() is not monotonic with IN
test <- function(a, b, c)
{
    cat(abs(a-b) + abs(b-c), "\n");
    cat(abs(a-c), "\n\n");

    cat(N.to.Qp(a)$q, N.to.Qp(b)$q, N.to.Qp(c)$q, "\n\n")

    d_ab <- abs(N.to.Qp(a)$q - N.to.Qp(b)$q);
    d_bc <- abs(N.to.Qp(b)$q - N.to.Qp(c)$q);
    d_ac <- abs(N.to.Qp(a)$q - N.to.Qp(c)$q);

    cat(d_ab + d_bc, "\n");
    cat(d_ac, "\n\n");

}

# t <- nthPhylo(Z.to.N(-2));

#-----------------------------------------------------------------------------------------------
# Adds new trees on each tip)
# https://www.rdocumentation.org/packages/treespace/versions/1.1.1/topics/simulateIndTree
#-----------------------------------------------------------------------------------------------

t <- nthPhylo(2);
save_tree_plot(t, "Mona_Before");
save_tree_plot(simulateIndTree(t, itips=2), "Mona_After_1");
save_tree_plot(simulateIndTree(t, itips=3), "Mona_After_2");

# random matrix
#M <- matrix(rbinom(num_row * num_col, 1, 0.5), ncol = num_col, nrow = num_row)
#Ma <- matrix(1, ncol = num_col, nrow = num_row)

#------------------------------------------------------------------------
# I.
# Find the transformation that maps ta to tb
# Case of one column vector
#------------------------------------------------------------------------

find_transformation <- function(a, b)
{
	a_ = matrix(1/a, ncol = n, nrow = 1)
	b_ = matrix(b, ncol = 1, nrow = n)
	R = 1/n * (b_ %*% a_)
	return (R)
}

k = 5 # number of leaves

ta = rtree(k);
tb = rtree(k);
tc = rtree(k);

treelabels(ta);
treelabels(tb);
treelabels(tc);

distunlab(ta, tb)
distunlab(ta, tc)
distunlab(tb, tc)

save_tree_plot(ta, "ta", flder);
save_tree_plot(tb, "tb", flder);
save_tree_plot(tc, "tc", flder);

a = treelabels(ta);
b = treelabels(tb);
c = treelabels(tc);
n = length(a)

R_ab = find_transformation(a, b)
R_bc = find_transformation(b, c)
R_ac = find_transformation(a, c)

# Checking composition

# R_ab a -> b
b_ <- setNames(c(R_ab %*% a), names(b))
all.equal(b, b_)

#z R_bc b -> c
c_ <- setNames(c(R_bc %*% b), names(c))
all.equal(c, c_)

# R_ac a -> c
c__ <- setNames(c(R_ac %*% a), names(c))
all.equal(c, c__)


# Dyadic product (dyadic tensor: https://en.wikipedia.org/wiki/Dyadics)
a = matrix(c(3, 4, 1), ncol = 1, nrow = 3);
b = matrix(c(2, 3, 0), ncol = 3, nrow = 1);
a %*% b

num_row = 3
num_col = 4

a <- round(matrix(runif(num_col * num_row, 1, 8), ncol=num_row));
b <- round(matrix(runif(num_col * num_row, 1, 8), ncol=num_row));
a_ <- solve(t(a) %*% a) %*% t(a);  # left inverse
a_a = round(a_ %*% a) # should be = I
ab = a %*% t(b) # dyadic product
R = a_ %*% ab

#--------------------------------------------------------------------------

N = 3
a <- round(matrix(runif(N * num_row, 1, 8), ncol=N));
a_inv_r = t(a) %*% solve(a %*% t(a));
all.equal(a %*% a_inv_r, diag(3))

b <- round(matrix(runif(N * num_row, 1, 8), ncol=N));

A <- b %*% a_inv_r

A %*% a

# A = b %*% t(a) %*% solve (a %*% t(a))

#-----------------------------------------------------------------------------
# The following matrices are random,
# they need to be tested on the matrices coming from rtree() as in (I)
#-----------------------------------------------------------------------------
go <- function(num_row, num_col)
{
	# Condition: num_row > num_col
	a <- round(matrix(runif(num_col* num_row, 1, 8), ncol=num_col))
	b <- round(matrix(runif(num_col * num_row, 1, 8), ncol=num_col))
	A1 = b %*% ginv(a);

	print ('Matrices:')
	print (a)
	print (b)

	print (all.equal(b, A1 %*% a))	

	# Compute the singular-value decomposition of a rectangular matrix.
	# X = U D V'
	s = svd(a)

	# generalised inverse of a, based on singular value decomposition.
	# http://astrostatistics.psu.edu/su07/R/html/base/html/svd.html
	inva = s$v %*% solve(diag(s$d)) %*% t(s$u) 
	A2 = b %*% inva
	print (all.equal(b, A2 %*% a))
	
	return (list(a=a, b=b, A1=A1, A2=A2))
}


go2 <- function(num_row, num_col)
{
	# Condition: num_row > num_col
	a <- round(matrix(runif(num_col* num_row, 1, 8), ncol=num_col))
	b <- round(matrix(runif(num_col * num_row, 1, 8), ncol=num_col))
	c <- round(matrix(runif(num_col * num_row, 1, 8), ncol=num_col))

	s = svd(a)
	inva = s$v %*% solve(diag(s$d)) %*% t(s$u) 
	A = b %*% inva

	sb = svd(b)
	invb = sb$v %*% solve(diag(sb$d)) %*% t(sb$u) 
	B = c %*% invb
	
	print (a)
	print (b)
	print (c)
	print(all.equal(ret$A %*% ret$a, ret$b))
	print(all.equal(ret$B %*% ret$b, ret$c))

	return (list(a=a, b=b, c=c, A=A, B=B))
}

#------------------------------------------------------------------------------------------
# Generates 3 rtrees, builds their MRCA matrices, and checks composition.

go4 <- function(num_leaves = 4)
{

	# Creating the trees and their matrices	
	ta = rtree(num_leaves);
	tb = rtree(num_leaves);
	tc = rtree(num_leaves);

	treelabels(ta);
	treelabels(tb);
	treelabels(tc);

	print (distunlab(ta, tb))
	print (distunlab(tb, tc))
	print (distunlab(ta, tc))

	# Saving trees
	flder = "tree_algebra/R/Transformations/";
	save_tree_plot(ta, "ta", flder);
	save_tree_plot(tb, "tb", flder);
	save_tree_plot(tc, "tc", flder);

	# vectors
	a = treelabels(ta);
	b = treelabels(tb);
	c = treelabels(tc);
	n = length(a)

	print ('Initial vectors, or tree topology:')
	print (a)
	print (b)
	print (c)

	print ('New matrices, topology and features:')
	# Redefining vectors with random 'num_features' random features
	min = 0;
	max = 8;

	# MRCA vectors, based on topology only (lambda=0) and branch lengths (lambda=1)
	# two features
	num_features = 2
	ncols = num_leaves * (num_leaves - 1) / 2 + num_leaves

	Ma <- matrix(0, nrow = num_features, ncol = ncols)
	Ma[1,] <- treeVec(ta, lambda=0) 
	Ma[2,] <- treeVec(ta, lambda=1) 

	Mb <- matrix(0, nrow = num_features, ncol = ncols)
	Mb[1,] <- treeVec(tb, lambda=0) 
	Mb[2,] <- treeVec(tb, lambda=1) 

	Mc <- matrix(0, nrow = num_features, ncol = ncols)
	Mc[1,] <- treeVec(tc, lambda=0) 
	Mc[2,] <- treeVec(tc, lambda=1) 

	print (Ma)
	print (Mb)
	print (Mc)

	# Transpose to get num_row > num_col
	a = t(Ma)
	b = t(Mb)
	c = t(Mc)
	colnames(a) <- colnames(b) <- colnames(c)  <- c("topo", "length");
	
	#pdf(file= paste(flder, "Ma.pdf", sep = "");
	#boxplot(Ma)
	#dev.off()

	s = svd(a);
	inva = s$v %*% solve(diag(s$d)) %*% t(s$u);
	A = b %*% inva;

	sb = svd(b);
	invb = sb$v %*% solve(diag(sb$d)) %*% t(sb$u);
	B = c %*% invb;
	
	print (all.equal(ret$A %*% ret$a, ret$b));
	print (all.equal(ret$B %*% ret$b, ret$c));

	return (list(ta=ta, tb=tb, tc=tc, a=a, b=b, c=c, A=A, B=B));

}

#------------------------------------------------
# Extracting edges, nodes, length

k = 4
t1 = rtree(k)
save_tree_plot(t1, "T1", flder);
t1$tip.label

t1%edge
t1$edge.length

m = length(t1$edge.length) # number of edges
n = t1$Nnode + 1; # number of leaves

#  Creating t2 by changing t1's lengths uniformly
t2 = t1
min = 1;
max = 8;
t2$edge.length = runif(m, min, max)
save_tree_plot(t2, "T2", flder);


#----------------------------------------------------------------------------------------------
# Actual extraction of the MRCA vector of the tree
# topological vector of MRCA distances from root (https://thibautjombart.github.io/treespace/reference/treeVec.html)
#----------------------------------------------------------------------------------------------

k = 4
t1 = rtree(k)

t1$edge.length

# based on topology only (lambda=0) else branch lengths (lambda=1)
treeVec(t1, lambda=0) 
treeVec(t1, lambda=1) 

# be careful, tips are not ordered as they appear when plotting the tree,
# they are based on ordered labels, using this vector to transform indices, in metrics.R
match(1:k, order(t1$tip.label))

save_tree_plot(t1, "T1", flder);

#----------------------------------------------------------------------------------------------
# generates rtree(k), builds its matrix, and finds the shortest path from root to leaves.
#----------------------------------------------------------------------------------------------

solution_shortest_path_based_on_length <- function(tx, tree_name, dir)
{
	# Setting descriptive matrix, with lengths are values.
	num_edges <- nrow(tx$edge);
	numNodes <- max(tx$edge);
	mat <- matrix(0.0, nrow = numNodes, ncol = numNodes);
	for (j in 1: num_edges)
	{
		e = tx$edge[j,];
		mat[e[1], e[2]] = tx$edge.length[j];
	}
	
	# edge, length
	#el <- cbind(tx$edge, tx$edge.length)
	
	# paths are found with nodepath() from 'ape'
	
	paths <- nodepath(tx);
	minimum <- c(100, -1) # value, index
	for (i in 1:length(paths))
	{
		path <- paths[[i]];
		cat (paste('\nPath', i, ':\n'));
		cat (path, '\t s = ');
		
		k <- 1: (length(path)-1);
		s = sum(mat[path[k] , path[k+1]]);
	
		cat (s, '\n');
		if ( minimum[1] > s)	minimum <- c(s, i)
	}
	
	sol1 = paths[[minimum[2]]];
	
	print(paste(" min cost path: (format 1)"));
	print (sol1)

	print(paste(" min cost path: (format 2)"));
	j <- 1: (length(sol1)-1);	
	sol2 = cbind(sol1[j], sol1[j+1])
	print(sol2)

	print(paste(" min cost path: (format 3)"));
	z <- apply(tx$edge, 1, function(a) apply(sol2, 1, function(b) all(a==b)))
	
	if (dim(sol2)[1] == 1)
	{
		sol3 <- as.integer(as.logical(z));
	}
	else
	{
		sol3 <- as.integer(as.logical(z[1,] + z[2,]));
	}
	print(tx$edge)
	print(sol3)
	save_tree_plot(tx,
					tree_name,
					dir,
					paste("Shortest path based on lengths\nS =",
					toString(sol3)));
	return (list(tree=tx, M=mat, sol1=sol1, sol2=sol2, sol3=sol3))
}
	
#------------------------------------------------------------------------------------------------	

test_1 <- function()
{
	dir = "tree_algebra/R/Generated/";
	k = 6
	num_trees = 10
	for (i in 1:num_trees)
	{
			T = rtree(k)
			s <- solution_shortest_path_based_on_length(T, paste("T_", i, sep=""), dir)
	}
	# TODO plot clusters of trees and clusters of solutions, given different transformations:
	# monotomic, non monotonic, etc.
}

#------------------------------------------------------------------------------------------------	
monotonic_length_transform <- function(tree, mode='monotonic')
{
	# bounds
	a <- 0.0
	b <- 2.0
	tree1 <- tree;
	n <- length(tree1$edge.length);

	if (mode == 'monotonic') 
	{
		epsilon = runif(1, a, b);	
	}
	if (mode == 'partial') 
	{
		# Partial transformation, cut in half, epsilon= (0 ... 0, uniform) 
		epsilon = c(rep(0, n/2), runif(n/2, a, b))	
	}
	if (mode == 'random') 
	{
		epsilon = runif(n, a, b);		
	}
	
	print (cat(' epsilon = ', epsilon));
	print (cat(' tree1$edge.length = ', tree1$edge.length));
	
	tree1$edge.length <- tree1$edge.length + epsilon;
	return (tree=tree1);
}

Commdiag <- function()
{
	dir = "tree_algebra/R/Generated/";
	k = 8
	
	T1 = rtree(k)
	S1 <- solution_shortest_path_based_on_length(T1,"_T_1", dir)$sol3

	T2 = rtree(k)
	S2 <- solution_shortest_path_based_on_length(T2,"_T_2", dir)$sol3

	T3 = monotonic_length_transform(T2, 'random')
	S3 <- solution_shortest_path_based_on_length(T3,"_T_3", dir)$sol3

	T4 = monotonic_length_transform(T3, 'partial')
	S4 <- solution_shortest_path_based_on_length(T4,"_T_4", dir)$sol3

	T5 = monotonic_length_transform(T1, 'monotonic')
	S5 <- solution_shortest_path_based_on_length(T5,"_T_5", dir)$sol3
	
	T6 = monotonic_length_transform(T4, 'monotonic')
	S6 <- solution_shortest_path_based_on_length(T6,"_T_6", dir)$sol3

	T7 = rtree(k)
	S7 <- solution_shortest_path_based_on_length(T7,"_T_7", dir)$sol3

	T8 = rtree(k)
	S8 <- solution_shortest_path_based_on_length(T8,"_T_8", dir)$sol3

	T9 = rtree(k)
	S9 <- solution_shortest_path_based_on_length(T9,"_T_9", dir)$sol3

	T10 = rtree(k)
	S10 <- solution_shortest_path_based_on_length(T10,"_T_10", dir)$sol3

	# Building multiPhylo object of the trees
	S_ <- rbind(S1, S2, S3, S4, S5, S6, S7, S8, S9, S10)
	y <- c(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10)
	names(y) <- paste("T", 1:length(y), sep = "")

	y <- .compressTipLabel(y)

	mat <- matrix(0, ncol = length(S10), nrow = 10)
	mat <- S_

	rownames(mat) <- paste("S", 1:length(taby), sep = "")
	print (mat)	

	# Using 'blockcluster'
	#-- out <- coclusterBinary(mat, nbcocluster=c(2,3))
	#" Summarize the output results
	#--summary(out)
	# Plot the original and Co-clustered data
	#--plot(out)

	# kmeans
	kmeans(mat, centers=3)

	# hierar clustering
	d <- dist(mat, method = "binary")
	hc <- hclust(d)

    pdf(file = paste(dir, "Solutions.pdf"), width = 6, height = 8, family = "Helvetica");
	plot(hc, main='Solutions')
    dev.off();

	# using treespace, https://thibautjombart.github.io/treespace/index.html
	res <- treespace(y, nf=3)
	# table.image
	#table.image(res$D, nclass=30)

	# table.value with some customization
	#table.value(res$D, nclass=6, method="color",  symbol="circle", col=redpal(10))

    pdf(file = paste(dir, "Problems.pdf"), width = 6, height = 8, family = "Helvetica");
	plotGroves(res$pco, lab.show=TRUE, lab.cex=1.0)
    dev.off();


}


solution_shortest_path_based_on_length_and_prob <- function(k)
{
}


# End.

