#---------------------------------------------------------------
# http://schmitzlab.info/phylo.html
# cran.r-project.org/web/packages/ape/ape.pdf
# http://rpubs.com/ematsen/ape-traversal-sample
#---------------------------------------------------------------

library(ape)

# returns the path with max sum of label values
# particulraly, it returns index of the leaf node ending such path

solution <- function(t, num, k, view = TRUE)
{
	paths <- nodepath(t);
	maximum <- c(-1, -1) # value, index
	for (i in 1:length(paths))
	{
		s <- sum(paths[[i]])
		#cat ("path", i, ": ", paths[[i]], " sum: ", s, "\n")
		if ( maximum[1] < s)	maximum <- c(s, i)
	}

	mycol <- rep("gray", num)
	index <- tail(paths[[maximum[2]]], n=1)
	mycol[index] <- "red"
	if (view)
	{
		#cat ("Maximum for (", paths[[maximum[2]]], ") with", maximum[1], "\n")
		#cat ("index = ", index)
		# store plot
		fname <- paste("tree_algebra/R/Trees/T_", num, "_", k, ".pdf", sep = "")
		pdf(file = fname, width = 6, height = 8, family = "Helvetica")
		plot(t, show.tip.label=FALSE, font = 1, label.offset = 0.4)
		nodelabels(); # shows internal node numbers
		tiplabels(col="black", bg=mycol) 
		dev.off()		
		
	}
	return (index)

}

# Loop over different n values

check <- function()
{
	for (n in 8:15)
	{
        for (k in 1:1)
        {
            set.seed(k)
            t <- rtree(n, FALSE)
            s <- solution(t, n, k);
            cat (n, " -> ", s, "\n");
            
            if (k>1)
            {
                ta <- rtree(n, FALSE)
                tb <- rtree(n, FALSE)
            }
        }
	}
}




