
# Multiple trees
# https://thibautjombart.github.io/treespace/index.html

library("treespace")
library("adegenet")
library("adegraphics")
library("ggplot2")

# generate list of trees
set.seed(1)
x <- rmtree(10, 20)
names(x) <- paste("tree", 1:10, sep = "")

# use treespace
res <- treespace(x, nf=3)
names(res)

# table.image
table.image(res$D, nclass=30)

# table.value with some customization
table.value(res$D, nclass=5, method="color", symbol="circle", col=redpal(5))

plotGroves(res$pco, lab.show=TRUE, lab.cex=1.5)

