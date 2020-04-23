# Load Packages 
library(IRanges)

# Create IRange Objects 
ir1 <- IRanges(start = c(1, 3, 5), end = c(3, 5, 7))
ir2 <- IRanges(start = c(1, 3, 5), width = 3) 

# Equal 
all.equal(ir1, ir2)

# Basic Info 
start(ir1)
end(ir1)
width(ir1)

width(ir1) <- 1
ir1
names(ir1) <- paste("A", 1:3, sep = "")
ir1

dim(ir1)

length(ir1)