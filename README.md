#The impact of tree prior and molecular clock models on the estimate of diversificaton rates

Here, I provide code and methods to simulate trees, simulate sequence data, estimate trees 
using BEAST, and estimate diversification parameters.

##Dependencies
In addition to R, you will need:
* BEAST v1.7.5 *
This approach ought to work with BEAST v2.0+, though I expect there to be differences in the XML input
* Seq-Gen
* BEASTifier

You will also need several R libraries:
* parallel
* ape
* laser
* TreeSim
* geiger

You can install these, and their dependencies, from the command line i.e.:
```
install.packages("ape")
```

Alternatively, you can download the source code, compile, and install locally if you are 
working on a system where you lack privileges.


##Execution

