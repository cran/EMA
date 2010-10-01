plotInertia <-
function(acp, ncp=5, ...) {

## Goal
## Barplot of inertia percentage of components

## Input
## acp : result from PCA function
## ncp : number of components displayed

## Output
## Barplot of inertia percentage of components
  barplot(acp$eig[1:ncp, 2], main = "Inertia percentage of components", las=2, names.arg = paste("Dim", 1:ncp, sep = ""), ...)
}

