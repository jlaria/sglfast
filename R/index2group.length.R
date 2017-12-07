#' Workaround to convert index vector to group lengths.
#'
#' @seealso group.length2index
index2group.length = function(index){
  ord <- order(index)
  index <- index[ord]

  groups <- unique(index)
  num.groups <- length(groups)
  range.group.ind <- rep(0,(num.groups+1))
  for(i in 1:num.groups){
    range.group.ind[i] <- min(which(index == groups[i])) - 1
  }
  range.group.ind[num.groups + 1] <- length(index)

  group.length <- diff(range.group.ind)
  return( list(group.length = group.length, ord = ord) )
}
