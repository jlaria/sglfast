#' Workaround to convert group lengths to index vector for the sparse-group lasso.
#'
#' \code{group.length2index} converts group lengths, e.g. \code{[1,2,3]} to \code{[1,2,2,3,3,3]}.
#'
#' @param group.lenght Vector with group lengths
#' @return Index vector
#' @seealso index2group.length
group.length2index = function(group.length){
  return(rep(1:length(group.length), group.length))
}
