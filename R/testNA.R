testNA <- function(i){
  out <- .Call("testNA", i)
  return (out)
}
