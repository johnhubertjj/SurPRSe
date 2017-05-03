tryCatch.W.E <- function(expr)
{
  W <- NULL
  w.handler <- function(w) { # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                          warning = w.handler),
                   	             warning = W)
}

str( tryCatch.W.E( log( 2 ) ) )
str( tryCatch.W.E( log( -1) ) )
tryCatch.W.E( log( -1))

vari <- 1
tryCatch(print("passes"), error = function(e) print(vari), finally=print("finished")) 
tryCatch(stop("fails"), error = function(e) print(vari), finally=print("finished")) 

tryCatch(1, finally = print("Hello"))
e <- simpleError("test error")
## Not run: 
stop(e)
tryCatch(stop(e), finally = print("Hello"))
tryCatch(stop("fred"), finally = print("Hello"))

## End(Not run)
tryCatch(stop(e), error = function(e) e, finally = print("Hello"))
tryCatch(stop("fred"),  error = function(e) e, finally = print("Hello"))
withCallingHandlers({ warning("A"); 1+2 }, warning = function(w) {})
## Not run: 
{ withRestarts(stop("A"), abort = function() {}); 1 }

## End(Not run)
withRestarts(invokeRestart("foo", 1, 2), foo = function(x, y) {x + y})
