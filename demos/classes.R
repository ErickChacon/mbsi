
     x <- 10
     class(x) # "numeric"
     oldClass(x) # NULL
     inherits(x, "a") #FALSE
     class(x) <- c("a", "b")
     inherits(x,"a") #TRUE
     inherits(x, "a", TRUE) # 1
     inherits(x, c("a", "b", "c"), TRUE) # 1 2 0

print.b <- function (x) {
}

print(x)

x <- rnorm(100)
y <- rnorm(100)
lm1 <- lm(y~x)

bla.lm <- function (x) {
  "hola"
}

bla(lm1)

f <- function(x) UseMethod("f")
f.a <- function(x) "Class a"

a <- structure(list(), class = "a")
class(a)

mean.a <- function(x) "a"
mean(a)
