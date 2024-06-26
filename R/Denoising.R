

BWF <- function(nc, fc) {
    if (nc %% 2 == 0)
        stop("filter length (nc) must be odd !")
    m <- (nc - 1) / 2
    sm <- -m:m
    val <- sapply(sm, function(i, nc, fc, m) {
        1 * sin(2 * pi * fc * i) * (0.42 - 0.5 * cos(2 * pi * (i + m) / (2 * m)) +
                                        0.08 * cos(4 * pi * (i + m) / (2 * m))) / i
    }, nc = nc, fc = fc, m = m)
    val[m + 1] <- 2 * pi * fc
    val <- val / sum(val)
    val
}


extendSize <- function(x, fsize, type = c("padding", "reflection")) {
    nx <- length(x)
    if (nx != fsize) {
        to_add <- fsize - nx
        if (to_add == 0)
            return(0)
        type <- match.arg(type)
        if (type == "padding") {
            x <- c(x, rep(0, to_add))
        }
        if (type == "reflection") {
            x <- c(x, x[nx:(nx - to_add + 1)])
        }
    }
    x
}


##mov avg from pracma
movavg <- function (x, n, type = c("s", "t", "w", "m", "e", "r"))
{
    stopifnot(is.numeric(x), is.numeric(n), is.character(type))
    if (length(n) != 1 || ceiling(n != floor(n)) || n <= 1)
        stop("Window length 'n' must be a single integer greater 1.")
    nx <- length(x)
    if (n >= nx)
        stop("Window length 'n' must be greater then length of time series.")
    y <- numeric(nx)
    if (type == "s") {
        for (k in 1:(n - 1))
            y[k] <- mean(x[1:k])
        for (k in n:nx)
            y[k] <- mean(x[(k - n + 1):k])
    }
    else if (type == "t") {
        n <- ceiling((n + 1) / 2)
        s <- movavg(x, n, "s")
        y <- movavg(s, n, "s")
    }
    else if (type == "w") {
        for (k in 1:(n - 1))
            y[k] <- 2 * sum((k:1) * x[k:1]) / (k *
                                                   (k + 1))
        for (k in n:nx)
            y[k] <- 2 * sum((n:1) * x[k:(k - n +
                                             1)]) / (n * (n + 1))
    }
    else if (type == "m") {
        y[1] <- x[1]
        for (k in 2:nx)
            y[k] <- y[k - 1] + (x[k] - y[k - 1]) / n
    }
    else if (type == "e") {
        a <- 2 / (n + 1)
        y[1] <- x[1]
        for (k in 2:nx)
            y[k] <- a * x[k] + (1 - a) * y[k - 1]
    }
    else if (type == "r") {
        a <- 1 / n
        y[1] <- x[1]
        for (k in 2:nx)
            y[k] <- a * x[k] + (1 - a) * y[k - 1]
    }
    else
        stop("The type must be one of 's', 't', 'w', 'm', 'e', or 'r'.")
    return(y)
}


convolveNormalized <-
    function(x,
             f,
             type = c("circular", "open", "filter"),
             extending = c("padding", "reflection")) {
        type <- match.arg(type)
        extending <- match.arg(extending)
        nx <- length(x)
        nf <- length(f)
        ##First putting the size to a multiple of 2f
        nnx <- nextn(nx, 2)
        nnf <- nextn(nf, 2)
        sizemax <- max(nnx, nnf)
        #if(sizemax!=nnx) stop("Filter length is too long")
        if ((nx != sizemax | nf != sizemax) & type != "filter") {
            x <- extendSize(x, sizemax, extending)
            f <- extendSize(f, sizemax, "padding")
        }
        resseq <- convolve(x, f, type = type)
        #plot(x,type="l")
        #lines(resseq,col="blue")
        if (type == "filter") {
            res <- resseq
        } else{
            res <- c(resseq[(nnx - floor(nf / 2) + 1):nnx], resseq[1:(nnx - floor(nf /
                                                                                     2))])#(1+floor((nf+1)/2)):(nx+floor((nf+1)/2))]
            res <- res[1:nx]
        }
        res
    }

####savogol function from pracma
SavGol <- function (vt,
                    fl,
                    forder = 4,
                    dorder = 0)
{
    stopifnot(is.numeric(vt), is.numeric(fl))
    if (fl <= 1 || fl %% 2 == 0)
        stop("Argument 'fl' must be an odd integer greater than 1.")
    n <- length(vt)
    fc <- (fl - 1) / 2
    X <- outer(-fc:fc, 0:forder, FUN = "^")
    Y <- pinv(X)
    T2 <- convolve(vt, rev(Y[(dorder + 1),]), type = "o")
    T2 <- T2[(fc + 1):(length(T2) - fc)]
    Tsg <- (-1) ^ dorder * T2
    return(Tsg)
}

smooth.SavGol <- function(x, fl = 21, ...) {
    return(SavGol(x, fl, 2))
}

smooth.BWF <- function(x, l = length(x), freq = 0.3) {
    fil <- BWF(l + ((l + 1) %% 2), freq)
    return(convolveNormalized(x, fil))
}
