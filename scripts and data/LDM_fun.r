form.call = function( object, data=NULL ) {
    forml=object
    if (is.null(data)) data=environment(forml)
    else data=model.frame(data, na.action=na.pass)
    obj=toString(object)
    obj=gsub('\\s','',obj)
    prefix=' ~ + 0 +'
    loc.comma=gregexpr(',',obj)[[1]]
    start.terms=loc.comma[2]
    terms=substr(obj,start=start.terms+1, stop=nchar(obj))
    model=list()
    j=1
    continue=TRUE
    while (continue) {
        if (substr(terms,1,1)=='(') {
            stop=regexpr(')\\+',terms)[1]
        }
        else {
            stop=regexpr('\\+',terms)[1] - 1
        }          
        
        if (stop<=0) stop=nchar(terms) 
        m=substr(terms, start=1, stop=stop)
        #      model[[j]]=model.matrix( as.formula( paste(prefix,m) ) , data=environment(forml) )
        model[[j]]=model.matrix( as.formula( paste(prefix,m) ) , data=data)
        if (stop+2<=nchar(terms)) {
            terms=substr(terms, start=stop+2, stop=nchar(terms))
            j=j+1
        }
        else {
            continue=FALSE
        }             
    }    
    #
    #  check for and extract confounders
    #
    loc.bar=regexpr('\\|',obj)[1]
    loc.minus=regexpr('-',obj)[1]
    loc.delim=max( loc.bar, loc.minus)
    if (loc.delim>0) {
        end.confound=loc.comma[2]
        c=substr(obj,start=loc.delim+1, stop=end.confound-1)
        #      conf=model.matrix( as.formula( paste(prefix,c) ), data=environment(forml) ) 
        conf=model.matrix( as.formula( paste(prefix,c) ), data=data ) 
    }
    else {
        conf=NULL
    }     
    #
    #  extract OTU table
    #      
    if (is.null(conf)) loc.delim=loc.comma[2]
    otu.name=substr(obj, start=loc.comma[1]+1, stop=loc.delim-1)
    otu.table=get(otu.name, envir=environment(forml))
    res=list( model=model,conf=conf, otu.table=otu.table)
    return(res)
}    



biplotQuantities = function (decomp.results, use.model=1:2) 
{
    #------------------------ 
    # copy from decomp.results
    #------------------------ 
    
    x <- decomp.results$x
    b <- decomp.results$b
    v <- decomp.results$v
    d <- decomp.results$d
    x.hat.model <- t(t(b[, use.model]) * d[use.model])
    x.hat.model <- x.hat.model %*% t(v[, use.model])
    
    if (dim(x)[1] != dim(b)[1]) x <- t(x)
    
    #--------------------------------------------------------
    # Quantities 
    #--------------------------------------------------------
    
    w <- t(d * t(v))
    otu.coords <- w[, use.model, drop=FALSE]              
    
    otu.score <- rowSums(otu.coords^2)                    # sum_k (Vjk . Dkk)^2, the norm of the coordinates in otu.coords
    otu.order <- order(otu.score, decreasing = TRUE)      # order by otu.score
    otu.coords.ordered <- w[otu.order, use.model, drop=FALSE]
    top.otus <- colnames(x)[otu.order]
    
    component.score <- d^2                                
    component.order <- order(component.score)
    
    ss.total <- sum(x^2)                                 
    ss.model <- sum(x.hat.model^2)
    ss.resid <- sum((x - x.hat.model)^2)
    
    out <- list(otu.coords = otu.coords, 
                otu.coords.ordered = otu.coords.ordered, 
                otu.score = otu.score, 
                otu.order = otu.order, 
                top.otus = top.otus, 
                component.score = component.score, 
                component.order = component.order, 
                ss.total = ss.total, 
                ss.resid = ss.resid, 
                ss.model = ss.model, 
                x.hat.model = x.hat.model)
    
    return(out)
}


gower = function (d, square=TRUE, center=TRUE) 
{
    
    a <- as.matrix(d)
    
    #--------------------------------------------------------------
    # squaring the matrix of d
    #--------------------------------------------------------------
    
    if (square) a <- a^2
    
    #--------------------------------------------------------------
    # centering rows and columns of the (squared) matrix of d
    #--------------------------------------------------------------
    
    if (center) {
        rs <- rowMeans(a)
        cs <- colMeans(a)
        
        n <- dim(a)[1]
        
        a <- a - matrix(rep(rs,each=n), ncol=n, byrow=TRUE)
        a <- a - matrix(rep(cs,each=n), nrow=n)
        
        a <- (-0.5) * (a + mean(cs)) 
    }
    
    return(a)
    
}# gower


setup.model = function( vars, adjust.for.confounders=FALSE, center.vars=TRUE, 
                        is.clustered=FALSE, cluster.id=NULL, 
                        permute.within.clusters=FALSE, permute.between.clusters=FALSE) {
    
    n.vars = length(vars)
    n.obs = dim(vars[[1]])[1]

    if (!is.clustered) 
    {
        n.cluster = NULL
        cluster.size = NULL
        cluster.paradigm = NULL
        cluster.order = NULL
    } else 
    {
        cluster.id = as.character(cluster.id)
        unique.id = unique(cluster.id)
        n.cluster = length(unique.id)
        
        cluster.size = rep(0,n.cluster)
        cluster.paradigm = rep(0,n.cluster)
        cluster.order = NULL
        
        for (i in 1:n.cluster) 
        {
            this.cluster = which(cluster.id==unique.id[i]) 
            cluster.order = c(cluster.order, this.cluster)
            cluster.paradigm[i] = this.cluster[1]
            cluster.size[i] = length(this.cluster)
        }
    }       

    m1 = NULL
    m1.index = NULL
    main.effect = NULL
    
    if (adjust.for.confounders) 
    {
        index = rep(0, n.vars-1)
        
        m.i  =  vars[[1]]
        if (center.vars) m.i = scale( m.i, center=TRUE, scale=FALSE )
        
        m1 = m.i
        m1.index = dim(m1)[2]
        svd.m1 = svd(m1)
        tol.d=10^-8 
        use = (svd.m1$d>tol.d)      
        proj = svd.m1$u[, use] %*% t( svd.m1$u[, use] )      
        
        for (i in 2:n.vars) {
            m.i  =  vars[[i]]
            if (center.vars) m.i = scale( m.i, center=TRUE, scale=FALSE )
            
            if (i==2) {
                m = m.i
                index[i-1] = dim(m.i)[2] 
            } else {
                m = cbind(m, m.i)   
                index[i-1] = index[i-2] + dim(m.i)[2]    
            }
        }
        
        low = 1
        up = index[n.vars-1]
        
        main.effect = proj %*% m[, low:up]
        resid = m[, low:up] - main.effect
        m[, low:up] = resid
        
    } else {
        
        index = rep(0, n.vars)
        
        for (i in 1:n.vars) {
            m.i  =  vars[[i]]
            if (center.vars) m.i = scale( m.i, center=TRUE, scale=FALSE )
            
            if (i==1) {
                m = m.i
                index[i] = dim(m.i)[2] 
            } else {
                m = cbind(m, m.i)   
                index[i] = index[i-1] + dim(m.i)[2]    
            }
            
        }
    }
    
    out = list( m=m, index=index, m1=m1, m1.index=m1.index, 
                adjust.for.confounders=adjust.for.confounders, center.vars=center.vars, 
                n.obs=n.obs, n.vars=n.vars, 
                is.clustered=is.clustered, 
                permute.within.clusters=permute.within.clusters,
                permute.between.clusters=permute.between.clusters,
                n.cluster=n.cluster, cluster.size=cluster.size, 
                cluster.order=cluster.order, cluster.paradigm=cluster.paradigm)   
    return(out)
    
} # setup.model


fdr.Sandev = function(p.otu) {
    
    m = length(p.otu)
    
    p.otu.sort = sort(p.otu)
    n.otu.detected = seq(1, m)
    pi0 = min(1, 2/m*sum(p.otu))
    
    qval.sort = m * pi0 * p.otu.sort / n.otu.detected
    j.min.q = 1
    while (j.min.q < m) {
        min.q = min( qval.sort[j.min.q:m] )
        new.j.min.q = (j.min.q-1) + max( which(qval.sort[j.min.q:m]==min.q) )
        qval.sort[j.min.q:new.j.min.q] = qval.sort[new.j.min.q]
        j.min.q = new.j.min.q+1
    }
    mat = match(p.otu, p.otu.sort)   
    qval.orig = qval.sort[mat]
    results = qval.orig
    return(results)
    
} # fdr.Sandev


factor.to.hat.matrix = function( factor ) 
{
    if (is.factor(factor)==FALSE) factor=as.factor(factor)
    n.factor=length(unique(factor))
    n.obs=length(factor)
    x=mat.or.vec(n.obs,n.factor)
    levels=unique(factor)
    for (i in 1:n.obs) {
        j=match( factor[i], levels)
        x[i,j]=1
    }
    
    QRx=qr(x)
    Q=qr.Q(QRx)
    hat.matrix= Q %*% t(Q)   
    #  project off constant
    r=rowMeans(hat.matrix)
    c=colMeans(hat.matrix)
    tot=sum(hat.matrix)/n.obs^2
    hat.matrix=hat.matrix - outer(r,rep(1,n.obs)) - outer(rep(1,n.obs),c) + tot*matrix( rep(1,n.obs*n.obs),nrow=n.obs)   
    
    res=list( x, hat.matrix )
    names(res)=c( 'x.matrix', 'hat.matrix' )
    return(res)
    
} # factor.to.hat.matrix


fit.m1 = function( d.gower, x.freq, x.tran, model) {
    
    #----------------------------------------------------------
    # output: d.resid, x.tilde, x.svd 
    #----------------------------------------------------------
    
    tol.d=10^-8
    
    m1.fit = NULL
    
    n.obs = dim(d.gower)[1]
    if (dim(x.freq)[1] != n.obs) x.freq = t(x.freq)
    if (dim(x.freq)[1] != n.obs) stop( 'numbers of observations mismatch between x.freq and d' )
    if (dim(x.tran)[1] != n.obs) x.tran = t(x.tran)
    if (dim(x.tran)[1] != n.obs) stop( 'numbers of observations mismatch between x.tran and d' )
    
    d.resid = d.gower
    
    var = model$m1
    svd.var = svd(var)  
    use1 = (svd.var$d>tol.d)    
    
    hat.matrix = svd.var$u[, use1] %*% t( svd.var$u[, use1] )
    
    #---------------------
    # calculate direction
    #---------------------
    
    n.dim = dim( hat.matrix)[1]
    
    d.model = hat.matrix %*% d.resid
    d.model = d.model %*% hat.matrix
    
    es.model = eigen(d.model, symmetric=TRUE) 
    
    use = ( abs(es.model$values)>tol.d )
    ndf.model = sum( use )
    
    b.model = svd.var$u[, use1] # es.model$vectors[, use, drop=FALSE]
    
    hat.matrix.bar = diag(n.dim)  - hat.matrix
    d.resid = hat.matrix.bar %*% d.resid
    d.resid = d.resid %*% hat.matrix.bar
    
    #-----------------------------
    # end of calculating direction
    #-----------------------------    
    
    m1.fit$d.resid = d.resid
    
    x1.tilda.freq = t(b.model) %*% x.freq
    x1.tilda.freq = b.model %*% x1.tilda.freq
    m1.fit$x.freq.proj = x.freq - x1.tilda.freq
    m1.fit$x.freq.svd = svd(m1.fit$x.freq.proj)
    
    x1.tilda.tran = t(b.model) %*% x.tran
    x1.tilda.tran = b.model %*% x1.tilda.tran
    m1.fit$x.tran.proj = x.tran - x1.tilda.tran
    m1.fit$x.tran.svd = svd(m1.fit$x.tran.proj)
    
    return(m1.fit)
    
} # fit.m1


fit.ldm.internal = function( d.gower, x.freq, x.tran, model, x.freq.svd, x.tran.svd,
                    n.iter.max=1000, tol.conv=10^-8, decomposition="orthogonal", verbose=FALSE) {
    
    index = model$index
    m = model$m
    m1 = model$m1
    
    n.vars = length(index)
    ndf.nominal = rep(0, n.vars+1)
    
    tol.d = 10^-8
    
    #--------------------------------------------------------------------------
    # construct directions matrix b 
    # from each set of variables in the list vars
    #--------------------------------------------------------------------------
    
    d.resid = d.gower
    
    for (i in 1:n.vars) 
    {
        
        if (model$adjust.for.confounders==TRUE) { # ???
            var = cbind(m1, m[,1:index[i]])
        } else {
            var = m[,1:index[i]]
        }
        
        svd.var = svd(var)   # not meaningful to save when conditioning on B1
        use = (svd.var$d>tol.d)    
        
        hat.matrix = svd.var$u[, use] %*% t( svd.var$u[, use] )
        
        #---------------------
        # calculate direction
        #---------------------
        
        n.dim = dim( hat.matrix)[1]
        
        d.model = hat.matrix %*% d.resid
        d.model = d.model %*% hat.matrix
        
        es.model = eigen(d.model, symmetric=TRUE) # es: eigen system in Mathematica
        
        use = ( abs(es.model$values)>tol.d )
        ndf.model = sum( use )
        
        b.model = es.model$vectors[, use]
        e.model = es.model$values[use]
        
        hat.matrix.bar = diag(n.dim)  - hat.matrix
        d.resid = hat.matrix.bar %*% d.resid
        d.resid = d.resid %*% hat.matrix.bar
        
        #-----------------------------
        # end of calculating direction
        #-----------------------------    
        
        if (i==1) {
            b = b.model
            e = e.model
        } else {   
            b = cbind(b, b.model)
            e = c(e, e.model )
        }
        
        ndf.nominal[i] = ndf.model
        
    }
    
    es.resid = eigen(d.resid, symmetric=TRUE) 
    use = ( abs(es.resid$values)>tol.d )
    
    ndf.nominal[n.vars+1] = sum(use)
    b = cbind(b, es.resid$vectors[, use])
    e = c(e, es.resid$values[use])
    
    #---------------------------------------------------------
    # fit LDM to b
    #---------------------------------------------------------
    
    if (tolower(decomposition) == tolower("orthogonal")) {
        ortho.res.freq = ortho.decomp( x=x.freq, b=b, x.svd=x.freq.svd,
                                   tol.conv=tol.conv, 
                                   n.iter.max=n.iter.max, verbose=FALSE )
    
        ortho.res.tran = ortho.decomp( x=x.tran, b=b, x.svd=x.tran.svd,
                                   tol.conv=tol.conv, 
                                   n.iter.max=n.iter.max, verbose=FALSE )
        
        d.freq = ortho.res.freq$d
        v.freq = ortho.res.freq$v
        wt.freq = d.freq * t( v.freq )
        
        d.tran = ortho.res.tran$d
        v.tran = ortho.res.tran$v
        wt.tran = d.tran * t( v.tran )
        
    } else if (tolower(decomposition) == tolower("non-orthogonal")) {

        wt.freq = t(b) %*% x.freq
        d.freq = sqrt(apply(wt.freq^2, 1, sum))
        v.freq = t(wt.freq/d.freq)
        
        wt.tran = t(b) %*% x.tran
        d.tran = sqrt(apply(wt.tran^2, 1, sum))
        v.tran = t(wt.tran/d.tran)
    } 
    
    #-------------------------------------------------
    # calculate test statistics for each hypothesis
    #-------------------------------------------------
    
    global.freq = rep(0, n.vars)
    global.tran = rep(0, n.vars)
    
    n.otu = dim(x.freq)[2]
    otu.freq = mat.or.vec(n.vars, n.otu)
    otu.tran = mat.or.vec(n.vars, n.otu)    
    
    up = 0
    
    for (j in 1:n.vars) 
    {
        low = up + 1
        up = up + ndf.nominal[j]
        
        global.freq[j] = sum( d.freq[low:up]^2 )
        global.tran[j] = sum( d.tran[low:up]^2 )
        otu.freq[j, ] = colSums(wt.freq[low:up, , drop=FALSE]^2)
        otu.tran[j, ] = colSums(wt.tran[low:up, , drop=FALSE]^2)
    }
    
    details=list()
    
    if (verbose) {
        details$b=b
        details$x.freq=x.freq
        details$d.freq=d.freq
        details$v.freq=v.freq
        details$x.tran=x.tran
        details$d.tran=d.tran
        details$v.tran=v.tran
    }
    
    res = list( global.freq=global.freq, otu.freq=otu.freq, 
                global.tran=global.tran, otu.tran=otu.tran, 
                details=details)
    
    return(res)
    
} # fit.ldm.internal


fit.m1.permanova = function( d.gower, model) {
    
    tol.d=10^-8
    
    svd.var = svd(model$m1)  
    use = (svd.var$d>tol.d)    
    
    hat.matrix = svd.var$u[, use] %*% t( svd.var$u[, use] )
    hat.matrix.bar = diag(dim(hat.matrix)[1])  - hat.matrix
    
    d.resid = hat.matrix.bar %*% d.gower
    d.resid = d.resid %*% hat.matrix.bar
    
    m1.fit = list(d.resid=d.resid)
    return(m1.fit)
    
} # fit.m1.permanova


fit.permanova = function( d.gower, model) {
    
    #--------------------------------------------------------------------------
    # construct e 
    # from each set of variables in the list vars
    #--------------------------------------------------------------------------
    
    tol.d=10^-8
    
    n.vars = length(model$index)
    n.obs = model$n.obs
    
    d.resid = d.gower
    
    for (i in 1:n.vars) 
    {
        
        if (model$adjust.for.confounders==TRUE) {
            var = cbind(model$m1, model$m[,1:model$index[i]])
        } else {
            var = model$m[,1:model$index[i]]
        }
        
        svd.var = svd(var)
        
        use = (svd.var$d>tol.d)    
        hat.matrix = svd.var$u[, use] %*% t( svd.var$u[, use] )
        
        d.model = hat.matrix %*% d.resid
        d.model = d.model %*% hat.matrix
        
        es.model = eigen(d.model, symmetric=TRUE) 
        
        use = ( abs(es.model$values)>tol.d )
        e.model = es.model$values[use]
        ndf.model = sum( use )
        
        #---------------------
        # d.resid
        #---------------------
        
        hat.matrix.bar = diag(dim(hat.matrix)[1]) - hat.matrix
        d.resid = hat.matrix.bar %*% d.resid
        d.resid = d.resid %*% hat.matrix.bar
        
        #-----------------------------
        # save
        #-----------------------------    
        
        if (i==1) {
            e = e.model
            ndf.nominal = ndf.model
        } else {   
            e = c(e, e.model)
            ndf.nominal = c(ndf.nominal, ndf.model)
        }
        
    }
    
    es.resid = eigen(d.resid, symmetric=TRUE) 
    
    use = ( abs(es.resid$values)>tol.d )
    e = c(e, es.resid$values[use])
    ndf.nominal = c(ndf.nominal, sum(use))
    
    #-------------------------------------------------
    # calculate test statistics for each hypothesis
    #-------------------------------------------------
    
    permanova.denom = sum(es.resid$values[use])
    permanova = rep(0, n.vars)
    
    up = 0
    for (j in 1:n.vars) 
    {
        low = up + 1
        up = up + ndf.nominal[j]
        
        permanova[j] = sum( e[low:up] )/ndf.nominal[j]/permanova.denom*(n.obs-1-ndf.nominal[j])
    }
    
    res = list( permanova=permanova)
    return(res)
    
} # fit.permanova

calculate.dist <- function(otu.table, tree=NULL, dist.type="Bray-Curtis") {
    
    if (min(rowSums(otu.table))==0) {
        stop("samples with a total of zero reads should be removed")
    }
    
    freq.table <- t( scale( t(otu.table), center=FALSE, scale=rowSums(otu.table) ) )
    
    lower.dist.type = tolower(dist.type)
    
    if (lower.dist.type==tolower("Bray-Curtis")) {
        # Bray-Curtis (normalization by proportion)
        dist <- 0.5*as.matrix( dist( x=freq.table, method='manhattan') )
    } else if (lower.dist.type==tolower("wt-UniFrac")) {
        # weighted UniFrac dist
        set.seed(123) # ??? improve?
        otu.rff <- Rarefy(otu.table)$otu.tab.rff
        unifracs <- GUniFrac(otu.rff, tree, alpha=c(1))$unifrac
        dist <- unifracs[,,"d_1"]
    } else if (lower.dist.type==tolower("Hellinger")) {
        # Euclidean/Hellinger
        dist <- 0.5*as.matrix( dist( x=sqrt(freq.table), method='euclidean') )
    } else if (lower.dist.type==tolower("uw-UniFrac")) {
        # unweighted UniFrac dist
        set.seed(123)
        otu.rff <- Rarefy(otu.table)$otu.tab.rff
        unifracs <- GUniFrac(otu.rff, tree, alpha=c(1))$unifrac
        dist <- unifracs[,,"d_UW"]
    } else if (lower.dist.type==tolower("Jaccard")) {
        # unweighted Jaccard dist
        otu.rff <- Rarefy(otu.table)$otu.tab.rff
        dist <- vegdist(otu.rff, method="jaccard", binary=TRUE)
    }
    
    dist <- as.matrix(dist)
    
    return(dist)
    
} # calculate.dist

#' Orthogonal decomposition
#' 
#' This function performs the (approximate) orthogonal decomposition which gives
#' a duality between principle components of a data matrix \code{x} and linear
#' combinations of a set of orthogonal directions \code{b} (e.g., principle
#' components of a distance matrix) Specifically, it finds the matrix \code{v}
#' with orthogonal columns and the diagonal matrix \code{D} with positive
#' diagonal entries that minimize \code{||x - b D v^T||^2}.
#' 
#' 
#' @param x the \code{n.obs} by \code{n.otu} data matrix, possibly scaled by the
#'   library size (converting counts into frequencies) and/or centered in the
#'   columns (as if for a Principal Coordinates Analysis).
#' @param b the \code{n.obs} by \code{n.vec} matrix with orthonormal columns,
#'   e.g., the matrix having columns given by all eigenvectors of a distance
#'   matrix.
#' @param x.svd the result object of \code{svd(x)}. The default is NULL and then
#'   \code{svd(x)} is performed internally.
#' @param tol.conv the stopping criterion for the tandem algorithm. The
#'   algorithm is stopped when the largest change in any element of \code{d} or
#'   \code{v} is smaller than \code{tol.conv}. The default is 10e-08.
#' @param n.iter.max the maximum number of iterations to take before stopping.
#'   The default is 1000.
#' @param verbose a logical variable indicating whether to return detailed
#'   information on the decomposition. The default is FALSE.
#' @return a list consisting of \item{v}{an \code{n.otu} by \code{n.vec} matrix
#'   having orthonormal columns that give the directions (OTU loading) in
#'   feature space.} \item{d}{a vector of length \code{n.vec} with the diagonal
#'   elements of \code{D}.} \item{n.iter}{the number of iterations used to
#'   converge.} \item{fail}{equals 0 if the algorithm converged to the desired
#'   tolerance or 1 if not.} \item{b.hat}{an \code{n.obs} by \code{n.vec} matrix
#'   containing the estimated values of \code{b}. Returned when
#'   \code{verbose=TRUE}.} \item{r2}{R2 between columns of \code{b} and columns
#'   of \code{b.hat}. Returned when \code{verbose=TRUE}.} \item{converge}{A
#'   message informing whether the algorithm converged or not and the number of
#'   iterations used. Returned when \code{verbose=TRUE}.}
#' @keywords microbiome PCA ordination distance
#' @author Yi-Juan Hu <yijuan.hu@emory.edu>, Glen A. Satten <gas0@cdc.gov>
#' @export
#' @references Satten G.A., Tyx R., Rivera A.J., Stanfill S. (2017). Restoring
#'   the Duality between Principal Components of a Distance Matrix and Linear
#'   Combinations of Predictors, with Application to Studies of the Microbiome.
#'   PLoS ONE 12(1):e0168131. doi:10.1371/journal.pone.0168131.
#' @examples
#' res <- ortho.decomp( x=data.matrix, b=dist.vectors, verbose=TRUE )

ortho.decomp = function (x, b, x.svd=NULL, tol.conv=10^-8, 
                         n.iter.max=1000, verbose=FALSE) 
{
    tol.d=10^-8
    
    if (verbose) {
        if (dim(x)[1] != dim(b)[1]) 
            x <- t(x)
        if (dim(x)[1] != dim(b)[1]) 
            stop( 'numbers of observations mismatch between x and b' )
    }
    
    n.obs <- dim(x)[1]
    n.otu <- dim(x)[2]
    n.vec <- dim(b)[2]
    
    dim.reduction <- FALSE
    
    if (n.otu > n.obs) 
    {
        x.orig <- x
        if (is.null(x.svd)) x.svd <- svd(x)
        use.x <- ( abs(x.svd$d)>tol.d ) 
        
        x <- t( x.svd$d[use.x] * t(x.svd$u[,use.x]) ) #|| LE - BDQ ||
        
        n.otu <- dim(x)[2]
        dim.reduction <- TRUE
    }
    
    #------------------------------
    # Tandem algorithm: d, v
    #------------------------------
    
    xtb <- t(x) %*% b
    
    n.iter <- 0
    fail <- 0
    diff <- 1
    
    d <- rep(1, n.vec)
    v <- mat.or.vec(n.otu, n.vec)
    x.hat <- x
    
    while ((diff > tol.conv) & (n.iter <= n.iter.max))       
    {
        # Solving for v given d
        
        XBD <- t(t(xtb) * d)
        XBD.svd <- svd(XBD)  # X'BD=WSZ'
        use <- (abs(XBD.svd$d) > tol.d)
        v <- XBD.svd$u[, use] %*% t(XBD.svd$v[, use]) # V=WZ'
        
        # Solving for d given v
        
        d <- colSums(xtb * v)  # diagonal elements of V'X'B
        
        # flip signs if d < 0
        
        sign.flip <- (d < 0)                          
        v[, sign.flip] <- -v[, sign.flip]
        d <- abs(d)
        
        # diff
        
        new.x.hat <- t(t(b) * d)
        new.x.hat <- new.x.hat %*% t(v)
        diff <- max(abs(new.x.hat - x.hat))
        x.hat <- new.x.hat        
        
        n.iter <- n.iter + 1
        
    }# end of Tandem algorithm
    
    if (n.iter > n.iter.max) {
        fail <- 1
    }
    
    #-----------------------
    # recover orig dim
    #-----------------------
    
    if (dim.reduction) {
        v <- x.svd$v[,use.x] %*% v  # V=RQ
    }
    
    
    b.hat <- NULL
    r2 <- NULL
    msg <- NULL
    
    if (verbose) {
        
        if (dim.reduction) {
            x.hat <- x.hat %*% t(x.svd$v[,use.x]) # X.hat=BDQ'R'
            x <- x.orig
        }
        
        #-----------------------
        # b.hat
        #-----------------------
        
        use.d <- (abs(d)>tol.d)
        dinv <- rep(0, n.vec)
        dinv[use.d] <- 1.0/d[use.d]
        b.hat <- x %*% v
        b.hat <- t(t(b.hat) * dinv) # X=BDV' => B=XVDinv
        
        #-----------------------
        # r2
        #-----------------------
        
        r2 <- rep(0, n.vec)
        for (k in 1:n.vec) {
            r2[k] <- sum(b[, k] * b.hat[, k])/sum(b.hat[, k]^2)
        }
        
        #-----------------------
        # message
        #-----------------------
        
        msg=ifelse( fail==0, 
                    paste("ortho.decomp converged in ", n.iter, " iterations"), 
                    paste("ortho.decomp failed to converge; final difference is ", diff))
        
        cat(sprintf("%s\n", msg))
    }
    
    results <- list(v = v, d = d, n.iter = n.iter, fail = fail, 
                    b.hat = b.hat, r2 = r2, converge = msg)
    return(results)
    
}# ortho.decomp


#' Adjust the data by covariates
#' 
#' This function produces adjusted distance matrix and otu table (if provided) 
#' after removing the effects of covariates (e.g., confounders). 
#' 
#' @param formula a symbolic description of the covariate model with the form \code{ ~ model}, 
#' where \code{model} is specified the same as for \code{lm} or \code{glm}. For example, 
#' \code{~ a + b} specifies a model with the main effects of covariates \code{a} and \code{b}, and 
#' \code{~ a*b}, equivalently \code{~ a + b + a:b}, specifies a model with the main effects of 
#' \code{a} and \code{b} as well as their interaction.
#' @param data an optional data frame, list or environment (or object coercible 
#' by as.data.frame to a data frame) containing the covariates. 
#' If not found in \code{data}, the covariates are taken from environment(formula), 
#' typically the environment from which \code{adjust.data.by.covariates} is called. 
#' The default is NULL.
#' @param otu.table the \code{n.obs} by \code{n.otu} matrix of read counts. 
#' If provided, the adjusted (column-centered) otu.table at both the frequency scale 
#' and arcsin-root-transformed frequency scale are outputted. If provided, 
#' it is also used for calculating the ditance matrix unless the distance matrix is directly 
#' imported through \code{dist}.
#' The default is NULL.
#' @param tree a phylogenetic tree. Used for calculating a
#'   phylogenetic-tree-based distance matrix. Not needed if the calculation of
#'   requested distance does not require a phylogenetic tree, or if the distance
#'   matrix is directly imported through \code{dist}. The default is NULL.
#' @param dist.type a string chosen from "Bray-Curtis", "Hellinger", 
#'   "wt-UniFrac" (case insensitive). Used for calculating the distance matrix. 
#'   Not needed if the distance matrix is directly imported through \code{dist}.
#'   The default is "Bray-Curtis".
#' @param dist a distance matrix. It
#'   can be the original distance matrix and then squared/centered by the
#'   default \code{square.dist=TRUE}/\code{center.dist=TRUE}; it can also be
#'   squared/centered a priori and used with
#'   \code{square.dist=FALSE}/\code{center.dist=FALSE}. If not provided, the
#'   distance matrix is calculated using \code{otu.table}, \code{tree}, and
#'   \code{dist.type}. The default is NULL.
#' @param square.dist a logical variable indicating whether to square the 
#'   distance matrix. The default is TRUE.
#' @param center.dist a logical variable indicating whether to center the 
#'   distance matrix as described by Gower (1966). The default is TRUE.
#' @param scale.otu.table a logical variable indicating whether to scale (i.e., 
#'   divide) the rows of the OTU table by the library size, so the read counts 
#'   becomes relative frequencies. The default is TRUE. It can be set to FALSE 
#'   if the input \code{otu.table} has been scaled a priori.
#' @param center.otu.table a logical variable indicating whether to center the 
#'   columns of the OTU table. The OTU table should be centered if the distance 
#'   matrix has been centered. The default is TRUE.
#' @return a list consisting of 
#'   \item{adj.dist}{the (squared/centered) distance matrix
#'   after adjusting for covariates.}
#'   \item{x.freq}{the (column-centered) frequency scale data matrix after adjusting 
#'   for covariates.} 
#'   \item{x.tran}{the (column-centered) arcsin-root-transformed frequency scale 
#'   data matrix after adjusting for covariates.} 

#' @keywords microbiome PCA ordination distance
#' @author Yi-Juan Hu <yijuan.hu@emory.edu>, Glen A. Satten <gas0@cdc.gov>
#' @export
#' @examples
#' adj.data <- adjust.data.by.covariates(formula= ~ Sex + AntibioticUse, data=throat.meta,
#'                                       otu.table=throat.otu.tab, dist.type="Bray-Curtis")
#'
#' #-------------------------------------------------
#' # Use the adjusted distance matrix for ordination
#' #-------------------------------------------------
#'
#' PCs <- eigen(adj.data$adj.dist, symmetric=TRUE)
#' 
#' color = rep("blue", length(SmokingStatus))
#' w = which(SmokingStatus=="Smoker")
#' color[w] = "red"
#' 
#' plot(PCs$vectors[,1], PCs$vectors[,2], xlab="PC1", ylab="PC2", 
#'      col=color, main="Smokers vs. non-smokers")
#' legend(x="topleft", legend=c("smokers","non-smokers"), pch=c(21,21), 
#'        col=c("red","blue"), lty=0)



adjust.data.by.covariates = function(formula=NULL, data=NULL, 
                                     otu.table=NULL, tree=NULL, dist.type="Bray-Curtis", dist=NULL, 
                                     square.dist=TRUE, center.dist=TRUE, 
                                     scale.otu.table=TRUE, center.otu.table=TRUE) {
    
    #------------------------
    # check parameter values
    #------------------------ 
    
    if (! tolower(dist.type) %in% c(tolower("Bray-Curtis"), tolower("wt-UniFrac"), tolower("Hellinger"))) {
        stop("dist.type must be \"Bray-Curtis\", \"wt-UniFrac\", or \"Hellinger\"")
    }
    
    #------------------------
    # dist matrix
    #------------------------
    
    if (is.null(dist) & is.null(otu.table)) {
        stop( 'must specify one of dist and otu.table' )
    }
    
    if (!is.null(otu.table) & !is.null(dist)) {
        dist <- as.matrix(dist)
        if (dim(otu.table)[1] != dim(dist)[1]) stop( 'numbers of observations mismatch between otu.table and dist' )
    }
    
    if (is.null(dist)) {
        dist <- calculate.dist(dist.type=dist.type, otu.table=otu.table, tree=tree)
    }
    
    d.gower <- gower(d=dist, square=square.dist, center=center.dist)
    
    if (!is.null(otu.table)) {
        rs <- rowMeans(d.gower)
        cs <- colMeans(d.gower)
        if ( (any(abs(rs) > 10^-6) & center.otu.table==TRUE) | (all(abs(rs) <= 10^-6) & center.otu.table==FALSE) ) {
            stop( 'discordant centering of the OTU table and distance matrix' )
        }
    }
    
    #------------------------
    # covariates (e.g., confounders)
    #------------------------
    
    center.m1 = TRUE
    
    m1 = model.matrix(object=formula, data=data)
    
    if (center.m1) m1 = scale( m1, center=TRUE, scale=FALSE )
    
    if (dim(m1)[1] != dim(dist)[1]) stop( 'numbers of observations mismatch between covariates and dist' )
    
    
    #---------------------
    # calculate d.resid
    #---------------------
    
    tol.d=10^-8
    svd.m1 = svd(m1)
    use = (svd.m1$d>tol.d)
    
    hat.matrix = svd.m1$u[, use] %*% t( svd.m1$u[, use] )
    hat.matrix.bar = diag(dim(d.gower)[1]) - hat.matrix
    
    d.resid = hat.matrix.bar %*% d.gower
    d.resid = d.resid %*% hat.matrix.bar
    
    #---------------------
    # calculate adj.otu.table
    #---------------------
    
    if (min(rowSums(otu.table))==0) {
        stop("samples with a total of zero reads should be removed")
    }
    
    if (min(colSums(otu.table))==0) {
        stop("OTUs with zero reads in each sample should be removed")
    }
    
    # freq
    if (scale.otu.table) {
        freq.table <- t( scale( t(otu.table), center=FALSE, scale=rowSums(otu.table) ) )
    } else {
        freq.table <- otu.table
    }
    x.freq <- scale( freq.table, center=center.otu.table, scale=FALSE )
    
    # arcsine
    theta <- asin(sqrt(freq.table))
    x.tran <- scale( theta, center=center.otu.table, scale=FALSE)
    
    # b.model
    b.model = svd.m1$u[, use]
    
    # adj.otu.table
    x1.tilda.freq = t(b.model) %*% x.freq
    x1.tilda.freq = b.model %*% x1.tilda.freq
    x.freq = x.freq - x1.tilda.freq
    
    x1.tilda.tran = t(b.model) %*% x.tran
    x1.tilda.tran = b.model %*% x1.tilda.tran
    x.tran = x.tran - x1.tilda.tran
    
    res <- list(x.freq=x.freq,
                x.tran=x.tran,
                adj.dist=d.resid)
    
    return(res)
    
} # adjust.data.by.covariates


#' Testing hypotheses using a linear decomposition model (LDM)
#' 
#' This function allows you to simultaneously test the global hypothesis of any 
#' microbiome association and test individual OTU associations to give coherent 
#' results. It is capable of handling complex design features such as 
#' confounders, interactions, and clustered data.
#' 
#' 
#' The formula has the form \code{otu.table ~ sets of variables} or 
#' \code{otu.table | confounders ~ sets of variables}, where \code{otu.table} is
#' the OTU table with rows for samples and columns for OTUs and \code{sets of 
#' variables} can be multiple sets of variables (i.e., models) that we test 
#' individually and sequentially (i.e., after accounting for the previous sets 
#' of variables); each set of variables are enclosed in parentheses.
#' 
#' Suppose there is an otu table \code{y} and a data frame \code{metadata} that contains 4 covariates, 
#' \code{a}, \code{b}, \code{c} and \code{d}. 
#' Some valid formulas would be:
#' 
#' \code{form.call( y ~ a + b + c + d, data=metadata)} # no confounders, 4 models
#' 
#' \code{form.call( y ~ (a+b) + (c+d), data=metadata)} # no confounders, 2 models each having 
#' 2 variables; it is OKay even when \code{a} has values ‘A’ and ‘B’ while 
#' \code{b} has numeric values
#' 
#' \code{form.call( y | b ~ (a+c) + d, data=metadata)} # \code{b} is a confounder, model 1 is 
#' \code{(a+c)}, and model 2 is \code{d}
#' 
#' \code{form.call( y | b+c ~ a*d, data=metadata)}     # there are 2 confounder: \code{b} 
#' and \code{c}; there are 3 models: \code{a}, \code{d}, and \code{a x d = a:d}
#' 
#' \code{form.call( y | as.factor(b) ~ (a+d) + a:d, data=metadata)}       # now confounder 
#' \code{b} will be treated as a factor variable, model 1 will have the main 
#' effects \code{a} and \code{d}, and model 2 will have only the interaction 
#' between \code{a} and \code{d}
#' 
#' \code{form.call( y | as.factor(b) ~ (a) + (d) + (a:d), data=metadata)} # it is OKay to put 
#' paratheses around a single variable
#' 
#' LDM uses two sequential stopping criteria. For the global test, LDM uses the 
#' stopping rule of Besag and Clifford (1991), which stops permutation when a 
#' pre-specified minimum number (default 10) of rejections (i.e., the permutation 
#' statistic exceeded the observed test statistic) has been reached. For the 
#' OTU-specific tests, LDM uses the stopping rule of Sandve et al. (2011), 
#' which stops permutation when every OTU test has either reached the pre-specified 
#' number (default 10) of rejections or yielded a q-value that is below the 
#' nominal FDR level (default 0.1). As a convention, we call a test ``stopped" 
#' if the corresponding stopping criterion has been satisfied. Although all tests 
#' are always terminated if a pre-specified maximum number (default 5,000) of 
#' permutations have been generated, some tests may not have ``stopped" and 
#' hence their results have not achieved enough precision. 
#' 
#' 
#' @param formula a symbolic description of the 
#'   model to be fitted. The details of model specification are given under 
#'   `Details'.
#' @param data an optional data frame, list or environment (or object coercible 
#' by as.data.frame to a data frame) containing the variables of interest and 
#' confounding covariates (but not the OTU data, which are contained in a separate OTU table). 
#' If not found in \code{data}, the variables are taken from environment(formula), 
#' typically the environment from which \code{ldm} is called. The default is NULL.
#' @param tree a phylogenetic tree. Used for calculating a 
#'   phylogenetic-tree-based distance matrix. Not needed if the calculation of 
#'   the requested distance does not involve a phylogenetic tree, or if the 
#'   distance matrix is directly imported through \code{dist}.
#' @param dist.type a string chosen from "Bray-Curtis", "Hellinger", 
#'   "wt-UniFrac" (case insensitive). Used for calculating the distance matrix. 
#'   Not needed if the distance matrix is directly imported through \code{dist}.
#'   The default is "Bray-Curtis".
#' @param dist a distance matrix. It 
#'   can be the original distance matrix and then squared/centered by the 
#'   default \code{square.dist=TRUE}/\code{center.dist=TRUE}; it can also be 
#'   squared/centered a priori and used with 
#'   \code{square.dist=FALSE}/\code{center.dist=FALSE}. If not specified, the 
#'   distance matrix is calculated using \code{otu.table}, \code{tree}, and 
#'   \code{dist.type}.
#' @param is.clustered a logical value indicating whether the observations are 
#'   clustered. The default is FALSE.
#' @param cluster.id cluster identifiers. Unequal cluster size is allowed. The 
#'   default is NULL.
#' @param permute.within.clusters a logical value indicating whether to permute 
#'   data within clusters. The default is FALSE.
#' @param permute.between.clusters a logical value indicating whether to permute
#'   data between clusters. The default is FALSE.
#' @param n.perm.max the maximum number of permutations to run.
#'   The default is 5000.
#' @param n.rej.stop the minimum number of rejections (i.e., the permutation 
#'   statistic exceeds the observed statistic) to obtain before stopping. The 
#'   default is 10.
#' @param seed an integer seed for the random number generator in the 
#'   permutation procedure. The default is NULL (i.e., the seed will be 
#'   generated from the system clock, so it will appear random).
#' @param test.global a logical value indicating whether to perform the global 
#'   test. The default is TRUE.
#' @param test.otu a logical value indicating whether to perform the 
#'   OTU-specific tests. The default is TRUE.
#' @param fdr.nominal the nominal FDR value. The default is 0.1.
#' @param out.stat a logical value indicating whether to output intermediate 
#'   statistics to external files. The default is FALSE. This feature can be
#'   used to parallelize the permutation onto multiple machines, each starting
#'   with a different \code{seed}.
#' @param out.prefix a character string specifying the prefix (including the 
#'   path) of the external files that are used to save intermediate statistics. 
#' @param decomposition a character string that takes either the default value 
#'   \code{"orthogonal"} (corresponding to the LDM method) or 
#'   \code{"non-orthogonal"} (corresponding to the non-standard Redundancy 
#'   Analysis).
#' @param square.dist a logical variable indicating whether to square the 
#'   distance matrix. The default is TRUE.
#' @param center.dist a logical variable indicating whether to center the 
#'   distance matrix as described by Gower (1966). The default is TRUE.
#' @param scale.otu.table a logical variable indicating whether to scale (i.e., 
#'   divide) the rows of the OTU table by the library size, so the read counts 
#'   becomes relative frequencies. The default is TRUE. It can be set to FALSE 
#'   if the input \code{otu.table} has been scaled a priori.
#' @param center.otu.table a logical variable indicating whether to center the 
#'   columns of the OTU table. The OTU table should be centered if the distance 
#'   matrix has been centered. The default is TRUE.
#' @return a list consisting of 
#'   \item{b}{the matrix B as defined in Hu and Satten (2018)} 
#'   \item{adj.dist}{the (squared/centered) distance matrix
#'   after adjusting for confounders. If there are no confounders, it is simply
#'   the (squared/centered) distance matrix.} 
#'   \item{x.freq}{the (column-centered) frequency scale data matrix after adjusting 
#'   for confounders (if there are any).} 
#'   \item{d.freq}{a vector of the diagonal elements of \code{D} that minimizes 
#'   \code{||x.freq - b D v^T||^2}} 
#'   \item{v.freq}{the v matrix with orthonormal columns that minimizes 
#'   \code{||x.freq - b D v^T||^2}} 
#'   \item{x.tran}{the (column-centered) arcsin-root-transformed frequency scale 
#'   data matrix after adjusting for confounders (if there are any).} 
#'   \item{d.tran}{a vector of the diagonal elements of \code{D} that minimizes 
#'   \code{||x.tran - b D v^T||^2}} 
#'   \item{v.tran}{the v matrix with orthonormal columns that minimizes 
#'   \code{||x.tran - b D v^T||^2}}
#'   \item{VE.global.freq}{variance explained by each set of variables, based on
#'   the frequency-scale data} 
#'   \item{VE.global.tran}{variance explained by each set of variables, based on 
#'   the arcsine-root-transformed frequency scale data} 
#'   \item{VE.otu.freq}{contributions (i.e., factor loadings) from each OTU to 
#'   \code{VE.global.freq}} 
#'   \item{VE.otu.tran}{contributions (i.e., factor loadings) from each OTU to 
#'   \code{VE.global.tran}} 
#'   \item{p.global.freq}{p-values for the global test of each set of variables
#'   based on the frequency scale data} 
#'   \item{p.global.tran}{p-values for the global test of each set of variables
#'   based on the arcsin-root-transformed frequency scale data} 
#'   \item{p.global.omni}{p-values for the global test of each set of variables 
#'   based on the omnibus statistics, which are the minima of the p-values obtained 
#'   from the frequency scale and the arcsin-root-transformed frequency scale data 
#'   as the final test statistics, and use the corresponding minima from the 
#'   permuted data to simulate the null distributions} 
#'   \item{p.otu.freq}{p-values for the OTU-specific tests based on the 
#'   frequency scale data} 
#'   \item{p.otu.tran}{p-values for the OTU-specific tests based on the 
#'   arcsin-root-transformed frequency scale data} 
#'   \item{p.otu.omni}{p-values for the OTU-specific tests based on the 
#'   omnibus statistics} 
#'   \item{q.otu.freq}{q-values (i.e., FDR-adjusted p-values) 
#'   for the OTU-specific tests based on the frequency scale data} 
#'   \item{q.otu.tran}{q-values for the OTU-specific tests based on 
#'   the arcsin-root-transformed frequency scale data} 
#'   \item{q.otu.omni}{q-values for the OTU-specific tests based on the 
#'   omnibus statistics} 
#'   \item{n.perm.completed}{number of permutations completed} 
#'   \item{global.tests.stopped}{a logical value indicating whether the 
#'   stopping criterion has been met by all global tests} 
#'   \item{otu.tests.stopped}{a logical value indicating whether the 
#'   stopping criterion has been met by all OTU-specific tests}
#' @keywords microbiome
#' @author Yi-Juan Hu <yijuan.hu@emory.edu>, Glen A. Satten <gas0@cdc.gov>
#' @export
#' @references Hu YJ, Satten GA. (2017) Testing hypotheses about microbiome 
#'   using an ordination-based linear decomposition model. 
#'   bioRXiv:doi.org/10.1101/229831.
#' @examples
#' 
#'#-----------------------------------------------
#'# fit only
#'#-----------------------------------------------
#'fit <- ldm(formula=throat.otu.tab | (Sex+AntibioticUse) ~ (SmokingStatus+PackYears)+Age, 
#'           data=throat.meta, 
#'           dist.type="Bray-Curtis", 
#'           n.perm.max=0)
#'
#'#-----------------------------------------------
#'# test the global hypothese only
#'#-----------------------------------------------
#'res1.ldm <- ldm(formula=throat.otu.tab | (Sex+AntibioticUse) ~ (SmokingStatus+PackYears)+Age, 
#'                data=throat.meta, 
#'                dist.type="Bray-Curtis", 
#'                test.global=TRUE, test.otu=FALSE,
#'                n.perm.max=10000, seed=1)
#'                     
#'#----------------------------------------------------
#'# test both the global hypothese and individual OTUs
#'#----------------------------------------------------
#'res2.ldm <- ldm(formula=throat.otu.tab | (Sex+AntibioticUse) ~ (SmokingStatus+PackYears)+Age, 
#'                data=throat.meta, 
#'                dist.type="Bray-Curtis", 
#'                test.global=TRUE, test.otu=TRUE, fdr.nominal=0.1,
#'                n.perm.max=85600, seed=100)
#' 
#'#----------------------------------------------------
#'# test both the global hypothese and individual OTUs
#'# save the intermediate statistis to external files
#'# parallelize permutation
#'#----------------------------------------------------
#'dir.create("tmp_files") # create the directory for external files to store intermediate statistics
#'
#'ext.seed = 1:9
#'res.tmp <- ldm(formula=throat.otu.tab | (Sex+AntibioticUse) ~ (SmokingStatus+PackYears)+Age, 
#'               data=throat.meta, 
#'               dist.type="Bray-Curtis", 
#'               test.global=TRUE, test.otu=TRUE, fdr.nominal=0.1,
#'               out.stat=TRUE, out.prefix="tmp_files/tmp", 
#'               n.perm.max=10000, seed=ext.seed[1]) 
#'
#'res3.ldm <- test.ldm.fromfile(in.prefix="tmp_files/tmp", n.perm.available=40000,
#'                              test.global=TRUE, test.otu=TRUE, fdr.nominal=0.1)  
#'
#'#----------------------------------------------------
#'# clustered data
#'#----------------------------------------------------
#'res4.ldm <- ldm(formula=sim.otu.tab | (C1+C2) ~ Y, 
#'                data=sim.meta,
#'                dist.type="Bray-Curtis", 
#'                is.clustered=TRUE, cluster.id=ID, permute.between.clusters=TRUE,
#'                test.global=TRUE, test.otu=TRUE, fdr.nominal=0.1,
#'                n.perm.max=5000, seed=1)


ldm = function( formula, data=NULL, tree=NULL, dist.type="Bray-Curtis", dist=NULL, 
                     is.clustered=FALSE, cluster.id=NULL, 
                     permute.within.clusters=FALSE, permute.between.clusters=FALSE,
                     n.perm.max=5000, n.rej.stop=10, seed=NULL, 
                     test.global=TRUE, test.otu=TRUE, fdr.nominal=0.1,
                     out.stat=FALSE, out.prefix=NULL, 
                     decomposition="orthogonal", 
                     square.dist=TRUE, center.dist=TRUE, 
                     scale.otu.table=TRUE, center.otu.table=TRUE) {  
    
    n.perm.init=100
    n.perm.step=100
    n.signif = 4
    n.iter.max=1000
    n.stable.max=10 # n.perm.step*n.stable.max=1000
    
    #------------------------
    # check parameter values
    #------------------------ 
    
    if (! tolower(dist.type) %in% c(tolower("Bray-Curtis"), tolower("wt-UniFrac"), tolower("Hellinger"))) {
        stop("dist.type must be \"Bray-Curtis\", \"wt-UniFrac\", or \"Hellinger\"")
    }
    
    if (! tolower(decomposition) %in% c(tolower("orthogonal"), tolower("non-orthogonal"))) {
        stop("decomposition must be \"orthogonal\" or \"non-orthogonal\"")
    }
    
    #------------------------
    # vars
    #------------------------
    
    dm = form.call(formula, data=data)
    
    otu.table = dm$otu.table
    
    adjust.for.confounders = !is.null(dm$conf)
    
    if (adjust.for.confounders) {
        vars <- list(dm$conf)
        for (j in 1:length(dm$model)) {
            vars[[j+1]] <- dm$model[[j]]
        }
    } else {
        vars <- dm$model
    }
    
    #------------------------
    # dist matrix
    #------------------------
    
    if (!is.null(dist)) {
        dist <- as.matrix(dist)
        if (dim(otu.table)[1] != dim(dist)[1]) stop( 'numbers of observations mismatch between otu.table and dist' )
    }
    
    if (is.null(dist)) {
        dist <- calculate.dist(dist.type=dist.type, otu.table=otu.table, tree=tree)
    }
    
    d.gower <- gower(d=dist, square=square.dist, center=center.dist)
    
    rs <- rowMeans(d.gower)
    cs <- colMeans(d.gower)
    if ( (any(abs(rs) > 10^-6) & center.otu.table==TRUE) | (all(abs(rs) <= 10^-6) & center.otu.table==FALSE) ) {
        stop( 'discordant centering of the OTU table and distance matrix' )
    }
    
    #------------------------
    # data matrix X
    #------------------------
    
    if (min(rowSums(otu.table))==0) {
        stop("samples with a total of zero reads should be removed")
    }
    
    if (min(colSums(otu.table))==0) {
        stop("OTUs with zero reads in each sample should be removed")
    }
    
    # freq
    if (scale.otu.table) {
        freq.table <- t( scale( t(otu.table), center=FALSE, scale=rowSums(otu.table) ) )
    } else {
        freq.table <- otu.table
    }
    x.freq <- scale( freq.table, center=center.otu.table, scale=FALSE )
    
    # arcsine
    theta <- asin(sqrt(freq.table))
    x.tran <- scale( theta, center=center.otu.table, scale=FALSE)
    
    #------------------------
    # setup model
    #------------------------
    
    center.vars=TRUE
    
    model <- setup.model( vars=vars, adjust.for.confounders=adjust.for.confounders, 
                          center.vars=center.vars, 
                          is.clustered=is.clustered, cluster.id=cluster.id, 
                          permute.within.clusters=permute.within.clusters, 
                          permute.between.clusters=permute.between.clusters)
    
    
    if (model$adjust.for.confounders == TRUE) {
        m1.fit = fit.m1( d.gower=d.gower, x.freq=x.freq, x.tran=x.tran, model=model) 
      
        d.gower=m1.fit$d.resid
        x.freq=m1.fit$x.freq.proj
        x.freq.svd=m1.fit$x.freq.svd
        x.tran=m1.fit$x.tran.proj
        x.tran.svd=m1.fit$x.tran.svd
        
    } else {
        x.freq.svd <- svd(x.freq)
        x.tran.svd <- svd(x.tran)
    }
  
    data.test.statistic = fit.ldm.internal( d.gower=d.gower, 
                                            x.freq=x.freq, x.freq.svd=x.freq.svd,
                                            x.tran=x.tran, x.tran.svd=x.tran.svd, 
                                            model=model, 
                                            n.iter.max=n.iter.max, tol.conv=10^-6,
                                            decomposition=decomposition,
                                            verbose=TRUE) 
    
    
    otu.names <- colnames(otu.table)
    colnames(data.test.statistic$otu.freq) <- otu.names
    colnames(data.test.statistic$otu.tran) <- otu.names
    
    
    p.global.freq = NULL
    p.global.tran = NULL
    p.global.omni = NULL
    
    p.otu.freq = NULL
    p.otu.tran = NULL
    p.otu.omni = NULL
    
    q.otu.freq = NULL
    q.otu.tran = NULL
    q.otu.omni = NULL
    
    n.perm.completed = NULL
    global.tests.stopped = NULL
    otu.tests.stopped    = NULL
    
    
    if (n.perm.max > 0) {
        
        n.otus = dim(x.freq)[2]
        n.vars = length(model$index)
      
        if (test.global) 
        {
            global.freq = data.test.statistic$global.freq
            global.tran = data.test.statistic$global.tran
            
            if (out.stat) {
                file.global.freq = paste(out.prefix,"_global.freq.txt",sep="")
                file.global.tran = paste(out.prefix,"_global.tran.txt",sep="")
            
                if (!file.exists(file.global.freq)) write.table(t(global.freq), file.global.freq, row.names=FALSE, col.names=FALSE, quote=FALSE)
                if (!file.exists(file.global.tran)) write.table(t(global.tran), file.global.tran, row.names=FALSE, col.names=FALSE, quote=FALSE)
    
                file.global.freq.perm = paste(out.prefix,"_global.freq.perm_seed",seed,".txt",sep="")
                file.global.tran.perm = paste(out.prefix,"_global.tran.perm_seed",seed,".txt",sep="")
                
            } else {
            
                global.freq.perm = matrix(NA, n.vars, n.perm.max)
                global.tran.perm = matrix(NA, n.vars, n.perm.max)
                
                n.global.freq = 0
                n.global.tran = 0
            }
        }
      
        if (test.otu)
        {
            otu.freq = data.test.statistic$otu.freq[, , drop=FALSE] 
            otu.tran = data.test.statistic$otu.tran[, , drop=FALSE]
            
            if (out.stat) {
                file.otu.freq = paste(out.prefix,"_otu.freq.txt",sep="")
                file.otu.tran = paste(out.prefix,"_otu.tran.txt",sep="")
                
                if (!file.exists(file.otu.freq)) write.table(otu.freq, file.otu.freq, row.names=FALSE, col.names=TRUE, quote=FALSE)
                if (!file.exists(file.otu.tran)) write.table(otu.tran, file.otu.tran, row.names=FALSE, col.names=TRUE, quote=FALSE)
                
                file.otu.freq.perm = paste(out.prefix,"_otu.freq.perm_var",1:n.vars,"_seed",seed,".txt",sep="")
                file.otu.tran.perm = paste(out.prefix,"_otu.tran.perm_var",1:n.vars,"_seed",seed,".txt",sep="")
                
            } else {
            
                otu.freq.perm = array(NA, c(n.vars, n.otus, n.perm.max))
                otu.tran.perm = array(NA, c(n.vars, n.otus, n.perm.max))
                
                n.otu.freq = 0
                n.otu.tran = 0
    
                Aset.freq = matrix(TRUE, n.vars, n.otus)
                Aset.tran = matrix(TRUE, n.vars, n.otus)
                Aset.omni = matrix(TRUE, n.vars, n.otus)
                
                Bset.freq = matrix(FALSE, n.vars, n.otus)
                Bset.tran = matrix(FALSE, n.vars, n.otus)
                Bset.omni = matrix(FALSE, n.vars, n.otus)
                
                p.otu.freq = matrix(NA, n.vars, n.otus)
                p.otu.tran = matrix(NA, n.vars, n.otus)
                p.otu.omni = matrix(NA, n.vars, n.otus)
                
            }
        }
      
        tol.eq = 10^-16
      
        perm.model = model
        
        n.perm.completed = 0
        global.tests.stopped = FALSE
        otu.tests.stopped    = FALSE
        
        n.stable.freq = 0
        n.stable.tran = 0
        n.stable.omni = 0
      
        set.seed(seed)
        
        for (i.sim in 1:n.perm.max) {
        
            if (i.sim%%100==0) {
                cat("permutations:", i.sim, "\n")
            }
        
            n.perm.completed = n.perm.completed + 1
        
        
            if (!model$is.clustered) {
                perm = sample(model$n.obs, replace=FALSE)
            } else {
                
                if (model$permute.between.clusters) {
                    perm.cluster = sample(model$n.cluster, replace=FALSE)
                }
          
                low = 0
                up = 0
                perm = rep(0, model$n.obs)
          
                for (i in 1:model$n.cluster) {
                    low = up+1
                    up = up+model$cluster.size[i]
            
                    if (model$permute.between.clusters) {
                        # assume that the covariate matrix is exactly the same for each member of the cluster
                        perm[model$cluster.order[low:up]] = rep( model$cluster.paradigm[ perm.cluster[i] ], model$cluster.size[i] ) 
                    }
                    if (model$permute.within.clusters) {
                        perm[model$cluster.order[low:up]] = sample(model$cluster.order[low:up], replace=FALSE)
                    }
                    
                }
            }           
        
        
            perm.model$m = model$m[perm, , drop=FALSE]
            
            perm.res = fit.ldm.internal( d.gower=d.gower, 
                                    x.freq=x.freq, x.freq.svd=x.freq.svd, 
                                    x.tran=x.tran, x.tran.svd=x.tran.svd, 
                                    model=perm.model, 
                                    n.iter.max=n.iter.max, tol.conv=10^-6,
                                    decomposition=decomposition,
                                    verbose=FALSE)
        
        
            if (test.global & !global.tests.stopped) 
            {
                if (out.stat) {
                    write.table(t(perm.res$global.freq), file.global.freq.perm, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
                    write.table(t(perm.res$global.tran), file.global.tran.perm, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
                } else {
                    global.freq.perm[,i.sim] = perm.res$global.freq
                    global.tran.perm[,i.sim] = perm.res$global.tran
                    
                    n.global.freq = n.global.freq + (perm.res$global.freq > global.freq + tol.eq) + (perm.res$global.freq > global.freq - tol.eq)
                    n.global.tran = n.global.tran + (perm.res$global.tran > global.tran + tol.eq) + (perm.res$global.tran > global.tran - tol.eq)
                    
                }
            }
        
            if (test.otu & !otu.tests.stopped) 
            {
                if (out.stat) {
                    for (f in 1:n.vars) {
                        write.table(signif(t(perm.res$otu.freq[f,]),n.signif), file.otu.freq.perm[f], append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
                        write.table(signif(t(perm.res$otu.tran[f,]),n.signif), file.otu.tran.perm[f], append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
                    }
                } else {
                    otu.freq.perm[,,i.sim] = perm.res$otu.freq
                    otu.tran.perm[,,i.sim] = perm.res$otu.tran
                    
                    n.otu.freq = n.otu.freq + (perm.res$otu.freq>otu.freq+tol.eq) + (perm.res$otu.freq>otu.freq-tol.eq)
                    n.otu.tran = n.otu.tran + (perm.res$otu.tran>otu.tran+tol.eq) + (perm.res$otu.tran>otu.tran-tol.eq)
                }
            }
        
        
            if (!out.stat & ((i.sim >= n.perm.init & i.sim %% n.perm.step == 0) | (i.sim==n.perm.max))) {
          
                if (test.global & !global.tests.stopped) 
                {
                    ################
                    # test global
                    ################
    
                    if ((all(n.global.freq >= n.rej.stop*2) & all(n.global.tran >= n.rej.stop*2)) | (i.sim==n.perm.max) )
                    {
                        
                        p.global.freq = ifelse((n.global.freq >= n.rej.stop*2), 0.5*n.global.freq/n.perm.completed, (0.5*n.global.freq+1)/(n.perm.completed+1))
                        p.global.tran = ifelse((n.global.tran >= n.rej.stop*2), 0.5*n.global.tran/n.perm.completed, (0.5*n.global.tran+1)/(n.perm.completed+1))
                        
                        p.global.freq.tmp <- 0.5*n.global.freq/n.perm.completed
                        p.global.tran.tmp <- 0.5*n.global.tran/n.perm.completed
              
                        pmin.global.omni <- pmin(p.global.freq.tmp, p.global.tran.tmp)
                        pnull.global.freq <- (n.perm.completed + 0.5 - apply(global.freq.perm[,1:n.perm.completed,drop=FALSE], 1 ,rank))/n.perm.completed
                        pnull.global.tran <- (n.perm.completed + 0.5 - apply(global.tran.perm[,1:n.perm.completed,drop=FALSE], 1 ,rank))/n.perm.completed
                        pnullmin.global.omni <- pmin(pnull.global.freq, pnull.global.tran)
              
                        n.global.omni <- apply((t(pnullmin.global.omni) < pmin.global.omni + tol.eq), 1, sum) + 0.5*apply((abs(t(pnullmin.global.omni) - pmin.global.omni) < tol.eq), 1, sum)  
                        p.global.omni = ifelse((n.global.omni >= n.rej.stop), n.global.omni/n.perm.completed, (n.global.omni+1)/(n.perm.completed+1))
              
                        if (all(n.global.omni >= n.rej.stop)) global.tests.stopped = TRUE
                    } 
                } 
          
                if (test.otu & !otu.tests.stopped) 
                {
            
                    ################
                    # test otu
                    ################
            
                    if (any(Aset.freq)) {
                        AtoB.freq <- Aset.freq & (n.otu.freq >= n.rej.stop*2)
                        Bset.freq <- Bset.freq | AtoB.freq
                        Aset.freq <- Aset.freq & !AtoB.freq
                        p.otu.freq[AtoB.freq] <- 0.5*n.otu.freq[AtoB.freq]/n.perm.completed
                        p.otu.freq[Aset.freq] <- (0.5*n.otu.freq[Aset.freq]+1)/(n.perm.completed+1) 
                    
                        q.otu.freq <- t(apply(p.otu.freq, 1, fdr.Sandev))
                    
                        Aset.freq.meet.criteria <- apply(((q.otu.freq < fdr.nominal) & Aset.freq) | (!Aset.freq), 1, all)
                        n.stable.freq <- ifelse(Aset.freq.meet.criteria, n.stable.freq + Aset.freq.meet.criteria, 0)
                        Aset.freq.rm.row <- (n.stable.freq >= n.stable.max)
                        Aset.freq[Aset.freq.rm.row,] = FALSE
                    }
            
                    if (any(Aset.tran)) {
                        AtoB.tran <- Aset.tran & (n.otu.tran >= n.rej.stop*2)
                        Bset.tran <- Bset.tran | AtoB.tran
                        Aset.tran <- Aset.tran & !AtoB.tran
                        p.otu.tran[AtoB.tran] <- 0.5*n.otu.tran[AtoB.tran]/n.perm.completed
                        p.otu.tran[Aset.tran] <- (0.5*n.otu.tran[Aset.tran]+1)/(n.perm.completed+1) 
                  
                        q.otu.tran <- t(apply(p.otu.tran, 1, fdr.Sandev))
                  
                        Aset.tran.meet.criteria <- apply(((q.otu.tran < fdr.nominal) & Aset.tran) | (!Aset.tran), 1, all)
                        n.stable.tran <- ifelse(Aset.tran.meet.criteria, n.stable.tran + Aset.tran.meet.criteria, 0)
                        Aset.tran.rm.row <- (n.stable.tran >= n.stable.max)
                        Aset.tran[Aset.tran.rm.row,] = FALSE
                    }
    
                    p.otu.freq.tmp <- 0.5*n.otu.freq/n.perm.completed
                    p.otu.tran.tmp <- 0.5*n.otu.tran/n.perm.completed
                
                    pmin.otu.omni <- pmin(p.otu.freq.tmp, p.otu.tran.tmp)
            
                    pnull.otu.freq <- (n.perm.completed + 0.5 - apply(otu.freq.perm[,,1:n.perm.completed,drop=FALSE], c(1,2), rank))/n.perm.completed
                    pnull.otu.tran <- (n.perm.completed + 0.5 - apply(otu.tran.perm[,,1:n.perm.completed,drop=FALSE], c(1,2), rank))/n.perm.completed
                    pnullmin.otu.omni <- pmin(pnull.otu.freq, pnull.otu.tran)
                    pnullmin.otu.omni <- aperm(pnullmin.otu.omni, c(2,3,1))
            
                    n.otu.omni <- apply( (pnullmin.otu.omni < rep(pmin.otu.omni, n.perm.completed) - tol.eq), c(1,2), sum) + 0.5 * apply( (abs(pnullmin.otu.omni - rep(pmin.otu.omni, n.perm.completed)) < tol.eq), c(1,2), sum)
            
                    if (any(Aset.omni)) {
                        AtoB.omni <- Aset.omni & (n.otu.omni >= n.rej.stop)
                        Bset.omni <- Bset.omni | AtoB.omni
                        Aset.omni <- Aset.omni & !AtoB.omni
                        p.otu.omni[AtoB.omni] <- n.otu.omni[AtoB.omni]/n.perm.completed
                        p.otu.omni[Aset.omni] <- (n.otu.omni[Aset.omni]+1)/(n.perm.completed+1) 
                  
                        q.otu.omni <- t(apply(p.otu.omni, 1, fdr.Sandev))
                  
                        Aset.omni.meet.criteria <- apply(((q.otu.omni < fdr.nominal) & Aset.omni) | (!Aset.omni), 1, all)
                        n.stable.omni <- ifelse(Aset.omni.meet.criteria, n.stable.omni + Aset.omni.meet.criteria, 0)
                        Aset.omni.rm.row <- (n.stable.omni >= n.stable.max)
                        Aset.omni[Aset.omni.rm.row,] = FALSE
                    }
                    
                    if (!any(Aset.freq) & !any(Aset.tran) & !any(Aset.omni)) otu.tests.stopped = TRUE 
            
                }# if test.otu
          
                if ( (test.global + test.otu == 2) & (global.tests.stopped + otu.tests.stopped) == 2) break
          
                if ( (test.global + test.otu == 1) & (global.tests.stopped + otu.tests.stopped) == 1) break
          
            }# check if stop early 
        
        }# permutation
        
    }# if (n.perm.max > 0)
    
    if (!is.null(p.otu.freq)) colnames(p.otu.freq) <- otu.names
    if (!is.null(p.otu.tran)) colnames(p.otu.tran) <- otu.names
    if (!is.null(p.otu.omni)) colnames(p.otu.omni) <- otu.names
    if (!is.null(q.otu.freq)) colnames(q.otu.freq) <- otu.names
    if (!is.null(q.otu.tran)) colnames(q.otu.tran) <- otu.names
    if (!is.null(q.otu.omni)) colnames(q.otu.omni) <- otu.names
    
    res = list( b=data.test.statistic$detail$b,
                    adj.dist=d.gower,
                    x.freq=x.freq,
                    v.freq=data.test.statistic$detail$v.freq,
                    d.freq=data.test.statistic$detail$d.freq,
                    x.tran=x.tran,
                    v.tran=data.test.statistic$detail$v.tran,
                    d.tran=data.test.statistic$detail$d.tran,
                    VE.global.freq=data.test.statistic$global.freq, 
                    VE.global.tran=data.test.statistic$global.tran, 
                    VE.otu.freq=data.test.statistic$otu.freq, 
                    VE.otu.tran=data.test.statistic$otu.tran,
                    p.global.freq=p.global.freq, 
                    p.global.tran=p.global.tran, 
                    p.global.omni=p.global.omni,
                    p.otu.freq=p.otu.freq,
                    p.otu.tran=p.otu.tran,
                    p.otu.omni=p.otu.omni,
                    q.otu.freq=q.otu.freq,
                    q.otu.tran=q.otu.tran,
                    q.otu.omni=q.otu.omni,
                    n.perm.completed=n.perm.completed,
                    global.tests.stopped=global.tests.stopped,
                    otu.tests.stopped=otu.tests.stopped)
      
    return(res)
    
} # ldm  


#' Testing hypotheses based on saved intermediate statistics
#'
#' This function summarizes the intermediate statistics saved in external files 
#' into p-values and q-values.
#' 
#' 
#' @param in.prefix a character string specifying the prefix (including the path) of 
#' the external files that save the intermediate statistics. 
#' @param n.perm.available a (approximate) estimate of the total number of permutations 
#' in the directory specified in \code{in.prefix}. 
#' Used for managing memory when reading large files.
#' @param n.rej.stop the minimum number of rejections (i.e., the permutation 
#'   statistic exceeds the observed statistic) to obtain before stopping. The 
#'   default is 10.
#' @param test.global a logical value indicating whether to perform the global 
#'   test. The default is TRUE.
#' @param test.otu a logical value indicating whether to perform the 
#'   OTU-specific tests. The default is TRUE.
#' @param fdr.nominal the nominal FDR value. 
#' The default is 0.1.
#' @return a list consisting of
#'   \item{p.global.freq}{p-values for the global test of each set of variables
#'   based on the frequency scale data} 
#'   \item{p.global.tran}{p-values for the global test of each set of variables
#'   based on the arcsin-root-transformed frequency scale data} 
#'   \item{p.global.omni}{p-values for the global test of each set of variables 
#'   based on the omnibus statistics, which are the minima of the p-values obtained 
#'   from the frequency scale and the arcsin-root-transformed frequency scale data 
#'   as the final test statistics, and use the corresponding minima from the permuted 
#'   data to simulate the null distributions} 
#'   \item{p.otu.freq}{p-values for the OTU-specific tests based on the 
#'   frequency scale data} 
#'   \item{p.otu.tran}{p-values for the OTU-specific tests based on the 
#'   arcsin-root-transformed frequency scale data} 
#'   \item{p.otu.omni}{p-values for the OTU-specific tests based on the 
#'   omnibus statistics} 
#'   \item{q.otu.freq}{q-values (i.e., FDR-adjusted p-values) 
#'   for the OTU-specific tests based on the frequency scale data} 
#'   \item{q.otu.tran}{q-values for the OTU-specific tests based on 
#'   the arcsin-root-transformed frequency scale data} 
#'   \item{q.otu.omni}{q-values for the OTU-specific tests based on the 
#'   omnibus statistics} 
#'   \item{n.perm.completed}{number of permutations completed} 
#'   \item{global.tests.stopped}{a logical value indicating whether the 
#'   stopping criterion has been met by all global tests} 
#'   \item{otu.tests.stopped}{a logical value indicating whether the 
#'   stopping criterion has been met by all OTU-specific tests}
#' @keywords microbiome
#' @author Yi-Juan Hu <yijuan.hu@emory.edu>, Glen A. Satten <gas0@cdc.gov>
#' @export
#' @examples
#'#----------------------------------------------------
#'# test both the global hypothese and individual OTUs
#'# save the intermediate statistis to external files
#'# parallelize permutation
#'#----------------------------------------------------
#'dir.create("tmp_files") # create the directory 
#'
#'ext.seed = 1:9
#'res.tmp <- ldm(formula=throat.otu.tab | (Sex+AntibioticUse) ~ (SmokingStatus+PackYears)+Age, 
#'                    data=throat.meta, dist.type="Bray-Curtis", 
#'                    test.global=TRUE, test.otu=TRUE, fdr.nominal=0.1,
#'                    out.stat=TRUE, out.prefix="tmp_files/tmp", 
#'                    n.perm.max=10000, seed=ext.seed[1]) 
#'
#'res3.ldm <- test.ldm.fromfile(in.prefix="tmp_files/tmp", n.perm.available=40000,
#'                              test.global=TRUE, test.otu=TRUE, fdr.nominal=0.1)  
#'in.prefix="tmp", n.perm.available=90000) 

# ??? remove n.perm.available
test.ldm.fromfile = function(in.prefix=NULL, n.perm.available, n.rej.stop=10, 
                     test.global=TRUE, test.otu=TRUE, fdr.nominal=0.1) { 
    
    tol.eq=10^-16
    
    p.global.freq = NULL
    p.global.tran = NULL
    p.global.omni = NULL
    
    p.otu.freq = NULL
    p.otu.tran = NULL
    p.otu.omni = NULL
    
    q.otu.freq = NULL
    q.otu.tran = NULL
    q.otu.omni = NULL
    
    global.tests.stopped = FALSE
    otu.tests.stopped = FALSE
    
    n.perm.completed = NULL

    # check existence of specified external files
    
    in.prefix.dir = dirname(in.prefix)
    in.prefix.file = basename(in.prefix)

    files <- list.files(path=in.prefix.dir, pattern=in.prefix.file, full.names=TRUE)
    if (length(files)==0) {
        stop("cannot find files with specified prefix")
    }
    
    ################################
    # test global
    ################################
    
    if (test.global) 
    {
        
        file.global.freq = paste(in.prefix,"_global.freq.txt",sep="")
        file.global.tran = paste(in.prefix,"_global.tran.txt",sep="")
        
        global.freq = as.numeric(unlist(read.table(file.global.freq, as.is=TRUE)))
        global.tran = as.numeric(unlist(read.table(file.global.tran, as.is=TRUE)))
        
        n.vars = length(global.freq)
        
        global.freq.perm = NULL
        global.tran.perm = NULL
        
        filelist.global.freq.perm = list.files(path=in.prefix.dir, pattern=paste(in.prefix.file, "_global.freq.perm_seed", sep=""), full.names=TRUE)
        filelist.global.tran.perm = list.files(path=in.prefix.dir, pattern=paste(in.prefix.file, "_global.tran.perm_seed", sep=""), full.names=TRUE)
        
        for (s in 1:length(filelist.global.freq.perm)) {
            global.freq.perm.tmp = read.table(filelist.global.freq.perm[s], colClasses = rep("numeric", n.vars))
            global.tran.perm.tmp = read.table(filelist.global.tran.perm[s], colClasses = rep("numeric", n.vars))

            global.freq.perm = rbind(global.freq.perm, global.freq.perm.tmp)
            global.tran.perm = rbind(global.tran.perm, global.tran.perm.tmp)
        }
        
        n.perm.completed = nrow(global.freq.perm) # must be the first file read in, because the length may have increased for the next file
        
        n.global.freq <- apply(t(global.freq.perm[1:n.perm.completed,]) > global.freq + tol.eq, 1, sum) + 0.5*apply(abs(t(global.freq.perm[1:n.perm.completed,]) - global.freq) < tol.eq, 1, sum)
        n.global.tran <- apply(t(global.tran.perm[1:n.perm.completed,]) > global.tran + tol.eq, 1, sum) + 0.5*apply(abs(t(global.tran.perm[1:n.perm.completed,]) - global.tran) < tol.eq, 1, sum)
        
        p.global.freq = ifelse((n.global.freq >= n.rej.stop), n.global.freq/n.perm.completed, (n.global.freq+1)/(n.perm.completed+1))
        p.global.tran = ifelse((n.global.tran >= n.rej.stop), n.global.tran/n.perm.completed, (n.global.tran+1)/(n.perm.completed+1))
        
        p.global.freq.tmp <- n.global.freq/n.perm.completed
        p.global.tran.tmp <- n.global.tran/n.perm.completed
        
        pmin.global.omni <- pmin(p.global.freq.tmp, p.global.tran.tmp)
        
        pnull.global.freq <- (n.perm.completed + 0.5 - apply(global.freq.perm[1:n.perm.completed,,drop=FALSE], 2 ,rank))/n.perm.completed
        pnull.global.tran <- (n.perm.completed + 0.5 - apply(global.tran.perm[1:n.perm.completed,,drop=FALSE], 2 ,rank))/n.perm.completed
        pnullmin.global.omni <- pmin(pnull.global.freq, pnull.global.tran)
        
        n.global.omni <- apply((t(pnullmin.global.omni) < pmin.global.omni - tol.eq), 1, sum) + 0.5*apply((abs(t(pnullmin.global.omni) - pmin.global.omni) < tol.eq), 1, sum)  
        p.global.omni = ifelse((n.global.omni >= n.rej.stop), n.global.omni/n.perm.completed, (n.global.omni+1)/(n.perm.completed+1))
        
        if (all(n.global.freq >= n.rej.stop) & all(n.global.tran >= n.rej.stop) & all(n.global.omni >= n.rej.stop)) global.tests.stopped = TRUE
        
    } # test.global
    
    
    ################################
    # test otu
    ################################
    
    if (test.otu) 
    {
        
        #n.perm.completed???
        
        otu.freq.full = read.table(list.files(path=in.prefix.dir, pattern=paste(in.prefix.file, "_otu.freq.txt", sep=""), full.names=TRUE), as.is=TRUE, header=TRUE, check.names=FALSE)
        otu.names = colnames(otu.freq.full)
        n.otus = ncol(otu.freq.full)
        n.vars = nrow(otu.freq.full)
        
        nttl = 5*10^8
        ncol.cur.block = floor(nttl / n.perm.available)
        ncol.cur.block = min(n.otus, ncol.cur.block)
        n.round = ceiling(n.otus/ncol.cur.block)
        
        n.otu.freq.full = NULL
        n.otu.tran.full = NULL
        n.otu.omni.full = NULL
        
        file.otu.freq = paste(in.prefix,"_otu.freq.txt",sep="")
        file.otu.tran = paste(in.prefix,"_otu.tran.txt",sep="")
        
        for (r in c(1:n.round)) {
            
            ncol.pre = ncol.cur.block*(r-1)
            ncol.cur = ifelse(r<n.round, ncol.cur.block, n.otus-ncol.pre)
            ncol.post = n.otus - ncol.pre - ncol.cur
            
            otu.freq = read.table(file.otu.freq, colClasses = c(rep("NULL", ncol.pre), rep("numeric", ncol.cur),rep("NULL", ncol.post)), header=TRUE, check.names=FALSE)
            otu.tran = read.table(file.otu.tran, colClasses = c(rep("NULL", ncol.pre), rep("numeric", ncol.cur),rep("NULL", ncol.post)), header=TRUE, check.names=FALSE)
            
            n.otu.freq = NULL
            n.otu.tran = NULL
            n.otu.omni = NULL
            
            for (f in 1:n.vars) {
            
                filelist.otu.freq.perm = list.files(path=in.prefix.dir, pattern=paste(in.prefix.file, "_otu.freq.perm_var",f,"_seed",sep=""), full.names=TRUE)
                filelist.otu.tran.perm = list.files(path=in.prefix.dir, pattern=paste(in.prefix.file, "_otu.tran.perm_var",f,"_seed",sep=""), full.names=TRUE)
                
                otu.freq.perm = NULL
                otu.tran.perm = NULL
                
                for (s in 1:length(filelist.otu.freq.perm)) {
                    
                    cat(paste("OTUs ", ncol.pre+1, "-", ncol.pre+ncol.cur, ", testing.vars[[", f, "]], reading from \"",filelist.otu.freq.perm[s], "\" and \"", filelist.otu.tran.perm[s], "\"\n", sep=""))
                    
                    otu.freq.perm.tmp = read.table(filelist.otu.freq.perm[s], colClasses = c(rep("NULL", ncol.pre), rep("numeric", ncol.cur),rep("NULL", ncol.post)))
                    otu.tran.perm.tmp = read.table(filelist.otu.tran.perm[s], colClasses = c(rep("NULL", ncol.pre), rep("numeric", ncol.cur),rep("NULL", ncol.post)))
                    
                    otu.freq.perm = rbind(otu.freq.perm, otu.freq.perm.tmp)
                    otu.tran.perm = rbind(otu.tran.perm, otu.tran.perm.tmp)
                }
                
                if (r==1 & f==1) {
                    n.perm.completed = nrow(otu.freq.perm)
                }
            
                # freq, tran
                n.otu.freq.var = apply( (t(otu.freq.perm[1:n.perm.completed,]) > t(otu.freq)[,f] + tol.eq), 1, sum) + 0.5 * apply( (abs(t(otu.freq.perm[1:n.perm.completed,]) - t(otu.freq)[,f]) < tol.eq), 1, sum)
                n.otu.tran.var = apply( (t(otu.tran.perm[1:n.perm.completed,]) > t(otu.tran)[,f] + tol.eq), 1, sum) + 0.5 * apply( (abs(t(otu.tran.perm[1:n.perm.completed,]) - t(otu.tran)[,f]) < tol.eq), 1, sum)
                
                n.otu.freq <- rbind(n.otu.freq, t(n.otu.freq.var))
                n.otu.tran <- rbind(n.otu.tran, t(n.otu.tran.var))
                
                # omni
                p.otu.freq.tmp <- n.otu.freq.var/n.perm.completed
                p.otu.tran.tmp <- n.otu.tran.var/n.perm.completed
                pmin.otu.omni <- pmin(p.otu.freq.tmp, p.otu.tran.tmp)
                
                pnull.otu.freq <- (n.perm.completed + 0.5 - apply(otu.freq.perm[1:n.perm.completed,], 2, rank))/n.perm.completed
                pnull.otu.tran <- (n.perm.completed + 0.5 - apply(otu.tran.perm[1:n.perm.completed,], 2, rank))/n.perm.completed
                pnullmin.otu.omni <- pmin(pnull.otu.freq, pnull.otu.tran)
                
                n.otu.omni.var <- apply( (t(pnullmin.otu.omni) < pmin.otu.omni - tol.eq), 1, sum) + 0.5 * apply( (abs(t(pnullmin.otu.omni) - pmin.otu.omni) < tol.eq), 1, sum)
                n.otu.omni <- rbind(n.otu.omni, t(n.otu.omni.var))
            } #f
            
            n.otu.freq.full = cbind(n.otu.freq.full, n.otu.freq)
            n.otu.tran.full = cbind(n.otu.tran.full, n.otu.tran)
            n.otu.omni.full = cbind(n.otu.omni.full, n.otu.omni)
            
        }#r
        
        Aset.freq = matrix(TRUE, n.vars, n.otus)
        Aset.tran = matrix(TRUE, n.vars, n.otus)
        Aset.omni = matrix(TRUE, n.vars, n.otus)
        
        Bset.freq = matrix(FALSE, n.vars, n.otus)
        Bset.tran = matrix(FALSE, n.vars, n.otus)
        Bset.omni = matrix(FALSE, n.vars, n.otus)
        
        p.otu.freq = matrix(NA, n.vars, n.otus)
        p.otu.tran = matrix(NA, n.vars, n.otus)
        p.otu.omni = matrix(NA, n.vars, n.otus)
        
        if (any(Aset.freq)) {
            AtoB.freq <- Aset.freq & (n.otu.freq.full >= n.rej.stop)
            Bset.freq <- Bset.freq | AtoB.freq
            Aset.freq <- Aset.freq & !AtoB.freq
            p.otu.freq[AtoB.freq] <- n.otu.freq.full[AtoB.freq]/n.perm.completed
            p.otu.freq[Aset.freq] <- (n.otu.freq.full[Aset.freq]+1)/(n.perm.completed+1) 
                
            q.otu.freq <- t(apply(p.otu.freq, 1, fdr.Sandev))
                
            Aset.freq.rm.row <- apply(((q.otu.freq < fdr.nominal) & Aset.freq) | (!Aset.freq), 1, all)
            Aset.freq[Aset.freq.rm.row,] = FALSE
        }
        
        if (any(Aset.tran)) {
            AtoB.tran <- Aset.tran & (n.otu.tran.full >= n.rej.stop)
            Bset.tran <- Bset.tran | AtoB.tran
            Aset.tran <- Aset.tran & !AtoB.tran
            p.otu.tran[AtoB.tran] <- n.otu.tran.full[AtoB.tran]/n.perm.completed
            p.otu.tran[Aset.tran] <- (n.otu.tran.full[Aset.tran]+1)/(n.perm.completed+1) 
            
            q.otu.tran <- t(apply(p.otu.tran, 1, fdr.Sandev))
            
            Aset.tran.rm.row <- apply(((q.otu.tran < fdr.nominal) & Aset.tran) | (!Aset.tran), 1, all)
            Aset.tran[Aset.tran.rm.row,] = FALSE
        }
        
        if (any(Aset.omni)) {
            AtoB.omni <- Aset.omni & (n.otu.omni.full >= n.rej.stop)
            Bset.omni <- Bset.omni | AtoB.omni
            Aset.omni <- Aset.omni & !AtoB.omni
            p.otu.omni[AtoB.omni] <- n.otu.omni.full[AtoB.omni]/n.perm.completed
            p.otu.omni[Aset.omni] <- (n.otu.omni.full[Aset.omni]+1)/(n.perm.completed+1) 
                
            q.otu.omni <- t(apply(p.otu.omni, 1, fdr.Sandev))
                
            Aset.omni.rm.row <- apply(((q.otu.omni < fdr.nominal) & Aset.omni) | (!Aset.omni), 1, all)
            Aset.omni[Aset.omni.rm.row,] = FALSE
        }
        
        if (!any(Aset.freq) & !any(Aset.tran) & !any(Aset.omni)) otu.tests.stopped = TRUE 
        
    }# if test.otu   
    
    if (!is.null(p.otu.freq)) colnames(p.otu.freq) <- otu.names
    if (!is.null(p.otu.tran)) colnames(p.otu.tran) <- otu.names
    if (!is.null(p.otu.omni)) colnames(p.otu.omni) <- otu.names
    if (!is.null(q.otu.freq)) colnames(q.otu.freq) <- otu.names
    if (!is.null(q.otu.tran)) colnames(q.otu.tran) <- otu.names
    if (!is.null(q.otu.omni)) colnames(q.otu.omni) <- otu.names
    
    res = list( p.global.freq=p.global.freq, 
                p.global.tran=p.global.tran, 
                p.global.omni=p.global.omni,
                p.otu.freq=p.otu.freq,
                p.otu.tran=p.otu.tran,
                p.otu.omni=p.otu.omni,
                q.otu.freq=q.otu.freq,
                q.otu.tran=q.otu.tran,
                q.otu.omni=q.otu.omni,
                global.tests.stopped=global.tests.stopped,
                otu.tests.stopped=otu.tests.stopped,
                n.perm.completed=n.perm.completed)
    
    return(res)
    
} # test.ldm.fromfile



#' PERMANOVA test of association
#' 
#' This function performs the PERMANOVA test that can allow adjustment of
#' confounders and clustered data.
#' 
#' @param formula a symbolic description of the model to be fitted in the form
#'   of \code{data.matrix ~ sets of variables} or \code{data.matrix |
#'   confounders ~ sets of variables}. The details of model specification are
#'   given in `Details' of \code{ldm}. The only difference between the
#'   model here and that for \code{ldm} is that, the \code{data.matrix}
#'   here is either an OTU table or a distance matrix. If it is an OTU table,
#'   the distance matrix will be calculated internally. If it is a distance
#'   matrix, it can be the original distance matrix and then squared/centered by
#'   the default \code{square.dist=TRUE}/\code{center.dist=TRUE}; it can also be
#'   squared/centered a priori and used with
#'   \code{square.dist=FALSE}/\code{center.dist=FALSE}.
#' @param data an optional data frame, list or environment (or object coercible 
#' by as.data.frame to a data frame) containing the variables of interest and 
#' confounding covariates (but not the OTU data, which are contained in a separate OTU table). 
#' If not found in \code{data}, the variables are taken from environment(formula), 
#' typically the environment from which \code{ldm} is called. The default is NULL.
#' @param tree a phylogenetic tree. Used for calculating a 
#'   phylogenetic-tree-based distance matrix. Not needed if the calculation of 
#'   the requested distance does not involve a phylogenetic tree, or if the 
#'   distance matrix is directly imported through \code{dist}.
#' @param dist.type a string chosen from "Bray-Curtis", "Hellinger", 
#'   "wt-UniFrac" (case insensitive). Used for calculating the distance matrix. 
#'   Not needed if the distance matrix is directly imported through \code{dist}.
#'   The default is "Bray-Curtis".
#' @param is.clustered a logical value indicating whether the observations are 
#'   clustered. The default is FALSE.
#' @param cluster.id cluster identifiers. Unequal cluster size is allowed. The 
#'   default is NULL.
#' @param permute.within.clusters a logical value indicating whether to permute 
#'   data within clusters. The default is FALSE.
#' @param permute.between.clusters a logical value indicating whether to permute
#'   data between clusters. The default is FALSE.
#' @param n.perm.max the maximum number of permutations to run.
#'   The default is 5000.
#' @param n.rej.stop the minimum number of rejections (i.e., the permutation 
#'   statistic exceeds the observed statistic) to obtain before stopping. The 
#'   default is 10.
#' @param seed an integer seed for the random number generator in the 
#'   permutation procedure. The default is NULL (i.e., the seed will be 
#'   generated from the system clock, so it will appear random).
#' @param square.dist a logical variable indicating whether to square the
#'   distance matrix. The default is TRUE.
#' @param center.dist a logical variable indicating whether to center the 
#'   distance matrix as described by Gower (1966). The default is TRUE.
#' @return  a list consisting of 
#' \item{F.statistics}{F statistics for the global test of each set of variables}
#' \item{p.permanova}{p-values for the global test
#'   of each set of variables} 
#' \item{n.perm.completed}{number of permutations completed}
#' \item{permanova.stopped}{a logical value indicating whether the 
#'   stopping criterion has been met by all global tests}
#' @keywords microbiome
#' @author Yi-Juan Hu <yijuan.hu@emory.edu>, Glen A. Satten <gas0@cdc.gov>
#' @export
#' @examples
#' res.permanova <- permanovaC(formula=throat.otu.tab | (Sex+AntibioticUse) ~ (SmokingStatus+PackYears)+Age, 
#'                             data=throat.meta, dist.type="Bray-Curtis", 
#'                             n.perm.max=85600, seed=1)


permanovaC = function(formula, data=NULL, tree=NULL, dist.type="Bray-Curtis",
                      is.clustered=FALSE, cluster.id=NULL, 
                      permute.within.clusters=FALSE, permute.between.clusters=FALSE,
                      n.perm.max=5000, n.rej.stop=10, seed=NULL,
                      square.dist=TRUE, center.dist=TRUE) {  
    
    #------------------------
    # check parameter values
    #------------------------ 
    
    if (! tolower(dist.type) %in% c(tolower("Bray-Curtis"), tolower("wt-UniFrac"), tolower("Hellinger"))) {
        stop("dist.type must be \"Bray-Curtis\", \"wt-UniFrac\", or \"Hellinger\"")
    }
    
    #------------------------
    # vars
    #------------------------
    
    dm = form.call(formula, data=data)
    
    adjust.for.confounders = !is.null(dm$conf)
    
    if (adjust.for.confounders) {
        vars <- list(dm$conf)
        for (j in 1:length(dm$model)) {
            vars[[j+1]] <- dm$model[[j]]
        }
    } else {
        vars <- dm$model
    }
    
    otu.or.dist <- as.matrix(dm$otu.table)
    
    otu.table <- NULL
    dist <- NULL
    if (isSymmetric(otu.or.dist)) {
        dist <- otu.or.dist
    } else {
        otu.table <- otu.or.dist
    }
    
    #------------------------
    # dist matrix
    #------------------------
    
    if (is.null(dist)) {
        dist <- calculate.dist(dist.type=dist.type, otu.table=otu.table, tree=tree)
    }
    
    d.gower <- gower(d=dist, square=square.dist, center=center.dist)
    
    #------------------------
    # setup model
    #------------------------
    
    center.vars=TRUE
    model <- setup.model( vars=vars, adjust.for.confounders=adjust.for.confounders, 
                          center.vars=center.vars, 
                          is.clustered=is.clustered, cluster.id=cluster.id, 
                          permute.within.clusters=permute.within.clusters, 
                          permute.between.clusters=permute.between.clusters)
    
    
    #---------------------
    # observed statistic
    #---------------------
    
    if (model$adjust.for.confounders == TRUE) 
    {
        m1.fit = fit.m1.permanova( d.gower=d.gower, model=model) 
        data.test.statistic = fit.permanova( d.gower=m1.fit$d.resid, model=model) 
        
    } else {
        data.test.statistic = fit.permanova( d.gower=d.gower, model=model) 
    }
    
    permanova = data.test.statistic$permanova
    
    #---------------------
    # permutation
    #---------------------
    
    tol.eq = 10^-16
    
    perm.model = model
    n.permanova = rep(0, length(model$index))
    n.perm.completed = 0
    
    set.seed(seed)
    
    for (i.sim in 1:n.perm.max) {
        
        if (i.sim%%1000==0) cat("perm:", i.sim, "\n")
        
        n.perm.completed = n.perm.completed + 1
        
        if (!model$is.clustered) {
            perm = sample(model$n.obs, replace=FALSE)
        } else {
            
            if (model$permute.between.clusters) {
                perm.cluster = sample(model$n.cluster, replace=FALSE)
            }
            
            low = 0
            up = 0
            perm = rep(0, model$n.obs)
            
            for (i in 1:model$n.cluster) {
                low = up+1
                up = up+model$cluster.size[i]
                
                if (model$permute.between.clusters) {
                    # assume that the covariate matrix is exactly the same for each member of the cluster
                    perm[model$cluster.order[low:up]] = rep( model$cluster.paradigm[ perm.cluster[i] ], model$cluster.size[i] ) 
                }
                if (model$permute.within.clusters) {
                    perm[model$cluster.order[low:up]] = sample(model$cluster.order[low:up], replace=FALSE)
                }
                
            }
        }           
        
        perm.model$m = model$m[perm, , drop=FALSE]
        
        if (model$adjust.for.confounders == TRUE) 
        {
            perm.res = fit.permanova( d.gower=m1.fit$d.resid, model=perm.model) 
        } else 
        {
            perm.res = fit.permanova( d.gower=d.gower, model=perm.model) 
        }
        
        n.permanova <- n.permanova + (perm.res$permanova > permanova + tol.eq) + 0.5*(abs(perm.res$permanova - permanova) < tol.eq)
        
        if (all(n.permanova >= n.rej.stop)) break
        
    }# permutation
    
    p.permanova <- ifelse((n.permanova >= n.rej.stop), n.permanova/n.perm.completed, (n.permanova+1)/(n.perm.completed+1))
    
    res = list( F.statistics=permanova,
                p.permanova=p.permanova, 
                n.perm.completed=n.perm.completed, 
                permanova.stopped=all(n.permanova >= n.rej.stop))
    return(res)
    
}# permanovaC


