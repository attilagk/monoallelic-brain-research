# add imprinted segments component 'seg' to data frame 'gs'
#
# Arguments
# gs: a "gene summary" data frame representing a single chromosome
#
# Value
# a new 'gene summary'-type data frame with an additional seg component.  The
# positive values of seg identify the imprinting clusters, the negative ones
# the intercluster segments and zero the segment before the first cluster on
# the first chromosome.
segs.one.chr <- function(gs, remove.str) {
    # return current segment given current index and previous segment
    #
    # Details
    #
    # Let's denote the imprinting cluster membership of a gene g with K
    # (kluster), the non-cluster membership with C (candidate) and the case
    # that g is the '0th' gene with 0. The notation X -> Y expresses that the
    # previos gene is in set X and the current gene in Y where both X and Y
    # may be one of K, C or 0. Then there are 6 cases:
    # case 1: 0 -> K
    # case 2: 0 -> C
    # case 3: K -> K
    # case 4: C -> K
    # case 5: K -> C
    # case 6: C -> C
    cur.seg <- function(i.cur, prev.seg) {
        is.known <- function(i)
            gs[rownames(gs)[i], "imprinting.status"] != remove.str
        i.prev <- i.cur - 1
        if(i.prev == 0) {
            if(is.known(i.cur)) 1 # 0 -> K
            else 0 # 0 -> C
        } else {
            if(is.known(i.cur)) {
                if(is.known(i.prev)) prev.seg # K -> K
                else 1 - prev.seg # C -> K
            } else {
                if(is.known(i.prev)) - prev.seg # K -> C
                else prev.seg # C -> C
            }
        }
    }
    # iteratively call 'cur.seg' and store results in 'gs$seg'
    for(i in seq_len(nrow(gs))) {
        prev.seg <- gs$seg[i - 1]
        gs$seg[i] <- cur.seg(i, prev.seg)
    }
    gs
}

# add imprinted segments component 'seg' to data frame 'gs'
#
# Arguments
# gs: a "gene summary" data frame representing all chromosomes
#
# Value
# a new data frame with the same components as 'gs' as well as a 'seg'
# component, an integer vector that describes a "step function" F with 2n + 1
# unique levels defined on the set of genes.  F has the following properties:
# 1. F(g) = 0 if g is the fist gene of the first chromosome and g is not in a cluster
#    a. if the first gene g is in a cluster (which must be the 1st one) then F(g) = 1
# 2. F(g) = k > 0 if g is in the k'th cluster
# 3. F(g) = -k < 0 if g is between the k'th and k+1'th imprinting cluster
#
# It follows that abs(F) is a non-decreasing function which steps 1 at the
# first gene of each cluster
make.impr.segs <- function(gs, remove.str = "candidate, >1MB") {
    # merge data frames for two consecutive chromosomes
    segs.two.chr <- function(df.prev, df.cur) {
        prev <- segs.one.chr(df.prev, remove.str)
        cur <- segs.one.chr(df.cur, remove.str)
        rbind(prev,
              segs.one.chr.adjust(cur, seg.0 = prev$seg[nrow(prev)]))
    }
    segs.one.chr.adjust <- function(df, seg.0 = 0) {
        is.seg <- df$seg > 0
        is.antiseg <- df$seg < 0
        is.neither <- df$seg == 0
        df$seg[is.seg] <- df$seg[is.seg] + abs(seg.0)
        df$seg[is.antiseg] <- df$seg[is.antiseg] - abs(seg.0)
        df$seg[is.neither] <- seg.0
        return(df)
    }
    # sort rows according to location
    gs.loc <- gs[with(gs, order(chr, start)), ]
    # break gs.loc into several data frames according to chromosomes and store them in list L
    L <- lapply(levels(gs.loc$chr), function(k) gs.loc[gs.loc$chr == k, ])
    Reduce(segs.two.chr, L)
}

do.beta.age <- function(l.M, conf.lev = 0.99) {
    B <- get.CI(l.M, coef.name = "Age", conf.lev = conf.lev)
    B$start <- gs[rownames(B), "start"]
    B$chr <- gs[rownames(B), "chr"]
    B <- B[with(B, order(chr, start)), ]
    B$cluster <- factor(x <- as.character(gs[rownames(B), "cluster"]), levels = unique(x), ordered = TRUE)
    B$Gene <- factor(B$Gene, levels = rev(B$Gene), ordered = TRUE)
    B$imprinting.status <- gs[rownames(B), "imprinting.status"]
    return(B)
}


do.beta <- function(l.M, coef.name = "Age", conf.lev = 0.99) {
    B <- get.CI(l.M, coef.name = coef.name, conf.lev = conf.lev)
    B$Coefficient <- coef.name
    B$start <- gs[rownames(B), "start"]
    B$chr <- gs[rownames(B), "chr"]
    B <- B[with(B, order(chr, start)), ]
    B$cluster <- factor(x <- as.character(gs[rownames(B), "cluster"]), levels = unique(x), ordered = TRUE)
    B$Gene <- factor(B$Gene, levels = rev(B$Gene), ordered = TRUE)
    B$imprinting.status <- gs[rownames(B), "imprinting.status"]
    return(B)
}


my.segplot2 <- function(data, ...) {
segplot(Gene ~ Lower.CL + Upper.CL | Coefficient, data = data,
                       groups = cluster,
               # groups = imprinting.status, # groups does not seem to work with segplot
               scales = list(x = list(draw = c(F, T)), y = list(relation = "free", rot = 0), alternating = FALSE, relation = "free"),
               panel = function(...) {
                   panel.grid(v = -1, h = -1)
                   panel.abline(v = 0, col = "red")
                   panel.segplot(...)
               },
               draw.bands = FALSE, centers = beta.hat,
               #as.table = FALSE,
               layout = c(1, length(levels(data$cluster))),
               par.settings = list(layout.heights = list(panel = (x <- table(data$cluster))[x > 0]),
                                   add.text = list(cex = 0.8)),
               xlab = expression(beta[age]),
               ...)
}
