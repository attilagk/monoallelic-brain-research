
cur.seg <- function(i.cur, prev.seg, df) {
    is.known <- function(i)
        df[rownames(df)[i], "imprinting.status"] != "candidate"
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

segs.one.chr <- function(df) {
    for(i in seq_len(nrow(df))) {
        prev.seg <- df$seg[i - 1]
        df$seg[i] <- cur.seg(i, prev.seg, df)
    }
    df
}

segs.one.chr.adjust <- function(df, seg.0 = 0) {
    is.seg <- df$seg > seg.0
    df$seg[is.seg] <- df$seg[is.seg] + abs(seg.0)
    is.antiseg <- df$seg < seg.0
    df$seg[is.antiseg] <- df$seg[is.antiseg] - abs(seg.0)
    return(df)
}

segs.two.chr <- function(df.prev, df.cur)
    rbind(df.prev, segs.one.chr.adjust(df.cur, seg.0 = df.prev$seg[nrow(df.prev)]))
