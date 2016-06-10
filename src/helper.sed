# delete rows above and below table
0,/^|=\+$/ d
/^|=\+$/, $ d
# replace 'ND's by the R friendly 'nastr'
s/\<ND\>/NA/g
# duplicate NA values for the L,H column
/^\(\(|[^|]*\)\{7\}\)\(|\s*NA\s*\)\(.*\)$/ s//\1|NA,NA\4/
# mark the start of each future row with an '\a'
/^|ANGEL1\s*\(.*\)$/ s//^\a\1/
