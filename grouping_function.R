
#  ----------------
#  ICD-10 structure
#  ================

gr <- read.csv(url("https://raw.githubusercontent.com/danlewer/group_tree/main/icd10_3l.csv"), stringsAsFactors = F)

#  -----------------
#  make example data
#  =================

set.seed(5)

codes <- c(outer(gr$icd3, 0:9, paste0))
codes <- sample(codes, 1e3)
probs <- sample(seq_len(100), length(codes), replace = T) ^ 4
sampleICD <- sample(codes, 1e5, replace = T, prob = probs)
hist(table(sampleICD), xlab = 'Number of appearances', ylab = 'Number of codes', col = 'grey', main = 'Distribution of codes')
length(unique(sampleICD)) # number of codes in sample

# assign ICD10 groups to sample data

sampleICDgr <- gr[match(substr(sampleICD, 0, 3), gr$icd3),]
rownames(sampleICDgr) <- NULL
sampleICDgr$icd4 <- sampleICD
sampleICDgr$chapter <- sub(" .*", "", sampleICDgr$chapter) # make shorter chapter names

# -----------------
# grouping function
# =================

group <- function(dat, levs, min.group, min.other = NA, headings = T) {
  
  nl <- length(levs)
  
  # make first summary table 'y'
  dat$n <- 1
  af <- formula(paste0('n~', paste0(levs, collapse = '+')))
  y <- aggregate(af, dat, length)
  y$meta <- 1
  levs2 <- c('meta', levs)
  
  # replace groups smaller than 'min.group' with 'zOth'
  for(i in seq_len(nl)) {
    ls <- levs2[1:(i + 1)]
    af2 <- formula(paste0('n~', paste0(ls, collapse = '+')))
    z <- aggregate(af2, y[, c(ls, 'n')], sum)
    names(z)[ncol(z)] <- 'ng'
    y <- merge(y, z)
    y[, levs[i]][y[, 'ng'] < min.group] <- 'zOth'
    y$ng <- NULL
  }
  
  # make second summary table q
  q <- aggregate(af, y, sum)
  
  # roll-up 'others'
  if (!is.na(min.other)) {
    for(i in seq_len(nl)) {
      q$meta <- 1
      ls <- levs2[1:(i + 1)]
      af2 <- formula(paste0('n~', paste0(ls, collapse = '+')))
      r <- aggregate(af2, q, sum)
      r <- r[order(r$n),]
      r <- split(r, f = r[, ls[-length(ls)]], drop = T)
      r <- lapply(r, function(x) {
        cs <- cumsum(x$n) < min.other
        if (any(cs)) {x[, tail(ls, 1)][seq_len(which.max(!cs))] <- 'zOth'}
        return (x)
      })
      r <- do.call(rbind, r)
      r <- aggregate(af2, r, sum)
      r$ng <- 'x'
      r$n <- NULL
      q <- merge(q, r, by = ls, all.x = T)
      h <- levs[i]
      q[, h] <- ifelse(is.na(q$ng), 'zOth', q[, h])
      q <- aggregate(af, q, sum)
    }
    q <- q[do.call(order, q), ]
  }
  
  # add headings
  
  if (headings == T) {
    heads <- lapply(seq_len(nl-1), function(i) {
      ls <- levs[seq_len(i)]
      h <- q[, c(ls, 'n')]
      h <- split(h, f = h[, ls], drop = T)
      h <- h[which(sapply(h, nrow) > 1)]
      if (length(h) == 0) return (NULL) else {
        h <- do.call(rbind, h)
        af2 <- formula(paste0('n~', paste0(ls, collapse = '+')))
        h <- aggregate(af2, h, sum)
        ec <- setdiff(levs, ls)
        h <- cbind(h, matrix(rep('aaTotal', nrow(h) * length(ec)), ncol = length(ec), dimnames = list(NULL, ec)))
        return (h[, names(q)])
      }
    })
    heads <- do.call(rbind, heads)
    q <- rbind(heads, q)
    q <- rbind(cbind(matrix(rep('aaTotal', length(levs)), ncol = length(levs), dimnames = list(NULL, levs)), data.frame(n = sum(y$n))), q)
  }
  
  q[, levs] <- lapply(q[, levs], as.character)
  q <- q[do.call(order, q), ]
  return (q)
  
}

# -----------
# example use
# ===========

icd_summary <- group(dat = sampleICDgr, levs = c('chapter', 'l1', 'l2', 'icd4'), min.group = 500, min.other = 500, headings = T)

