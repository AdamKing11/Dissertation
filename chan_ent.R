library(parallel)
library(pbapply)
library(tidyverse)
library(plyr)

library(lme4)
library(lmerTest)
library(car)
library(mgcv)

library(xtable)
library(ggrepel)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)

rm(list = ls())

nb_cores <- 1
nb_scramble <- 64
p_threshold <- .95
z. <- scale

fdir <- '~/rdy/'
lexdir <- paste0(fdir, 'lexicons_ce/')
noncedir <- paste0(fdir, 'new_forms/')
outdir <- paste0(fdir, 'ce/')
lex.files <- list.files(lexdir)
lex.files <- subset(lex.files, grepl('lex.txt$', lex.files))

filtered_H <- function(P, p_threshold = .95) {
  P <- sort(P, decreasing = T)
  for (i in 1:length(P)) {
    if (sum(P[1:i]) >= p_threshold) {
      break
    }
  }
  P <- P[1:i]
  P <- P / sum(P)
  sum(P * -log2(P))
}

filtered_H_length <- function(P, p_threshold = .95) {
  P <- sort(P, decreasing = T)
  for (i in 1:length(P)) {
    if (sum(P[1:i]) >= p_threshold) {
      break
    }
  }
  P <- P[1:i]
  i
}

lmer_to_latex <- function(m) {
  df <- as.data.frame(summary(m)$coefficients)
  # add p-vals
  if ('df' %in% names(df)) {
    df <- subset(df, select = -df)
  }
  fixed.effects <- print(xtable(df, digits = 3))
  # add extra columns to table
  # add title for "fixed effects
  fixed.effects <- sub('\\n  ', '\n {\\\\bf A. Fixed Effects:} \\\\\\\\\n', fixed.effects, perl = T)
  fixed.effects <- gsub('\\n\\\\end.tab.*', '', fixed.effects, perl = T)
  fixed.effects <- substr(fixed.effects, 1, nchar(fixed.effects)-10)
  fixed.effects
  
  random.effects <- c('\\hline \\hline',
                      '{\\bf B. Random Effects:} \\\\',
                      '\\hline',
                      '& Name & Variance & Std.Dev. \\\\', 
                      '\\hline')
  random.factors <- sapply(1:length(ranef(m)), function(i) {
    factor.name <- names(ranef(m))[i]
    factor.df <- as.data.frame(ranef(m)[i])
    names(factor.df) <- sapply(names(factor.df), function(s) { substr(s, nchar(factor.name)+2, nchar(s))})
    intercept.var <- substr(var(factor.df[,1]),1,5)
    intercept.sd <- substr(sd(factor.df[,1]),1,5)
    factor_str <- paste0(c(factor.name, ' & (Intercept) & ', intercept.var, 
                           ' & ', intercept.sd, ' \\\\'), collapse = '')
    if (length(names(factor.df))>1) { 
      for (j in 2:length(names(factor.df))) {
        slope.name <- names(factor.df)[j]
        slope.var <- substr(var(factor.df[,j]),1,5)
        slope.sd <- substr(sd(factor.df[,j]),1,5)
        factor_str <- paste0(c(factor_str, '\n & ', slope.name, ' & ', slope.var, 
                             ' & ', slope.sd, ' \\\\'), collapse = '')
      }
    }
    factor_str
  })
  random.effects <- c(random.effects, random.factors)
  latex.table <- paste0(c(fixed.effects, random.effects,
                          '\\end{tabular}\n\\end{table}'), collapse =  '\n')
  latex.table <- gsub('family', 'Family', latex.table)
  latex.table <- gsub('language', 'Language', latex.table)
  latex.table <- gsub('000', '001', latex.table)
  # comment out some stuff
  latex.table <- gsub('\\\\begin.table', '%\\begin{table', latex.table)
  latex.table <- gsub('\\\\centering', '%\\centering', latex.table)
  latex.table <- gsub('\\\\end.table', '%\\end{table', latex.table)
  latex.table
}

find_up <- function(d) {
  phones <- sort(d$phones,method = 'radix')
  ups <- sapply(1:(length(phones)), function(i) {
    w <- phones[i]
    for (j in nchar(w):0) {
      if (i < length(phones) && substr(w,1,j) == substr(phones[i+1],1,j)) {
        j <- j + 1
        break
      } else if (i > 1 && substr(w,1,j) == substr(phones[i-1],1,j)) {
        j <- j + 1
        break
      }
    }
    if (j == 0) {
      j <- nchar(w) + 1
    }
    j
  })
  df <- data.frame(phones = phones, up = ups)
  d$up <- df$up[match(d$phones, df$phones)] 
  d
}

shuffle_frequencies <- function(d, shuffle_same_tail = F, shuffle_same_nbrhood = F, keep_length = T) {
  total_freq <- sum(d$freq)
  if (keep_length & shuffle_same_tail) {
    for (i in 1:max(d$word.length)) {
      if (nrow(subset(d, word.length == i)) > 1) {
        for (j in 1:i) {
          if (nrow(subset(d, word.length == i & up == j)) > 1) {
            d$freq[d$word.length == i & d$up == j] <- 
              sample(d$freq[d$word.length == i & d$up == j])
          }
        }
      }
    }
  } else if (keep_length & shuffle_same_nbrhood) {
    for (i in 1:max(d$word.length)) {
      if (nrow(subset(d, word.length == i)) > 1) {
        for (j in 0:max(subset(d, word.length == i)$min.pairs)) {
          if (nrow(subset(d, word.length == i & min.pairs == j)) > 1) {
            d$freq[d$word.length == i & d$min.pairs == j] <-
              sample(d$freq[d$word.length == i & d$min.pairs == j])
          }
        }
      }
    }
  } else if (keep_length) {
    for (i in 1:max(d$word.length)) {
      if (nrow(subset(d, word.length == i)) > 1) {
        d$freq[d$word.length == i] <- sample(d$freq[d$word.length == i])
      }
    }
  } else {
    d$freq <- sample(d$freq)
  }
  d$unig.p <- d$freq / total_freq
  d$unig.h <- -log2(d$unig.p)
  d$log.prob <- -d$unig.h
  
  d
}

lf <- function(f, do_up = T, min_freq_rank = 10000) {
  d <- read.csv(f, sep = '\t', header = F)
  names(d) <- c('word', 'phones', 'freq') 
  d$word <- as.character(d$word)
  d$phones <- as.character(d$phones)
  d$phones <- gsub(' ', '', d$phones)
  # for tone languages, exclude tone marking from length
  if (grepl('Vietnamese|Cantonese', f)) {
    d$word.length <- nchar(gsub('[0-9]', '', d$phones))
  } else {
    d$word.length <- nchar(d$phones)
  }
  d <- subset(d, word.length <= 12)
  d <- subset(d, !(nchar(word) == 1 & nchar(phones) > 1))
  
  d$freq <- as.double(d$freq)
  d$freq.rank <- rank(-d$freq, ties.method = 'random')
  
  if (min_freq_rank > 0) {
    d <- subset(d, freq.rank <= min_freq_rank)
  }
  
  total_freq <- sum(d$freq)
  d$unig.p <- d$freq / total_freq
  d$unig.h <- -log2(d$unig.p)
  d$log.prob <- -d$unig.h
  
  if (do_up) {
    d <- find_up(d)
  }
  
  d
}

load_nonce <- function(f) {
  dn <- read.csv(f, header = F)
  names(dn)[1] <- 'phones'
  dn$phones <- enc2utf8(as.character(dn$phones))
  if (grepl('Cantonese|Vietnamese', f)) {
      dn$word.length <- nchar(gsub('[0-9]', '', dn$phones))
  } else {
    dn$word.length <- nchar(dn$phones)
  }
  
  dn
}

nonce_lexicon <- function(d, lang_name, keep_length = T, do_up = T) {
  #d <- d[!duplicated(d$phones), ]
  if (is.data.frame(lang_name)) {
    dn <- lang_name
    lang_name <- d$language
  } else {
    dn <- load_nonce(paste0(noncedir, lang_name, '.new.txt'))
    lang_name <- d$language
  }
  if (keep_length) {
    for (j in min(d$word.length):max(d$word.length)) {
      should_replace <- length(d$phones[d$word.length == j]) > length(dn$phones[dn$word.length == j])
      d$phones[d$word.length == j] <-
        sample(dn$phones[dn$word.length == j], 
               size = length(d$phones[d$word.length == j]),
               replace = should_replace)
    }
  } else{
    d$phones <- sample(dn$phones, 
                       size = length(d$phones),
                       replace = F)
    if (is.character(lang_name) && grepl('Cantonese|Vietnamese', lang_name)) {
      d$word.length <- nchar(gsub('[0-9]', '', d$phones))
    } else {
      d$word.length <- nchar(d$phones)
    }
  }
  
  if(do_up) {
    d <- find_up(d)
  }
  
  d
}

add_language_family <- function(d) {
  d$language <- as.character(d$language)
  d$family <- d$language
  d$family[d$language %in% c('Armenian', 'Bengali', 'Dutch', 
                             'English', 'French', 'German', 'Italian', 
                             'Russian', 'Slovak', 'Spanish')] <- 'Indo-European'
  d$family[d$language %in% c('Arabic', 'Hausa', 'Hebrew')] <- 'Afroasiatic'
  d$family[d$language %in% c('Malay', 'Tagalog')] <- 'Austronesean'
  d$family[d$language %in% c('Finnish', 'Hungarian')] <- 'Uralic'
  d$family[d$language %in% c('Cantonese', 'Vietnamese')] <- 'Sino-Tibetan'
  d$language <- as.factor(d$language)

  d
}

build_cohorts <- function(d, shuffle = F, ignore_phones = '') {
  total_freq <- sum(d$freq)
  total_count <- length(d$freq)
  
  if (shuffle) {
    d <- shuffle_frequencies(d, F, F)
  }
  
  cohorts <- lapply(1:max(nchar(d$phones)), function(i) {
    len.i <- subset(d, nchar(phones) >= i)
    len.i$prefix <- substr(len.i$phones,1,i)
    len.i.aggr <- aggregate(freq ~ prefix, len.i, FUN = sum)
    len.i.aggr.counts <- aggregate(freq ~ prefix, len.i, FUN = length)
    names(len.i.aggr.counts)[2] <- 'cohort.count'
    len.i.aggr <- left_join(len.i.aggr, len.i.aggr.counts, by = 'prefix')
    len.i.aggr$cohort.count.scaled <- scale(len.i.aggr$cohort.count)
    
    len.i.aggr
  })
  cohorts <- ldply(cohorts, data.frame)
  
  if (nchar(ignore_phones) > 0) {
    cohorts <- subset(cohorts, !grepl(paste0('[', ignore_phones, ']'), 
                                      substr(prefix, nchar(prefix), nchar(prefix))))
  }
  cohorts$last.one <- substr(cohorts$prefix, 1, nchar(cohorts$prefix) - 1)
  cohorts$last.one[nchar(cohorts$prefix) == 1] <- '0'
  
  cohorts$last.freq <- cohorts$freq[match(cohorts$last.one, cohorts$prefix)]
  cohorts$last.freq[nchar(cohorts$prefix) == 1] <- total_freq
  
  cohorts$last.count <- cohorts$cohort.count[match(cohorts$last.one, cohorts$prefix)]
  cohorts$last.count[nchar(cohorts$prefix) == 1] <- total_count
  
  cohorts
}

lex2si <- function(d, shuffle = F, shuffle_same_tail = F, shuffle_same_nbrhood = F) {
  total_freq <- sum(d$freq)
  total_count <- length(d$freq)
  if (shuffle) {
    d <- shuffle_frequencies(d, shuffle_same_tail = shuffle_same_tail,
                             shuffle_same_nbrhood = shuffle_same_nbrhood)
  }
  
  cohorts <- build_cohorts(d, shuffle = F)
  
  si <- lapply(1:max(nchar(d$phones)), function(i) {
    len.i <- subset(d, nchar(phones) >= i)
    len.i$position <- i  
    len.i$phone <- substr(len.i$phones, i, i)
    len.i$prefix.plus.seg <- substr(len.i$phones, 1, i)
    if (i > 1) {
      len.i$prefix <- substr(len.i$phones, 1, i-1)
      len.i$prefix.freq <- cohorts$freq[match(len.i$prefix, cohorts$prefix)]
      len.i$phone.freq <- cohorts$freq[match(len.i$prefix.plus.seg, cohorts$prefix)]      
    } else {
      len.i$prefix.freq <- total_freq
      len.i$phone.freq <- cohorts$freq[match(len.i$phone, cohorts$prefix)]  
    }
    len.i
  })
  si <- ldply(si, data.frame)
  si$cohort.size <- cohorts$cohort.count[match(si$prefix.plus.seg, cohorts$prefix)]
  si$prefix.cohort.size <- cohorts$cohort.count[match(si$prefix, cohorts$prefix)]
  si$prefix.cohort.size[si$position == 1] <- total_count
  
  si$seg.p.type <- si$cohort.size / si$prefix.cohort.size
  
  si$seg.p <- si$phone.freq / si$prefix.freq
  si$seg.p.excluded <- (si$phone.freq - si$freq) / (si$prefix.freq - si$freq)
  
  si$seg.p.excluded[si$cohort.size == 1] <- 
    si$seg.p[si$cohort.size == 1] 
  
  si <- subset(si, select =-c(prefix, prefix.plus.seg, prefix.cohort.size))
  
  si$seg.h <- -log2(si$seg.p)
  si$seg.h.type <- -log2(si$seg.p.type)
  
  si
}

med <- function(d, min_len = 2, max_len = 8) {
  
  d <- subset(d, word.length >= min_len - 1 & word.length <= max_len + 1)
  mean_edit_dist <- lapply(min_len:max_len, function(i) {
    lex.i <- subset(d, word.length == i)
    print(paste0(i, ' - ', nrow(lex.i)))
    
    if (nrow(lex.i) >= 2) {
      lex.i$scaled.unig.h <- scale(lex.i$unig.h)
      lex.i.plusminus <- subset(d, word.length >= i - 1 & word.length <= i + 1)
      mean_edit_dist <- pblapply(unique(lex.i$word), function(w){
        p <- lex.i$phones[match(w, lex.i$word)]
        edit_dists_all <- adist(p, unique(lex.i.plusminus$phones))
        ###
        edit_dists <- edit_dists_all[edit_dists_all <= 3]
        edit_dists_nonmp <- edit_dists_all[edit_dists_all > 1]
        ###
        data.frame(word = w, phones = p, 
                   min.pairs = sum(edit_dists == 1),
                   min.pairs.2 = sum(edit_dists == 2),
                   min.pairs.3 = sum(edit_dists == 3),
                   mean.dist = mean(edit_dists),
                   mean.dist.nonmp = mean(edit_dists_nonmp),
                   mean.dist.all = mean(edit_dists_all))
      }, cl = nb_cores)
      mean_edit_dist <- ldply(mean_edit_dist, data.frame)
      mean_edit_dist$scaled.unig.h <- 
        lex.i$scaled.unig.h[match(mean_edit_dist$word, lex.i$word)]
      mean_edit_dist$network <- mean_edit_dist$min.pairs + mean_edit_dist$min.pairs.2
      
      mean_edit_dist
    }
  })
  mean_edit_dist <- ldply(mean_edit_dist, data.frame)
  mean_edit_dist[is.na(mean_edit_dist)] <- 0
  
  d <- join(d, mean_edit_dist, by = c('word', 'phones'))
  d <- subset(d, word.length >= min_len & word.length <= max_len)
  d
}

diphone_h <- function(lex, use_freq = T) {
  print('calculating diphone informativity')
  df <- pblapply(unique(subset(lex, nchar(phones) >= 2)$word), function(w) {
    phones <- lex$phones[match(w, lex$word)]
    diphones <- sapply(2:nchar(phones)-1, function(i) { 
      substr(phones, i, i + 1)
    })
    if (use_freq) {
      freq <- lex$freq[match(w, lex$word)]
    } else {
      freq <- 1
    }
    data.frame(diphone = diphones, freq = freq)
  }, cl = nb_cores)
  a <- aggregate(freq ~ diphone, ldply(df, data.frame), FUN = sum)
  a$first.seg <- substr(a$diphone, 1, 1)
  first.seg.a <- aggregate(freq ~ first.seg, a, FUN = sum)
  names(first.seg.a)[2] <- 'first.seg.freq'
  a <- left_join(a, first.seg.a, by = 'first.seg')
  a$p <- a$freq / a$first.seg.freq
  a$h <- -log2(a$p)
  print('assigning information values to words')
  diphone <- pblapply(unique(subset(lex, nchar(phones) >= 2)$word), function(w) {
    phones <- lex$phones[match(w, lex$word)]
    diphones <- sapply(2:nchar(phones)-1, function(i) { 
      substr(phones, i, i + 1)
    })
    data.frame(word = w, phones = phones, diphone = diphones, position = 1 + 2:nchar(phones)-1)
  }, cl = nb_cores)
  diphone <- ldply(diphone, data.frame)
  
  diphone$phones <- as.character(diphone$phones)
  diphone$diphone.h <- a$h[match(diphone$diphone, a$diphone)]
  diphone$unig.h <- lex$unig.h[match(diphone$phones, lex$phones)]
  diphone$word.length <- nchar(diphone$phones)
  diphone$freq <- lex$freq[match(diphone$word, lex$word)]
  diphone$up <- lex$up[match(diphone$phones, lex$phones)]
  diphone$post_up <- diphone$position > diphone$up
  diphone
}

diphone_redun <- function(diphone) {
  diphone.all <- aggregate(diphone.h ~ word + phones + freq + unig.h + word.length + up, diphone, FUN = sum)
  names(diphone.all)[ncol(diphone.all)] <- 'sum.diphone.h'
  diphone.postup <- aggregate(diphone.h ~ word, subset(diphone, post_up), FUN = sum)
  names(diphone.postup)[2] <- 'postup.sum.diphone.h'
  diphone.a <- left_join(diphone.all, diphone.postup, by = 'word')
  diphone.a$postup.sum.diphone.h[is.na(diphone.a$postup.sum.diphone.h)] <- 0
  diphone.a$diphone.h.redundant_proportion <- diphone.a$postup.sum.diphone.h / diphone.a$sum.diphone.h
  
  diphone.a$segments.redundant <- diphone.a$word.length - diphone.a$up
  diphone.a$redundancy.prop <- diphone.a$segments.redundant / diphone.a$word.length
  diphone.a
}

entropy_by_cohort <- function(d, include_post_up = F, min_per_cohort = 1, min_segs = 1, ignore_phones = '') {
  if (!include_post_up) {
    d$phones <- substr(d$phones, 1, d$up)
  }
  
  total_freq <- sum(d$freq)
  cohorts <- build_cohorts(d, ignore_phones = ignore_phones)
  cohorts$cohort.p <- cohorts$last.freq / total_freq
  cohorts$seg.p <- cohorts$freq / cohorts$last.freq
  
  cohorts.H <- aggregate(seg.p ~ last.one + cohort.p, cohorts, filtered_H)
  names(cohorts.H)[3] <- 'entropy'
  cohorts.nb_segs <- aggregate(seg.p ~ last.one, cohorts, filtered_H_length)
  names(cohorts.nb_segs)[2] <- 'nb.cont.segs'
  cohorts.cohort_count <- aggregate(cohort.count ~ last.one, cohorts, sum)
  
  chrt <- join(cohorts.H, cohorts.nb_segs, by = 'last.one')
  chrt <- join(chrt, cohorts.cohort_count, by = 'last.one')
  names(chrt)[1] <- 'cohort'
  
  chrt$position <- nchar(chrt$cohort) + 1
  chrt$position[chrt$cohort == 0] <- 1
  
  chrt$cohort.h <- -log2(chrt$cohort.p) 
  chrt$log.prob <- -chrt$cohort.h
  
  ###
  chrt <- subset(chrt, cohort.count >= min_per_cohort & nb.cont.segs >= min_segs)
  chrt$log.cont.segs <- log2(chrt$nb.cont.segs)
  chrt$log.cohort.count <- log2(chrt$cohort.count)
  ###
  
  chrt$cohort.h.by_length <- 0
  chrt$nb.cont.segs.by_length <- 0
  chrt$log.cont.segs.by_length <- 0
  chrt$entropy.by_length <- 0
  chrt$cohort.count.by_length <- 0
  chrt$log.cohort.count.by_length <- 0
  for (i in 1:max(chrt$position)) {
    chrt$cohort.h.by_length[chrt$position == i] <- scale(chrt$cohort.h[chrt$position == i])
    chrt$nb.cont.segs.by_length[chrt$position == i] <- scale(chrt$nb.cont.segs[chrt$position == i])
    chrt$log.cont.segs.by_length[chrt$position == i] <- scale(chrt$log.cont.segs[chrt$position == i])
    chrt$entropy.by_length[chrt$position == i] <- scale(chrt$entropy[chrt$position == i])
    chrt$cohort.count.by_length[chrt$position == i] <- scale(chrt$cohort.count[chrt$position == i])
    chrt$log.cohort.count.by_length[chrt$position == i] <- scale(chrt$log.cohort.count[chrt$position == i])
  }
  chrt$log.prob.by_length <- -chrt$cohort.h.by_length
  
  chrt$entropy.by_segs <- 0
  chrt$log.prob.by_segs <- 0
  for (i in 1:max(chrt$nb.cont.segs)) {
    chrt$entropy.by_segs[chrt$nb.cont.segs == i] <-
      scale(chrt$entropy[chrt$nb.cont.segs == i])
    chrt$log.prob.by_segs[chrt$nb.cont.segs == i] <-
      scale(chrt$log.prob[chrt$nb.cont.segs == i])
  }
  chrt[is.na(chrt)] <- 0
  
  chrt
}

phones_by_position <- function(d, N = 1) {
  d <- subset(d, word.length > (N-1))
  for (i in 1:(max(d$word.length) - (N-1))) {
    len_i <- subset(d, word.length - (N-1) >= i)
    if (i == 1) {
      df <- data.frame(position = i, phone = substr(len_i$phones, i, i + (N-1)), freq = len_i$freq)
    } else{
      df <- rbind(df, 
                  data.frame(position = i, phone = substr(len_i$phones, i, i + (N-1)), freq = len_i$freq))
    }
  }
  df
}

lex_ent <- function(l) {
  #  l <- l[!duplicated(paste(l$phones, l$position)),]
  
  sum(l$unig.p * l$seg.p * l$seg.h)
}

bag_ent <- function(l) {
  l.a <- aggregate(freq ~ phone, l, FUN = sum)
  l.a$seg.p <- l.a$freq / sum(l.a$freq)
  l.a$seg.h <- -log2(l.a$seg.p)
  sum(l.a$seg.p * l.a$seg.h)
}

phone_ent <- function(l) {
  l$total_f <- sum(l$freq)
  l$phone_p <- l$freq / l$total_f
  sum(l$phone_p * -log2(l$phone_p))
}

first_dist <- function(d, N = 1, p_threshold = .95, shuffle = F, keep_running = F) {
  total_freq <- sum(d$freq)
  total_count <- length(d$freq)
  if (shuffle) {
    d <- shuffle_frequencies(d)
  }
  d$first.phone <- substr(d$phones, 1, N)
  d$log.freq <- log2(d$freq)
  
  fp.token <- aggregate(freq ~ first.phone, d, sum)
  fp.log.token <- aggregate(log.freq ~ first.phone, d, sum)
  fp.type <- aggregate(freq ~ first.phone, d, length)
  names(fp.type)[2] <- 'count'
  first.phones <- join(fp.token, fp.type, by = 'first.phone')
  first.phones <- join(first.phones, fp.log.token, by = 'first.phone')
  names(first.phones)[1] <- 'phone'
  first.phones$seg.p <- first.phones$freq / total_freq
  first.phones$seg.p.type <- first.phones$count / total_count
  first.phones$mean <- first.phones$freq / first.phones$count
  first.phones$log.mean <- first.phones$log.freq / first.phones$count
  first.phones$p.rank <- rank(-first.phones$seg.p, ties.method = 'random')
  
  if (p_threshold <= 1) {
  # figure the running probability totals for segments
    
    running_token_p <- 0
    first.phones$running.token.p <- 0
    for (phone in first.phones$phone[order(first.phones$seg.p)]) {
      phone_p <- first.phones$seg.p[match(phone, first.phones$phone)]
      running_token_p <- running_token_p + phone_p
      first.phones$running.token.p[match(phone, first.phones$phone)] <- running_token_p
    }
    first.phones <- subset(first.phones, running.token.p >= (1 - p_threshold))
    if (!keep_running) {
      first.phones <- subset(first.phones, select = -c(running.token.p))
      first.phones <- subset(first.phones, select = -c(p.rank))
    }
  } else if (p_threshold > 1) {
    first.phones <- subset(first.phones, p.rank <= p_threshold)
    if (!keep_running) {
      first.phones <- subset(first.phones, select = -c(p.rank))
    }
  }
  first.phones
}

first_vs_later <- function(lex, p_threshold = .95, pos = 1) {
  # split into Nth (default first) segment and all segments after that
  lex.1 <- aggregate(freq ~ phone, subset(lex, position == pos), sum)
  lex.2 <- aggregate(freq ~ phone, subset(lex, position > pos), sum)
  
  # figure out the relative prob for segments at position N
  lex.1$p <- lex.1$freq / sum(lex.1$freq)
  
  if (p_threshold < 1) {
    # b/c strong Zipfian dist., subset out segments that together are less than
    # 10%
    running_p <- 0
    lex.1$running.p <- 0
    for (phone in lex.1$phone[order(lex.1$p)]) {
      p <- lex.1$p[match(phone, lex.1$phone)]
      running_p <- running_p + p
      lex.1$running.p[match(phone, lex.1$phone)] <- running_p
    }
    lex.1 <- subset(lex.1, running.p >= (1 - p_threshold))
    lex.1 <- subset(lex.1, select = -c(running.p))
  } else if (p_threshold > 1) {
    # if passed an intiger greater than 1, just select N most probable
    lex.1$p.rank <- rank(-lex.1$p, ties.method = 'random')
    lex.1 <- subset(lex.1, p.rank <= p_threshold)
    lex.1 <- subset(lex.1, select = -c(p.rank))
  }
  lex.2 <- subset(lex.2, phone %in% unique(lex.1$phone))
  lex.1$p <- lex.1$freq / sum(lex.1$freq)
  lex.2$p <- lex.2$freq / sum(lex.2$freq)
  
  lex.1$position <- 'first'
  lex.2$position <- 'later'
  rbind(lex.1, lex.2)
}

cohort_H_m <- function(chrt, p_threshold = .95, nb_segs = F) {
  running_p <- 0
  chrt$running.p <- 0
  chrt$running.p <- 0
  for (cohort in chrt$cohort[order(chrt$cohort.p)]) {
    p <- chrt$cohort.p[match(cohort, chrt$cohort)]
    running_p <- running_p + p
    chrt$running.p[match(cohort, chrt$cohort)] <- running_p
  }
  chrt <- subset(chrt, running.p >= (1 - p_threshold))
  
  if (nb_segs) {
    m <- lm(data = chrt,
          entropy ~ log.prob + position)
  } else {
    m <- lm(data = chrt,
          entropy ~ log.prob + cohort.count + position)
  }
  m
}

pval_str <- function(p) {
  if (p < .001) {
    pval.str <- '***'
  } else if (p < .01) {
    pval.str <- '**'
  } else if (p <= .05) {
    pval.str <- '*'
  } else if (p <= .1) {
    pval.str <- '-'
  } else {
    pval.str <- 'n.s.'
  }
  pval.str
}

####
# test functions
# chapter 3
# neighborhood
tst_redundancy_1_a <- function(load_files = F) {
  if (!load_files) {
    all_langs <- lapply(lex.files, function(f) {
      d <- lf(paste0(lexdir, f))
      f <- sub('.lex.txt', '', f)
      print(f)
      
      if (grepl('Cantonese|Vietnamese', f)) {
        d$phones <- gsub('[0-9]', '', d$phones)
      }
      d <- med(d, min(d$word.length), max(d$word.length))
      #d <- minpair_location(d, min(d$word.length), max(d$word.length))
      d$language <- f
    
      write_tsv(x = d, path = paste0(outdir, 'rdn1/', f, '.txt'))
      1
    })
  }
  
  all_langs <- pblapply(lex.files, function(f) {
    f <- sub('.lex.txt', '', f)
    d <- read_tsv(paste0(outdir, 'rdn1/', f, '.txt'), col_types = cols())
    
    d$network.prop <- d$min.pairs / d$network
    
    d$scaled.mean.dist <- 0
    d$scaled.mean.dist.all <- 0
    d$scaled.mean.dist.nonmp <- 0
    d$scaled.min.pairs <- 0
    d$scaled.min.pairs.2 <- 0
    d$scaled.min.pairs.3 <- 0
    d$scaled.network <- 0
    d$scaled.network.prop <- 0
    for (i in min(d$word.length):max(d$word.length)) {
      d$scaled.mean.dist[d$word.length == i] <- 
          scale(d$mean.dist[d$word.length == i])
      
      d$scaled.mean.dist.all[d$word.length == i] <- 
        scale(d$mean.dist.all[d$word.length == i])
      
      d$scaled.mean.dist.nonmp[d$word.length == i] <-
          scale(d$mean.dist.nonmp[d$word.length == i])
      
      d$scaled.min.pairs[d$word.length == i] <- 
          scale(d$min.pairs[d$word.length == i])
      
      d$scaled.min.pairs.2[d$word.length == i] <- 
          scale(d$min.pairs.2[d$word.length == i])
      
      d$scaled.min.pairs.3[d$word.length == i] <- 
          scale(d$min.pairs.3[d$word.length == i])
      
      d$scaled.network[d$word.length == i] <- 
          scale(d$network[d$word.length == i])
      
      d$scaled.network.prop[d$word.length == i] <-
          scale(d$network.prop)
    }
    d[is.na(d)] <- 0
    
    d
  })
  all_langs <- ldply(all_langs, data.frame) 
  all_langs <- add_language_family(all_langs)
  all_langs <- subset(all_langs, word.length >= 2 & word.length <= 8)
  all_langs$scaled.log.prob <- -all_langs$scaled.unig.h
  
  plots.1 <- lapply(unique(all_langs$language), function(f) {
    d <- subset(all_langs, language == f)
    
    ct1 <- cor.test(d$scaled.min.pairs, d$scaled.log.prob)
    rho1 <- substr(as.character(ct1$estimate), 1, 5)
    print(paste0(f, ' - min pairs:'))
    print(paste0(ct1$estimate, ', p = ', ct1$p.value))
    
    x_range <- c(max(-2.5, min(d$scaled.log.prob)), min(2.5, max(d$scaled.log.prob)))
    y_range <- c(max(-2.5, min(d$scaled.min.pairs)), min(2.5, max(d$scaled.min.pairs)))
    pl1 <- ggplot(data = d, aes(x = scaled.log.prob, y = scaled.min.pairs)) +
      xlim(x_range) + ylim(y_range) +
      theme(text = element_text(size = 20)) +
      geom_point(data = subset(d, scaled.min.pairs > -4 & scaled.min.pairs < 4), alpha = .05) + 
      geom_smooth(method = 'lm', size = 2) + 
      ggtitle(f) +
      xlab('log Word Probability\n(z-scored)') + ylab('Neighborhood size\n(z-scored)')
    xadj <- layer_scales(pl1)$x$range$range[2] - .15 * (layer_scales(pl1)$x$range$range[2] - layer_scales(pl1)$x$range$range[1])
    yadj <- layer_scales(pl1)$y$range$range[2] - .15 * (layer_scales(pl1)$y$range$range[2] - layer_scales(pl1)$y$range$range[1])
    pl1 <- pl1 + geom_text(x = xadj, y = yadj, size = 7, 
                           label = paste0('r = ', rho1, '\n', pval_str(ct1$p.value)))
    pl1
  }) 
  
  plots.2 <- lapply(unique(all_langs$language), function(f) {
    d <- subset(all_langs, language == f)
    
    ct2 <- cor.test(d$scaled.mean.dist.nonmp, d$scaled.log.prob)
    rho2 <- substr(as.character(ct2$estimate), 1, 5)
      
    print(paste0(f, ' - mean dist:'))
    print(paste0(ct2$estimate, ', p = ', ct2$p.value))
    
    x_range <- c(max(-2.5, min(d$scaled.log.prob)), min(2.5, max(d$scaled.log.prob)))
    y_range <- c(max(-2.5, min(d$scaled.mean.dist.nonmp)), min(2.5, max(d$scaled.mean.dist.nonmp)))
    pl2 <- ggplot(data = d, aes(x = scaled.log.prob, y = scaled.mean.dist.nonmp)) +
      xlim(x_range) + ylim(y_range) +
      theme(text = element_text(size = 20)) +
      geom_point(data = subset(d, scaled.mean.dist.nonmp > -4 & scaled.mean.dist.nonmp < 4), alpha = .1) + 
      geom_smooth(method = 'lm', size = 2) + 
      ggtitle(f) +
      xlab('log Word Probability\n(z-scored)') + ylab('Mean Edit Dist.\n(z-scored)')
    xadj <- layer_scales(pl2)$x$range$range[2] - .15 * (layer_scales(pl2)$x$range$range[2] - layer_scales(pl2)$x$range$range[1])
    yadj <- layer_scales(pl2)$y$range$range[2] - .15 * (layer_scales(pl2)$y$range$range[2] - layer_scales(pl2)$y$range$range[1])
    
    pl2 <- pl2 + geom_text(x = xadj, y = yadj, size = 7, 
                           label = paste0('r = ', rho2, '\n', pval_str(ct2$p.value)))
    pl2
  })  
  
  plots.3 <- lapply(unique(all_langs$language), function(f) {
    d <- subset(all_langs, language == f)
    
    ct3 <- cor.test(d$scaled.network.prop, d$scaled.log.prob)
    rho3 <- substr(as.character(ct3$estimate), 1, 5)
    
    print(paste0(f, ' - network proportion:'))
    print(paste0(ct3$estimate, ', p = ', ct3$p.value))
  
    x_range <- c(max(-2.5, min(d$scaled.log.prob)), min(2.5, max(d$scaled.log.prob)))
    y_range <- c(max(-2.5, min(d$scaled.network.prop)), min(2.5, max(d$scaled.network.prop)))
    pl3 <- ggplot(data = d, aes(x = scaled.log.prob, y = scaled.network.prop)) +
      xlim(x_range) + ylim(y_range) +
      theme(text = element_text(size = 20)) +
      geom_point(data = subset(d, scaled.network.prop > -4 & scaled.network.prop < 4), alpha = .1) + 
      geom_smooth(method = 'lm', size = 2) + 
      ggtitle(f) +
      xlab('log Word Probability\n(z-scored)') + ylab('Rel. Nbrhood Size\n(z-scored)')
    xadj <- layer_scales(pl3)$x$range$range[2] - .15 * (layer_scales(pl3)$x$range$range[2] - layer_scales(pl3)$x$range$range[1])
    yadj <- layer_scales(pl3)$y$range$range[2] - .15 * (layer_scales(pl3)$y$range$range[2] - layer_scales(pl3)$y$range$range[1])
    pl3 <- pl3 + geom_text(x = xadj, y = yadj, size = 7, 
                           label = paste0('r = ', rho3, '\n', pval_str(ct3$p.value)))
    pl3
  })
  
  for (i in seq(1,length(unique(all_langs$language)), 5)) {
    print(paste0(i, ': min_pair'))
    pl <- ggarrange(plots.1[[i]], plots.1[[i + 1]], 
                    plots.1[[i + 2]], plots.1[[i + 3]],
                    plots.1[[i + 4]],
                      common.legend = TRUE, legend="right")

    ggsave(pl, filename = paste0(outdir, 'rdn1/minpair/', i, '.png'), width = 12, height = 7)
  }
  
  for (i in seq(1,length(unique(all_langs$language)), 5)) {
    print(paste0(i, ': mean_dist'))
    pl <- ggarrange(plots.2[[i]], plots.2[[i + 1]], 
                    plots.2[[i + 2]], plots.2[[i + 3]],
                    plots.2[[i + 4]],
                      common.legend = TRUE, legend="right")  
    ggsave(pl, filename = paste0(outdir, 'rdn1/meandist/', i, '.png'), width = 12, height = 7)
  }
  
  for (i in seq(1,length(unique(all_langs$language)), 5)) {
    print(paste0(i, ': network'))
    pl <- ggarrange(plots.3[[i]], plots.3[[i + 1]], 
                    plots.3[[i + 2]], plots.3[[i + 3]],
                    plots.3[[i + 4]],
                    common.legend = TRUE, legend="right")  
    ggsave(pl, filename = paste0(outdir, 'rdn1/network/', i, '.png'), width = 12, height = 7)
  }
  
  # LME's
  all_langs$log.prob <- -all_langs$unig.h
  all_langs$network.prop <- all_langs$min.pairs / all_langs$network
  
  all_langs$log.prob.all <- 0
  all_langs$min.pairs.all <- 0
  all_langs$network.prop.all <- 0
  all_langs$mean.dist.all <- 0
  for (i in 2:8) {
    all_langs$log.prob.all[all_langs$word.length == i] <-
      scale(all_langs$log.prob[all_langs$word.length == i])
    all_langs$min.pairs.all[all_langs$word.length == i] <- 
      scale(all_langs$min.pairs[all_langs$word.length == i])
    all_langs$network.prop.all[all_langs$word.length == i] <- 
      scale(all_langs$network.prop[all_langs$word.length == i])
    all_langs$mean.dist.all[all_langs$word.length == i] <- 
      scale(all_langs$mean.dist.nonmp[all_langs$word.length == i])
  }
  
  m.minpair <- lmer(data = all_langs,
                    min.pairs.all ~ log.prob.all +
                      (1|family/language),
                    REML = F)
  
  lmer_str <- lmer_to_latex(m.minpair)
  lmer_str <- gsub('log.prob.all', 'log Word Prob.', lmer_str)
  write(x = lmer_str, file = paste0(outdir, 'tex_tables/rdn1/sig/minpair.tex'))
  
  m.np <- lmer(data = all_langs,
               network.prop.all ~ log.prob.all + 
                 (1|family/language),
               REML = F)
  lmer_str <- lmer_to_latex(m.np)
  lmer_str <- gsub('log.prob.all', 'log Word Prob.', lmer_str)
  write(x = lmer_str, file = paste0(outdir, 'tex_tables/rdn1/sig/network.tex'))
  
  m.med <- lmer(data = all_langs,
                mean.dist.all ~ log.prob.all + 
                  (1|family/language),
                REML = F)
  lmer_str <- lmer_to_latex(m.med)
  lmer_str <- gsub('log.prob.all', 'log Word Prob.', lmer_str)
  write(x = lmer_str, file = paste0(outdir, 'tex_tables/rdn1/sig/meandist.tex'))
  
  m.all <- lmer(data = all_langs,
                log.prob.all ~ min.pairs.all + network.prop.all + mean.dist.all +
                  (1|family/language),
                REML = F)
  lmer_str <- lmer_to_latex(m.all)
  lmer_str <- gsub('min.pairs.all', 'Neighborhood Size', lmer_str)
  lmer_str <- gsub('network.prop.all', 'Rel. Nbrhood', lmer_str)
  lmer_str <- gsub('mean.dist.all', 'Mean Edit Dist.', lmer_str)
  write(x = lmer_str, file = paste0(outdir, 'tex_tables/rdn1/sig/all.tex'))
 
  1 
}

tst_redundancy_1_b <- function(load_files = F) {
  if (!load_files) {
    all_langs <- lapply(lex.files, function(f) {
      f <- sub('.lex.txt', '', f)
      d <- read_tsv(paste0(outdir, 'rdn1/', f, '.txt'), col_types = cols())
      d <- subset(d, word.length >= 2 & word.length <= 8)
      print(f)
    
      d$scaled.mean.dist <- 0
      d$scaled.mean.dist.all <- 0
      d$scaled.mean.dist.nonmp <- 0
      d$scaled.min.pairs <- 0
      d$scaled.min.pairs.2 <- 0
      d$scaled.min.pairs.3 <- 0
      d$scaled.network <- 0
      d$scaled.network.prop <- 0
      for (i in min(d$word.length):max(d$word.length)) {
        d$scaled.mean.dist[d$word.length == i] <- 
          scale(d$mean.dist[d$word.length == i])
      
        d$scaled.mean.dist.all[d$word.length == i] <- 
          scale(d$mean.dist.all[d$word.length == i])
      
        d$scaled.mean.dist.nonmp[d$word.length == i] <-
          scale(d$mean.dist.nonmp[d$word.length == i])
      
        d$scaled.min.pairs[d$word.length == i] <- 
          scale(d$min.pairs[d$word.length == i])
      
        d$scaled.min.pairs.2[d$word.length == i] <- 
          scale(d$min.pairs.2[d$word.length == i])
      
        d$scaled.min.pairs.3[d$word.length == i] <- 
          scale(d$min.pairs.3[d$word.length == i])
      
        d$scaled.network[d$word.length == i] <- 
          scale(d$network[d$word.length == i])
      
        d$scaled.network.prop[d$word.length == i] <-
          scale(d$min.pairs[d$word.length == i] / d$network[d$word.length == i])
      }
      d[is.na(d)] <- 0
      d$scaled.log.prob <- -d$scaled.unig.h
      
      ct.minpairs <- cor.test(d$scaled.log.prob, d$scaled.min.pairs)
      ct.networkprop <- cor.test(d$scaled.log.prob, d$scaled.network.prop)
      ct.med <- cor.test(d$scaled.log.prob, d$scaled.mean.dist.nonmp)
      
      rhos <- pblapply(1:nb_scramble, function(i) {
        d.scr <- shuffle_frequencies(d)
        for (j in 2:8) {
          d.scr$scaled.unig.h[d.scr$word.length == j] <-
            scale(d.scr$unig.h[d.scr$word.length == j])
        }
        d.scr$scaled.log.prob <- -d.scr$scaled.unig.h
        
        ct.scr.minpairs <- cor.test(d.scr$scaled.log.prob, d.scr$scaled.min.pairs)
        ct.scr.networkprop <- cor.test(d.scr$scaled.log.prob, d.scr$scaled.network.prop)
        ct.scr.med <- cor.test(d.scr$scaled.log.prob, d.scr$scaled.mean.dist.nonmp)
        
        rbind(
          data.frame(i = i + 1, lex.type = 'scr', test = 'min.pairs', 
                   b = ct.scr.minpairs$estimate, b0 = 0,
                   r = ct.scr.minpairs$estimate),
          data.frame(i = i + 1, lex.type = 'scr', test = 'network.prop', 
                   b = ct.scr.networkprop$estimate, b0 = 0,
                   r = ct.scr.networkprop$estimate),
          data.frame(i = i + 1, lex.type = 'scr', test = 'mean.edit', 
                   b = ct.scr.med$estimate, b0 = 0,
                   r = ct.scr.med$estimate)
          )
      }, cl = nb_cores)
      rhos <- rbind(
        rbind(
          data.frame(i = 1, lex.type = 'og', test = 'min.pairs', 
                   b = ct.minpairs$estimate, b0 = 0,
                   r = ct.minpairs$estimate),
          data.frame(i = 1, lex.type = 'og', test = 'network.prop', 
                   b = ct.networkprop$estimate, b0 = 0,
                   r = ct.networkprop$estimate),
          data.frame(i = 1, lex.type = 'og', test = 'mean.edit', 
                   b = ct.med$estimate, b0 = 0,
                   r = ct.med$estimate)
              ),
        ldply(rhos, data.frame)
      )
      rhos$language <- f
      write_tsv(x = rhos, path = paste0(outdir, 'rdn1/', f, '.novel.txt'))
      1
    })
  }
  all_langs <- lapply(lex.files, function(f) {
    f <- sub('.lex.txt', '', f)
    d <- read_tsv(paste0(outdir, 'rdn1/', f, '.novel.txt'), col_types = cols())
    d$r.z
    d$r.p <- 0
    for (i in c('min.pairs', 'network.prop')) {
      d$r.z[d$test == i] <- scale(d$r[d$test == i])
      d$r.p[d$test == i] <- max(.0001, 
                                sum(d$b[d$test == i & d$lex.type == 'og'] < d$b[d$test == i]) / length(d$b[d$test == i]))
    }
    d$r.z[d$test == 'mean.edit'] <- scale(d$r[d$test == 'mean.edit'])
    d$r.p[d$test == 'mean.edit'] <- max(.0001,
                                        sum(d$b[d$test == 'mean.edit' & d$lex.type == 'og'] > d$b[d$test == 'mean.edit']) / length(d$b[d$test == 'mean.edit']))
    
    d
  })
  all_langs <- ldply(all_langs, data.frame) 
  all_langs <- add_language_family(all_langs)
  
  unique_langs <- unique(all_langs$language)
  for (i in seq(1,length(unique_langs), 5)) {
    print(i)
    langs <- unique_langs[i:(i+4)]
    
    #### min pairs
    all_langs.sub <- subset(all_langs, language %in% langs & test == 'min.pairs')
    pl1 <- ggplot(data = all_langs.sub, aes(color = lex.type)) +
      geom_text(data = subset(all_langs.sub, lex.type == 'og'),
                x = Inf, y = Inf, 
                hjust   = 'right',
                vjust   = 'top',
                size = 6,
                color = 'black',
                aes(label = paste0('p <= ', substr(r.p, 1, 5)))) +
      geom_abline(data = sample_n(subset(all_langs.sub, lex.type == 'scr'), 
                                     max(50, nrow(subset(all_langs.sub, lex.type == 'scr')))),
                  size = .5, alpha = .3, linetype = 'dashed',
                  aes(slope = b, intercept = b0, color = lex.type)) +
      geom_abline(data = subset(all_langs.sub, lex.type == 'og'),
                  size = 1.5, alpha = 1, 
                  aes(slope = b, intercept = b0, color = lex.type)) +
      facet_wrap(~language, ncol = 3) +
      scale_x_continuous(limits = c(-2,2)) + xlab('log Word Probablity\n(z-scored)') +
      scale_y_continuous(limits = c(-1,1)) + ylab('Neighborhood Size\n(z-scored)') +
      scale_color_manual(name = 'Lexicon Type',
                         labels = c('Real-world', 'Probability-shuffled'),
                         values = c('#00BFC4', '#F8766D')) +
      theme(text = element_text(size=18), 
            legend.position = c(.95, .15), 
            legend.justification = c(1, 0))
    
    ggsave(paste0(outdir,'rdn1/shuffled/minpair/', i, '.png'), pl1,
           width = 9, height = 7)
    
    #### relative min pair
    all_langs.sub <- subset(all_langs, language %in% langs & test == 'network.prop')
    pl2 <- ggplot(data = all_langs.sub, aes(color = lex.type)) +
      geom_text(data = subset(all_langs.sub, lex.type == 'og'),
                x = Inf, y = Inf, 
                hjust   = 'right',
                vjust   = 'top',
                size = 6,
                color = 'black',
                aes(label = paste0('p <= ', substr(r.p, 1, 5)))) +
      geom_abline(data = sample_n(subset(all_langs.sub, lex.type == 'scr'), 
                                     max(50, nrow(subset(all_langs.sub, lex.type == 'scr')))),
                  size = .5, alpha = .3, linetype = 'dashed',
                  aes(slope = b, intercept = b0, color = lex.type)) +
      geom_abline(data = subset(all_langs.sub, lex.type == 'og'),
                  size = 1, alpha = 1, 
                  aes(slope = b, intercept = b0, color = lex.type)) +
      facet_wrap(~language, ncol = 3) +
      scale_x_continuous(limits = c(-2,2)) + xlab('log Word Probability\n(z-scored)') +
      scale_y_continuous(limits = c(-.5,.5)) + ylab('Relative Neighborhood Size\n(z-scored)') +
      scale_color_manual(name = 'Lexicon Type',
                         labels = c('Real-world', 'Probability-shuffled'),
                         values = c('#00BFC4', '#F8766D')) +
      theme(text = element_text(size=18), 
            legend.position = c(1, .15), 
            legend.justification = c(1, 0))
    
    ggsave(paste0(outdir,'rdn1/shuffled/network/', i, '.png'), pl2,
           width = 9, height = 7)
    
    #### mean edit dist
    all_langs.sub <- subset(all_langs, language %in% langs & test == 'mean.edit')
    pl3 <- ggplot(data = all_langs.sub, aes(color = lex.type)) +
      geom_text(data = subset(all_langs.sub, lex.type == 'og'),
                x = Inf, y = Inf, 
                hjust   = 'right',
                vjust   = 'top',
                size = 6,
                color = 'black',
                aes(label = paste0('p <= ', substr(r.p, 1, 5)))) +
      geom_abline(data = sample_n(subset(all_langs.sub, lex.type == 'scr'), 
                                     max(50, nrow(subset(all_langs.sub, lex.type == 'scr')))),
                  size = .5, alpha = .3, linetype = 'dashed',
                  aes(slope = b, intercept = b0, color = lex.type)) +
      geom_abline(data = subset(all_langs.sub, lex.type == 'og'),
                  size = 1.5, alpha = 1, 
                  aes(slope = b, intercept = b0, color = lex.type)) +
      facet_wrap(~language, ncol = 3) +
      scale_x_continuous(limits = c(-2,2)) + xlab('log Word Probability\n(z-scored)') +
      scale_y_continuous(limits = c(-.5,.5)) + ylab('Mean Edit Distance (excl. min. pairs)\n(z-scored)') +
      scale_color_manual(name = 'Lexicon Type',
                         labels = c('Real-world', 'Probability-shuffled'),
                         values = c('#00BFC4', '#F8766D')) +
      theme(text = element_text(size=18), 
            legend.position = c(1, .15), 
            legend.justification = c(1, 0))
    
    ggsave(paste0(outdir,'rdn1/shuffled/meandist/', i, '.png'), pl3,
           width = 9, height = 7)
  }
  
  all_langs$lex.type <- factor(all_langs$lex.type, levels = c('scr', 'og'))
  m.mp <- glmer(data = subset(all_langs, test == 'min.pairs'),
                lex.type ~ r + (1|family/language), family = binomial)
  lmer_str <- lmer_to_latex(m.mp)
  lmer_str <- gsub(' r ', ' $r$ ', lmer_str)
  write(x = lmer_str, file = paste0(outdir, 'tex_tables/rdn1/shuffled/minpair.tex'))
  
  m.np <- glmer(data = subset(all_langs, test == 'network.prop'),
                lex.type ~ r + (1| family/language), family = binomial)
  lmer_str <- lmer_to_latex(m.np)
  lmer_str <- gsub(' r ', ' $r$ ', lmer_str)
  write(x = lmer_str, file = paste0(outdir, 'tex_tables/rdn1/shuffled/network.tex'))
  
  m.med <- glmer(data = subset(all_langs, test == 'mean.edit'),
                lex.type ~ r + (1| family/language), family = binomial)
  lmer_str <- lmer_to_latex(m.med)
  lmer_str <- gsub(' r ', ' $r$ ', lmer_str)
  write(x = lmer_str, file = paste0(outdir, 'tex_tables/rdn1/shuffled/meandist.tex'))
  
  1
}

# post-UP
tst_redundancy_2_a <- function(load_files = F) {
  if (!load_files) {
    all_langs <- lapply(lex.files, function(f) {
      d <- lf(paste0(lexdir, f), do_up = T)
      f <- sub('.lex.txt', '', f)
      print(f)
    
      lex <- lex2si(d)
      diphone <- diphone_h(lex, use_freq = F)
      # word redundancy with diphones
      wrwd <- diphone_redun(diphone)
      #wrwd <- subset(wrwd, up <= word.length)
    
      wrwd$mean.redun.seg.h <- wrwd$postup.sum.diphone.h / wrwd$segments.redundant 
      wrwd$seg.h.redun.prop <- wrwd$postup.sum.diphone.h / wrwd$sum.diphone.h
    
      wrwd$scaled.unig.h <- 0
      wrwd$scaled.redundancy.prop <- 0
      wrwd$scaled.mean.redun.seg.h <- 0
      wrwd$scaled.seg.h.redun.prop <- 0
      for (i in 1:max(wrwd$word.length)) {
        wrwd$scaled.unig.h[wrwd$word.length == i] <- 
          scale(wrwd$unig.h[wrwd$word.length == i])
        wrwd$scaled.redundancy.prop[wrwd$word.length == i] <- 
          scale(wrwd$redundancy.prop[wrwd$word.length == i])
        wrwd$scaled.mean.redun.seg.h[wrwd$word.length == i] <- 
          scale(wrwd$mean.redun.seg.h[wrwd$word.length == i])
        wrwd$scaled.seg.h.redun.prop[wrwd$word.length == i] <- 
          scale(wrwd$seg.h.redun.prop[wrwd$word.length == i])
      }
      wrwd$language <- f
      wrwd[is.na(wrwd)] <- 0.
      
      if (grepl('Vietnamese|Cantonese', f)) {
        wrwd$word.length <- nchar(gsub('[0-9]', '', wrwd$phones))
      }
      write_tsv(x = wrwd, path = paste0(outdir, 'rdn2/', f, '.txt'))
      1
    })
  }
  all_langs <- pblapply(lex.files, function(f) {
      f <- sub('.lex.txt', '', f)
      d <- read_tsv(paste0(outdir, 'rdn2/', f, '.txt'), col_types = cols())
      d
  })
  all_langs <- ldply(all_langs, data.frame)
  all_langs <- add_language_family(all_langs)
  all_langs$log.prob <- -all_langs$unig.h
  all_langs$scaled.log.prob <- -all_langs$scaled.unig.h
  all_langs$bigram.prob <- -all_langs$mean.redun.seg.h
  all_langs$scaled.bigram.prob <- -all_langs$scaled.mean.redun.seg.h
  
  # all together
  pl <- ggplot(data = all_langs, aes(x = scaled.log.prob, y = scaled.redundancy.prop)) +
    scale_x_continuous(breaks = -3:5) +
    theme(text = element_text(size = 20)) +
    geom_smooth(method = 'lm', size = 1.5) + 
    xlab('log Word Probability\n(z-scored)') + ylab('Redundant Segments\n(z-scored)') +
    facet_wrap(~language)
  ggsave(pl, filename = paste0(outdir, 'rdn2/redun_prop/all.png'), width = 11, height = 9)
  
  plots.1 <- lapply(unique(all_langs$language), function(f) {
    wrwd <- subset(all_langs, language == f)
    print(paste0(f, ' - redundant segments:'))
    
    ct1 <- cor.test(wrwd$scaled.log.prob, wrwd$scaled.redundancy.prop)
    rho1 <- substr(as.character(ct1$estimate), 1, 5)
    print(paste0(ct1$estimate, ', p = ', ct1$p.value))
      
    
    x_range <- c(max(-2.5, min(wrwd$scaled.log.prob)), min(2.5, max(wrwd$scaled.log.prob)))
    y_range <- c(max(-2.5, min(wrwd$scaled.redundancy.prop)), min(2.5, max(wrwd$scaled.redundancy.prop)))
    pl1 <- ggplot(data = wrwd, aes(x = scaled.log.prob, y = scaled.redundancy.prop)) +
      xlim(x_range) + ylim(y_range) +
      theme(text = element_text(size = 20)) +
      geom_point(alpha = .05) + geom_smooth(method = 'lm', size = 2) + 
      ggtitle(f) +
      xlab('log Word Probability\n(z-scored)') + ylab('Redundant Segments\n(z-scored)')
    xadj <- layer_scales(pl1)$x$range$range[2] - .15 * (layer_scales(pl1)$x$range$range[2] - layer_scales(pl1)$x$range$range[1])
    yadj <- layer_scales(pl1)$y$range$range[2] - .15 * (layer_scales(pl1)$y$range$range[2] - layer_scales(pl1)$y$range$range[1])
    pl1 <- pl1 + geom_text(x = xadj, y = yadj, size = 7, 
                           label = paste0('r = ', rho1, '\n', pval_str(ct1$p.value)))
    pl1
  }) 
  
  plots.2 <- lapply(unique(all_langs$language), function(f) {
    wrwd <- subset(all_langs, language == f & up < word.length)
    # re-scale after subset
    for (i in min(wrwd$word.length):max(wrwd$word.length)) {
      if (nrow(subset(wrwd, word.length == i)) > 1) {
        wrwd$scaled.log.prob[wrwd$word.length == i] <-
          scale(wrwd$log.prob[wrwd$word.length == i])
        wrwd$scaled.bigram.prob[wrwd$word.length == i] <-
          scale(wrwd$bigram.prob[wrwd$word.length == i])
      }
    }
    
    ct2 <- cor.test(wrwd$scaled.log.prob, wrwd$scaled.bigram.prob)
    rho2 <- substr(as.character(ct2$estimate), 1, 5)
    print(paste0(f, ' - mean log bigram prob.'))
    print(paste0(ct2$estimate, ', p = ', ct2$p.value))
    
    x_range <- c(max(-2.5, min(wrwd$scaled.log.prob)), min(2.5, max(wrwd$scaled.log.prob)))
    y_range <- c(max(-2.5, min(wrwd$scaled.bigram.prob)), min(2.5, max(wrwd$scaled.bigram.prob)))
    pl2 <- ggplot(data = wrwd, aes(x = scaled.log.prob, y = scaled.bigram.prob)) +
      xlim(x_range) + ylim(y_range) +
      theme(text = element_text(size = 20)) +
      geom_point(alpha = .05) + geom_smooth(method = 'lm', size = 2) + 
      ggtitle(f) +
      xlab('log Word Probability\n(z-scored)') + ylab('Mean log Bigram Prob.\n(z-scored)')
    xadj <- layer_scales(pl2)$x$range$range[2] - .15 * (layer_scales(pl2)$x$range$range[2] - layer_scales(pl2)$x$range$range[1])
    yadj <- layer_scales(pl2)$y$range$range[2] - .15 * (layer_scales(pl2)$y$range$range[2] - layer_scales(pl2)$y$range$range[1])
    pl2 <- pl2 + geom_text(x = xadj, y = yadj, size = 7, 
                           label = paste0('r = ', rho2, '\n', pval_str(ct2$p.value)))
    pl2
  }) 
  
  plots.3 <- lapply(unique(all_langs$language), function(f) {
    wrwd <- subset(all_langs, language == f & up < word.length)
    # re-scale after subset
    for (i in min(wrwd$word.length):max(wrwd$word.length)) {
      if (nrow(subset(wrwd, word.length == i)) > 1) {
        wrwd$scaled.log.prob[wrwd$word.length == i] <-
          scale(wrwd$log.prob[wrwd$word.length == i])
        wrwd$scaled.seg.h.redun.prop[wrwd$word.length == i] <-
          scale(wrwd$seg.h.redun.prop[wrwd$word.length == i])
      }
    }
    ct3 <- cor.test(wrwd$scaled.log.prob, wrwd$scaled.seg.h.redun.prop)
    rho3 <- substr(as.character(ct3$estimate), 1, 5)
    print(paste0(f, ' - ratio of total seg info:'))
    print(paste0(ct3$estimate, ', p = ', ct3$p.value))
    
    x_range <- c(max(-2.5, min(wrwd$scaled.log.prob)), min(2.5, max(wrwd$scaled.log.prob)))
    y_range <- c(max(-2.5, min(wrwd$scaled.seg.h.redun.prop)), min(2.5, max(wrwd$scaled.seg.h.redun.prop)))
    pl3 <- ggplot(data = wrwd, aes(x = scaled.log.prob, y = scaled.seg.h.redun.prop)) +
      xlim(x_range) + ylim(y_range) +
      theme(text = element_text(size = 20)) +
      geom_point(alpha = .05) + geom_smooth(method = 'lm', size = 2) + 
      ggtitle(f) +
      xlab('log Word Probability\n(z-scored)') + ylab('Rel. Phonotactic Prob.\n(z-scored)')
    xadj <- layer_scales(pl3)$x$range$range[2] - .15 * (layer_scales(pl3)$x$range$range[2] - layer_scales(pl3)$x$range$range[1])
    yadj <- layer_scales(pl3)$y$range$range[2] - .15 * (layer_scales(pl3)$y$range$range[2] - layer_scales(pl3)$y$range$range[1])
    pl3 <- pl3 + geom_text(x = xadj, y = yadj, size = 7, 
                           label = paste0('r = ', rho3, '\n', pval_str(ct3$p.value)))
    pl3
  }) 
  
  for (i in seq(1,length(unique(all_langs$language)), 5)) {
    print(paste0(i, ': redun_prop'))
    pl <- ggarrange(plots.1[[i]], plots.1[[i + 1]], 
                    plots.1[[i + 2]], plots.1[[i + 3]], 
                    plots.1[[i + 4]], 
                      common.legend = TRUE, legend="right")
    
    ggsave(pl, filename = paste0(outdir, 'rdn2/redun_prop/', i, '.png'), width = 12, height = 7)
  }
  
  for (i in seq(1,length(unique(all_langs$language)), 5)) {
    print(paste0(i, ': mean_segh'))
    pl <- ggarrange(plots.2[[i]], plots.2[[i + 1]], 
                    plots.2[[i + 2]], plots.2[[i + 3]], 
                    plots.2[[i + 4]], 
                      common.legend = TRUE, legend="right")
    
    ggsave(pl, filename = paste0(outdir, 'rdn2/bigram_prob/', i, '.png'), width = 12, height = 7)
  }
  
  for (i in seq(1,length(unique(all_langs$language)), 5)) {
    print(paste0(i, ': segh_prop'))
    pl <- ggarrange(plots.3[[i]], plots.3[[i + 1]], 
                    plots.3[[i + 2]], plots.3[[i + 3]], 
                    plots.3[[i + 4]], 
                      common.legend = TRUE, legend="right")
      
    ggsave(pl, filename = paste0(outdir, 'rdn2/segh_prop/', i, '.png'), width = 12, height = 7)
  }
  
  all_langs$log.prob.all <- 0
  all_langs$redundancy.prop.all <- 0
  all_langs$bigram.prob.all <- 0
  all_langs$seg.h.redun.prop.all <- 0
  for (i in min(all_langs$word.length):max(all_langs$word.length)) {
    all_langs$log.prob.all[all_langs$word.length == i] <- 
      scale(all_langs$log.prob[all_langs$word.length == i])
    all_langs$redundancy.prop.all[all_langs$word.length == i] <- 
      scale(all_langs$redundancy.prop[all_langs$word.length == i])
    all_langs$bigram.prob.all[all_langs$word.length == i] <- 
      scale(all_langs$bigram.prob[all_langs$word.length == i])
    all_langs$seg.h.redun.prop.all[all_langs$word.length == i] <- 
      scale(all_langs$seg.h.redun.prop[all_langs$word.length == i])
  }
  
  m.redun_prop <- lmer(data = all_langs,
                    redundancy.prop.all ~ log.prob.all + (1|family/language),
                    REML = F)
  lmer_str <- lmer_to_latex(m.redun_prop)
  lmer_str <- gsub('log.prob.all', 'log Word Prob.', lmer_str)
  write(x = lmer_str, file = paste0(outdir, 'tex_tables/rdn2/sig/redun_prop.tex'))
  
  all_langs.sub <- subset(all_langs, up < word.length)
  for (i in min(all_langs$word.length):max(all_langs$word.length)) {
    all_langs.sub$log.prob.all[all_langs.sub$word.length == i] <- 
      scale(all_langs.sub$log.prob[all_langs.sub$word.length == i])
    all_langs.sub$redundancy.prop.all[all_langs.sub$word.length == i] <- 
      scale(all_langs.sub$redundancy.prop[all_langs.sub$word.length == i])
    all_langs.sub$bigram.prob.all[all_langs.sub$word.length == i] <- 
      scale(all_langs.sub$bigram.prob[all_langs.sub$word.length == i])
    all_langs.sub$seg.h.redun.prop.all[all_langs.sub$word.length == i] <- 
      scale(all_langs.sub$seg.h.redun.prop[all_langs.sub$word.length == i])
  }
  
  m.bigram_prob <- lmer(data = all_langs.sub,
                       bigram.prob.all ~ log.prob.all + (1|family/language),
                       REML = F)
  lmer_str <- lmer_to_latex(m.bigram_prob)
  lmer_str <- gsub('log.prob.all', 'log Word Prob.', lmer_str)
  write(x = lmer_str, file = paste0(outdir, 'tex_tables/rdn2/sig/bigram_prob.tex'))

  m.segh_prop <- lmer(data = all_langs.sub,
                      seg.h.redun.prop.all ~ log.prob.all + (1|family/language),
                      REML = F)
  lmer_str <- lmer_to_latex(m.segh_prop)
  lmer_str <- gsub('log.prob.all', 'log Word Prob.', lmer_str)
  write(x = lmer_str, file = paste0(outdir, 'tex_tables/rdn2/sig/segh_prop.tex'))
  
  m.all <- lmer(data = all_langs.sub,
                log.prob.all ~ redundancy.prop.all + bigram.prob.all + seg.h.redun.prop.all +
                  (1|family/language),
                REML = F)
  lmer_str <- lmer_to_latex(m.all)
  lmer_str <- gsub('redundancy.prop.all', 'Num. Redundant Segs.', lmer_str)
  lmer_str <- gsub('bigram.prob.all', 'Mean Phonotactic Prob. in Tail', lmer_str)
  lmer_str <- gsub('seg.h.redun.prop.all', 'Rel. Phonotactic Prob. in Tail', lmer_str)
  write(x = lmer_str, file = paste0(outdir, 'tex_tables/rdn2/sig/all.tex'))

  1
}

tst_redundancy_2_b <- function(load_files = F) {
  if (!load_files) {
    all_langs <- lapply(lex.files, function(f) {
      f <- sub('.lex.txt', '', f)
      print(f)
      
      d.all <- read_tsv(paste0(outdir, 'rdn2/', f, '.txt'), col_types = cols())
      d.all$bigram.prob <- -d.all$mean.redun.seg.h
      d.all$scaled.bigram.prob <- -d.all$scaled.mean.redun.seg.h
      d.all$log.prob <- -d.all$unig.h
      d.all$scaled.log.prob <- -d.all$scaled.unig.h
      d <- subset(d.all, up < word.length)
      
      # have to re-scale for subset of words with at least 1 redundant seg
      for (j in 1:max(d.all$word.length)) {
        if (nrow(subset(d.all, word.length == j)) > 1) {
          d$scaled.unig.h[d$word.length == j] <-
            scale(d$unig.h[d$word.length == j])
          d$scaled.redundancy.prop[d$word.length == i] <-
            scale(d$redundancy.prop[d$word.length == i])
        }
      }
      d.all$scaled.log.prob <- -d.all$scaled.unig.h
      
      for (j in 1:max(d$word.length)) {
        if (nrow(subset(d, word.length == j)) > 1) {
          d$scaled.unig.h[d$word.length == j] <-
            scale(d$unig.h[d$word.length == j])
          d$scaled.bigram.prob[d$word.length == i] <-
            scale(d$bigram.prob[d$word.length == i])
          d$scaled.seg.h.redun.prop[d$word.length == i] <-
            scale(d$seg.h.redun.prop[d$word.length == i])
        }
      }
      d$scaled.log.prob <- -d$scaled.unig.h
      
      ct1 <- cor.test(d.all$scaled.log.prob, d.all$scaled.redundancy.prop)
      ct2 <- cor.test(d$scaled.log.prob, d$scaled.bigram.prob)
      ct3 <- cor.test(d$scaled.log.prob, d$scaled.seg.h.redun.prop)
      
      rhos <- pblapply(1:nb_scramble, function(i) {
        d.all.scr <- shuffle_frequencies(d.all)
        d.scr <- subset(d.all.scr, up < word.length)
        
        for (j in 1:max(d.all.scr$word.length)) {
          if (nrow(subset(d.all.scr, word.length == j)) > 1) {
            d.all.scr$scaled.unig.h[d.all.scr$word.length == j] <-
              scale(d.all.scr$unig.h[d.all.scr$word.length == j])
          }
        }
        d.all.scr$scaled.log.prob <- -d.all.scr$scaled.unig.h
        # need to re-scale for subset of words w/ UP's 
        for (j in 1:max(d.scr$word.length)) {
          if (nrow(subset(d.scr, word.length == j)) > 1) {
            d.scr$scaled.unig.h[d.scr$word.length == j] <-
              scale(d.scr$unig.h[d.scr$word.length == j])
            d.scr$scaled.bigram.prob[d.scr$word.length == i] <-
              scale(d.scr$bigram.prob[d.scr$word.length == i])
            d.scr$scaled.seg.h.redun.prop[d.scr$word.length == i] <-
              scale(d.scr$seg.h.redun.prop[d.scr$word.length == i])
          }
        }
        d.scr$scaled.log.prob <- -d.scr$scaled.unig.h
        
        ct1.scr <- cor.test(d.all.scr$scaled.log.prob, d.all.scr$scaled.redundancy.prop)
        ct2.scr <- cor.test(d.scr$scaled.log.prob, d.scr$scaled.bigram.prob)
        ct3.scr <- cor.test(d.scr$scaled.log.prob, d.scr$scaled.seg.h.redun.prop)
        
        rbind(data.frame(r = ct1.scr$estimate,
                         b = ct1.scr$estimate, b0 = 0,
                         lex.type = 'scr', test = 'redun.prop', i = i + 1),
              data.frame(r = ct2.scr$estimate,
                         b = ct2.scr$estimate, b0 = 0,
                         lex.type = 'scr', test = 'bigram.prob', i = i + 1),
              data.frame(r = ct3.scr$estimate,
                         b = ct3.scr$estimate, b0 = 0,
                         lex.type = 'scr', test = 'seg.h.prop', i = i + 1)
              )
      }, cl = nb_cores)
      rhos <- rbind(data.frame(r = ct1$estimate,
                               b = ct1$estimate, b0 = 0,
                               #b = m1$coefficients[2], b0 = m1$coefficients[1],
                                     lex.type = 'og', test = 'redun.prop', i = 1),
                    data.frame(r = ct2$estimate,
                               b = ct2$estimate, b0 = 0,
                               #b = m2$coefficients[2], b0 = m2$coefficients[1],
                               lex.type = 'og', test = 'bigram.prob', i = 1),
                    data.frame(r = ct3$estimate, 
                               b = ct3$estimate, b0 = 0,
                               #b = m3$coefficients[2], b0 = m3$coefficients[1],
                               lex.type = 'og', test = 'seg.h.prop', i = 1),
                    ldply(rhos, data.frame))
      
      rhos$language <- f
      write_tsv(x = rhos, path = paste0(outdir, 'rdn2/', f, '.novel.txt'))
      1
    })
  }
  all_langs <- lapply(lex.files, function(f) {
    f <- sub('.lex.txt', '', f)
    d <- read_tsv(paste0(outdir, 'rdn2/', f, '.novel.txt'), col_types = cols())
    d$r.p <- 0
    for (i in c('redun.prop', 'seg.h.prop')) {
      d$r.p[d$test == i] <- max(.0001, 
                                sum(d$b[d$test == i & d$lex.type == 'og'] > d$b[d$test == i]) / length(d$b[d$test == i]))
    }
    d$r.p[d$test == 'bigram.prob'] <- max(.0001,
                                        sum(d$b[d$test == 'bigram.prob' & d$lex.type == 'og'] < d$b[d$test == 'bigram.prob']) / length(d$r[d$test == 'bigram.prob']))
    d
  })
  all_langs <- ldply(all_langs, data.frame)
  all_langs <- add_language_family(all_langs)

  unique_langs <- unique(all_langs$language)
  for (i in seq(1,length(unique_langs), 5)) {
    print(i)
    langs <- unique_langs[i:(i+4)]
    
    # relative up positions
    all_langs.sub <- subset(all_langs, language %in% langs & test == 'redun.prop')
    pl1 <- ggplot(data = all_langs.sub, aes(color = lex.type)) +
      geom_text(data = subset(all_langs.sub, lex.type == 'og'),
                x = Inf, y = Inf, 
                hjust   = 'right',
                vjust   = 'top',
                size = 6,
                color = 'black',
                aes(label = paste0('p <= ', substr(r.p, 1, 5)))) +
      geom_abline(data = sample_n(subset(all_langs.sub, lex.type == 'scr'),
                                  max(50, nrow(subset(all_langs.sub, lex.type == 'scr')))),
                  size = .5, alpha = .3, linetype = 'dashed',
                  aes(slope = b, intercept = b0, color = lex.type)) +
      geom_abline(data = subset(all_langs.sub, lex.type == 'og'),
                  size = 1.5, alpha = 1, 
                  aes(slope = b, intercept = b0, color = lex.type)) +
      facet_wrap(~language, ncol = 3) +
      scale_x_continuous(limits = c(-2,2)) + xlab('log Word Probability\n(z-scored)') +
      scale_y_continuous(limits = c(-.6,.6)) + ylab('Redundant Segments\n(z-scored)') +
      scale_color_manual(name = 'Lexicon Type',
                         labels = c('Real-world', 'Probability-shuffled'), 
                         values = c('#00BFC4', '#F8766D')) +
      theme(text = element_text(size=18), 
            legend.position = c(1, .15), 
            legend.justification = c(1, 0))
    
    ggsave(paste0(outdir, 'rdn2/shuffled/redun_prop/', i, '.png'), pl1,
           width = 9, height = 7)
    
    # mean seg h in tail
    all_langs.sub <- subset(all_langs, language %in% langs & test == 'bigram.prob')
    pl2 <- ggplot(data = all_langs.sub, aes(color = lex.type)) +
      geom_text(data = subset(all_langs.sub, lex.type == 'og'),
                x = Inf, y = Inf, 
                hjust   = 'right',
                vjust   = 'top',
                size = 6,
                color = 'black',
                aes(label = paste0('p <= ', substr(r.p, 1, 5)))) +
      geom_abline(data = sample_n(subset(all_langs.sub, lex.type == 'scr'),
                                  max(50, nrow(subset(all_langs.sub, lex.type == 'scr')))),
                  size = .5, alpha = .3, linetype = 'dashed',
                  aes(slope = b, intercept = b0, color = lex.type)) +
      geom_abline(data = subset(all_langs.sub, lex.type == 'og'),
                  size = 1.5, alpha = 1, 
                  aes(slope = b, intercept = b0, color = lex.type)) +
      facet_wrap(~language, ncol = 3) +
      scale_x_continuous(limits = c(-2,2)) + xlab('log Word Probability\n(z-scored)') +
      scale_y_continuous(limits = c(-.6,.6)) + ylab('Mean log Bigram Prob.\n(z-scored)') +
      scale_color_manual(name = 'Lexicon Type',
                         labels = c('Real-world', 'Probability-shuffled'),
                         values = c('#00BFC4', '#F8766D')) +
      theme(text = element_text(size=18), 
            legend.position = c(1, .15), 
            legend.justification = c(1, 0))
    
    ggsave(paste0(outdir, 'rdn2/shuffled/bigram_prob/', i, '.png'), pl2,
           width = 9, height = 7)
    
    # seg h in tail proportion 
    all_langs.sub <- subset(all_langs, language %in% langs & test == 'seg.h.prop')
    pl3 <- ggplot(data = all_langs.sub, aes(color = lex.type)) +
      geom_text(data = subset(all_langs.sub, lex.type == 'og'),
                x = Inf, y = Inf, 
                hjust   = 'right',
                vjust   = 'top',
                size = 6,
                color = 'black',
                aes(label = paste0('p <= ', substr(r.p, 1, 5)))) +
      geom_abline(data = sample_n(subset(all_langs.sub, lex.type == 'scr'),
                                  max(50, nrow(subset(all_langs.sub, lex.type == 'scr')))),
                  size = .5, alpha = .3, linetype = 'dashed',
                  aes(slope = b, intercept = b0, color = lex.type)) +
      geom_abline(data = subset(all_langs.sub, lex.type == 'og'),
                  size = 1.5, alpha = 1, 
                  aes(slope = b, intercept = b0, color = lex.type)) +
      facet_wrap(~language, ncol = 3) +
      scale_x_continuous(limits = c(-2,2)) + xlab('log Word Probability\n(z-scored)') +
      scale_y_continuous(limits = c(-.6,.6)) + ylab('Rel. Phonotactic Prob.\n(z-scored)') +
      scale_color_manual(name = 'Lexicon Type',
                         labels = c('Real-world', 'Probability-shuffled'),
                         values = c('#00BFC4', '#F8766D')) +
      theme(text = element_text(size=18), 
            legend.position = c(1, .15), 
            legend.justification = c(1, 0))
    
    ggsave(paste0(outdir, 'rdn2/shuffled/segh_prop/', i, '.png'), pl3,
           width = 9, height = 7)
  }
  
  all_langs$lex.type <- factor(all_langs$lex.type, levels = c('scr', 'og'))
  m.redun_prop <- glmer(data = subset(all_langs, test == 'redun.prop'),
                        lex.type ~ r  + (1|family/language), family = binomial)
  lmer_str <- lmer_to_latex(m.redun_prop)
  lmer_str <- gsub(' r ', ' $r$ ', lmer_str)
  write(x = lmer_str, file = paste0(outdir, 'tex_tables/rdn2/shuffled/redun_prop.tex'))
  
  m.bigram_prob <- glmer(data = subset(all_langs, test == 'bigram.prob'),
                        lex.type ~ r + (1|family/language), family = binomial)
  lmer_str <- lmer_to_latex(m.bigram_prob)
  lmer_str <- gsub(' r ', ' $r$ ', lmer_str)
  write(x = lmer_str, file = paste0(outdir, 'tex_tables/rdn2/shuffled/bigram_prob.tex'))
  
  m.seg_h_prop <- glmer(data = subset(all_langs, test == 'seg.h.prop'),
                         lex.type ~ r + (1|family/language), family = binomial)
  lmer_str <- lmer_to_latex(m.seg_h_prop)
  lmer_str <- gsub(' r ', ' $r$ ', lmer_str)
  write(x = lmer_str, file = paste0(outdir, 'tex_tables/rdn2/shuffled/seg_h_prop.tex'))
  

  1
}

# chapter 4
# biphone balance
tst_balanced_sig <- function(p_threshold = .95) {
  all_langs <- pblapply(lex.files, function(f) {
    d <- lf(paste0(lexdir, f), do_up = F)
    d <- subset(d, word.length > 1)
    f <- sub('.lex.txt', '', f)
    
    first.phones <- first_dist(d, p_threshold = p_threshold, N = 2)
    first.phones$log.mean.z <- z.(log(first.phones$mean))
    first.phones$log.count.z <- z.(log(first.phones$count))
    first.phones$language <- f
    first.phones 
  }, cl = nb_cores)
  all_langs <- ldply(all_langs, data.frame)
  
  pl <- ggplot(data = all_langs, aes(x = log.mean.z, y = log.count.z)) + 
    theme(text = element_text(size=18)) +
    geom_smooth(method = 'lm', color = 'blue') + 
    xlab('log Mean Token Freq.\n(z-scored)') + ylab('log Type Count\n(z-scored)') +
    scale_x_continuous(breaks = c(-2,0,2,4), labels = c(-2,0,2,4)) +
    facet_wrap(~language)
  ggsave(pl, filename = paste0(outdir, 'bal/sig/all.png'), 
         width = 9, height = 8)
  
  plots <- lapply(unique(all_langs$language), function(f) {
    print(paste0(f, ' - balanced (sig)'))
    first.phones <- subset(all_langs, language == f)
    ct <- cor.test(log(first.phones$mean), log(first.phones$count))
    print(paste0(ct$estimate, ', p = ', ct$p.value))
    rho <- sub('0.','.', substr(as.character(ct$estimate), 1, 6))
    
    pl <- ggplot(data = first.phones, aes(x = log(mean), y = log(count), label = phone)) + 
      theme(text = element_text(size=18)) +
      geom_smooth(data = subset(first.phones), 
                  method = 'lm', color = 'blue') + 
      geom_text(alpha = .6) + 
      ggtitle(f) +
      xlab('log Mean Token Freq.') + ylab('log Type Count')
    xadj <- layer_scales(pl)$x$range$range[2] - .15 * (layer_scales(pl)$x$range$range[2] - layer_scales(pl)$x$range$range[1])
    yadj <- layer_scales(pl)$y$range$range[2] - .15 * (layer_scales(pl)$y$range$range[2] - layer_scales(pl)$y$range$range[1])
    pl <- pl + geom_text(x = xadj, y = yadj, size = 6, 
                         label = paste0('r = ', rho, '\n', pval_str(ct$p.value)))
    pl
  })
  
  for (i in seq(1,length(lex.files), 5)) {
    print(i)
    pl <- grid.arrange(plots[[i]], plots[[i + 1]], plots[[i + 2]], 
                       plots[[i + 3]], plots[[i + 4]], 
                       nrow = 2)  
    ggsave(pl, filename = paste0(outdir, 'bal/sig/', i, '.png'), 
           width = 10, height = 7)
  }
 
  all_langs <- add_language_family(all_langs)
  m <- lmer(data = all_langs, 
            log(count) ~ log(mean) + (1|family/language),
            REML = F)
  summary(m)
  lmer_str <- lmer_to_latex(m)
  lmer_str <- sub('log.mean.', 'log2 Mean Token Freq.', lmer_str)
  write(x = lmer_str, file = paste0(outdir, 'tex_tables/bal/sig.tex'))
  
  1 
}

tst_balanced_nonce <- function(load_files = F, p_threshold = .95) {
  if(!load_files) {
    lapply(lex.files, function(f) {
      d <- lf(paste0(lexdir, f), do_up = F)
      d <- subset(d, word.length > 1)
      f <- sub('.lex.txt', '', f)
      print(paste0(f, ' - balanced (nonce)'))
      
      fp <- first_dist(d, p_threshold = p_threshold, N = 2)
      ct <- cor.test(log(fp$mean), log(fp$count))
      m <- lm(data = fp, scale(log2(fp$count)) ~ scale(log2(fp$mean)))
      dn <- load_nonce(paste0(noncedir, f, '.new.txt'))

      nonce_rhos <- pblapply(1:nb_scramble, function(i) {
        
        # no need to keep length the same, we're only looking at first segments
        d.scr <- nonce_lexicon(d, dn, keep_length = F)
        d.scr <- subset(d.scr, word.length > 1)
        fp.scr <- first_dist(d.scr, p_threshold = p_threshold, N = 2)
        j <- 0
        # make sure the novel lexicon has more or less equal data points
        while (j < (nb_scramble/4) & nrow(fp) > nrow(fp.scr) + 2 & nrow(fp) < nrow(fp.scr) - 2) {
          d.scr <- nonce_lexicon(d, dn, keep_length = kl)
          fp.scr <- first_dist(d.scr, p_threshold = p_threshold, N = 2)
          j <- j + 1
        }
        
        ct.scr <- cor.test(log(fp.scr$mean), log(fp.scr$count))
        m.scr <- lm(data = fp.scr, scale(log(fp.scr$count)) ~ scale(log(fp.scr$mean)))
        
        data.frame(b0 = m.scr$coefficients[1], b1 = m.scr$coefficients[2],
                   r = ct.scr$estimate, R2 = summary(m.scr)$r.squared,
                 lex.type = 'nonce', i = i + 1)
      }, cl = nb_cores)
      nonce_rhos <- rbind(data.frame(b0 = m$coefficients[1], b1 = m$coefficients[2],
                                      r = ct$estimate, R2 = summary(m)$r.squared,
                                      lex.type = 'og', i = 1), 
                        ldply(nonce_rhos, data.frame))
      nonce_rhos$r.z <- scale(nonce_rhos$r)
      nonce_rhos$R2.z <- scale(nonce_rhos$R2)
      nonce_rhos$og.greater.r <- max(.0001, sum(nonce_rhos$r[1] > nonce_rhos$r) / nb_scramble) 
      nonce_rhos$og.greater.R2 <- sum(nonce_rhos$R2[1] > nonce_rhos$R2) / nb_scramble
      nonce_rhos$language <- f
      
      print(paste0(nonce_rhos$r[1], ' : ', mean(nonce_rhos$r), 
                                                '. p < ', sum(nonce_rhos$r[1] > nonce_rhos$r) / nb_scramble))
      write_tsv(nonce_rhos, paste0(outdir, out_folder, '/', f, '.novel.txt'))
      1  
    })
  }
  all_langs <- lapply(lex.files, function(f) {
    f <- sub('.lex.txt', '', f)
    d <- read_tsv(paste0(outdir, 'bal/', f, '.novel.txt'), col_types = cols())
    d
  })
  all_langs <- ldply(all_langs, data.frame)
  all_langs$lex.type <- factor(all_langs$lex.type, levels = c('og', 'nonce'))
  
  ## all langs
  pl <- ggplot(data = all_langs, aes(color = lex.type)) +
    xlim(c(-2,2)) + xlab('log Mean Frequency\n(z-scored)') +
    ylim(c(-2,2)) + ylab('log Type Count\n(z-scored)') +
    geom_abline(data = subset(all_langs, lex.type == 'nonce'), 
                size = .25, alpha = .15, linetype = 'dashed',
                aes(slope = b1, intercept = b0, color = lex.type)) +
    geom_abline(data = subset(all_langs, lex.type == 'og'),
                size = 1.75, alpha = 1,
                aes(slope = b1, intercept = b0, color = lex.type)) +
    facet_wrap(~language) +
    scale_color_manual(name = 'Lexicon Type',
                       labels = c('Nonce Forms', 'Real-world'),
                       values = c('#00BFC4', '#F8766D')) +
    theme(text = element_text(size=20),
          legend.position='bottom') 
  
  ggsave(paste0(outdir, 'bal/nonce/all.png'), pl,
         width = 9, height = 8)
  
  unique_langs <- unique(all_langs$language)
  for (i in seq(1,length(unique_langs), 5)) {
    print(i)
    langs <- unique_langs[i:(i+4)]
    
    all_langs.sub <- subset(all_langs, language %in% langs)
    pl <- ggplot(data = all_langs.sub, aes(color = lex.type)) +
       geom_text(data = subset(all_langs.sub, lex.type == 'og'),
                x = Inf, y = Inf, 
                hjust   = 'right',
                vjust   = 'top',
                size = 6,
                color = 'black',
                aes(label = paste0('p < ', substr(og.greater.r, 1, 5)))) +
                                   #'\nR2 < ', substr(og.greater.R2, 1, 5)))) +
        xlim(c(-2,2)) + xlab('log Mean Frequency\n(z-scored)') +
        ylim(c(-2,2)) + ylab('log Type Count\n(z-scored)') +
        geom_abline(data = subset(all_langs.sub, lex.type == 'nonce'), 
                size = .5, alpha = .3, linetype = 'dashed',
                aes(slope = b1, intercept = b0, color = lex.type)) +
        geom_abline(data = subset(all_langs.sub, lex.type == 'og'),
                size = 1.5, alpha = 1,
                aes(slope = b1, intercept = b0, color = lex.type)) +
        facet_wrap(~language) +
        scale_color_manual(name = 'Lexicon Type',
                         labels = c('Nonce Forms', 'Real-world'),
                         values = c('#00BFC4', '#F8766D')) +
        theme(text = element_text(size=20), 
            legend.position = c(1, .15), 
            legend.justification = c(1, 0))
    
    ggsave(paste0(outdir, 'bal/nonce/', i, '.png'), pl,
           width = 9, height = 7)
  }
  
  
  all_langs <- add_language_family(all_langs)
  all_langs$lex.type <- factor(all_langs$lex.type, levels = c('nonce', 'og'))
  m <- glmer(data = all_langs,
            lex.type ~ r + (1|family/language), 
            family = binomial, verbose = 2)
  lmer_str <- lmer_to_latex(m)
  lmer_str <- gsub(' r ', ' $r$ ', lmer_str)
  write(x = lmer_str, file = paste0(outdir, 'tex_tables/bal/nonce.tex'))
  1
}

tst_first_later_ent <- function(load_files = F, p_threshold = .95) {
  if (!load_files) {
    lapply(lex.files, function(f) {
      d <- lf(paste0(lexdir, f), do_up = F)
      d <- subset(d, word.length > 1)
      f <- sub('.lex.txt', '', f)
      print(f)
      
      nb_to_remove <- floor(nrow(d) / 10)
      bootstrap <- pblapply(1:nb_scramble, function(i) {
        #  # phone by position
        pbp <- phones_by_position(d[-sample(1:nrow(d), nb_to_remove, replace = F), ], N = 1)
        fvl <- first_vs_later(pbp, p_threshold = p_threshold, pos = 1)
        H.1 <- sum(fvl$p[fvl$position == 'first'] * -log2(fvl$p[fvl$position == 'first']))
        H.2 <- sum(fvl$p[fvl$position != 'first'] * -log2(fvl$p[fvl$position != 'first']))
        rbind(data.frame(i = i, position = 'first', H = H.1),
              data.frame(i = i, position = 'later', H = H.2))
      }, cl = nb_cores)
      bootstrap <- ldply(bootstrap, data.frame)
      print(paste0(mean(bootstrap$H[bootstrap$position == 'first']), 
                   ' -- ', mean(bootstrap$H[bootstrap$position != 'first'])))
      bootstrap$language <- f
      write_tsv(bootstrap, paste0(outdir, 'first_v_later/', f, '.novel.txt'))
      1
    })
  }
  all_langs <- lapply(lex.files, function(f) {
    f <- sub('.lex.txt', '', f)
    d <- read_tsv(paste0(outdir, 'first_v_later/', f, '.novel.txt'), col_types = cols())
    d$H.z <- scale(d$H)
    H.1 <- d$H[d$position == 'first']
    H.2 <- d$H[d$position != 'first']
    t.t <- t.test(H.1, H.2, paired = T, alternative = 'greater')
    d$t.dist <- t.t$estimate
    d$p.val <- t.t$p.value
    d
  })
  all_langs <- ldply(all_langs, data.frame)
  all_langs$t.str <- substr(sub('0.', '.', all_langs$t.dist), 1, 4)
  all_langs$p.str <- sapply(all_langs$p.val, pval_str)
  
  pl <- ggplot(data = all_langs, aes(H, fill = position)) + 
    geom_histogram(data = subset(all_langs, position == 'first'), 
                   bins = max(30, nb_scramble/20), alpha = .7, color = 'gray') +
    geom_histogram(data = subset(all_langs, position != 'first'),
                   bins = max(30, nb_scramble/20), alpha = .7, color = 'gray') +
    xlab('Entropy (H)') + ylab('Count') +
    scale_x_continuous(breaks = c(2,2.5,3,3.5,4,4.5,5), labels = c(2,2.5,3,3.5,4,4.5,5)) +
    scale_fill_manual(name = 'Contrast Position',
                      labels = c('Word-initial', 'Non-initial'),
                      values = c('#00BFC4', '#F8766D')) +
    theme(text = element_text(size = 20),
          legend.position = 'bottom') +
    facet_wrap(~language)
  ggsave(paste0(outdir, 'first_v_later/pngs/all.png'), pl,
         width = 9, height = 8)
  
  unique_langs <- unique(all_langs$language)
  for (i in seq(1,length(unique_langs), 5)) {
    print(i)
    langs <- unique_langs[i:(i+4)]
    
    all_langs.sub <- subset(all_langs, language %in% langs)
    pl <- ggplot(data = all_langs.sub, aes(H, fill = position)) + 
      geom_histogram(data = subset(all_langs.sub, position == 'first'), 
                     bins = max(30, nb_scramble/20), alpha = .5, color = 'gray') +
      geom_histogram(data = subset(all_langs.sub, position != 'first'),
                     bins = max(30, nb_scramble/20), alpha = .5, color = 'gray') +
      xlab('Entropy (H)') + ylab('Count') +
      scale_x_continuous(breaks = c(2,2.5,3,3.5,4,4.5,5), labels = c(2,2.5,3,3.5,4,4.5,5)) +
      scale_fill_manual(name = 'Contrast Position',
                        labels = c('Word-initial', 'Non-initial'),
                        values = c('#00BFC4', '#F8766D')) +
      theme(text = element_text(size = 20)) +
      facet_wrap(~language)
    
    xadj <- mean(layer_scales(pl)$x$range$range)
    yadj <- layer_scales(pl)$y$range$range[1] + (.65 * (layer_scales(pl)$y$range$range[2] - layer_scales(pl)$y$range$range[1]))
    pl <- pl + geom_text(data = subset(all_langs.sub, i == 1 & position == 'first'),
                         x = xadj, y = yadj, size = 7, fontface = 'bold',
                         aes(label = paste0(t.str, '\n', p.str)))
    
    ggsave(plot = pl, filename = paste0(outdir, 'first_v_later/pngs/', i, '.png'),
           width = 10, height = 7)
  }
  
  all_langs <- add_language_family(all_langs)
  all_langs$position <- factor(all_langs$position, levels = c('later', 'first'))
  m <- glmer(data = all_langs,
             position ~ H + (1|family/language), family = binomial)
  lmer_str <- lmer_to_latex(m)
  write(x = lmer_str, file = paste0(outdir, 'tex_tables/first_v_later/lme.tex'))
  
  1
  all_langs$diff <- all_langs$H[match(paste(all_langs$language,
                                            all_langs$i, 
                                            'first'), 
                                      paste(all_langs$language,
                                            all_langs$i,
                                            all_langs$position))] - 
    all_langs$H[match(paste(all_langs$language,
                            all_langs$i, 
                            'later'), 
                      paste(all_langs$language,
                            all_langs$i,
                            all_langs$position))]
  
  all_langs$not.sig.in.prev <- F
  all_langs$not.sig.in.prev[all_langs$language %in% c('Armenian', 'Bengali', 'Cantonese',
                                                      'Georgian', 'Hausa', 'Hungarian',
                                                      'Kaqchikel', 'Korean', 'Swahili',
                                                      'Turkish', 'Vietnamese')] <- T
  
  all_langs$not.sig.in.prev <- factor(all_langs$not.sig.in.prev, 
                                      levels = c(F, T))
  
  m <- glmer(data = subset(all_langs, position == 'later'),
             not.sig.in.prev ~ diff +(1|family/language),
             family = binomial)
  summary(m)
  lmer_str <- lmer_to_latex(m)
  write(x = lmer_str, file = paste0(outdir, 'tex_tables/first_v_later/diff-lme.tex'))
  1
}

# cohortp
tst_cohort_entropy_sig <- function(load_files = F) {
  if (!load_files) {
    print('prepping cohort entropy data....')
    all_chrts <- pblapply(lex.files, function(f) {
      d <- lf(paste0(lexdir, f), do_up = T)
      f <- sub('.lex.txt', '', f)
      
      # remove all cohorts post UP here, will do the pruning in the LM function
      chrt <- entropy_by_cohort(d, include_post_up = F)
      chrt$language <- f
      chrt
    }, cl = nb_cores)
    all_chrts <- ldply(all_chrts, data.frame) 
    write_tsv(all_chrts, paste0(outdir, 'cohort_p/all_chrts.txt'))
    1
  }
  all_chrts <- read_tsv(paste0(outdir, 'cohort_p/all_chrts.txt'), col_types = cols())
  
  all_ogs <- pblapply(unique(all_chrts$language), function(f) {
    chrt <- subset(all_chrts, language == f)
    
    m <- cohort_H_m(chrt)
    sm <- summary(m)
    
    rbind(data.frame(language = f, b = sm$coefficients[2,1], t = sm$coefficients[2,3], factor = 'log.prob'),
          data.frame(language = f, b = sm$coefficients[3,1], t = sm$coefficients[3,3], factor = 'cohort.count'),
          data.frame(language = f, b = sm$coefficients[4,1], t = sm$coefficients[4,3], factor = 'position'))
          
          
  }, cl = nb_cores)
  all_ogs <- ldply(all_ogs, data.frame)
  all_ogs$greater.12 <- abs(all_ogs$t) > 12
  all_ogs$t <- pmin(12, pmax(-12, all_ogs$t))
  all_ogs$t[all_ogs$factor == 'log.prob'] <- 
    pmin(11, pmax(-11, all_ogs$t[all_ogs$factor == 'log.prob']))
  
  pl <- ggplot(data = all_ogs, aes(x = language, y = t, color = factor, shape = greater.12)) +
    xlab('Language') + ylab('t value') +
    scale_color_discrete(name = 'Model Factor',
                         labels = c('log Cohort Prob', 'Cohort Size (type)', 'Position')) +
    scale_shape_discrete(name = 't value',
                         labels = c('|t| <= 12', '|t| > 12')) +
    theme(text = element_text(size = 18),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_jitter(size = 4, width = 0) +
    geom_hline(yintercept = 2, linetype = 'dashed') +
    geom_hline(yintercept = -2, linetype = 'dashed')
  ggsave(plot = pl, filename = paste0(outdir, 'cohort_p/sig/all.png'),
         width = 9, height = 8)
  
  
  unique_langs <- unique(all_ogs$language)
  for (i in seq(1, length(unique_langs), 8)[1:3]) {
    print(i)
    all_ogs.sub <- subset(all_ogs, language %in% unique_langs[i:(i+7+(i==17))])
    pl <- ggplot(data = all_ogs.sub, aes(x = language, y = t, color = factor, shape = greater.12)) +
      xlab('Language') + ylab('t value') +
      scale_color_discrete(name = 'Model Factor',
                           labels = c('log Cohort Prob', 'Cohort Size (type)', 'Position')) +
      scale_shape_discrete(name = 't value',
                           labels = c('|t| <= 12', '|t| > 12')) +
      theme(text = element_text(size = 18),
            axis.text.x = element_text(angle = 90, hjust = 1)) +
      geom_jitter(size = 4, width = 0) +
      geom_hline(yintercept = 2, linetype = 'dashed') +
      geom_hline(yintercept = -2, linetype = 'dashed')
    
    ggsave(plot = pl, filename = paste0(outdir, 'cohort_p/sig/',i,'.png'),
           width = 9, height = 7)
  }
  
  # LME's
  all_chrts <- add_language_family(all_chrts)
  
  print('cohort p - LME')
  m <- lmer(data = all_chrts,
            z.(entropy) ~ z.(log.prob) + z.(cohort.count) + position +
              (1|family/language), 
            REML = F)
  lmer_str <- lmer_to_latex(m)
  lmer_str <- gsub('z..log.prob.', 'log Cohort Prob.', lmer_str)
  lmer_str <- gsub('z..cohort.count.', 'Cohort Size (type)', lmer_str)
  lmer_str <- gsub('position', 'Cohort Position', lmer_str)
  write(lmer_str, file = paste0(outdir, 'tex_tables/cohort_p/lme.tex'))
  
  1
}

tst_cohort_entropy_nonce <- function(load_files = F) {
  if (!load_files) {
    all_chrts <- read_tsv(paste0(outdir, 'cohort_p/all_chrts.txt'), col_types = cols())
    all_langs <- lapply(unique(all_chrts$language), function(f) {
      print(paste0(f, ' - cohort ent by prob'))
      
      d <- lf(paste0(lexdir, f, '.lex.txt'), do_up = T)
      dn <- load_nonce(paste0(noncedir, f, '.new.txt'))
      chrt <- subset(all_chrts, language == f)
      
      m <- cohort_H_m(chrt)
      sm <- summary(m)
      print(sm)
      print(vif(m))
    
      print(paste0('nonce lexicons for: ', f))
      nonce_bs <- pblapply(1:nb_scramble, function(i) {
        d.nonce <- nonce_lexicon(d, dn, keep_length = T, do_up = T)
        chrt.nonce <- entropy_by_cohort(d.nonce, include_post_up = F)
      
        m.nonce <- cohort_H_m(chrt.nonce)
        sm.nonce <- summary(m.nonce) 
      
        rbind(data.frame(b = sm.nonce$coefficients[2], t = sm.nonce$coefficients[2,3],
                       p = sm.nonce$coefficients[2,4], lex.type = 'nonce'))
      }, cl = nb_cores)
      df <- rbind(data.frame(b = sm$coefficients[2], t = sm$coefficients[2,3],
                           p = sm$coefficients[2,4], lex.type = 'og'),
                ldply(nonce_bs, data.frame))
      print(paste0(df$b[1], ' -- ', mean(df$b), ' : ',
                 sum(df$b > df$b[1]) / nb_scramble))
    
      df$language <- f
      df$b.z <- 0
      df$b.z <- scale(df$b)
      df$b.less <- sum(df$b > df$b[1]) / nb_scramble
      write_tsv(df, paste0(outdir, 'cohort_p/', f, '.txt'))
      df
    })
  }
  all_langs <- lapply(lex.files, function(f) {
    f <- sub('.lex.txt', '', f)
    df <- read_tsv(file = paste0(outdir, 'cohort_p/', f, '.txt'), col_types = cols())
    df
  })
  all_langs <- ldply(all_langs, data.frame)
  all_langs$lex.type <- factor(all_langs$lex.type, levels = c('og', 'nonce'))
  
  pl <- ggplot(data = all_langs, aes(b.z, fill = lex.type)) +
    xlab('Model Estimate: log Cohort Prob.') + ylab('Count') +
    scale_x_continuous(breaks = c(-4,-2,0,2,4), labels = c(-4,2,0,2,4)) +
    geom_histogram(alpha = .3, bins = max(30, nb_scramble / 100)) +
    geom_vline(xintercept = -2, size = 1, linetype = 'dotted') +
    geom_vline(xintercept = 2, size = 1, linetype = 'dotted') +
    geom_vline(data = subset(all_langs, lex.type == 'og'), 
               aes(xintercept = b.z), size = 1.5, color = '#F8766D') +
    scale_fill_manual(name = 'Lexicon Type',
                      labels = c('Real-world', 'Nonce Forms'),
                      values = c('#F8766D', '#00BFC4')) +
    facet_wrap(~language) +
    theme(text = element_text(size=20), 
          legend.position = 'bottom')
  ggsave(filename = paste0(outdir, 'cohort_p/nonce/all.png'), pl,
         width = 9, height = 8)
  
  print('cohort p (nonce lexicons)')
  unique_langs <- unique(all_langs$language)
  for (i in seq(1,length(unique_langs), 5)) {
    print(i)
    langs <- unique_langs[i:(i+4)]
    all_langs.sub <- subset(all_langs, language %in% langs)
    
    pl <- ggplot(data = all_langs.sub, aes(b.z, fill = lex.type)) +
      xlab('Model Estimate: log Cohort Prob.') + ylab('Count') +
      scale_x_continuous(breaks = c(-4,-2,0,2,4), labels = c(-4,2,0,2,4)) +
      geom_histogram(alpha = .3, bins = max(30, nb_scramble / 100)) +
      geom_vline(xintercept = -2, size = 1, linetype = 'dotted') +
      geom_vline(xintercept = 2, size = 1, linetype = 'dotted') +
      geom_vline(data = subset(all_langs.sub, lex.type == 'og'), 
                 aes(xintercept = b.z), size = 1.5, color = '#F8766D') +
      geom_text(mapping = aes(x = Inf, y = Inf, 
                              label = paste0('p < ', substr(pmax(b.less, .001), 1, 5))),
                hjust   = 'right',
                vjust   = 'top',
                size = 5) +
      scale_fill_manual(name = 'Lexicon Type',
                        labels = c('Real-world', 'Nonce Forms'),
                        values = c('#F8766D', '#00BFC4')) +
      facet_wrap(~language) +
      theme(text = element_text(size=20), 
            legend.position = c(.95, .15), 
            legend.justification = c(1, 0))
    
    
    ggsave(filename = paste0(outdir, 'cohort_p/nonce/',i ,'.png'), pl,
           width = 9, height = 7)
  }
  
  all_langs$lex.type <- factor(all_langs$lex.type, levels = c('nonce', 'og'))
  all_langs <- add_language_family(all_langs)
  
  m <- glmer(data = all_langs,
             lex.type ~ b + (1|family/language), family = binomial)
  lmer_str <- lmer_to_latex(m)
  lmer_str <- gsub(' b ', ' $b$ ', lmer_str)
  write(lmer_str, paste0(outdir, 'tex_tables/cohort_p/nonce.tex'))
  
  1
}

# chapter 5
#done
tst_synergy_same_nbrhood <- function(load_files = F) {
  if (!load_files) {
    # load neighborhood size from part 1
    mps <- lapply(lex.files, function(f) {
      f <- sub('.lex.txt', '', f)
      d <- read_tsv(paste0(outdir, 'rdn1/', f, '.txt'), col_types = cols())
      d
    })
    mps <- ldply(mps, data.frame) 
    mps <- subset(mps, select = c(word, language, min.pairs))
    
    all_langs <- lapply(lex.files, function(f) {
      d <- lf(paste0(lexdir, f), do_up = F)
      f <- sub('.lex.txt', '', f)
      d$language <- f
      print(f)
      
      d.mp <- subset(mps, language == f)
      d$min.pairs <- d.mp$min.pairs[match(d$word, d.mp$word)]
      d <- subset(d, min.pairs >= 0)
      
      print(nrow(d))
      lex <- lex2si(d, shuffle = F, shuffle_same_nbrhood = F)
      # remove words of length 1
      lex <- subset(lex, word.length >= 2)
      
      df <- data.frame(language = f, i = 1, lex.type = 'original',
                       seq.H = lex_ent(lex))
      
      df.shuffled <- pblapply(1:nb_scramble + 1, function(i) {
        lex.scr <- lex2si(d, shuffle = T, shuffle_same_nbrhood = T)
        lex.scr <- subset(lex.scr, word.length >= 2)
        
        data.frame(language = f, i = i, lex.type = 'shuffled',
                   seq.H = lex_ent(lex.scr))
      }, cl = nb_cores)
      df.shuffled <- ldply(df.shuffled, data.frame)
      df <- rbind(df, df.shuffled)
      
      df$seq.H.scaled <- scale(df$seq.H)
      
      print(sum(df$seq.H > df$seq.H[1]) / (length(df$seq.H) - 1))
      print(paste0(df$seq.H[1], '  -  ', mean(df$seq.H[2:nrow(df)])))
      
      write_tsv(x = df, path = paste0(outdir, 'total_const_redun/nbrhood/', f, '.txt'))
      1
    })
  }
  all_langs <- lapply(lex.files, function(f) {
    f <- sub('.lex.txt', '', f)
    df <- read_tsv(file = paste0(outdir, 'total_const_redun/nbrhood/', f, '.txt'), col_types = cols())
    df$diff.from.mean <- df$seq.H - mean(df$seq.H)
    df$seq.greater <- sum(df$seq.H > df$seq.H[1]) / (nrow(df) - 1)
    df
  })
  all_langs <- ldply(all_langs, data.frame)
  
  print('synergy - nbrs')
  
  pl <- ggplot(data = all_langs, aes(seq.H.scaled, fill = lex.type)) +
    geom_histogram(alpha = .3, bins = max(30, nb_scramble / 10)) +
    scale_fill_manual(name = 'Lexicon Type',
                      labels = c('Real-world', 'Probability-shuffled'),
                      values = c('#F8766D', '#00BFC4')) +
    xlab('Total Lexical Entropy, Neighborhood Size Held Contstant\n(z-scored)') + ylab('Count') +
    scale_x_continuous(breaks = c(-4,-2,0,2,4), labels = c(-4,-2,0,2,4)) +
    geom_vline(data = subset(all_langs, lex.type == 'original'), 
               aes(xintercept = seq.H.scaled), color = '#F8766D', size = 1.5) +
    geom_text(mapping = aes(x = Inf, y = Inf, 
                            label = paste0('  p < ', substr(pmax(seq.greater, .001), 1, 5))),
              hjust   = 'right',
              vjust   = 'top',
              size = 6) +
    theme(text = element_text(size=22), legend.position = 'none') +
    facet_wrap(~language)
  ggsave(filename = paste0(outdir, 'total_const_redun/nbrhood/pngs/all.png'), pl,
         width = 10, height = 9)
  
  li <- c(0, 9, 18, 25)
  for (i in 1:(length(li)-1)) {
    print(i)
    sub_langs <- unique(all_langs$language)[(li[i]+1):li[i+1]]
    all_langs_sub <- subset(all_langs, language %in% sub_langs)
    
    pl <- ggplot(data = all_langs_sub, aes(seq.H.scaled, fill = lex.type)) +
      geom_histogram(alpha = .3, bins = max(30, nb_scramble / 10)) +
      scale_fill_manual(name = 'Lexicon Type',
                        labels = c('Real-world', 'Probability-shuffled'),
                        values = c('#F8766D', '#00BFC4')) +
      xlab('Total Lexical Entropy, Neighborhood Size Held Contstant\n(z-scored)') + ylab('Count') +
      scale_x_continuous(breaks = c(-4,-2,0,2,4), labels = c(-4,2,0,2,4)) +
      geom_vline(data = subset(all_langs_sub, lex.type == 'original'), 
                 aes(xintercept = seq.H.scaled), color = '#F8766D', size = 1.5) +
      geom_text(mapping = aes(x = -Inf, y = Inf, 
                              label = paste0('  p < ', substr(pmax(seq.greater, .001), 1, 5))),
                hjust   = 'left',
                vjust   = 'top',
                size = 6) +
      theme(text = element_text(size=22)) +
      facet_wrap(~language)
    
    ggsave(filename = paste0(outdir, 'total_const_redun/nbrhood/pngs/',i,'.png'), pl,
           width = 12, height = 10)
    1
  }
  all_langs <- add_language_family(all_langs)
  all_langs$lex.type <- factor(all_langs$lex.type, levels = c('shuffled', 'original'))
  
  m <- glmer(data = all_langs,
             lex.type ~ z.(diff.from.mean) + (1|family/language), family = binomial)
  summary(m)
  lmer_str <- lmer_to_latex(m)
  lmer_str <- sub('z..diff.from.mean.', '$H^*_L$; difference from mean', lmer_str)
  write(lmer_str, paste0(outdir, 'tex_tables/synergy/neighbors.tex'))
  1
}

tst_synergy_same_tail <- function(load_files = F) {
  if (!load_files) {
    all_langs <- lapply(lex.files, function(f) {
      d <- lf(paste0(lexdir, f), do_up = T)
      f <- sub('.lex.txt', '', f)
      d$language <- f
      
      print(f)
      print(nrow(d))
    
      lex <- lex2si(d, shuffle = F, shuffle_same_tail = F)

      df <- data.frame(language = f, i = 1, lex.type = 'original',
                   seq.H = lex_ent(lex))
    
      df.shuffled <- pblapply(1:nb_scramble + 1, function(i) {
        lex.scr <- lex2si(d, shuffle = T, shuffle_same_tail = T)

        data.frame(language = f, i = i, lex.type = 'shuffled',
                   seq.H = lex_ent(lex.scr))
      }, cl = nb_cores)
      df.shuffled <- ldply(df.shuffled, data.frame)
      df <- rbind(df, df.shuffled)
    
      df$seq.H.scaled <- scale(df$seq.H)
  
      print(sum(df$seq.H > df$seq.H[1]) / (length(df$seq.H) - 1))
      print(paste0(df$seq.H[1], '  -  ', mean(df$seq.H[2:nrow(df)])))
    
      write_tsv(x = df, path = paste0(outdir, 'total_const_redun/tail/', f, '.txt'))
      1
    })
  }
  all_langs <- lapply(lex.files, function(f) {
    f <- sub('.lex.txt', '', f)
    df <- read_tsv(file = paste0(outdir, 'total_const_redun/tail/', f, '.txt'), col_types = cols())
    df$diff.from.mean <- df$seq.H - mean(df$seq.H)
    df$seq.greater <- sum(df$seq.H > df$seq.H[1]) / (nrow(df) - 1)
    #df <- df[1:100,]
    df
  })
  all_langs <- ldply(all_langs, data.frame)
  
  print('synergy - tail')
  
  pl <- ggplot(data = all_langs, aes(seq.H.scaled, fill = lex.type)) +
    geom_histogram(alpha = .3, bins = max(30, nb_scramble / 10)) +
    scale_fill_manual(name = 'Lexicon Type',
                      labels = c('Real-world', 'Probability-shuffled'),
                      values = c('#F8766D', '#00BFC4')) +
    xlab('Total Lexical Entropy, Redundant Tails Held Constant\n(z-scored)') + ylab('Count') +
    scale_x_continuous(breaks = c(-4,-2,0,2,4), labels = c(-4,-2,0,2,4)) +
    geom_vline(data = subset(all_langs, lex.type == 'original'), 
               aes(xintercept = seq.H.scaled), color = '#F8766D', size = 1.5) +
    geom_text(mapping = aes(x = Inf, y = Inf, 
                            label = paste0('  p < ', substr(pmax(seq.greater, .001), 1, 5))),
              hjust   = 'right',
              vjust   = 'top',
              size = 6) +
    theme(text = element_text(size=22), legend.position = 'none') +
    facet_wrap(~language)
  ggsave(filename = paste0(outdir, 'total_const_redun/tail/pngs/all.png'), pl,
         width = 10, height = 9)
  
  li <- c(0, 9, 18, 25)
  for (i in 1:(length(li)-1)) {
    print(i)
    sub_langs <- unique(all_langs$language)[(li[i]+1):li[i+1]]
    all_langs_sub <- subset(all_langs, language %in% sub_langs)
    
    pl <- ggplot(data = all_langs_sub, aes(seq.H.scaled, fill = lex.type)) +
      scale_fill_manual(name = 'Lexicon Type',
                         labels = c('Real-world', 'Probability-shuffled'),
                         values = c('#F8766D', '#00BFC4')) +
      xlab('Total Lexical Entropy, Redundant Tails Held Constant\n(z-scored)') + ylab('Count') +
      scale_x_continuous(breaks = c(-4,-2,0,2,4), labels = c(-4,2,0,2,4)) +
      geom_histogram(alpha = .3, bins = max(30, nb_scramble / 10)) +
      geom_vline(data = subset(all_langs_sub, lex.type == 'original'), 
                 aes(xintercept = seq.H.scaled), color = '#F8766D', size = 1.5) +
      geom_text(mapping = aes(x = -Inf, y = Inf, 
                              label = paste0('  p < ', substr(pmax(seq.greater, .001), 1, 5))),
                hjust   = 'left',
                vjust   = 'top',
                size = 6) +
      theme(text = element_text(size=22)) +
      facet_wrap(~language)
    
    ggsave(filename = paste0(outdir, 'total_const_redun/tail/pngs/',i,'.png'), pl,
           width = 12, height = 10)
    1
  }
  
  all_langs <- add_language_family(all_langs)
  all_langs$lex.type <- factor(all_langs$lex.type, levels = c('shuffled', 'original'))
  
  m <- glmer(data = all_langs,
            lex.type ~ z.(diff.from.mean) + (1|family/language), family = binomial)
  summary(m)
  lmer_str <- lmer_to_latex(m)
  lmer_str <- sub('z..diff.from.mean.', '$H^*_L$; difference from mean', lmer_str)
  write(lmer_str, paste0(outdir, 'tex_tables/synergy/tail.tex'))
  1
}

# intro stuff
intro_graphs <- function(do_all_langs = T) {
  
  # p(b|b) in a noisy channel
  print('p(b|b)')
  df <- data.frame(log.b = seq(-7, -.8, .1))
  df$b <- 10 ** df$log.b
  df$a <- 1- df$b
  df$ratio <- df$b / df$a
  
  df.9 <- df
  df.9$bb <- .9 * df.9$b
  df.9$ab <- .1 * df.9$a
  df.9$b.given.b <- '.9'
  
  df.99 <- df
  df.99$bb <- .99 * df.99$b
  df.99$ab <- .01 * df.99$a
  df.99$b.given.b <- '.99'
  
  df <- rbind(df.9, df.99)
  df$given.ratio <- df$bb / df$ab
  
  pl <- ggplot(df, aes(x = log.b, y = given.ratio, color = b.given.b)) + 
    theme(text = element_text(size=16)) +
    geom_line() + 
    geom_hline(linetype = 'dashed', yintercept = 1) +
    scale_y_continuous(breaks = seq(0, 1.3, .2), limits = c(0, 1.3)) +
    scale_color_discrete(name = 'p(b|b)') +
    xlab('log p(b)') + ylab('p(b|b) / p(b|a)')
  ggsave(pl, filename = paste0(outdir, 'intro/p-ab.png'), width = 8, height = 8)
  
  
  # Armenian stops
  d <- lf(paste0(lexdir, 'Armenian.lex.txt'))
  pbp <- phones_by_position(d)
  fvl <- first_vs_later(pbp, p_threshold = 1)
  fvl.stops <- subset(fvl, phone %in% c('p', 'b', 'P', 't', 'd', 'T', 'k', 'g', 'K'))
  
  fvl.stops$place <- ''
  fvl.stops$place[fvl.stops$phone %in% c('p', 'b', 'P')] <- 'Labial'
  fvl.stops$place[fvl.stops$phone %in% c('t', 'd', 'T')] <- 'Coronal'
  fvl.stops$place[fvl.stops$phone %in% c('k', 'g', 'K')] <- 'Velar'
  fvl.stops$place <- as.factor(fvl.stops$place)
  fvl.stops$place <- factor(fvl.stops$place, levels = c('Labial', 'Coronal', 'Velar'))
  
  fvl.stops$voice <- ''
  fvl.stops$voice[fvl.stops$phone %in% c('P', 'T', 'K')] <- 'Voiceless Aspirated'
  fvl.stops$voice[fvl.stops$phone %in% c('p', 't', 'k')] <- 'Voiceless'
  fvl.stops$voice[fvl.stops$phone %in% c('b', 'd', 'g')] <- 'Voiced'
  fvl.stops$voice <- as.factor(fvl.stops$voice)
  fvl.stops$voice <- factor(fvl.stops$voice, levels = c('Voiceless Aspirated', 'Voiceless', 'Voiced'))
  
  H.1 <- sum(fvl.stops$p[fvl.stops$position == 'first'] * -log2(fvl.stops$p[fvl.stops$position == 'first']))
  H.2 <- sum(fvl.stops$p[fvl.stops$position != 'first'] * -log2(fvl.stops$p[fvl.stops$position != 'first']))
  
  pl <- ggplot(fvl.stops, aes(fill = voice, density = place, x = position, y = p)) + 
    theme(text = element_text(size=20)) +
    scale_fill_discrete(name = 'Voicing') +
    scale_x_discrete(name = 'Position', labels = c('Word-initial', 'Later')) + 
    ylab('Total Segment Probability') +
    geom_col(position = 'dodge') +     
    geom_text(aes(label=phone), 
              position=position_dodge(width=0.9), vjust=-0.25) +
    geom_text(data = fvl.stops[1,], label = paste0('H = ', substr(H.1, 1, 4)), x = 1, y = .01, size = 6) +
    geom_text(data = fvl.stops[1,], label = paste0('H = ', substr(H.2, 1, 4)), x = 2, y = .01, size = 6)
  
  ggsave(plot = pl, filename = paste0(outdir, 'intro/armenian-stops.png'),
         width = 11, height = 8)
  
  # Zipf's LoA for all langs
  print('Zipf LoA')
  if (T) {
    all_langs <- pblapply(lex.files, function(f) {
      d <- lf(paste0(lexdir, f))
      f <- sub('.lex.txt', '', f)
      d <- subset(d, select = c(freq, unig.h, word.length))
      d$Language <- f
      d$wl.z <- scale(d$word.length)
      d$u.z <- scale(d$unig.h)
      d
    })
    all_langs <- ldply(all_langs, data.frame)
  
    pl <- ggplot(data = all_langs, aes(x = log(freq), y = wl.z)) + 
      theme(text = element_text(size=16)) +
      geom_smooth(method = 'lm') +
      xlab('log Word Frequency') + ylab('Word length (segments)\n(z-scored)') +
      facet_wrap(~Language)
    rm(all_langs)
  } else {
    d <- lf(paste0(lexdir, 'English.lex.txt'))
    pl <- ggplot(data = d, aes(x = log(freq), y = word.length, label = word)) +
      xlab('log Word Frequency') + ylab('Word length (segments)') +
      scale_y_continuous(breaks = 1:10) +
      theme(text = element_text(size=20)) +
      geom_smooth(method = 'lm', alpha = .75, linetype = 'dashed', se = 0) +
      geom_text(size = 8, data = sample_n(subset(d, log(freq) < 8 & log(freq) > 1), 22), check_overlap = T) +
      geom_text(size = 8, data = sample_n(subset(d, log(freq) > 8), 4), check_overlap = T)
  }
  ggsave(pl, filename = paste0(outdir, 'intro/zipf-loa-all.png'), width = 10, height = 8)
  
  if (do_all_langs) {
    mclapply(lex.files, function(f) {
      d <- lf(paste0(lexdir, f))
      f <- sub('.lex.txt', '', f)
    
      print(f)
      # Zipf's freq
      pl1 <- ggplot(data = subset(d, freq.rank <= 4000), aes(x = freq.rank, y = freq)) + 
        theme(legend.position = "none", text = element_text(size=16)) +
        geom_point() +
        xlab('Frequency Rank') + ylab('Frequency per million')
    
      # Zipf's Law of Abbr.
      pl2 <-ggplot(data = d, aes(x = log(freq), y = word.length)) +
        geom_point(alpha = .05) + geom_smooth(method = 'lm', size = 1.5) +
        theme(legend.position = "none", text = element_text(size=16)) +
        xlab('Log Frequency') + ylab('Word Length')
    
      # type freq by segment
      lex <- lex2si(d)
      lex <- lex[!duplicated(paste(lex$word, lex$phone)),]
      if (grepl('Cantonese|Vietnamese', f)) {
        lex <- subset(lex, !grepl('[0-9]', lex$phone))
      }
      lex.c <- aggregate(freq ~ phone, lex, length)
      lex.f <- aggregate(freq ~ phone, lex, sum)
      lex.a <- join(lex.f, lex.c, by = 'phone')
      names(lex.a)[3] <- 'count'
    
      pl3 <- ggplot(data = lex.a, aes(y = log(freq), x = count, label = phone)) +
        theme(legend.position = "none", text = element_text(size=16)) +
        geom_text_repel(segment.alpha = .25, 
                      arrow = arrow(length = unit(0.02, "npc"))) +
        ylab('Log Frequency') + xlab('Num. Word Types')
        
      # word length
      pl4 <- ggplot(data = d, aes(word.length)) + 
        geom_histogram(binwidth = 1) +
        theme(legend.position = "none", text = element_text(size=16)) +
        scale_x_continuous(breaks = 1:max(d$word.length)) +
        xlab('Word Length') + ylab('Count')
        
      pl <- grid.arrange(pl1, pl2, pl3, pl4, 
                       nrow = 2, 
                       top = textGrob(paste0(f, ' (', nrow(d), ' lemmas)'), gp=gpar(fontsize=16,font=3)))
      ggsave(pl, filename = paste0(outdir, 'intro/', f, '.png'), width = 12, height = 7)  
    }, mc.cores = nb_cores)
  }
  
  # shuffled, length not held constant
  d <- lf(paste0(lexdir, 'English.lex.txt'))
  ttl_length <- data.frame(i = 1, lex.type = 'og', total_len =sum(d$freq * d$word.length))
  ttl_length_scr <- pblapply(2:1001, function(i) {
    d.scr <- shuffle_frequencies(d, keep_length = F)
    data.frame(i = i, lex.type = 'scr', total_len = sum(d.scr$freq * d.scr$word.length))
  })
  ttl_length <- rbind(ttl_length, ldply(ttl_length_scr, data.frame))
  
  pl <- ggplot(data = ttl_length, aes(total_len, fill = lex.type)) +
    xlab('Total Corpus Length (segments)') + ylab('Count') +
    geom_histogram(bins = floor(nrow(ttl_length) / 30)) + 
    geom_vline(data = subset(ttl_length, lex.type == 'og'), 
                          aes(xintercept = total_len), size = 1.5, color = '#F8766D') +
    theme(legend.position = 'none', text = element_text(size = 20))
  ggsave(plot = pl, filename = paste0(outdir, 'intro/shuffled-nolength.png'), 
         width = 10, height = 8)
  
  pl <- ggplot(data = all_langs.sub, aes(b.z, fill = lex.type)) +
    xlab('Model Estimate: log Cohort Prob.') + ylab('Count') +
    geom_histogram(alpha = .3, bins = max(30, nb_scramble / 100)) +
    geom_vline(xintercept = -2, size = 1, linetype = 'dotted') +
    geom_vline(xintercept = 2, size = 1, linetype = 'dotted') +
    geom_vline(data = subset(all_langs.sub, lex.type == 'og'), 
               aes(xintercept = b.z), size = 1.5, color = '#F8766D') +
    geom_text(mapping = aes(x = Inf, y = Inf, 
                            label = paste0('p < ', substr(pmax(b.less, .001), 1, 5))),
              hjust   = 'right',
              vjust   = 'top',
              size = 5) +
    scale_fill_manual(name = 'Lexicon Type',
                      labels = c('Real-world', 'Nonce Forms'),
                      values = c('#F8766D', '#00BFC4')) +
}

redun_methods_graphs <- function() {
  d <- read_tsv(paste0(outdir, 'rdn1/', 'English', '.txt'), col_types = cols())
  d <- subset(d, word.length >= 2 & word.length <= 8)
  
  d$scaled.min.pairs <- 0
  d$scaled.mean.dist.nonmp <- 0
  for (i in 2:8) {
    d$scaled.min.pairs[d$word.length == i] <- scale(d$min.pairs[d$word.length == i])  
    d$scaled.mean.dist.nonmp[d$word.length == i] <- scale(d$mean.dist.nonmp[d$word.length == i])
  }
  
  pl.unscaled <- ggplot(data = d, aes(x = word.length, y = min.pairs, label = word)) +
    theme(text = element_text(size = 20)) +
    ggtitle('English') + scale_x_continuous(breaks = 2:8) +  
    xlab('Word Length') + ylab('Neighborhood size') +
    geom_point(alpha = .2) + geom_smooth(method = 'loess')
  
  ggsave(pl.unscaled, filename = paste0(outdir, 'redun_methods/minpair-unscaled.png'), 
         width = 10, height = 8)   
  
  pl.scaled <- ggplot(data = d, aes(x = word.length, y = scaled.min.pairs, label = word)) +
    theme(text = element_text(size = 20)) +
    ggtitle('English') + scale_x_continuous(breaks = 2:8) +  
    xlab('Word Length') + ylab('Neighborhood size\n(z-scored by length)') +
    geom_point(alpha = .2) + geom_smooth(method = 'loess')
  
  ggsave(pl.scaled, filename = paste0(outdir, 'redun_methods/minpair-scaled.png'), 
         width = 10, height = 8)
  
  pl.unscaled <- ggplot(data = d, aes(x = word.length, y = mean.dist.nonmp, label = word)) +
    theme(text = element_text(size = 20)) +
    ggtitle('Raw Values') + scale_x_continuous(breaks = 2:8) +  
    xlab('Word Length') + ylab('Mean Edit Distance (excl. min pairs)') +
    geom_point(alpha = .2) + geom_smooth(method = 'loess')
  
  pl.scaled <- ggplot(data = d, aes(x = word.length, y = scaled.mean.dist.nonmp, label = word)) +
    theme(text = element_text(size = 20)) +
    ggtitle('z-scored by length') + scale_x_continuous(breaks = 2:8) +  
    xlab('Word Length') + ylab('Mean Edit Distance (excl. min pairs)\n(z-scored by length)') +
    geom_point(alpha = .2) + geom_smooth(method = 'loess')
  
  pl.grid <- grid.arrange(pl.unscaled, pl.scaled, ncol = 2)
  
  ggsave(pl.grid, filename = paste0(outdir, 'redun_methods/med.png'), 
         width = 11, height = 7)

}

efficiency_methods_graphs <- function() {
  # truncated Zipf graph for segment distribution
  d <- lf(paste0(lexdir, 'English.lex.txt'))
  d <- subset(d, word.length > 1)
  fd <- first_dist(d, p_threshold = 1, N = 2, keep_running = T)
  fd$closest <- rank(abs(fd$running.token.p - .05))
  fd$Included <- 'Not Included'
  fd$Included[fd$running.token.p >= .05] <- 'Included'
  pl <- ggplot(data = fd, aes(x = p.rank, y = seg.p, label = phone, color = Included)) + 
    xlab('Probability Rank') + ylab('Total Probability') + 
    geom_text(data = subset(fd, p.rank > 20)) + 
    geom_text_repel(data = subset(fd, p.rank <= 20),
                    segment.alpha = 1, 
                    arrow = arrow(length = unit(0.02, "npc"))) +
    scale_color_manual(name = 'Included in Analysis', values = c('#00BFC4', '#F8766D')) +
    theme(text = element_text(size=18)) +
    geom_vline(data = subset(fd, closest == 1), 
               aes(xintercept = I(p.rank)),
               linetype = 'dotted')
  ggsave(paste0(outdir, 'efficiency_methods/trunc-zipf-segs.png'), plot = pl,
         width = 9, height = 7)
  
  pl <- ggplot(data = fd, aes(x = log(mean), y= log(count), label = phone)) + 
    xlab('log Mean Token Frequency') + ylab('log Type Count') +
    theme(text = element_text(size=18)) +
    scale_color_manual(name = 'Included in Analysis', values = c('#00BFC4', '#F8766D')) +
    geom_point(data = subset(fd, Included == 'Not Included'), aes(color = Included),
               alpha = .5) +
    geom_point(data = subset(fd, Included == 'Included'), aes(color = Included),
               alpha = .5) +
    geom_smooth(data = subset(fd, Included == 'Included'), aes(color = Included),
                method = 'lm', size = 1.5) +
    geom_smooth(method = 'lm', linetype = 'dashed', color = 'black', size = 1.5)
  ggsave(paste0(outdir, 'efficiency_methods/trunc-zipf-lines.png'), plot = pl,
         width = 9, height = 7)
  
  pl <- ggplot(data = fd, aes(y = log(freq), x = log(count), label = phone)) + 
    ylab('log Total Token Frequency') + xlab('log Type Count') +
    theme(text = element_text(size=18)) +
    scale_color_manual(name = 'Included in Analysis', values = c('#00BFC4', '#F8766D')) +
    geom_point(data = subset(fd, Included == 'Not Included'), aes(color = Included),
               alpha = .5) +
    geom_point(data = subset(fd, Included == 'Included'), aes(color = Included),
               alpha = .5) +
    geom_smooth(data = subset(fd, Included == 'Included'), aes(color = Included),
                method = 'lm', size = 1.5) +
    geom_smooth(method = 'lm', linetype = 'dashed', color = 'black', size = 1.5)
  ggsave(paste0(outdir, 'efficiency_methods/trunc-zipf-lines-total.png'), plot = pl,
         width = 9, height = 7)
  
  # all langs
  all_langs <- lapply(lex.files, function(f) {
    d <- lf(paste0(lexdir, f))
    f <- gsub('.lex.txt', '', f)
    print(f)
    
    d <- subset(d, word.length > 1)
    fd <- first_dist(d, p_threshold = 1, N = 2, keep_running = T)
    fd$closest <- rank(abs(fd$running.token.p - .05))
    fd$Included <- 'Not Included'
    fd$Included[fd$running.token.p >= .05] <- 'Included'
    f <- gsub('.lex.txt', '', f)
    fd$language <- f
    fd$mean.z <- z.(log(fd$mean))
    fd$count.z <- z.(log(fd$count))
    fd$freq.z <- z.(log(fd$freq))
    fd
  })
  all_langs <- ldply(all_langs, data.frame)
  
  pl <- ggplot(data = all_langs, aes(x = mean.z, y = count.z, label = phone)) + 
    xlab('log Mean Token Frequency (z-scored)') + ylab('log Type Count (z-scored)') +
    xlim(c(-1.5,4)) + ylim(c(-2.5,3)) +
    theme(text = element_text(size=16)) +
    scale_color_manual(name = 'Included in Analysis', values = c('#00BFC4', '#F8766D')) +
    geom_point(data = subset(all_langs, Included == 'Not Included'), 
               aes(color = Included), alpha = .1) +
    geom_point(data = subset(all_langs, Included == 'Included'), 
               aes(color = Included), alpha = .1) +
    geom_smooth(data = subset(all_langs, Included == 'Included'), aes(color = Included),
                method = 'lm', size = 1.5) +
    geom_smooth(linetype = 'dotted', color = 'black',
                method = 'lm', size = 1.25) +
    facet_wrap(~language)
  ggsave(paste0(outdir, 'efficiency_methods/trunc-zipf-lines-all.png'), plot = pl,
         width = 9, height = 7)
  
  pl <- ggplot(data = all_langs, aes(y = freq.z, x = count.z, label = phone)) + 
    ylab('log Total Token Frequency (z-scored)') + xlab('log Type Count (z-scored)') +
    #ylim(c(-1.5,4)) + xlim(c(-2.5,3)) +
    theme(text = element_text(size=16)) +
    scale_color_manual(name = 'Included in Analysis', values = c('#00BFC4', '#F8766D')) +
    geom_point(data = subset(all_langs, Included == 'Not Included'), 
               aes(color = Included), alpha = .1) +
    geom_point(data = subset(all_langs, Included == 'Included'), 
               aes(color = Included), alpha = .1) +
    geom_smooth(data = subset(all_langs, Included == 'Included'), aes(color = Included),
                method = 'lm', size = 1.5) +
    geom_smooth(linetype = 'dotted', color = 'black',
                method = 'lm', size = 1.25) +
    facet_wrap(~language)
  ggsave(paste0(outdir, 'efficiency_methods/trunc-zipf-lines-all-total.png'), plot = pl,
         width = 9, height = 7)
    

  fd <- subset(fd, Included == 'Included')
  token.pl <- ggplot(data = fd, aes(mean)) + 
    theme(text = element_text(size=16)) +
    ylab('Count') + xlab('Mean Token Freq.') +
    geom_histogram(bins = 30)
  log.token.pl <- ggplot(data = fd, aes(log(mean))) + 
    theme(text = element_text(size=16)) +
    ylab('Count') + xlab('log Mean Token Freq.') +
    geom_histogram(bins = 30)
  type.pl <- ggplot(data = fd, aes(count)) + 
    theme(text = element_text(size=16)) +
    ylab('Count') + xlab('Type Count') +
    geom_histogram(bins = 30)
  log.type.pl <- ggplot(data = fd, aes(log(count))) + 
    theme(text = element_text(size=16)) +
    ylab('Count') + xlab('log Type Count') +
    geom_histogram(bins = 30)
  
  pl <- grid.arrange(token.pl, type.pl,
               log.token.pl, log.type.pl)
  ggsave(paste0(outdir, 'efficiency_methods/trunc-zipf-log.png'), plot = pl,
         width = 9, height = 7)
  
  
  all_langs <- pblapply(lex.files, function(f) {
    d <- lf(paste0(lexdir, f))
    d$language <- sub('.lex.txt', '', f)
    d
  })
  all_langs <- ldply(all_langs, data.frame)  
  all_langs <- subset(all_langs, word.length > 1)
  
  all_langs$first.seg <- substr(all_langs$phones,1,1)
  all_langs$first.biphone <- substr(all_langs$phones,1,2)
  
  first.segs <- aggregate(freq ~ first.seg + language, subset(all_langs, language == 'English'), length)
  first.biphones <- aggregate(freq ~ first.biphone + language, subset(all_langs, language == 'English'), length)
  
  ggplot(data = first.segs, aes(freq)) +
    geom_histogram() + facet_wrap(~language)
  
  ggplot(data = first.biphones, aes(freq)) +
    geom_histogram() + facet_wrap(~language)
  
  #######
  all_chrts <- read_tsv(paste0(outdir, 'cohort_p/all_chrts.txt'), col_types = cols())
  df <- pblapply(unique(all_chrts$language), function(f) {
    chrt <- subset(all_chrts, language == f)
    
    running_p <- 0
    chrt$running.p <- 0
    chrt$running.p <- 0
    for (cohort in chrt$cohort[order(chrt$cohort.p)]) {
      p <- chrt$cohort.p[match(cohort, chrt$cohort)]
      running_p <- running_p + p
      chrt$running.p[match(cohort, chrt$cohort)] <- running_p
    }
    chrt <- subset(chrt, running.p >= (1 - .95))
    
    rbind(data.frame(language = f, var = 'cohort',
                     r = cor.test(chrt$log.prob, chrt$cohort.count)$estimate),
          data.frame(language = f, var = 'segs',
                     r = cor.test(chrt$log.prob, chrt$nb.cont.segs)$estimate))
  }, cl = nb_cores)
  df <- ldply(df, data.frame)
  df$b0 <- 0
  df$r.str <- paste('r=', substr(df$r, 2, 4))
  
  pl <- ggplot(data = df, aes(color = var, label = r.str)) +
    theme(text = element_text(size = 20),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank()) +
    xlim(c(0,1)) + ylim(c(0,.6)) +
    xlab('Factor (see legend)') + ylab('log Cohort Prob.') + 
    scale_color_discrete(name = 'Corr. w/ log Cohort Prob.',
                         labels = c('Cohort size (type)',
                                    'Nb. Cont. Segs.')) +
    geom_text(data = subset(df, var == 'segs'), 
              x = 0, y = .5, size = 6,
              hjust   = 'left') +
    geom_text(data = subset(df, var == 'cohort'), 
              x = 0, y = .35, size = 6,
              hjust   = 'left') +
    geom_abline(aes(slope = r, intercept = b0, color = var), size = 1.35) + 
    facet_wrap(~language)
  ggsave(filename = paste0(outdir, 'efficiency_methods/seg-cohort-r.png'), plot = pl,
         width = 10, height = 8)
  1
}

syn_methods_graphs <- function() {
  d <- lf(paste0(lexdir, 'English.lex.txt'))
  d.minpairs <- read_tsv(paste0(outdir, 'rdn1/English.txt'), col_types = cols())
  d <- subset(d, word.length > 1)
  d$min.pairs <- d.minpairs$min.pairs[match(d$word, d.minpairs$word)]
  
  shuf_groups_up <- subset(aggregate(freq ~ word.length + up, d, length), freq > 1)
  shuf_groups_mp <- subset(aggregate(freq ~ word.length + min.pairs, d, length), freq > 1)
  
  pl1 <- ggplot(data = shuf_groups_up, aes(freq)) + 
    geom_histogram(bins = 10, fill = '#00BFC4') +
    theme(text = element_text(size = 18)) +
    ggtitle('post-Uniqueness Point Segs.') +
    xlab('Length ~ Redun. Segs. group size') + ylab('Count') 
  pl2 <- ggplot(data = shuf_groups_mp, aes(freq)) + 
    geom_histogram(bins = 10, fill = '#F8766D') +
    theme(text = element_text(size = 18)) +
    ggtitle('Number of Neighbors') +
    xlab('Length ~ Neighbors group size') + ylab('Count') 
  
  pl <- ggarrange(pl1, pl2) 
  ggsave(plot = pl, filename = paste0(outdir, 'synergy_methods/group_size.png'),
         width = 10, height = 7)
  }

# conclusion graphs
conclusion_graphs <- function() {
  all_langs <- lapply(lex.files, function(f) {
    d <- lf(paste0(lexdir, f))
    f <- sub('.lex.txt', '', f)
    print(f)
    d$language <- f
    d$mean.seg.h <- d$unig.h / d$word.length
    d$mean.mean <- mean(d$mean.seg.h)
    d$mean.seg.h.z <- z.(d$mean.seg.h)
    d
  })
  all_langs <- ldply(all_langs, data.frame)
  
  pl <- ggplot(data = subset(all_langs, mean.seg.h.z > -3 & mean.seg.h.z < 3),
               aes(mean.seg.h)) + geom_histogram(bins = 20) +
    geom_vline(aes(xintercept = mean.mean), size = 1.5, linetype = 'dashed', color = '#F8766D') +
    theme(text = element_text(size=18)) +
    xlab('Mean Token-based Segment Information per Word') + ylab('Count') +
    facet_wrap(~language)
  ggsave(paste0(outdir, 'conclusion/mean-seg-h.png'), pl,
         width = 10, height = 8)
  
  all_langs <- lapply(lex.files, function(f) {
    d <- lf(paste0(lexdir, f))
    f <- sub('.lex.txt', '', f)
    print(f)
    
    d$mean.seg.h <- d$unig.h / d$word.length
    lex <- lex2si(d)
    
    rbind(data.frame(language = f, 
                     H = lex_ent(lex), 
                     lex.type = 'og'),
          data.frame(language = f,
                     H = sum(d$unig.p * (d$mean.seg.h ** -2) * d$mean.seg.h * nchar(d$phones)),
                     lex.type = 'optimal')
    )
  })
  all_langs <- ldply(all_langs, data.frame)
  
  pl <- ggplot(data = all_langs, aes(x = language, y = H, color = lex.type, shape = lex.type)) +
    geom_point(size = 3) +
    xlab('Language') + ylab('Modified Lexical Entropy') +
    theme(axis.text.x = element_text(angle = 90),
          text = element_text(size=18)) +
    scale_color_manual(name = 'Lexicon Type',
                       labels = c('Real-world', "'Optimal' lexicon"),
                       values = c('#F8766D', '#00BFC4')) +
    scale_shape_discrete(name = 'Lexicon Type',
                         labels = c('Real-world', "'Optimal' lexicon"))
  ggsave(paste0(outdir, 'conclusion/real-v-optimal.png'), pl,
         width = 10, height = 8)
  
  # graph for slopes of 'thwart' and 'story'
  d <- lf(paste0(lexdir, 'English.lex.txt'))
  lex <- lex2si(d)
  pl <- ggplot(subset(lex, word %in% c('thwart', 'story')), aes(x = position, y = seg.h.type, color = word)) +
    geom_point(size = 2) + geom_line(size = 1.25) +
    xlab('Segment Position') + ylab('Segment Information') +
    theme(text = element_text(size=20)) +
    scale_color_discrete(name = 'Word')
  ggsave(paste0(outdir, 'conclusion/segh-slopes.png'), pl,
         width = 9, height = 8)
  
  d <- lf(paste0(lexdir, 'English.lex.txt'))
  d <- subset(d, word.length > 1)
  fd <- first_dist(d, p_threshold = 1, N = 2, keep_running = T)
  
  fd$closest <- rank(abs(fd$running.token.p - .05))
  fd$Included <- 'Not Included'
  fd$Included[fd$running.token.p >= .05] <- 'Included'
  fd_included <- sum(fd$Included == 'Included')
  pl.og <- ggplot(data = fd, aes(x = p.rank, y = seg.p, label = phone, color = Included)) + 
    xlab('Probability Rank') + ylab('Total Probability') +
    ggtitle('Real-world lexicon') +
    geom_text(data = subset(fd, p.rank > 20)) + 
    geom_text_repel(data = subset(fd, p.rank <= 20),
                    segment.alpha = 1, 
                    arrow = arrow(length = unit(0.02, "npc"))) +
    scale_color_manual(name = 'Included in Analysis', values = c('#00BFC4', '#F8766D')) +
    theme(text = element_text(size=18)) +
    geom_vline(data = subset(fd, closest == 1), 
               aes(xintercept = I(p.rank)),
               linetype = 'dotted')
  
  d.scr <- shuffle_frequencies(d, keep_length = F)
  fd.scr <- first_dist(d.scr, p_threshold = 1, N = 2, keep_running = T)
  
  fd.scr$closest <- rank(abs(fd.scr$running.token.p - .1))
  fd.scr$Included <- 'Not Included'
  fd.scr$Included[fd.scr$running.token.p >= .1] <- 'Included'
  fd.scr_included <- sum(fd.scr$Included == 'Included')
  pl.scr <- ggplot(data = fd.scr, aes(x = p.rank, y = seg.p, label = phone, color = Included)) + 
    xlab('Probability Rank') + ylab('Total Probability') + 
    ggtitle('Probability-shuffled lexicon') +
    geom_text(data = subset(fd.scr, p.rank > 20)) + 
    geom_text_repel(data = subset(fd.scr, p.rank <= 20),
                    segment.alpha = 1, 
                    arrow = arrow(length = unit(0.02, "npc"))) +
    scale_color_manual(name = 'Included in Analysis', values = c('#00BFC4', '#F8766D')) +
    theme(text = element_text(size=18)) +
    geom_vline(data = subset(fd.scr, closest == 1), 
               aes(xintercept = I(p.rank)),
               linetype = 'dotted')
  
  pl <- ggarrange(pl.og, pl.scr,
                  nrow = 2,
                  common.legend = TRUE, legend="right") 
  ggsave(paste0(outdir, 'conclusion/removed-segs.png'), pl,
         width = 10, height = 8)
  
  
}

# post-hoc
posthoc_rdn_1a <- function() {
  all_langs <- pblapply(lex.files, function(f) {
    f <- sub('.lex.txt', '', f)
    d <- read_tsv(paste0(outdir, 'rdn1/', f, '.txt'), col_types = cols())
    
    d$network.prop <- d$min.pairs / d$network
    
    d$scaled.mean.dist <- 0
    d$scaled.mean.dist.all <- 0
    d$scaled.mean.dist.nonmp <- 0
    d$scaled.min.pairs <- 0
    d$scaled.min.pairs.2 <- 0
    d$scaled.min.pairs.3 <- 0
    d$scaled.network <- 0
    d$scaled.network.prop <- 0
    for (i in min(d$word.length):max(d$word.length)) {
      d$scaled.mean.dist[d$word.length == i] <- 
        scale(d$mean.dist[d$word.length == i])
      
      d$scaled.mean.dist.all[d$word.length == i] <- 
        scale(d$mean.dist.all[d$word.length == i])
      
      d$scaled.mean.dist.nonmp[d$word.length == i] <-
        scale(d$mean.dist.nonmp[d$word.length == i])
      
      d$scaled.min.pairs[d$word.length == i] <- 
        scale(d$min.pairs[d$word.length == i])
      
      d$scaled.min.pairs.2[d$word.length == i] <- 
        scale(d$min.pairs.2[d$word.length == i])
      
      d$scaled.min.pairs.3[d$word.length == i] <- 
        scale(d$min.pairs.3[d$word.length == i])
      
      d$scaled.network[d$word.length == i] <- 
        scale(d$network[d$word.length == i])
      
      d$scaled.network.prop[d$word.length == i] <-
        scale(d$network.prop)
    }
    d[is.na(d)] <- 0
    
    d
  })
  all_langs <- ldply(all_langs, data.frame) 
  all_langs <- add_language_family(all_langs)
  all_langs <- subset(all_langs, word.length >= 2 & word.length <= 8)
  all_langs$scaled.log.prob <- -all_langs$scaled.unig.h
  all_langs$is.semitic <- 'Other'
  all_langs$is.semitic[all_langs$language %in% c('Arabic', 'Hebrew')] <- 'Semitic'
  all_langs$is.semitic <- factor(all_langs$is.semitic, levels = c('Semitic', 'Other'))
  
  pl <- ggplot(all_langs, aes(scaled.network.prop, fill = is.semitic)) +
    theme(text = element_text(size = 20)) +
    xlim(c(-2,6)) + xlab('Rel. Nbrhood Size (z-scored by length per language)') + ylab('Count') +
    scale_fill_discrete(name = 'Language Family', guide = guide_legend(reverse=TRUE)) +
    geom_histogram(data = subset(all_langs, is.semitic == 'Other'), bins = 30, alpha = .5) +
    geom_histogram(data = subset(all_langs, is.semitic == 'Semitic'), bins = 30)
    
  ggsave(plot = pl, filename = paste0(outdir, 'posthoc/rdn1/np-hist.png'),
         width = 10, height = 7)
  
  df <- aggregate(min.pairs ~ language + word.length + is.semitic, all_langs, mean)
  df.2 <- aggregate(min.pairs.2 ~ language + word.length + is.semitic, all_langs, mean)
  df <- gather(join(df, df.2), 'dist.type', 'mean.word.types', -language, -word.length, -is.semitic)
  
  pl <- ggplot(data = df, 
         aes(x = word.length, y = mean.word.types, linetype = dist.type, color = is.semitic)) +
    theme(text = element_text(size = 20)) +
    xlab('Word Length') + ylab('Mean Word Forms') +
    scale_color_manual(name = 'Language Family', values = c('#00BFC4', '#F8766D')) +
    scale_linetype_discrete(name = 'Distance Type', labels = c('Distance 1', 'Distance 2')) +
    geom_line(data = aggregate(mean.word.types ~ word.length + dist.type + is.semitic, df, mean),
              size = 1.5)
  
  ggsave(plot = pl, filename = paste0(outdir, 'posthoc/rdn1/np-line.png'),
         width = 11, height = 7)
  
}

posthoc_rdn_2a <- function() {
  ### Turksih
  d <- lf(paste0(outdir, 'posthoc/dfs/Turkish_wiki.lex.txt'), T)
  d$redun.prop <- (d$word.length - d$up) / d$word.length
  d$scaled.unig.h <- 0
  d$scaled.redun.prop <- 0
  for (i in 1:max(d$word.length)) {
    d$scaled.unig.h[d$word.length == i] <-
      scale(d$unig.h[d$word.length == i])
    d$scaled.redun.prop[d$word.length == i] <-
      scale(d$redun.prop[d$word.length == i])
  }
  d$scaled.log.prob <- -d$scaled.unig.h
  
  ct1 <- cor.test(d$scaled.log.prob, d$scaled.redun.prop)
  rho1 <- substr(ct1$estimate, 1, 5)
  x_range <- c(max(-2.5, min(d$scaled.log.prob)), min(2.5, max(d$scaled.log.prob)))
  y_range <- c(max(-2.5, min(d$scaled.redun.prop)), min(2.5, max(d$scaled.redun.prop)))
  pl1 <- ggplot(data = d, aes(x = scaled.log.prob, y = scaled.redun.prop)) +
    xlim(x_range) + ylim(y_range) +
    theme(text = element_text(size = 20)) +
    geom_point(alpha = .1) + geom_smooth(method = 'lm', size = 2) + 
    ggtitle('Turkish (Wikipedia)') +
    xlab('log Word Probability\n(z-scored)') + ylab('Redundant Segments\n(z-scored)')
  xadj <- layer_scales(pl1)$x$range$range[2] - .15 * (layer_scales(pl1)$x$range$range[2] - layer_scales(pl1)$x$range$range[1])
  yadj <- layer_scales(pl1)$y$range$range[2] - .15 * (layer_scales(pl1)$y$range$range[2] - layer_scales(pl1)$y$range$range[1])
  pl1 <- pl1 + geom_text(x = xadj, y = yadj, size = 7, 
                         label = paste0('r = ', rho1, '\n', pval_str(ct1$p.value)))
  
  ggsave(plot = pl1, filename = paste0(outdir, 'posthoc/rdn2/turkish-1.png'),
         width = 9, height = 7)
  
  d <- subset(d, up < word.length)
  lex <- lex2si(d)
  diphone <- diphone_h(lex, use_freq = F)
  # word redundancy with diphones
  wrwd <- diphone_redun(diphone)
  wrwd$mean.redun.seg.h <- wrwd$postup.sum.diphone.h / wrwd$segments.redundant 
  wrwd$seg.h.redun.prop <- wrwd$postup.sum.diphone.h / wrwd$sum.diphone.h
  wrwd$scaled.unig.h <- 0
  wrwd$scaled.redundancy.prop <- 0
  wrwd$scaled.mean.redun.seg.h <- 0
  wrwd$scaled.seg.h.redun.prop <- 0
  for (i in 1:max(wrwd$word.length)) {
    wrwd$scaled.unig.h[wrwd$word.length == i] <- 
      scale(wrwd$unig.h[wrwd$word.length == i])
    wrwd$scaled.redundancy.prop[wrwd$word.length == i] <- 
      scale(wrwd$redundancy.prop[wrwd$word.length == i])
    wrwd$scaled.mean.redun.seg.h[wrwd$word.length == i] <- 
      scale(wrwd$mean.redun.seg.h[wrwd$word.length == i])
    wrwd$scaled.seg.h.redun.prop[wrwd$word.length == i] <- 
      scale(wrwd$seg.h.redun.prop[wrwd$word.length == i])
  }
  wrwd[is.na(wrwd)] <- 0.
  wrwd$scaled.bigram.prob <- -wrwd$scaled.mean.redun.seg.h
  wrwd$scaled.log.prob <- -wrwd$scaled.unig.h
  
  ct2 <- cor.test(wrwd$scaled.log.prob, wrwd$scaled.bigram.prob)
  rho2 <- substr(as.character(ct2$estimate), 1, 5)
  x_range <- c(max(-2.5, min(wrwd$scaled.log.prob)), min(2.5, max(wrwd$scaled.log.prob)))
  y_range <- c(max(-2.5, min(wrwd$scaled.bigram.prob)), min(2.5, max(wrwd$scaled.bigram.prob)))
  pl2 <- ggplot(data = wrwd, aes(x = scaled.log.prob, y = scaled.bigram.prob)) +
    xlim(x_range) + ylim(y_range) +
    theme(text = element_text(size = 20)) +
    geom_point(alpha = .1) + geom_smooth(method = 'lm', size = 2) + 
    ggtitle('Turkish (Wikipedia)') +
    xlab('log Word Probability\n(z-scored)') + ylab('Mean log Bigram Prob.\n(z-scored)')
  xadj <- layer_scales(pl2)$x$range$range[2] - .15 * (layer_scales(pl2)$x$range$range[2] - layer_scales(pl2)$x$range$range[1])
  yadj <- layer_scales(pl2)$y$range$range[2] - .15 * (layer_scales(pl2)$y$range$range[2] - layer_scales(pl2)$y$range$range[1])
  pl2 <- pl2 + geom_text(x = xadj, y = yadj, size = 7, 
                         label = paste0('r = ', rho2, '\n', pval_str(ct2$p.value)))
  
  ggsave(plot = pl2, filename = paste0(outdir, 'posthoc/rdn2/turkish-2.png'),
         width = 9, height = 7)
  
  ct3 <- cor.test(wrwd$scaled.log.prob, wrwd$scaled.seg.h.redun.prop)
  rho3 <- substr(as.character(ct3$estimate), 1, 5)
  x_range <- c(max(-2.5, min(wrwd$scaled.log.prob)), min(2.5, max(wrwd$scaled.log.prob)))
  y_range <- c(max(-2.5, min(wrwd$scaled.seg.h.redun.prop)), min(2.5, max(wrwd$scaled.seg.h.redun.prop)))
  pl3 <- ggplot(data = wrwd, aes(x = scaled.log.prob, y = scaled.seg.h.redun.prop)) +
    xlim(x_range) + ylim(y_range) +
    theme(text = element_text(size = 20)) +
    geom_point(alpha = .1) + geom_smooth(method = 'lm', size = 2) + 
    ggtitle('Turkish (Wikipedia)') +
    xlab('log Word Probability\n(z-scored)') + ylab('Rel. Phonotactic Prob.\n(z-scored)')
  xadj <- layer_scales(pl3)$x$range$range[2] - .15 * (layer_scales(pl3)$x$range$range[2] - layer_scales(pl3)$x$range$range[1])
  yadj <- layer_scales(pl3)$y$range$range[2] - .15 * (layer_scales(pl3)$y$range$range[2] - layer_scales(pl3)$y$range$range[1])
  pl3 <- pl3 + geom_text(x = xadj, y = yadj, size = 7, 
                         label = paste0('r = ', rho3, '\n', pval_str(ct3$p.value)))
  
  ggsave(plot = pl3, filename = paste0(outdir, 'posthoc/rdn2/turkish-3.png'),
         width = 9, height = 7)
  
  pl <- ggarrange(pl1, pl2, pl3, rows = 2)
  ggsave(plot = pl, filename = paste0(outdir, 'posthoc/rdn2/grid.png'),
         width = 9, height = 7)
  1
}

posthoc_cohort_entropy_a <- function(load_files = F) {
  if (!load_files) {
    all_chrts <- read_tsv(paste0(outdir, 'cohort_p/all_chrts.txt'), col_types = cols())
    
    langs_to_do <- subset(unique(all_chrts$language), grepl('Georgian|Hebrew|Russian|Slovak|Spanish', unique(all_chrts$language)))
    all_langs <- lapply(langs_to_do, function(f) {
      print(paste0(f, ' - cohort ent by prob'))
      chrt <- subset(all_chrts, language == f)
      
      d <- lf(paste0(lexdir, f, '.lex.txt'), do_up = T)
      dn <- load_nonce(paste0(noncedir, f, '.new.txt'))
      
      chrt <- entropy_by_cohort(d, ignore_phones = 'aeiouAEIOU')
      m <- cohort_H_m(chrt, nb_segs = F)
      sm <- summary(m)
      print(sm)
      print(vif(m))
      
      print(paste0('nonce lexicons for: ', f, ' (nb. segs)'))
      nonce_bs <- pblapply(1:nb_scramble, function(i) {
        d.nonce <- nonce_lexicon(d, dn, keep_length = T, do_up = T)
        d.nonce$phones <- substr(d.nonce$phones, 1, d.nonce$up)
        chrt.nonce <- entropy_by_cohort(d.nonce, min_per_cohort = 1)
        
        m.nonce <- cohort_H_m(chrt.nonce, nb_segs = T)
        sm.nonce <- summary(m.nonce) 
        
        rbind(data.frame(b = sm.nonce$coefficients[2], t = sm.nonce$coefficients[2,3],
                         p = sm.nonce$coefficients[2,4], lex.type = 'nonce'))
      }, cl = nb_cores)
      df <- rbind(data.frame(b = sm$coefficients[2], t = sm$coefficients[2,3],
                             p = sm$coefficients[2,4], lex.type = 'og'),
                  ldply(nonce_bs, data.frame))
      print(paste0(df$b[1], ' -- ', mean(df$b), ' : ',
                   sum(df$b > df$b[1]) / nb_scramble))
      
      df$language <- f
      df$b.z <- 0
      df$b.z <- scale(df$b)
      df$b.less <- sum(df$b > df$b[1]) / nb_scramble
      write_tsv(df, paste0(outdir, 'posthoc/cohort_p/', f, '.segs.txt'))
      df
    })
  }
  
  langs_to_do <- subset(lex.files, grepl('Georgian|Hebrew|Russian|Slovak|Spanish', lex.files))
  all_langs <- lapply(langs_to_do, function(f) {
    f <- sub('.lex.txt', '', f)
    df <- read_tsv(file = paste0(outdir, 'posthoc/cohort_p/', f, '.segs.txt'), col_types = cols())
    df
  })
  all_langs <- ldply(all_langs, data.frame)
  all_langs$lex.type <- factor(all_langs$lex.type, levels = c('og', 'nonce'))
  
  
  print('cohort p (nonce lexicons; nb. segs)')
  unique_langs <- unique(all_langs$language)
  for (i in seq(1,length(unique_langs), 5)) {
    print(i)
    langs <- unique_langs[i:(i+4)]
    all_langs.sub <- subset(all_langs, language %in% langs)
    
    pl <- ggplot(data = all_langs.sub, aes(b.z, fill = lex.type)) +
      xlab('Model Estimate: log Cohort Prob.') + ylab('Count') +
      geom_histogram(alpha = .3, bins = max(30, nb_scramble / 100)) +
      geom_vline(xintercept = -2, size = 1, linetype = 'dotted') +
      geom_vline(xintercept = 2, size = 1, linetype = 'dotted') +
      geom_vline(data = subset(all_langs.sub, lex.type == 'og'), 
                 aes(xintercept = b.z), size = 1.5, color = '#F8766D') +
      geom_text(mapping = aes(x = Inf, y = Inf, 
                              label = paste0('p < ', substr(pmax(b.less, .001), 1, 5))),
                hjust   = 'right',
                vjust   = 'top',
                size = 5) +
      scale_fill_manual(name = 'Lexicon Type',
                        labels = c('Real-world', 'Nonce Forms'),
                        values = c('#F8766D', '#00BFC4')) +
      facet_wrap(~language) +
      theme(text = element_text(size=20), 
            legend.position = c(.95, .15), 
            legend.justification = c(1, 0))
    
    
    ggsave(filename = paste0(outdir, 'posthoc/cohort_p/nonce/',i ,'.png'), pl,
           width = 9, height = 7)
    1
  }
  1
}

posthoc_cohort_entropy_b <- function(load_files = F) {
  if (!load_files) {
    d <- lf(paste0(outdir, 'posthoc/dfs/SP_wiki.lex.txt'))
    dn <- load_nonce(paste0(outdir, 'posthoc/dfs/SP_wiki.new.txt'))
    
    chrt <- entropy_by_cohort(d)
    m <- cohort_H_m(chrt)
    sm <- summary(m)
    print(sm)
    print(vif(m))
    
    f <- 'Spanish_wiki'
    print(paste0('nonce lexicons for: ', f, ' (nb. segs)'))
    nonce_bs <- pblapply(1:nb_scramble, function(i) {
      d.nonce <- nonce_lexicon(d, dn, keep_length = T, do_up = T)
      d.nonce$phones <- substr(d.nonce$phones, 1, d.nonce$up)
      chrt.nonce <- entropy_by_cohort(d.nonce, min_per_cohort = 1)
      
      m.nonce <- cohort_H_m(chrt.nonce)
      sm.nonce <- summary(m.nonce) 
      
      rbind(data.frame(b = sm.nonce$coefficients[2], t = sm.nonce$coefficients[2,3],
                       p = sm.nonce$coefficients[2,4], lex.type = 'nonce'))
    }, cl = nb_cores)
    df <- rbind(data.frame(b = sm$coefficients[2], t = sm$coefficients[2,3],
                           p = sm$coefficients[2,4], lex.type = 'og'),
                ldply(nonce_bs, data.frame))
    print(paste0(df$b[1], ' -- ', mean(df$b), ' : ',
                 sum(df$b > df$b[1]) / nb_scramble))
    
    df$language <- f
    df$b.z <- 0
    df$b.z <- scale(df$b)
    df$b.less <- sum(df$b > df$b[1]) / nb_scramble  
    write_tsv(x = df, paste0(outdir, 'posthoc/cohort_p/Spanish_wiki.txt'))
    1
  }
  df <- read_tsv(paste0(outdir, 'posthoc/cohort_p/Spanish_wiki.txt'), col_types = cols())
  
  df$lex.type <- factor(df$lex.type, levels = c('og', 'nonce'))
  pl <- ggplot(data = df, aes(b.z, fill = lex.type)) +
    ggtitle('Spanish (Wikipedia data)') + ylab('Count') +
    xlab('Model Estimate: log Cohort Prob.') + 
    geom_histogram(alpha = .3, bins = max(30, nb_scramble / 100)) +
    geom_vline(xintercept = -2, size = 1, linetype = 'dotted') +
    geom_vline(xintercept = 2, size = 1, linetype = 'dotted') +
    geom_vline(data = subset(df, lex.type == 'og'), 
               aes(xintercept = b.z), size = 1.5, color = '#F8766D', alpha = .7) +
    geom_text(mapping = aes(x = Inf, y = Inf, 
                            label = paste0('p < ', substr(pmax(b.less, .001), 1, 5))),
              hjust   = 'right',
              vjust   = 'top',
              size = 8) +
    scale_fill_manual(name = 'Lexicon Type',
                      labels = c('Real-world', 'Nonce Forms'),
                      values = c('#F8766D', '#00BFC4')) +
    theme(text = element_text(size=18))
  
  ggsave(paste0(outdir, 'posthoc/cohort_p/nonce/Spanish_wiki.png'), pl,
         width = 9, height = 7)

  1
}

posthoc_cohort_entropy_c <- function(load_files = F) {
  if (!load_files) {
    d <- lf(paste0(lexdir, 'Georgian.lex.txt'), do_up = T)
    dn <- load_nonce(paste0(noncedir, 'Georgian.new.txt'))
    
    chrt <- entropy_by_cohort(d, include_post_up = F, ignore_phones = '^ptkPTKbdgBDGqQ')
    m <- cohort_H_m(chrt)
    sm <- summary(m)
    print(sm)
    print(vif(m))
    
    f <- 'Georgian_stops'
    print(paste0('nonce lexicons for: ', f, ' (nb. segs)'))
    nonce_bs <- pblapply(1:nb_scramble, function(i) {
      d.nonce <- nonce_lexicon(d, dn, keep_length = T, do_up = T)
      chrt.nonce <- entropy_by_cohort(d.nonce, , include_post_up = F, ignore_phones = '^ptkPTKbdgBDGqQ')
      
      m.nonce <- cohort_H_m(chrt.nonce)
      sm.nonce <- summary(m.nonce) 
      
      rbind(data.frame(b = sm.nonce$coefficients[2], t = sm.nonce$coefficients[2,3],
                       p = sm.nonce$coefficients[2,4], lex.type = 'nonce'))
    }, cl = nb_cores)
    df <- rbind(data.frame(b = sm$coefficients[2], t = sm$coefficients[2,3],
                           p = sm$coefficients[2,4], lex.type = 'og'),
                ldply(nonce_bs, data.frame))
    print(paste0(df$b[1], ' -- ', mean(df$b), ' : ',
                 sum(df$b > df$b[1]) / nb_scramble))
    
    df$language <- f
    df$b.z <- 0
    df$b.z <- scale(df$b)
    df$b.less <- sum(df$b > df$b[1]) / nb_scramble  
    write_tsv(x = df, paste0(outdir, 'posthoc/cohort_p/Georgian_stops.txt'))
    1
  }
  df <- read_tsv(paste0(outdir, 'posthoc/cohort_p/Georgian_stops.txt'), col_types = cols())
  
  df$lex.type <- factor(df$lex.type, levels = c('og', 'nonce'))
  pl <- ggplot(data = df, aes(b.z, fill = lex.type)) +
    xlab('Model Estimate: log Cohort Prob.') + ylab('Count') +
    ggtitle('Georgian Stops') +
    geom_histogram(alpha = .3, bins = max(30, nb_scramble / 100)) +
    geom_vline(xintercept = -2, size = 1, linetype = 'dotted') +
    geom_vline(xintercept = 2, size = 1, linetype = 'dotted') +
    geom_vline(data = subset(df, lex.type == 'og'), 
               aes(xintercept = b.z), size = 1.5, color = '#F8766D', alpha = .7) +
    geom_text(mapping = aes(x = Inf, y = Inf, 
                            label = paste0('p < ', substr(pmax(b.less, .001), 1, 5))),
              hjust   = 'right',
              vjust   = 'top',
              size = 8) +
    scale_fill_manual(name = 'Lexicon Type',
                      labels = c('Real-world', 'Nonce Forms'),
                      values = c('#F8766D', '#00BFC4')) +
    theme(text = element_text(size=18))
  
  ggsave(paste0(outdir, 'posthoc/cohort_p/nonce/Georgian_stops.png'), pl,
         width = 9, height = 7)
  
  1
}

posthoc_synergy <- function() {
  d <- lf(paste0(outdir, 'posthoc/dfs/Dutch_wiki.lex.txt'))
  mps <- med(d, min(d$word.length), max(d$word.length))
  mps <- subset(mps, select = c(word, min.pairs))
  
  d$min.pairs <- mps$min.pairs[match(d$word, mps$word)]
  d <- subset(d, min.pairs >= 0)
  
  lex <- lex2si(d, shuffle = F, shuffle_same_nbrhood = F)
  # remove words of length 1
  lex <- subset(lex, word.length >= 2)
  
  df <- data.frame(i = 1, lex.type = 'original',
                   seq.H = lex_ent(lex))
  
  # neighbors
  df.nbr.shuffled <- pblapply(1:nb_scramble + 1, function(i) {
    lex.scr <- lex2si(d, shuffle = T, shuffle_same_nbrhood = T)
    lex.scr <- subset(lex.scr, word.length >= 2)
    
    data.frame(i = i, lex.type = 'shuffled',
               seq.H = lex_ent(lex.scr))
  }, cl = nb_cores)
  df.nbr.shuffled <- ldply(df.nbr.shuffled, data.frame)
  
  # tail
  df.tail.shuffled <- pblapply(1:nb_scramble + 1, function(i) {
    lex.scr <- lex2si(d, shuffle = T, shuffle_same_tail = T)
    
    data.frame(i = i, lex.type = 'shuffled',
               seq.H = lex_ent(lex.scr))
  }, cl = nb_cores)
  df.tail.shuffled <- ldply(df.tail.shuffled, data.frame)
  
  df.tail <- rbind(df, df.tail.shuffled)
  df.tail$seq.H.scaled <- scale(df.tail$seq.H)
  df.tail$diff.from.mean <- df.tail$seq.H - mean(df.tail$seq.H)
  df.tail$seq.greater <- sum(df.tail$seq.H > df.tail$seq.H[1]) / (nrow(df.tail) - 1)
  
  df.nbr <- rbind(df, df.nbr.shuffled)
  df.nbr$seq.H.scaled <- scale(df.nbr$seq.H)
  df.nbr$diff.from.mean <- df.nbr$seq.H - mean(df.nbr$seq.H)
  df.nbr$seq.greater <- sum(df.nbr$seq.H > df.nbr$seq.H[1]) / (nrow(df.nbr) - 1)
  
  pl.tail <- ggplot(data = df.tail, aes(seq.H.scaled, fill = lex.type)) +
    scale_fill_manual(name = 'Lexicon Type',
                      labels = c('Real-world', 'Probability-shuffled'),
                      values = c('#F8766D', '#00BFC4')) +
    xlab('Modified Lexical Entropy\n(z-scored)') + ylab('Count') +
    geom_histogram(alpha = .3, bins = max(30, nb_scramble / 10)) +
    geom_vline(data = subset(df.tail, lex.type == 'original'), 
               aes(xintercept = seq.H.scaled), color = '#F8766D', size = 1.5) +
    geom_text(mapping = aes(x = -Inf, y = Inf, 
                            label = paste0('  p < ', substr(pmax(seq.greater, .001), 1, 5))),
              hjust   = 'left',
              vjust   = 'top',
              size = 6) +
    theme(text = element_text(size=16)) +
    ggtitle('Segments in Redundant Tails')
  
  pl.nbr <- ggplot(data = df.nbr, aes(seq.H.scaled, fill = lex.type)) +
    scale_fill_manual(name = 'Lexicon Type',
                      labels = c('Real-world', 'Probability-shuffled'),
                      values = c('#F8766D', '#00BFC4')) +
    xlab('Modified Lexical Entropy\n(z-scored)') + ylab('Count') +
    geom_histogram(alpha = .3, bins = max(30, nb_scramble / 10)) +
    geom_vline(data = subset(df.nbr, lex.type == 'original'), 
               aes(xintercept = seq.H.scaled), color = '#F8766D', size = 1.5) +
    geom_text(mapping = aes(x = -Inf, y = Inf, 
                            label = paste0('  p < ', substr(pmax(seq.greater, .001), 1, 5))),
              hjust   = 'left',
              vjust   = 'top',
              size = 6) +
    theme(text = element_text(size=16)) +
    ggtitle('Neighborhood Size')
    
  pl <- ggarrange(pl.nbr, pl.tail, common.legend = TRUE, legend="right")
  ggsave(paste0(outdir, 'posthoc/synergy/Dutch_wiki.png'), pl,
         width = 10, height = 7)
  
  1
}

############################################
lex.files <- list.files(lexdir)
lex.files <- subset(lex.files, grepl('lex.txt$', lex.files))

nb_cores <- 7
nb_scramble <- 1000

if (T) {
##### 
# redundancy
# part 1 - minimal pairs and mead edit distance
tst_redundancy_1_a()
tst_redundancy_1_b()
  
tst_redundancy_2_a()
tst_redundancy_2_b()

# part 2 - efficiency
tst_balanced_sig(load_files = T)
tst_balanced_nonce(load_files = T)

tst_first_later_ent()
tst_cohort_entropy_by_prob()
tst_total_entropy()

# part 3 - efficiency and redundancy
tst_synergy_same_tail()
tst_synergy_same_nbrhood()
}



defense_graphs <- function() {
  d <- lf(paste0(lexdir, 'English.lex.txt'))

  lex <- lex2si(d, shuffle = F)
  D.s <- subset(lex, phone %in% c('D', 's') & position == 1)
  D.s$log.freq <- log(D.s$freq)

  colors <- c('#00BFC4', '#F8766D')
  a <- aggregate(freq ~ phone, D.s, sum)
  a$log.freq <- log(a$freq)
  pl.D <- ggplot(data = D.s, aes(log.freq, fill = phone)) +
    scale_x_continuous(breaks = 2:12, limits = c(2, 12)) + theme(text = element_text(size = 20)) +
    xlab('log Word Frequency') + ylab('Count') +
    geom_histogram(data = subset(D.s, phone == 'D'), bins = 8, fill =  colors[1]) +
    ggtitle('[ð]-initial words')

  pl.s <- ggplot(data = D.s, aes(log.freq, fill = phone)) +
    scale_x_continuous(breaks = 2:12, limits = c(2, 12)) + theme(text = element_text(size = 20)) +
    xlab('log Word Frequency') + ylab('Count') +
    geom_histogram(data = subset(D.s, phone == 's'), bins = 16, fill =  colors[2]) +
    ggtitle('[s]-initial words')

  pl <- ggarrange(pl.D, pl.s, nrow = 2)
  ggsave(plot = pl, filename = paste0(outdir, 'defense/D-s.1.png'), width = 9, height = 8)

  pl.D <- pl.D + geom_vline(data = subset(a, phone == 'D'), aes(xintercept = log.freq), color = colors[1], size = 2)
  pl.s <- pl.s + geom_vline(data = subset(a, phone == 's'), aes(xintercept = log.freq), color = colors[2], size = 2)

  pl <- ggarrange(pl.D, pl.s, nrow = 2)
  ggsave(plot = pl, filename = paste0(outdir, 'defense/D-s.2.png'), width = 9, height = 8)


  ### shuffled
  lex <- lex2si(d, shuffle = T)
  D.s <- subset(lex, phone %in% c('D', 's') & position == 1)
  D.s$log.freq <- log(D.s$freq)

  colors <- c('#00BFC4', '#F8766D')
  a <- aggregate(freq ~ phone, D.s, sum)
  a$log.freq <- log(a$freq)
  pl.D <- ggplot(data = D.s, aes(log.freq, fill = phone)) +
    scale_x_continuous(breaks = 2:12, limits = c(2, 12)) + theme(text = element_text(size = 20)) +
    xlab('log Word Frequency') + ylab('Count') +
    geom_histogram(data = subset(D.s, phone == 'D'), bins = 8, fill =  colors[1]) +
    ggtitle('[ð]-initial words')

  pl.s <- ggplot(data = D.s, aes(log.freq, fill = phone)) +
    scale_x_continuous(breaks = 2:12, limits = c(2, 12)) + theme(text = element_text(size = 20)) +
    xlab('log Word Frequency') + ylab('Count') +
    geom_histogram(data = subset(D.s, phone == 's'), bins = 16, fill =  colors[2]) +
    ggtitle('[s]-intitial words')

  pl <- ggarrange(pl.D, pl.s, nrow = 2)
  ggsave(plot = pl, filename = paste0(outdir, 'defense/shuf_D-s.1.png'), width = 9, height = 8)

  pl.D <- pl.D + geom_vline(data = subset(a, phone == 'D'), aes(xintercept = log.freq), color = colors[1], size = 2)
  pl.s <- pl.s + geom_vline(data = subset(a, phone == 's'), aes(xintercept = log.freq), color = colors[2], size = 2)

  pl <- ggarrange(pl.D, pl.s, nrow = 2)
  ggsave(plot = pl, filename = paste0(outdir, 'defense/shuf_D-s.2.png'), width = 9, height = 8)

}