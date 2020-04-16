#' ---
#' title: "Window plots"
#' author: "John Flournoy"
#' date: "`r format(Sys.Date())`"
#' output:
#'   html_document:
#'     toc: true
#'     toc_float: true
#'     code_folding: hide
#' ---

knitr::opts_chunk$set(results = 'asis', warning = FALSE, message = FALSE)

windowed_mean <- function(d, winmiddle, column, wincol, window = .5){
  se <- function(x) sqrt(var(x, na.rm = TRUE)/sum(!is.na(x)))
  win_min <- winmiddle - window
  win_max <- winmiddle + window
  d_ <- d[d[[wincol]] > win_min & d[[wincol]] < win_max, ]
  amean <- mean(d_[[column]], na.rm = TRUE)
  ase <- se(d_[[column]])
  wm <- list(mean = amean,
             se = ase,
             u = amean + qnorm(.975)*ase,
             l = amean + qnorm(.025)*ase)
  return(wm)
}
fit.gompertz <- function(data, time){
  d <- data.frame(y=data, t=time)
  
  # Must have at least 3 datapoints at different times
  if (length(unique(d$t)) < 3) stop("too few data points to fit curve")
  
  # Pick starting values ###
  i <- which.max(diff(d$y))
  starting.values <- c(a=max(d$y), 
                       mu=max(diff(d$y))/(d[i+1,"t"]-d[i, "t"]), 
                       lambda=d$t[i])
  print("Starting Values for Optimization: ")
  print(starting.values)
  ##########################
  
  formula.gompertz <- "y~a*exp(-exp(mu*exp(1)/a*(lambda-t)+1))"
  nls(formula.gompertz, d, starting.values)
}
window_mean_df <- function(d, winmiddle, column, wincol, window, othercols = c(), excludei = TRUE){
  df_winmean <- dplyr::bind_rows(lapply(1:dim(d)[1], function(i){
    winmiddle <- d[i, wincol][[1]]
    if(excludei){
      d <- d[-i,]
    }
    windowed_mean(d, winmiddle = winmiddle, wincol = wincol, column = column, window = window)
  }))
  
  newd <- dplyr::bind_cols(d[, c(column, wincol, othercols)], df_winmean)
  return(newd)
}
plot_windowed <- function(d, wincol, column, window, titletext = NULL, excludei = TRUE, showpoints = TRUE, wavecol = NULL, gomp = FALSE){
  library(ggplot2)
  if(gomp){
    if(inherits(d, 'list')){
      if(is.null(wavecol)){
        stop('Lists of data frames are assumed to be a list of data per wave. Nothing to do if they\'re not')
      }
      plotds <- lapply(d, function(ad){
        ad.g <- na.omit(ad[,c(column, wincol, wavecol)])
        mod <- fit.gompertz(ad.g[[column]], ad.g[[wincol]])
        ad.g$mean <- predict(mod)
        return(ad.g)
      })
      plotd <- dplyr::bind_rows(plotds)
      aplot <- ggplot(plotd, aes_string(x = wincol, y = 'mean', group = wavecol)) 
      if(showpoints){
        aplot <- aplot + 
          geom_segment(aes_string(x = wincol, xend = wincol, y = column, yend = 'mean'), alpha = .2, size = .5) + 
          geom_point(aes_string(y = column), size = 1, shape = 1)
      }
      aplot <- aplot +   
        geom_line() + 
        theme_minimal() + 
        labs(x = wincol, y = column, title = titletext)
    } else {
      d.g <- na.omit(d[,c(column, wincol, wavecol)])
      mod <- fit.gompertz(d.g[[column]], d.g[[wincol]])
      d.g$mean <- predict(mod)
      plotd <- d.g
      aplot <- ggplot(plotd, aes_string(x = wincol, y = 'mean'))
      if(showpoints){
        aplot <- aplot + 
          geom_segment(aes_string(x = wincol, xend = wincol, y = column, yend = 'mean'), alpha = .4, size = .5) + 
          geom_point(aes_string(y = column), size = 1, shape = 1)
      }
      aplot <- aplot +   
        geom_line() + 
        theme_minimal() + 
        labs(x = wincol, y = column, title = titletext)
    }
  } else {
    if(inherits(d, 'list')){
      if(is.null(wavecol)){
        stop('Lists of data frames are assumed to be a list of data per wave. Nothing to do if they\'re not')
      }
      plotds <- lapply(d, window_mean_df, wincol = wincol, column = column, window = window, excludei = excludei, othercols = wavecol)
      plotd <- dplyr::bind_rows(plotds)
      aplot <- ggplot(plotd, aes_string(x = wincol, y = 'mean', group = wavecol)) + 
        geom_ribbon(aes(ymin = l, ymax = u), alpha = .2)
      if(showpoints){
        aplot <- aplot + 
          geom_segment(aes_string(x = wincol, xend = wincol, y = column, yend = 'mean'), alpha = .2, size = .5) + 
          geom_point(aes_string(y = column), size = 1, shape = 1)
      }
      aplot <- aplot +   
        geom_line() + 
        theme_minimal() + 
        labs(x = wincol, y = column, title = titletext)
    } else {
      plotd <- window_mean_df(d, wincol = wincol, column = column, window = window, excludei = excludei) 
      aplot <- ggplot(plotd, aes_string(x = wincol, y = 'mean')) + 
        geom_ribbon(aes(ymin = l, ymax = u), alpha = .4)
      if(showpoints){
        aplot <- aplot + 
          geom_segment(aes_string(x = wincol, xend = wincol, y = column, yend = 'mean'), alpha = .4, size = .5) + 
          geom_point(aes_string(y = column), size = 1, shape = 1)
      }
      aplot <- aplot +   
        geom_line() + 
        theme_minimal() + 
        labs(x = wincol, y = column, title = titletext)
    }
  }
  return(aplot)
}

p <- readr::read_csv('Z:/dsnlab/TAG/projects/W1_W2_pubertal_timing/w1w2_all_sepres_long.csv')
w1 <- p[p$wave==1,]
w2 <- p[p$wave==2,]

cols <- names(p)[5:16]

#'
#' # Wave 1
#' 

invisible(lapply(cols, function(col){
  cat(paste0('

## ', col, '{.tabset}

### Exclude *i*th observation

'))
    print(plot_windowed(w1, wincol = 'age', column = col, window = .5, excludei = TRUE, 
                      titletext = paste0('Wave 1 ', col), showpoints = TRUE))
    cat('

### Include *i*th observation

')
    print(plot_windowed(w1, wincol = 'age', column = col, window = .5, excludei = FALSE, 
                        titletext = paste0('Wave 1 ', col), showpoints = TRUE))
    
    cat('

### Gompertz

')
    try(print(plot_windowed(w1, wincol = 'age', column = col, window = .5, excludei = FALSE, 
                        titletext = paste0('Wave 1 ', col), showpoints = TRUE, gomp = TRUE)))
}))

#'
#' # Wave 2
#' 

invisible(lapply(cols, function(col){
  cat(paste0('

## ', col, '{.tabset}

### Exclude *i*th observation

'))
  print(plot_windowed(w2, wincol = 'age', column = col, window = .5, excludei = TRUE, 
                      titletext = paste0('Wave 2 ', col), showpoints = TRUE))
  cat('

### Include *i*th observation

')
  print(plot_windowed(w2, wincol = 'age', column = col, window = .5, excludei = FALSE, 
                      titletext = paste0('Wave 2 ', col), showpoints = TRUE))
  
  cat('

### Gompertz

')
  try(print(plot_windowed(w2, wincol = 'age', column = col, window = .5, excludei = FALSE, 
                          titletext = paste0('Wave 1 ', col), showpoints = TRUE, gomp = TRUE)))
}))

#'
#' # Both waves
#' 

invisible(lapply(cols, function(col){
  cat(paste0('

## ', col, '{.tabset}

### Exclude *i*th observation

'))
  print(plot_windowed(list(w1, w2), wincol = 'age', column = col, window = .5, excludei = TRUE, 
                      titletext = paste0('Wave 1 ', col), showpoints = TRUE, wavecol = 'wave'))
  cat('

### Include *i*th observation

')
  print(plot_windowed(list(w1, w2), wincol = 'age', column = col, window = .5, excludei = FALSE, 
                      titletext = paste0('Wave 1 ', col), showpoints = TRUE, wavecol = 'wave'))

}))

#'
#' # Both waves TOGETHER
#' 

invisible(lapply(cols, function(col){
  cat(paste0('

## ', col, '{.tabset}

### Exclude *i*th observation

'))
  print(plot_windowed(p, wincol = 'age', column = col, window = .5, excludei = TRUE, 
                      titletext = paste0('both together', col), showpoints = TRUE))
  cat('

### Include *i*th observation

')
  print(plot_windowed(p, wincol = 'age', column = col, window = .5, excludei = FALSE, 
                      titletext = paste0('both together ', col), showpoints = TRUE))
  
  cat('

### Gompertz

')
  try(print(plot_windowed(p, wincol = 'age', column = col, window = .5, excludei = FALSE, 
                          titletext = paste0('both together ', col), showpoints = TRUE, gomp = TRUE)))
}))