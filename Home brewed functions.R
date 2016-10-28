  #### FUNCTIONS FILE ####

## Spectral analysis functions ##

trascosine<-function(data, q){
  #generates a matrix with the harmonics matching with the dimensions of the univariate dataset
  size<-dim(data)
  tempo<-c(1:length(data))
  peri<-(tempo-1/2)/length(tempo)
  psyt<-matrix(ncol=q, nrow=length(data))
  for (i in 1:q){
    psyt[,i]<-sqrt(2)*cos(pi*i*peri) 
  }
  psyt<-as.data.frame(psyt)
  headers<-as.character(c(1:q))
  return(psyt)
}
#generates a matrix with the harmonics matching with the dimensions of the univariate dataset


cosine<-function(span, period, phase, ampl){
  time<-rep(0:span)
  osc<-ampl*cos(2*pi*(period^(-1))*time+phase*pi)
  return(osc)  
}
# generates a matrix of waves with given timespan, period, phase and amplitude


randcosine<-function(ser, len){
  # generates random time series of given length (len variable), by summing up a given number of 
  # harmonics (as given by ser variable)
  ts<-vector(length=len+1)
  top<-round(len/3)
  
  for (i in 1:ser) {
    
    period<-max(ceiling(rbeta(1, 1, 15, ncp = 1)*1000), round(top))#should be a beta here
    phase<-ceiling(runif(1, min=-1, max=10)/10)
    ampl<-ceiling(runif(1,min=-1, max=5))
    ts<-ts+cosine(len, period, phase, ampl)
    
  }
  
  return(ts)
}
# generates random time series of given length (len variable), by summing up a given number of 
# harmonics (as given by ser variable)


periodogram<-function(data, oth){
  
  #case for even length of data
  lung<-length(data)
  lung.h<-lung/2
  lung.p<-lung.h+1
  fftdata<-fft(data)
  nonscal<-abs(fftdata/sqrt(lung))^2
  freq<-(0:lung.h)/lung
  P<-(oth/lung)*nonscal[1:lung.p]
  plot(freq, P, type='l', col='blue', xlab='Frequency', ylab='FFT value')
  title(main='Periodogram')
  mass<-which.max(nonscal)
  peak<-1/freq[mass]
  #print(max(nonscal))
  print(peak)
}
# returns the periodogram of a given time series, with the scaling factor precised by 'oth' variable


spergram<-function(serie, k){
  perio<-spec.pgram(serie, k, taper=0, log='no')
  mas<-which.max(perio$spec)
  peak<-1/perio$freq[mas]
  return(peak)
}
#applies spec.pgram with some standard parametrisation and returns the peak frequency


getYearQuarter <- function(x, 
                           firstMonth=7, 
                           fy.prefix='FY', 
                           quarter.prefix='Q',
                           sep='-',
                           level.range=c(min(x), max(x)) ) {
  if(level.range[1] > min(x) | level.range[2] < max(x)) {
    warning(paste0('The range of x is greater than level.range. Values ',
                   'outside level.range will be returned as NA.'))
  }
  quarterString <- function(d) {
    year <- as.integer(format(d, format='%Y'))
    month <- as.integer(format(d, format='%m'))
    y <- ifelse(firstMonth > 1 & month >= firstMonth, year+1, year)  
    q <- cut( (month - firstMonth) %% 12, breaks=c(-Inf,2,5,8,Inf), 
              labels=paste0(quarter.prefix, 1:4))
    return(paste0(fy.prefix, y, sep, q))
  }
  vals <- quarterString(x)
  levels <- unique(quarterString(seq(
    as.Date(format(level.range[1], '%Y-%m-01')), 
    as.Date(format(level.range[2], '%Y-%m-28')), by='month')))
  return(factor(vals, levels=levels, ordered=TRUE))
}
# to get and manipulate quarters in time series and generate dates


perdev<-function(indspf2){
  media<-apply(indspf2, 1, mean, na.rm=T)
  relativa<-as.data.frame(matrix(ncol=ncol(indspf2), nrow=nrow(indspf2)))
  for(i in 1:nrow(indspf2)) {
    relativa[i,]<-(indspf2[i,]-media[i])/media[i]
  }
  return(relativa)
}

# to transform a matrix in its deviations from the mean, row-wise -- TIME


perdevj<-function(indspf2){
  med<-apply(indspf2, 1, mean, na.rm=T)
  relativa<-(indspf2/med)-1
  return(relativa)
}

# to transform a matrix in its deviations from the mean, row-wise -- JAIME


trendev<-function(mat){
  matdat<-mat[,2:ncol(mat)]
  temp<-1:nrow(mat)
  temp2<-temp^2
  regr<-function(x){
    dta<-data.frame(x, temp, temp2)
    names(dta)<-c('x', 'temp', 'temp2')
    model<-lm(x~temp+temp2, data=dta)
    GAPS<-(model$residuals/(x-model$residuals))
    gaps<-as.matrix(na.omit(GAPS))
    gap<-gaps[nrow(gaps)]
    return(gap)
  }
  outcome<-apply(matdat, 2, regr)
  outcome<-as.matrix(outcome)
  return(outcome*100)
}

# given an array of different TS variables in the columns, this regress each coulmn
# against a quadratic time trend and recovers the percentage deviations of each column's
# last observation from the trend -- rather time consuming


oner<-function(x){
  x<-as.matrix(x)
  rifer<-sd(x, na.rm=T)*2+mean(x, na.rm=T)
  index<-matrix(nrow=nrow(x), ncol=1)
  for (i in 1:nrow(x)){
    
    if (is.na(x[i,1])){
      index[i,1]=NA
    }else if (abs(x[i,1])<rifer){
      index[i,1]=0
    } else #if (abs(x[i,1])>rifer){
      index[i,1]=1#}
  }
  return(as.data.frame(index))
}

# takes a variable vector and creates a dummy vector, taking value 1 when the corresponding observation
# is 2sd's above the mean



shift<-function(x,shift_by){
  stopifnot(is.numeric(shift_by))
  stopifnot(is.numeric(x))
  
  if (length(shift_by)>1)
    return(sapply(shift_by,shift, x=x))
  
  out<-NULL
  abs_shift_by=abs(shift_by)
  if (shift_by > 0 )
    out<-c(tail(x,-abs_shift_by),rep(NA,abs_shift_by))
  else if (shift_by < 0 )
    out<-c(rep(NA,abs_shift_by), head(x,-abs_shift_by))
  else 
    out<-x
  out
}

# provides a re-usable txt list of the installed packages
pack_list <- function(){
  pax <- as.data.frame(installed.packages())
  pax <- as.character(pax$Package)
  write.table(pax, file.path(getwd(),'pakketty.txt'), row.names=F, col.names=F, sep=';')
}


# creates a dummy vector: if there is increase in the input data it gives one,
# otherwise it assigns 0 for decrease and inaction.
dummifier <- function(x){
  y <- as.vector(x)
  n <- length(y)
  index <- matrix(nrow=n, ncol=1, rep(99, n))
  index[1,1] <- NA
  for (i in 2:n){
    if (is.na(y[i]) | is.na(y[i-1])){
      index[i] <- NA
    } else if (y[i]-y[i-1]>0){
      index[i]<-1
    } else {
      index[i] <- 0
    }
  }
  return(index)
}

# generating the indexes matrix: an upper triangle matrix of ones,
# as used in the theory from Hendry & Co.
# this is particularly slow

index_gen <- function(datafram){
  
  # generating the indexes matrix: an upper triangle matrix of ones
  # this is particularly slow
  
  nobs=nrow(datafram)
  indexes <- matrix(0, ncol=nobs, nrow=nobs)
  
  for (i in 1:nobs){
    for (l in 1:nobs){
      if (l>=i){
        indexes[i,l]=1
      }
    }
  }
  colnames(indexes) <- paste('obs', (1:nobs), sep='_')
  return(indexes[,-nobs])
}

# takes the output of pack_list to install back all the 
# previously installed packages
old_pax <- function(){
  pax <- read.table(file.choose())
  apply(pax,1,install.packages)
}
