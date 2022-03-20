
##' @title Convert genotype coded in A/T/C/G to 0/1
##' @description Convert base (A/T/C/G) coded genotype to numeric (0/1) coded 
##' @param gmt a dataframe of genotype data of gamete cells
##' @param ref a character represents reference allele
##' @param alt a character represents alternative allele
##' @return a dataframe containing converted genotype
##' @export
##' @author Ruidong Li
##' @examples 
##' ref <- sample(c('A','T'),500, replace=TRUE)
##' alt <- sample(c('C','G'),500, replace=TRUE)
##' 
##' gmt <- data.frame(chr=rep(1,500), pos=seq_len(500),
##'     ref=ref, alt=alt, gmt1=ref, gmt2=alt, gmt3=ref,
##'     gmt4=ref, gmt5=c(alt[1:250], ref[251:500]),
##'     stringsAsFactors = FALSE)
##'     
##' gmtDa <- base2num(gmt=gmt[5:9], ref=ref, alt=alt)
base2num <- function(gmt, ref, alt) {
    newGMT <- apply(gmt, 2, function(x) x != ref)
    newGMT[newGMT==TRUE] <- 1
    newGMT[newGMT==FALSE] <- 0
    return (data.frame(newGMT))
}


##' @title Convert genotype coded in 0/1 to A/T/C/G
##' @description Convert numeric (0/1) coded genotype to base (A/T/C/G) coded 
##' @param hap a dataframe of consensus haplotypes
##' @param ref a character represents reference allele
##' @param alt a character represents alternative allele
##' @return a dataframe containing converted haplotypes
##' @export
##' @author Ruidong Li
##' @examples 
##' ref <- sample(c('A','T'),500, replace=TRUE)
##' alt <- sample(c('C','G'),500, replace=TRUE)
##' 
##' consensusHap <- data.frame(hap1=rep(0,500),hap2=rep(1,500),
##'     total=rep(5,500),rate=rep(1,500),
##'     confidence=rep('F',500),
##'     stringsAsFactors = FALSE)
##' rownames(consensusHap) <- seq_len(500)
##' 
##' hap <- num2base(hap=consensusHap, ref=ref, alt=alt)
num2base <- function(hap, ref, alt) {
    
    hap$hap1[hap$hap1==0] <- ref[hap$hap1==0]
    hap$hap1[hap$hap1=='1'] <- alt[hap$hap1=='1']
    hap$hap1[hap$hap1=='7'] <- NA
    
    hap$hap2[hap$hap2==0] <- ref[hap$hap2==0]
    hap$hap2[hap$hap2=='1'] <- alt[hap$hap2=='1']
    hap$hap2[hap$hap2=='7'] <- NA
    
    return (hap)
}



### Flip one allele to the althernative allele ###
flipFun <- function(v){
    v2 <- ifelse(v==7,7,ifelse(v==0, 1, 0))
    return (v2)
}

### HMM functions (from Dr. Lin Himmelmann's HMM package)
### https://rdrr.io/cran/HMM/src/R/HMM.r
### https://cran.r-project.org/web/packages/HMM/index.html
initHMM = function(States, Symbols, startProbs=NULL, transProbs=NULL,
    emissionProbs=NULL)
{
  nStates    = length(States)
  nSymbols   = length(Symbols)
  S          = rep(1/nStates,nStates)
  T          = 0.5*diag(nStates) + array(0.5/(nStates),c(nStates,nStates))
  E          = array(1/(nSymbols),c(nStates,nSymbols))
  names(S)   = States
  dimnames(T)= list(from=States,to=States)
  dimnames(E)= list(states=States,symbols=Symbols)
  if(!is.null(startProbs)){S[]  = startProbs[]}
  if(!is.null(transProbs)){T[,] = transProbs[,]}
  if(!is.null(emissionProbs)){E[,] = emissionProbs[,]}
  return(list(States=States,Symbols=Symbols,startProbs=S,transProbs=T,
      emissionProbs=E))
}

simHMM = function(hmm, length)
{
  hmm$transProbs[is.na(hmm$transProbs)]       = 0
  hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
  states   = c()
  emission = c()
  states   = c(states, sample(hmm$States,1,prob=hmm$startProbs))
  for(i in 2:length)
  {
    state  = sample(hmm$States, 1, prob=hmm$transProbs[states[i-1],])
  	states = c(states, state)
  }
  for(i in 1:length)
  {
    emi      = sample(hmm$Symbols, 1, prob=hmm$emissionProbs[states[i],])
  	emission = c(emission, emi)
  }
  return(list(states=states,observation=emission))
}

viterbi = function(hmm, observation)
{
  hmm$transProbs[is.na(hmm$transProbs)]       = 0
  hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
  nObservations  = length(observation)
  nStates    = length(hmm$States)
  v          = array(NA,c(nStates,nObservations))
  dimnames(v)= list(states=hmm$States,index=1:nObservations)
  # Init
  for(state in hmm$States)
  {
    v[state,1] = log(hmm$startProbs[state] * hmm$emissionProbs[state,observation[1]])
  }
  # Iteration
  for(k in 2:nObservations)
  {
    for(state in hmm$States)
    {
      maxi = NULL
      for(previousState in hmm$States)
      {
        temp = v[previousState,k-1] + log(hmm$transProbs[previousState,state]) 
        maxi = max(maxi, temp)
      }
      v[state,k] = log(hmm$emissionProbs[state,observation[k]]) + maxi
    }
  }
  # Traceback
  viterbiPath = rep(NA,nObservations)
  for(state in hmm$States)
  {
    if(max(v[,nObservations])==v[state,nObservations])
    {
      viterbiPath[nObservations] = state
      break
    }
  }
  for(k in (nObservations-1):1)
  {
    for(state in hmm$States)
    {
      if(max(v[,k]+log(hmm$transProbs[,viterbiPath[k+1]]))
          ==v[state,k]+log(hmm$transProbs[state,viterbiPath[k+1]]))
      {
        viterbiPath[k] = state
        break
      }
    }
  }
  return(viterbiPath)
}

forward = function(hmm, observation)
{
  hmm$transProbs[is.na(hmm$transProbs)]       = 0
  hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
  nObservations  = length(observation)
  nStates    = length(hmm$States)
  f          = array(NA,c(nStates,nObservations))
  dimnames(f)= list(states=hmm$States,index=1:nObservations)
  # Init
  for(state in hmm$States)
  {
    f[state,1] = log(hmm$startProbs[state] * hmm$emissionProbs[state,observation[1]])
  }
  # Iteration
  for(k in 2:nObservations)
  {
    for(state in hmm$States)
    {
      logsum = -Inf
      for(previousState in hmm$States)
      {
        temp   = f[previousState,k-1] + log(hmm$transProbs[previousState,state])
		if(temp > - Inf)
		{
			logsum = temp + log(1 + exp(logsum - temp ))
		}
      }
      f[state,k] = log(hmm$emissionProbs[state,observation[k]]) + logsum
    }
  }
  return(f)
}

backward = function(hmm, observation)
{
  hmm$transProbs[is.na(hmm$transProbs)]       = 0
  hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
  nObservations  = length(observation)
  nStates    = length(hmm$States)
  b          = array(NA,c(nStates,nObservations))
  dimnames(b)= list(states=hmm$States,index=1:nObservations)
  # Init
  for(state in hmm$States)
  {
    b[state,nObservations] = log(1)
  }
  # Iteration
  for(k in (nObservations-1):1)
  {
    for(state in hmm$States)
    {
      logsum = -Inf
      for(nextState in hmm$States)
      {
        temp   = b[nextState,k+1] + log(hmm$transProbs[state,nextState]*hmm$emissionProbs[nextState,observation[k+1]])
		if(temp > - Inf)
		{
        	logsum = temp + log(1 + exp(logsum-temp))
		}
      }
      b[state,k] = logsum
    }
  }
  return(b)
}
