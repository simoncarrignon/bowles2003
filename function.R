createGroup <- function(N,baseline=10){
  pop=matrix(nrow=N,ncol=2)
  colnames(pop)=c("type","payoff")
  pop[,1]=0
  pop[,2]=baseline
  pop
}

getAgentPairs <- function(group,s=0,log=F){
  pairs=c()
  #pairing:
  if(s==0){
      pairs <- sample(1:nrow(group))
      midp <- floor(length(pairs)/2)
      pairs <- sapply(1:midp,function(i)pairs[c(i,i+midp)])
  }
  else{
      topair=sample(nrow(group))
      while(length(topair)>1){ ##this won't work yet as need to remove agent already paired
          i <- topair[1]
          topair <- topair[-1]
          i.type  <-  group[i,"type"]
          a.group <<-group
          a.topair <<-topair
          candidates <- c()
          if(runif(1)<s) candidates <- topair[group[topair,"type"]==i.type]
          else candidates  <- topair
          if(length(candidates)==0) {topair=c(i,topair);warnings("no candidate found")}
          else{
              j <- ifelse(length(candidates)==1,candidates,sample(candidates,size=1))
              pairs <- cbind(pairs,c(i,j))
              topair <- topair[-which(topair == j)]
              if(log)print(paste("pairing",i,"and",j,"left to pair=",length(topair),"(",paste0(topair,collapse=","),")"))
          }
          i=j=NULL;
      }
  }
  return(pairs)
}


PD <- function(i,j,c=1,b=2){
    if(i==0 && j==0)return(c(0,0))
    if(i==1 && j==1)return(c(b-c,b-c))
    if(i==1 && j==0)return(c(-c,b))
    if(i==0 && j==1)return(c(b,-c))
}


mutateInstitutions <- function(init.value,min.value=0,max.value=1,increment=0.1,mut.rate=.1)
{
    #institution mutaiton
    if(mut.rate>0){
        mutations <- runif(length(init.value))<mut.rate
        new.values=sapply(init.value[mutations],function(i)i+ifelse(runif(1)<.5,-increment,increment))
        new.values=ifelse(new.values<=min.value,min.value,new.values)
        new.values=ifelse(new.values>=max.value,max.value,new.values)
        init.value[mutations]=new.values
    }
    return(init.value)
}
