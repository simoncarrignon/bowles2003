model <- function(allgroups,s,t,steps,m,k,e,mu_s=0.1,mu_t=0.1,increment=.1,logs=c(),baseline=10,c=1,b=2){
	store.q=c()
	store.s=c()
	store.t=c()
	store.n=c()
	print("starting simulation")
	for(ts in 0:steps){
		if(ts %%100 ==0)print(paste("steps ",ts))
		for( ig in seq_along(allgroups)){
			if("pairing" %in% logs)print(paste("pairing, dilema and mutation gp",ig))
			#PAIRING and PD
			group <- allgroups[[ig]]
			oldpayoff <- sum(group[,"payoff"])
			group[,"payoff"]=baseline
			#print(paste("group",ig,"starts with",oldpayoff))
			if(nrow(group)>2){
				pairs <- getAgentPairs(group,s=s[ig])
				stopifnot(all(table(pairs)==1))
				for(i in 1:ncol(pairs)){
					a=pairs[1,i]
					b=pairs[2,i]
					results <- PD(group[a,"type"],group[b,"type"])
					group[c(a,b),"payoff"] <- group[c(a,b),"payoff"]+results
					if("PD" %in% logs) print(paste0("pairing ", a, " (type:", group[a,"type"], ") with ", b, " (type:", group[b,"type"], "), out:", paste0(results,collapse=",")))
				}
				if(t[ig]!=0){
					taxtotal=sum(group[,"payoff"]*t[ig])/nrow(group)
					group[,"payoff"] = group[,"payoff"]*(1-t[ig]) + taxtotal
				}
				newpayoff <- sum(group[,"payoff"])
				if("PD" %in% logs)print(paste("group",ig," payoff change: ",newpayoff-oldpayoff))
			}

			#Reproduction
			newind <- sample(1:nrow(group),size=nrow(group),prob=group[,"payoff"]/newpayoff,replace=TRUE)
			group <- group[newind,]

			#Mutation

			mutate <- runif(nrow(group))<e
			if(sum(mutate)>0) group[mutate,"type"] <- sample(c(0,1),size=sum(mutate),replace=TRUE)
			if("mutation" %in% logs) print(paste(sum(mutate),"agent mutate"))
			allgroups[[ig]] <- group
		}
		#institution mutaiton
        s  <- mutateInstitutions(s,min.value = 0,max.value = 0.5,increment = 0.1,mut.rate = mu_s)
		t  <- mutateInstitutions(t,min.value = 0,max.value = 1 ,increment = 0.1,mut.rate = mu_t)

		########
		#Miration 

		#count migrant one way:
		# migrants <- lapply(allgroups,function(g)sample(x=1:nrow(g),size=nrow(g)*m))
		#       another way:
		migrants <- lapply(allgroups,function(g)which(runif(nrow(g))<m))
		if("migration" %in% logs)print(paste(sum(lengths(migrants)),"will be moving"))
		aftermigration <- allgroups
		if(sum(lengths(migrants))>0){
			## removing migrant from original groups
			for(i in 1:length(migrants)){
				cur.migrants <- migrants[[i]]
				if(length(cur.migrants)>0){
					aftermigration[[i]] <- allgroups[[i]][-cur.migrants,]
					if("migration" %in% logs)print(paste(" ",length(cur.migrants),"are moving from",i,"(pop will decrease from",nrow(allgroups[[i]]),"to",nrow(aftermigration[[i]]),")"))
				}
			}
			## moving them to their new groups
			for(i in 1:length(migrants)){
				cur.migrants <- migrants[[i]]
				potential.new <- seq_along(allgroups)[-i] #group where the migrand can move
				if(length(potential.new)==1) #if one group, no choice they all move there
					new.groups <- rep(potential.new,length(migrants[[i]]))
				else #else the randomly choose one
					new.groups <- sample(potential.new,size=length(migrants[[i]]),replace=T)
				#actual move' people (think I could do that cleaner
				for(ng in unique(new.groups)){
					movers <- allgroups[[i]][cur.migrants[new.groups==ng],,drop=F]
					aftermigration[[ng]] <- rbind(aftermigration[[ng]],movers)
					if("migration" %in% logs)print(paste(" ",nrow(movers),"are moving to",ng,"from",i,"(thus pop will go from",nrow(allgroups[[ng]]),"to",nrow(aftermigration[[ng]]),")"))
				}
			}
		}
		allgroups <- aftermigration
		########Group competition
		# attakers <- sample(x=1:length(allgroups),size=length(allgroups)*k)
		#       another way:
		stopifnot(length(allgroups)>1)

		attakers <- which(runif(length(allgroups))<k) #in the paper, it looks like all group in 'attakers' are involved, and if uneven then another one is drawn
		if(length(attakers)>0){
			nonaattaker=seq_along(allgroups)[-attakers]
			if(length(attakers)%%2 !=0){
				forced=ifelse(length(nonaattaker)==1,nonaattaker,sample(nonaattaker,size=1))
				attakers=c(attakers,forced)
			}

			#potentially one group can attack another, and the group attacked can then again be attacked?
			stopifnot(length(attakers)>1)
			attakers=sample(attakers) #randomize order of attack so it's not always the same attacking first
			for(ai in 1:(length(attakers)/2)){
				at <- attakers[ai]
				def <- attakers[ai+(length(attakers)/2)]
				atp <- sum(allgroups[[at]][,"payoff"])-0.5*(s[at]^2+t[at]^2) 
				defp <- sum(allgroups[[def]][,"payoff"])-0.5*(s[at]^2+t[at]^2) 
				win=loose=NULL
				if(atp>defp){win=at;loose=def}
				if(atp<defp){win=def;loose=at}
				if(atp==defp){draw=sample(c(at,def));win=draw[1];loose=draw[2]}
				w=allgroups[[win]]
				l=allgroups[[loose]]
				l=w[sample(x=1:nrow(w),size=nrow(l),replace=T),] #replace loosing group by sample of winner
				new=rbind(w,l) #merge group
				splitorder <- sample(nrow(new)) #the fusionned individuals are evenly redistributed in two groups
				new.win <- splitorder[1:floor(length(splitorder)/2)] 
				new.loose <- splitorder[floor((length(splitorder)/2)+1):length(splitorder)]
				stopifnot(all(c(new.win,new.loose)==splitorder))
				allgroups[[win]] <- new[new.win,]
				allgroups[[loose]] <- new[new.loose,]

				#group level transmission
				s[loose] <- s[win]
				t[loose] <- t[win]
				# Print additional information
				if ("fusion" %in% logs) {
					print(paste("Attack between groups", at, "(payoff:", atp, ") and", def, "(payoff:", defp, ")"))
					print(paste("Group", win, "wins against", loose))
					print(paste("Before split: win group size:", nrow(w), "loose group size:", nrow(l)))
					print(paste("After split: win group size:", length(new.win), "loose group size:", length(new.loose)))
				}
			}

		}
		store.q=c(store.q,sum(unlist(sapply(allgroups,function(i)(i[,"type"]))))/sum(sapply(allgroups,nrow)))
		store.s=c(store.s,mean(s))
		store.t=c(store.t,mean(t))
		#store.n=c(store.n,sum(sapply(allgroups,nrow)))
	}

	return(cbind(store.q,store.t,store.s))
}
