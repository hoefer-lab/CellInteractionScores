



#######
# 'interactionmatrix' function below retrieves expression value of ligands and receptors in the interaction database from a count matrix

#
####arguments###
#-counts' is a count matrix with genes in row, if our interaction database is to be used, the rownames in the count matrix should be ENSEMBL mouse genes identifiers.
#-interactions': the interaction database
#-datasheet: description of the samples in the count matrix, the rows of the datasheet must be ordered in the same way as the columns in the count matrix. the first column should contains  sample names, the second column should specify to what cell_type or condition belongs the sample, the script generates interactomes for pairs of cell populations. The third column must indicate either "receptor" or "ligd", this determines whether receptors or ligands expression values are extracted from the count matrix for a given sample.

interactionmatrix=function(counts,interactions,datasheet)
{
    interactionmatPlus=matrix(ncol=ncol(counts),nrow=nrow(interactions))
    interactionmatPlus=data.frame(interactionmatPlus)
    
    for (i in 1:length(interactions[,1]))
    {
        
        if (interactions[i,2]%in%rownames(counts)&&interactions[i,1]%in%rownames(counts))#we exclude genes which are not in the the count matrix
        {
            interactionmatPlus[i,datasheet[,3]%in%"ligd"]=counts[rownames(counts)==interactions[i,2],datasheet[,3]%in%"ligd"]
            interactionmatPlus[i,datasheet[,3]%in%"receptor"]=counts[rownames(counts)==interactions[i,1],datasheet[,3]%in%"receptor"]
        }
        else
        {
            interactionmatPlus[i,]=rep(NA,ncol(interactionmatPlus))
        }
    }
    interactionmatPlus[interactionmatPlus<0]=0 # if there are negative values in the count matrix(can happen with log normalization for instance, they will be counted as 0)
    return(interactionmatPlus)
}





#######
#'interactionranking' generates sets of Interaction Scores for pairs of interacting cells

#you will need to install the library dplyr for dense ranking. The ranking is based, for each cell population (the datasheet determines to what cell type belongs a sample), on the mean expression values of the receptors or ligands.
#### arguments###
#datasheet (same as the one from 'interactionmatrix')
# interactionmatrix is the result of the 'interactionmatrix' function

interactionranking=function(datasheet,interactionmatrix)
{
    library(dplyr)
    l=1
    name=vector()
    datasheetlgd=datasheet[datasheet[,3]%in%"ligd",]
    datasheetreceptor=datasheet[datasheet[,3]%in%"receptor",]
    #first we create matrices to contain the interaction scores, the ranks of the ligands and the ranks of the receptors, each column is a pair of interacting cells, and each row is an interaction from the interaction database
    Score=matrix(ncol=length(unique(datasheetlgd[,2]))*length(unique(datasheetreceptor[,2])),nrow=nrow(interactionmatrix))
    rankLigs=matrix(ncol=length(unique(datasheetlgd[,2]))*length(unique(datasheetreceptor[,2])),nrow=nrow(interactionmatrix))
    rankRecs=matrix(ncol=length(unique(datasheetlgd[,2]))*length(unique(datasheetreceptor[,2])),nrow=nrow(interactionmatrix))
    
    
    #we now compute the ranks of receptors and ligands and then sum these ranks
    for(i in 1:length(unique(datasheetlgd[,2])))
    {
        for (j in 1:length(unique(datasheetreceptor[,2])))
        {
            rankLig=dense_rank(rowMeans(interactionmatrix[,which(datasheet[,2]%in%unique(datasheetlgd[,2])[i])]))
            
            rankLigs[,l]=rankLig
            rankRec=dense_rank(rowMeans(interactionmatrix[,which(datasheet[,2]%in%unique(datasheetreceptor[,2])[j])]))
            rankRecs[,l]=rankRec
            Score[,l]=(rankLig+rankRec)
            name[l]=paste(unique(datasheetreceptor[,2])[j],"/",unique(datasheetlgd[,2])[i],sep="")
            l=l+1
        }
    }
    #for each pair of cells divides the scores to the highest score in the interactome so that each interaction is scored from 0 to 1
    for (i in 1:ncol(Score))
    {
        Score[,i]=Score[,i]/max(Score[,i],na.rm=TRUE)
    }
    colnames(Score)=name
    Score=cbind(Score,1:nrow(Score))
    
    #stores the Scores, the ranks of the receptors and the ranks of the ligands in a list
    EverythingInItsRightPlace=list()
    EverythingInItsRightPlace[[1]]=Score
    EverythingInItsRightPlace[[2]]=rankRecs
    EverythingInItsRightPlace[[3]]=rankLigs

    return(EverythingInItsRightPlace)
}

##################################################

#### the functions below can be used to look at the specificity of the interactions. Here I try to reply to the following question: for a given pair of interacting cells and a given interaction how high is the Interactionscore compared to interactionscores of other pairs of cells####

#####
#  the function 'norm' converts values into z-scores, it is used in the next function 'divergent'
norm=function(values)
{
    normal=((values-mean(values))/sd(values))
}

#####
# for each interaction in the DB, 'divergent' computes the z-scores for every pairs of cells. It returns a list of matrices (one per pair of cells) containing Interaction Scores and z-scores

####Arguments:###
#- interactions is the the RL interaction database
# - Scores contains selected sets of Interaction Scores (first of the 3 matrices produced by the the "interactionranking" function),


divergent=function(interactions,Scores)
{
    normScores=apply(Scores[,1:ncol(Scores)-1],1,norm)
    colnames(normScores)=1:(nrow(Scores))

    informativeResult=list()
    for(i in 1:nrow(normScores))
    {
        SortedScores=sort(normScores[i,],decreasing=TRUE)
        
        SortedScores2=interactions[as.numeric(names(SortedScores)),]
        informativeResult[[i]]=cbind(SortedScores2,Scores[SortedScores2[,5],i],SortedScores)
        colnames(informativeResult[[i]])[c(5:7)]=c("Interaction.ID","InteractionScore","Z-score")
        informativeResult[[i]][,c(6:7)]=round(informativeResult[[i]][,c(6:7)],digits=3)
    }
    return(informativeResult)
}




