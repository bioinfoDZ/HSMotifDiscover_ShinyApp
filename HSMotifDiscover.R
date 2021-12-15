HSMotifDiscover <- function(input_HSseq_file, motifLenVec,  charGrpFile, seq_weight_file, numCores=1,  affinity_threshold=0,  itr=40000, if_revCompStrand=FALSE)
{
  
  
  
  
  isRStudio <- Sys.getenv("RSTUDIO") == "1"
  
  if((Sys.info()[['sysname']]== 'Linux') | (Sys.info()[['sysname']]== 'Darwin')){
    os_Flag= 'good'
  }
  
  #motifLenVec=seq(motif_range[1], motif_range[2])
  
  if(missing(input_HSseq_file) & missing(motifLenVec)){
    geterrmessage("input sequence file and motif length/ range can't be left empty")
  }
  if (missing(seq_weight_file) & missing(charGrpFile))
  {
    print('>>>C1')
    if((os_Flag=='good') & (numCores > 1) & (isRStudio== FALSE)){
      motif_range_data_results=mclapply(motifLenVec , function(x) CompileResults(input_HSseq_file=input_HSseq_file, motif_length=x, affinity_threshold=affinity_threshold,  itr=itr, if_revCompStrand=if_revCompStrand), mc.cores = getOption("mc.cores", 2L))
    }else(
      motif_range_data_results=lapply(motifLenVec , function(x) CompileResults(input_HSseq_file=input_HSseq_file, motif_length=x,  affinity_threshold=affinity_threshold,  itr=itr, if_revCompStrand=if_revCompStrand))
    )
  }
  if(missing(seq_weight_file) & !missing(charGrpFile))
  {
    print('>>>C2')
    if((os_Flag=='good') & (numCores > 1) & (isRStudio== FALSE)){
      motif_range_data_results=mclapply(motifLenVec , function(x) CompileResults(input_HSseq_file=input_HSseq_file, motif_length=x, charGrpFile=charGrpFile, affinity_threshold=affinity_threshold, itr=itr, if_revCompStrand=if_revCompStrand),mc.cores = getOption("mc.cores", 2L) )
    }else(
      motif_range_data_results=lapply(motifLenVec , function(x) CompileResults(input_HSseq_file=input_HSseq_file, motif_length=x,  charGrpFile=charGrpFile, affinity_threshold=affinity_threshold,  itr=itr, if_revCompStrand=if_revCompStrand))
    )
  }
  if(!missing(seq_weight_file) & missing(charGrpFile))
  {
    print('>>>C3')
    if((os_Flag=='good') & (numCores > 1) & (isRStudio== FALSE)){
      motif_range_data_results=mclapply(motifLenVec , function(x) CompileResults(input_HSseq_file=input_HSseq_file, motif_length=x,  seq_weight_file=seq_weight_file, affinity_threshold=affinity_threshold,  itr=itr, if_revCompStrand=if_revCompStrand), mc.cores = getOption("mc.cores", 2L))
    }else(
      motif_range_data_results=lapply(motifLenVec , function(x) CompileResults(input_HSseq_file=input_HSseq_file, motif_length=x, seq_weight_file=seq_weight_file, affinity_threshold=affinity_threshold,  itr=itr, if_revCompStrand=if_revCompStrand))
    )
  }
  if(!missing(seq_weight_file) & !missing(charGrpFile))
  {
    print('>>>C4')
    if((os_Flag=='good') & (numCores > 1) & (isRStudio== FALSE)){
      motif_range_data_results=mclapply(motifLenVec , function(x) CompileResults(input_HSseq_file=input_HSseq_file, motif_length=x,  charGrpFile=charGrpFile,  seq_weight_file=seq_weight_file, affinity_threshold=affinity_threshold, itr=itr, if_revCompStrand=if_revCompStrand),  mc.cores = getOption("mc.cores", 2L))
    }else(
      motif_range_data_results=lapply(motifLenVec , function(x) CompileResults(input_HSseq_file=input_HSseq_file, motif_length=x, charGrpFile=charGrpFile,  seq_weight_file=seq_weight_file, affinity_threshold=affinity_threshold,  itr=itr, if_revCompStrand=if_revCompStrand))
    )
  }
  
  print(' out > HSMotifDiscover')
  
 
  
  for(j in 1:length(motifLenVec))
  {
    motif_range_data_results[[j]]=motif_range_data_results[[j]]$motifData
    print(paste0('Motif Length: ',motifLenVec[j] ))
    cat('\n=================\n')
    cat('Position Specific Scoring Matrix (PSSM): \n')
    print(motif_range_data_results[[j]]$PSSM)
    cat('Reverse complementry PSSM: \n')
    print(motif_range_data_results[[j]]$ReverseComplement_PSSM)
    cat('\nMotif statistics: \n')
    cat(paste0('Motif Entropy: ', motif_range_data_results[[j]]$MotifEntropy,', Information Content: ',motif_range_data_results[[j]]$IC,', Motif P-value: ',format(motif_range_data_results[[j]]$motif_Pval, digits = 4, scientific=T )))
    cat('\n\n')
    cat('\n')
  }
  
  names(motif_range_data_results)=paste0('MotifLength_',motifLenVec)
  


  
  return(list(motifData=motif_range_data_results[[1]]))
  
  #return(list(motifData=motif_range_data_results$motifData))
}
