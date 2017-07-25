rm(list=ls())

library(ggplot2)
library(reshape2)
library(gridExtra)
library (scales)

##  For debugging, otherwise pass in file root as parameter
## base = "Y:\\LiBiNorm\\Wold\\run intersection nonempty\\raw\\results_0"

if (!exists("base")) {
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  stop("A file root must be supplied",call.=FALSE)
}  
base = args[1]
}

ppi <- 300

filename <- paste(base,"_results.txt",sep="")
if (!file.exists(filename)) {
  cat("Unable to open ",filename,"\n")
} else {
  con <- file(filename, "r")
  params <- as.data.frame(t(read.table(con,skip=1,nrows=10,row.names=1,skipNul = TRUE,sep="\t")))
  close(con)
  for (i in c(3:10)) {params[,i] <- as.numeric(as.character(params[,i]))}
  nulls <- c(0,0,0,0,0,0,0,0)
  params <-  rbind(params[1:2,],c("t1","A",nulls),c("t2","A",nulls),c("a","A",nulls),params[-(1:2),])
  params <-  rbind(params[1:11,],c("a","B",nulls),params[-(1:11),])
  params <-  rbind(params[1:16,],c("t1","C",nulls),params[-(1:16),])
  params <-  rbind(params[1:18,],c("a","C",nulls),params[-(1:18),])
  params <-  rbind(params[1:25,],c("a","D",nulls),params[-(1:25),])
  params <-  rbind(params[1:32,],c("a","E",nulls),params[-(1:32),])

  for (i in c(3:10)) {params[,i] <- as.numeric(as.character(params[,i]))}

  # So that the models are plotted in the right order
  params$Model <- factor(params$Model, levels = c("A","B","C","D","E","BD"))

  # Extract out data for each parameter
  dParams = params[c(1,8,15,22,29,36),]
  hParams = params[c(2,9,16,23,30,37),]
  tParams = params[c(3,4,10,11,17,18,24,25,31,32,38,39),]
  aParams = params[c(5,12,19,26,33,40),]
  LLParams = params[c(6,13,20,27,34,41),]

  # t1 and t2 are premultipled  
  tParams$"Abs Opt" <- tParams$"Abs Opt"*10000
  tParams$"Abs Opt"[c(9,10)] <- tParams$"Abs Opt"[c(9,10)]/100
  tParams$"Spread Abs" <- tParams$"Spread Abs"*10000
  tParams$"Spread Abs"[c(9,10)] <- tParams$"Spread Abs"[c(9,10)]/100

  errSize <- 0.5

  p1 <- ggplot(dParams, aes(x=Name, y=dParams$"Abs Opt", fill=Model)) +  theme_bw() +
    geom_bar(position=position_dodge(width=0.9), stat="identity") +
    geom_errorbar(aes(ymin=dParams$"Abs Opt" - dParams$"Spread Abs",
                      ymax=dParams$"Abs Opt" + dParams$"Spread Abs"), 
                width=.5,position=position_dodge(width=0.9),size=errSize) +
    scale_y_continuous() + guides(fill=FALSE) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
        axis.title.y=element_blank())  

  p2 <- ggplot(hParams, aes(x=Name, y=hParams$"Abs Opt", fill=Model)) +  theme_bw() +
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=hParams$"Abs Opt" - hParams$"Spread Abs",
                      ymax=hParams$"Abs Opt" + hParams$"Spread Abs"), 
                width=.5,position=position_dodge(width=0.9),size=errSize) +
    scale_y_continuous() + guides(fill=FALSE) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
        axis.title.y=element_blank())

  p3 <- ggplot(tParams, aes(x=Name, y=tParams$"Abs Opt", fill=Model)) +  theme_bw() +
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=tParams$"Abs Opt" - tParams[7], ymax=tParams$"Abs Opt" + tParams[7]), 
                width=.5,position=position_dodge(width=0.9),size=errSize) +  
    ylab(bquote('x10'^-4~' Model E: x10'^-2~'')) +
    scale_y_continuous() +guides(fill=FALSE) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank())

  p4 <- ggplot(aParams, aes(x=Name, y=aParams$"Abs Opt", fill=Model)) +  theme_bw() +
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=aParams$"Abs Opt" - aParams$"Spread Abs",
                      ymax=aParams$"Abs Opt" + aParams$"Spread Abs"), 
                width=.5,position=position_dodge(width=0.9),size=errSize) +  
    scale_y_continuous() +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
        axis.title.y=element_blank())

  png(paste(base,"_resultsParam.png",sep=""),units = "in",height=3,width=6,res=ppi)
  grid.arrange(p1,p2,p3,p4, ncol=4,widths = c(8,8,16,15))
# invisible to hide unwanted output to the terminal
  invisible(dev.off())

  png(paste(base,"_resultsLL.png",sep=""),units = "in",height=3,width=4,res=ppi)
  
  LLParams$"Abs Opt" <- LLParams$"Abs Opt"/1000 
  LLParams$"Spread Abs" <- LLParams$"Spread Abs"/1000 
  ymin = floor(min(LLParams$"Abs Opt"))
    # The wonders of R means that ggplot must be inside a print statement if it 
  # is inside an if statement.   Why???
  print(ggplot(LLParams, aes(x=Name, y=LLParams$"Abs Opt", fill=Model)) +  theme_bw() +
    geom_bar(width = 0.6,position=position_dodge(width=1.0), stat="identity") +
    geom_errorbar(aes(ymin=LLParams$"Abs Opt" - LLParams$"Spread Abs",
                      ymax=LLParams$"Abs Opt" + LLParams$"Spread Abs"), 
                  width=.4,position=position_dodge(width=1.0),size=errSize) +  
    scale_y_continuous(limits= c(ymin,NA),oob=rescale_none) +
    ylab(bquote('Negative Log likelihood: x10'^3~'')) +
    theme(axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank()))  
  invisible(dev.off())    
}

filename <-paste(base,"_norm.txt",sep="") 
if (!file.exists(filename)) {
  cat("Unable to open ",filename,"\n")
} else {
  
  #  Read in each of the sets of data for each line, add an identifier and stack them
  #  into a single dataframe using rbind
  con = file(filename, "r")
  models <- as.data.frame(t(read.table(con,skip=2,nrows=2,row.names=1)))
  models <- cbind (models,"A")
  colnames(models)[3] <- "Model";
  model <- as.data.frame(t(read.table(con,skip=3,nrows=2,row.names=1)))
  model <- cbind (model,"B")
  colnames(model)[3] <- "Model";
  models <- rbind(models,model)
  model <- as.data.frame(t(read.table(con,skip=3,nrows=2,row.names=1)))
  model <- cbind (model,"C")
  colnames(model)[3] <- "Model";
  models <- rbind(models,model)
  model <- as.data.frame(t(read.table(con,skip=3,nrows=2,row.names=1)))
  model <- cbind (model,"D")
  colnames(model)[3] <- "Model";
  models <- rbind(models,model)
  model <- as.data.frame(t(read.table(con,skip=3,nrows=2,row.names=1)))
  model <- cbind (model,"E")
  colnames(model)[3] <- "Model";
  models <- rbind(models,model)
  model <- as.data.frame(t(read.table(con,skip=3,nrows=2,row.names=1)))
  model <- cbind (model,"BD")
  colnames(model)[3] <- "Model";
  models <- rbind(models,model)
  close(con)

  png(paste(base,"_norm.png",sep=""),units = "in",height=5,width=6,res=ppi)
  print(qplot(Length,Bias,data=models,colour=Model,geom="line") +theme_bw() +
    scale_x_log10(breaks=c(100,200,400,700,1000,2000,4000,10000,20000)) +
    scale_y_continuous(breaks=c(0:20)/10) +
    xlab ("Length (bp)") + ylab ("Abundance relative to 1 kb"))
    
  invisible(dev.off())
}

filename <- paste(base,"_bias.txt",sep="")

if (!file.exists(filename)) {
  cat("Unable to open ",filename,"\n")
} else {
  geneMatrix <- as.matrix( read.table(filename,sep = "\t"))
  geneMeltData <- melt(geneMatrix)
  my_palette <- colorRampPalette(c("white", "orange","red","black"))(n = 30)
  png(paste(base,"_bias.png",sep=""),units = "in",height=10,width=3,res=ppi)
  print(ggplot(geneMeltData, aes(x = Var2, y = Var1, fill = value)) + 
    coord_fixed(ratio = 0.7) +
    geom_tile() +
    scale_fill_gradientn(colours = my_palette) +
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background = element_blank()) +
    guides(fill=FALSE)) 
  invisible(dev.off())
}







