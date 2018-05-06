rm(list=ls())

library(ggplot2)
library(reshape2)
library(gridExtra)
library (scales)

##  For debugging, otherwise pass in file root as parameter
## base = "Y:\\LiBiNorm validate\\Combs\\testModel2\\SRR1743145"
if (!exists("base")) {
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  stop("A file root must be supplied",call.=FALSE)
}  
base = args[1]
}

## Change this for better, or worse resolution
ppi <- 300


filename <- paste(base,"_distribution.txt",sep="")

if (!file.exists(filename)) {
  cat("Unable to open ",filename,"\n")
} else {

  con = file(filename, "r")
  info <- read.table(con,skip=0,nrows=1,row.names=1,sep="\t")
  xlab <- "Position along transcript"
  ylab <- "Relative read density"

  if (info[1,1] == 0) 
  {
    dist <- as.data.frame(t(read.table(con,skip=0,nrows=3,row.names=1,sep="\t")))
    dist[, 'Length'] <- as.factor(dist[, 'Length'])
    
    ymin = min(dist$Reads)
    png(paste(base,"_distribution.png",sep=""),units = "in",height=5,width=6,res=ppi)
    print(ggplot(dist,aes(Position,Reads,colour=Length)) + theme_bw() + geom_line() + 
            theme(axis.text = element_text(colour = "black")) +
            labs(x = xlab,y=ylab) +
            annotate("text",0.5,ymin+0.05,label = "Read distribution")) 
    invisible(dev.off())
  }  else   {
    dist <- as.data.frame(t(read.table(con,skip=0,nrows=4,row.names=1,sep="\t")))
    dist[, 'Length'] <- as.factor(dist[, 'Length'])
  
    ymin = min(min(dist$Reads),min(dist[4]))
    png(paste(base,"_distribution.png",sep=""),units = "in",height=5,width=6,res=ppi)
    print(ggplot(dist,aes(Position,Reads,colour=Length)) + theme_bw() + geom_line() + 
            theme(axis.text = element_text(colour = "black")) +
            geom_line(aes(Position,dist[4],colour=Length),linetype="dashed") +
            labs(x = xlab,y=ylab) +
            annotate("text",0.5,ymin+0.05,label = colnames(dist)[4]))

    invisible(dev.off())

    if (info[1,1] > 1) 
    {
      dist2 <- as.data.frame(t(read.table(con,skip=0,nrows=6,row.names=1,sep="\t")))
      plots <- list()

      for (i in 1:6)
      local({
        i <- i
        ymin = min(min(dist$Reads),min(dist2[i]))
        p1 <- ggplot(dist,aes(Position,Reads,colour=Length)) + theme_bw() + geom_line(linetype="dashed") + 
        geom_line(aes(Position,dist2[i],colour=Length))  +scale_colour_discrete(guide = FALSE) + 
          annotate("text",0.5,ymin+0.05,label = colnames(dist2)[i])
        plots[[i]] <<- p1
      })
      png(paste(base,"_distAll.png",sep=""),units = "in",height=6,width=9,res=ppi)
      grid.arrange(grobs = plots, ncol=3,nrow=2)

      invisible(dev.off())
      }
    }
    close(con)
 }



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
    theme(axis.text = element_text(colour = "black")) +
    geom_bar(position=position_dodge(width=0.9), stat="identity") +
    geom_errorbar(aes(ymin=dParams$"Abs Opt" - dParams$"Spread Abs",
                      ymax=dParams$"Abs Opt" + dParams$"Spread Abs"), 
                width=.5,position=position_dodge(width=0.9),size=errSize) +
    scale_y_continuous() + guides(fill=FALSE) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
        axis.title.y=element_blank())  

  p2 <- ggplot(hParams, aes(x=Name, y=hParams$"Abs Opt", fill=Model)) +  theme_bw() +
    theme(axis.text = element_text(colour = "black")) +
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=hParams$"Abs Opt" - hParams$"Spread Abs",
                      ymax=hParams$"Abs Opt" + hParams$"Spread Abs"), 
                width=.5,position=position_dodge(width=0.9),size=errSize) +
    scale_y_continuous() + guides(fill=FALSE) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
        axis.title.y=element_blank())

  p3 <- ggplot(tParams, aes(x=Name, y=tParams$"Abs Opt", fill=Model)) +  theme_bw() +
    theme(axis.text = element_text(colour = "black")) +
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
        axis.title.y=element_blank(),
        axis.text = element_text(colour = "black"))

  png(paste(base,"_resultsParam.png",sep=""),units = "in",height=3,width=6,res=ppi)
  grid.arrange(p1,p2,p3,p4, ncol=4,widths = c(8,8,16,15))
# invisible to hide unwanted output to the terminal
  invisible(dev.off())

  png(paste(base,"_resultsLL.png",sep=""),units = "in",height=3,width=4,res=ppi)
  
  LLParams$"Abs Opt" <- LLParams$"Abs Opt"/1000 
  LLParams$"Spread Abs" <- LLParams$"Spread Abs"/1000 
  ymin = floor(min(LLParams$"Abs Opt")*10)/10
# The wonders of R means that ggplot must be inside a print statement if it 
# is inside an if statement.   Why???
  print(ggplot(LLParams, aes(x=Name, y=LLParams$"Abs Opt", fill=Model)) +  theme_bw() +
    theme(axis.text = element_text(colour = "black")) +
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
  
  models$Length <- models$Length/1000

  png(paste(base,"_norm.png",sep=""),units = "in",height=5,width=6,res=ppi)
  print(qplot(Length,Bias,data=models,colour=Model,geom="line") +theme_bw() +
  scale_x_log10(breaks=c(0.1,0.2,0.4,1,2,4,10,20)) + 
    theme(axis.text = element_text(colour = "black")) +
    scale_y_continuous(breaks=c(0:20)/10) +
    xlab ("Length (kbp)") + ylab ("Bias relative to 1 kbp"))

  invisible(dev.off())
}

filename <- paste(base,"_bias.txt",sep="")

if (!file.exists(filename)) {
  cat("Unable to open ",filename,"\n")
} else {
  geneMatrix <- as.matrix( read.table(filename,sep = "\t"))
  geneMatrix <- geneMatrix[,-1]
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







