wd <- "/Users/SJC/Documents/practice/internship"
datafile <- "ova_db.csv" 
varfile <- "ova_variable.csv" 

setwd(wd) 


ova.db.csv=as.data.frame(read.csv(datafile,skip=2,stringsAsFactors=F,check.names=F)) 
NAval <- sapply(strsplit(colnames(ova.db.csv),"\\n"), 
                function(v) { q  <- trimws(v[-1]) 
                qq <- grep("unknown", q, ignore.case=T) 
                ifelse(length(qq),strsplit(q[qq[1]], "\\.")[[1]][1],NA) } ) 
for (i in 1:ncol(ova.db.csv)) if (!is.na(NAval[i])) ova.db.csv[which(NAval[i]==ova.db.csv[,i]),i] <- NA 
colnames(ova.db.csv)=gsub(" ", "_", sapply(strsplit(colnames(ova.db.csv),"\\n"),function(v)trimws(v[1]))) 


idx <- which(apply(ova.db.csv, 2, function(v)length(table(v)) == 1)) 
ova.db.csv <- ova.db.csv[,-idx] 
idx <- which(apply(ova.db.csv, 2, function(v)sum(is.na(v)))/nrow(ova.db.csv) > 0.2) 
ova.db.csv <- ova.db.csv[,-idx] 
options(warn=2); 
idx <- which(apply(ova.db.csv, 2, function(v)class(try(as.numeric(v), T))=='try-error')) 
options(warn=1); 
ova.db.csv <- ova.db.csv[,-idx] 


ova.var.csv=read.csv(varfile) 
vars <- apply(ova.var.csv[seq(2,18,4),], 1, function(v)sapply(strsplit(v,"\\n"),function(vv)trimws(vv[1])))[-1,] 
rownames(vars) <- NULL 
ova.db.csv <- as.data.frame(apply(ova.db.csv, 2, as.numeric)) 




ova.db.csv$Recurrence[ova.db.csv$Recurrence==2] <- 1 






cols <- list() 
library(MASS) 
library(ROCR) 
for (i in 1:ncol(vars)) { 
  cols[[i]] <- na.omit(colnames(ova.db.csv)[match(vars[,i], colnames(ova.db.csv))]) 
  v <- cols[[i]] 
  ova.db.csv[,v] 
} 
vars[which(is.na(match(tolower(vars[,1]), tolower(colnames(ova.db.csv))))),1] 


tbl <- apply(ova.db.csv[,cols[[1]]], 2, table) 
newcols <- names(which(unlist(lapply(tbl, function(v) { 
  length(v) != 2 | min(v[1]/v[2], v[2]/v[1]) >= 0.05 
})))) 
newcols2 <- names(which((colSums(is.na(ova.db.csv[,newcols])) / nrow(ova.db.csv)) == 0)) 


if (0) { 
  vv <- cor(ova.db.csv[,newcols2], use="pair") 
  vv[abs(vv)<0.5] <- 0 
  vx <- vv[colSums(abs(vv))>1,colSums(abs(vv))>1] 
  
  
  image(vx, xaxt='n', yaxt='n', zlim=c(-1,1)) 
  box() 
  axis(1, seq(0, 1, length.out=nrow(vx)), labels=rep("", nrow(vx))) 
  text(seq(0, 1, length.out=nrow(vx)), par("usr")[1]-0.05, srt=-45, adj=0, labels=colnames(vx), xpd=T) 
  
  
} 


#colSums(is.na(ova.db.csv[,newcols2])) / nrow(ova.db.csv) 






XX <- c("Recurrence", "Platinum_resistance_6mo") 
X <- 1 


ova.db.csv.sub <- ova.db.csv[!is.na(ova.db.csv[XX[X]]),] 
newcols1 <- na.omit(colnames(ova.db.csv.sub)[match(tolower(vars[,X]), 
                                                   tolower(colnames(ova.db.csv.sub)))]) 


mod <- glm(paste(XX[X], "~",paste(newcols1,collapse="+")), data=ova.db.csv.sub, family="binomial") 
nul <- glm(paste(XX[X], "~1"), data=ova.db.csv.sub) 
step.0 <- stepAIC(nul, scope=list(lower=nul,upper=mod), direction="both") 
step.0x <- step(nul, scope=list(lower=nul,upper=mod), direction="both") 
step.f <- stepAIC(mod, direction="both") 
step.fx <- step(mod, direction="both") 
step.0$anova 
step.f$anova 


# Picked 
var.full <- setdiff(newcols2, substr(as.character(step.f$anova$Step), 3, 1000)) 
mean(do.cv(var.full, 3)$auc) 
var.full <- substr(as.character(step.0$anova$Step), 3, 1000)[-1] 
mean(do.cv(var.full, 3)$auc) 


var.full <- setdiff(newcols2, substr(as.character(step.fx$anova$Step), 3, 1000)) 
mean(do.cv(var.full, 3)$auc) 
var.full <- substr(as.character(step.0x$anova$Step), 3, 1000)[-1] 
mean(do.cv(var.full, 3)$auc) 




# Evaluate 
n <- nrow(ova.db.csv.sub) 
do.cv <- function(var, ncv) { 
  ret <- c() 
  for (i in 1:100) { 
    idx <- sample(1:n)[1:round(n/ncv)] 
    curmod <- glm(paste(XX[X], "~",paste(var,collapse="+")), data=ova.db.csv.sub[-idx,], family="binomial") 
    curpr <- predict(curmod, ova.db.csv.sub[idx,], type="response") 
    curpr <- prediction(curpr, ova.db.csv.sub[idx,XX[X]]) 
    prof <- performance(curpr, measure="tpr", x.measure="fpr") 
    auc <- performance(curpr, measure = "auc") 
    auc <- auc@y.values[[1]] 
    ret <- c(ret, auc) 
  } 
  list( 
    v = prof, 
    auc = ret 
  ) 
} 


rex <- NULL 
for (i in newcols1) { 
  ret <- do.cv(i, 3)$auc 
  if (is.null(rex)) rex <- ret 
  else rex <- rbind(rex, ret) 
} 
rownames(rex) <- newcols1 
req <- t(rex)[,which(colMeans(t(rex))>0.7)] 
req <- req[,order(colMeans(req), decreasing=T)] 
xxx <- xx <- par("mai") 
xxx[1] <- par("mai")[1] + 0.5 
xxx[2] <- par("mai")[2] + 1 
par(mai=xxx) 
barplot(colMeans(req), ylim=c(0,1), xaxt='n', xpd=F, space=1, ylab="AUC", 
        main="Individual variable prediction performance w/ AUC<0.6");box(); 
grid(NA, 5, lwd=2) 
par(new=T) 
barplot(colMeans(req), ylim=c(0,1), xaxt='n', xpd=F, space=1, ylab="AUC", 
        main="Individual variable prediction performance w/ AUC<0.6");box(); 
end_point = 0.5 + ncol(req) + ncol(req) - 1 
text(seq(1.5,end_point,by=2), par("usr")[3]-0.03,  
     srt = 25, adj= 1, xpd = TRUE, 
     labels = paste(colnames(req)), cex=.8) 
for (i in seq(1.5,end_point,by=2)) { 
  J <- i/2+0.25 
  y <- c(mean(req[,J])+sd(req[,J]),mean(req[,J])-sd(req[,J])) 
  lines(x=c(i,i), y=y) 
  lines(x=c(i-0.2,i+0.2), y=c(y[1],y[1])) 
  lines(x=c(i-0.2,i+0.2), y=c(y[2],y[2])) 
} 
par(mai=xx) 




vf <- do.cv(var.full, 3)$auc 
vf 


plot(prof) 
text(0.8, 0.1, paste("AUC =", round(auc, 3))) 
} 
vars[,1] 
colnames(ova.db.csv) 