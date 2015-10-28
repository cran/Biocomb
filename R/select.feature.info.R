compute.aucs <- function(dattable){
  labs <- dattable[,ncol(dattable)]
  aucvals <- rep(0,ncol(dattable)-1)

  val=levels(labs)
  pos <- rep(val[2],ncol(dattable)-1)
  if(length(val)==2)
  {
    aucvals <- rep(0,ncol(dattable)-1)
    for (i in 1:(ncol(dattable)-1)){
      pred <- prediction(dattable[,i],labs)
      aucv <- performance(pred,"tpr", "fpr",measure="auc")
      aucval <- attr(aucv,"y.values")[[1]]
      if (aucval<0.5){
        aucval <- 1-aucval
        pos[i] <- val[1]  ####Positive Correlation ist richtig, AUC-Wert >=0.5, sonst AUC-Wert <0.5
      }
      aucvals[i] <- aucval
    }
    auctab<-data.frame(names(dattable)[1:c(ncol(dattable)-1)],aucvals,pos)
    names(auctab)<-c("Biomarker","AUC","Positive class")
  }
  else
  {
  for (i in 1:(ncol(dattable)-1)){
    aucval <- multiclass.roc(labs,dattable[,i])$auc
    aucval2 <- multiclass.roc(labs,dattable[,i])$auc
    if (aucval<aucval2){
      aucval <- aucval2
    }
    aucvals[i] <- aucval
  }
    auctab<-data.frame(names(dattable)[1:c(ncol(dattable)-1)],aucvals)
    names(auctab)<-c("Biomarker","AUC")
  }

  return(auctab)
}

fun1.chi<-function(data,class)
{
  int.list=NULL
  d1=dim(data)
  for(i in 1:d1[2])
  {
    alist=NULL
    out=sort(data[,i],index.return=TRUE)
    cl=class[out$ix]
    out.unique=unique(out$x)
    for(j in 1:length(out.unique))
    {
      index=which(out$x==out.unique[j])
      alist=c(alist,list(list(var=out$ix[index],cl=cl[index])))
    }
    int.list=c(int.list,list(alist))
  }

  int.list
}

fun2.chi<-function(int.list,mat.int)
{
  dm=dim(mat.int)
  till=length(int.list)

  chi.stat=NULL

  for(i in 1:till)
  {
    chi.vrem=numeric()

    vrem=int.list[[i]]
    for(j in 1:(length(vrem)-1))
    {
      for(ij in 1:2)
      {
        if((j>1)&&(ij==1))
        {
          mat.int[ij,]=mat.int[ij+1,]
        }
        else
        {
        res=table(vrem[[j+(ij-1)]]$cl)
        mat.int[ij,]=rep(0,dm[2])
        for(kj in 1:length(res))
        {
          mat.int[ij,names(res)[kj]]=res[kj]
        }
        }
      }
      chi.vrem=c(chi.vrem,check.stat(mat.int))
    }
    chi.stat=c(chi.stat,list(chi.vrem))
  }
  chi.stat
}

check.stat<-function(mat.int)
{
  dm=dim(mat.int)

  row.sum=rowSums(mat.int)
  col.sum=colSums(mat.int)
  res=row.sum%*%t(col.sum)
  res=res/sum(mat.int)

  index=which(res==0)
  res[index]=0.1

  out=((mat.int-res)**2)/res
  result=sum(out)
  result
}


check.incons<-function(data,vrem.nominal,class)
{
  data.check=cbind(data,vrem.nominal)
  data.check=t(data.check)
  input=dim(data)[1]

  incons=0
  while(input>1)
  {
    vrem=data.check-data.check[,1]
    out=sapply(1:dim(vrem)[2], function(z) all(vrem[,z]==0))
    if(length(out[out])>1)
    {
      e.max=max(table(class[out]))
      incons=incons+(length(out[out])-e.max)
    }
    class=class[!out]
    input=input-length(out[out])
    data.check=data.check[,!out,drop=FALSE]
  }
  incons=incons/dim(data)[1]

  incons
}

chi2.algorithm<- function(matrix,attrs.nominal,threshold)
{
  dd=dim(matrix)
  
  if(length(attrs.nominal)>0)
  {
  for(i in 1:length(attrs.nominal))
  {
    matrix[,attrs.nominal[i]]=as.factor(matrix[,attrs.nominal[i]])
  }
  }

  #for inconsistency for nominal
  vrem.nominal=matrix[,attrs.nominal,drop=FALSE]
  if(length(attrs.nominal)>0)
  {
  for(i in 1:length(attrs.nominal))
  {
    vrem.nominal[,i]=as.numeric(vrem.nominal[,i])
  }
  }
  #-------

  data=matrix[,-c(attrs.nominal,dd[2]),drop=FALSE]
  data.start=data
  class=matrix[,dd[2]]
  class=as.character(class)

  d1=dim(data)

  label=unique(class)
  mat.int=matrix(0,2,length(label))
  colnames(mat.int)=label


  int.list=fun1.chi(data,class)
  int.list.start=int.list

  #Phase 1
  sig.value=0.6
  df=length(label)-1
  chi.value=qchisq(1-sig.value, df=df)

  chi.stat=fun2.chi(int.list,mat.int)
  chi.stat.start=chi.stat

  #threshold=(d1[1]-max(table(class))-1)/d1[1]
  incons=0

  #flag.end=FALSE
  step=0.1
  delta=6
  shag=1

  while(incons<=threshold)
  {
    calc=0

    #if(flag.end) break


    sig.value0=sig.value

    if(shag==delta)
    {
      step=step*0.1
      delta=delta+9
    }

    sig.value=sig.value-step
    shag=shag+1

    if(sig.value<0.000000000011)
    {
      #browser()
    }
    chi.value=qchisq(1-sig.value, df=df)


    check=sapply(chi.stat,function(z) length(z))
    if(all(check==0))
    {
      break
    }


    for(irow in 1:d1[2])
    {
      while(TRUE)
      {
        #select the min element of chi.stat
        check=length(chi.stat[[irow]])
        if(check==0)
        {
          break
        }


    #cat(paste("Value",irow,"Value",icol,"\n"))
        icol=which.min(chi.stat[[irow]]) #which interval
        if(chi.stat[[irow]][icol]>chi.value)
        {
          break
        }
    vrem=int.list[[irow]]
    #cat(paste("Value",length(vrem),"\n"))
    vrem[[icol]]$var=c(vrem[[icol]]$var,vrem[[icol+1]]$var)
    vrem[[icol]]$cl=c(vrem[[icol]]$cl,vrem[[icol+1]]$cl)

    if(icol!=(length(vrem)-1))
    {
    for( i in (icol+1):(length(vrem)-1))
    {
      vrem[[i]]=vrem[[i+1]]
      if(i==(length(vrem)-1)) break
      chi.stat[[irow]][i]=chi.stat[[irow]][i+1]
    }
    }

    vrem[[length(vrem)]]=NULL
    chi.stat[[irow]]=chi.stat[[irow]][-length(vrem)]

    #new chi values

    for(j in (icol-1):icol)
    {
      if((j>0)&&(j<length(vrem)))
      {
        for(ij in 1:2)
        {
          res=table(vrem[[j+(ij-1)]]$cl)
          mat.int[ij,]=rep(0,length(label))
          for(kj in 1:length(res))
          {
            mat.int[ij,names(res)[kj]]=res[kj]
          }
        }
        chi.stat[[irow]][j]=check.stat(mat.int)
      }
    }


    data[vrem[[icol]]$var,irow]=data[vrem[[icol]]$var[1],irow]

    int.list[[irow]]=vrem

    calc=calc+1
    #cat(paste("Iteration",calc,"\n"))
      }
    }



    incons=check.incons(data,vrem.nominal,class)


  }

  #Phase 2
  #data=matrix[,-c(attrs.nominal,dd[2])]
  #class=matrix[,dd[2]]
  #d1=dim(data)
  data=data.start

  sig.attr=rep(sig.value0,d1[2])
  chi.value=qchisq(1-sig.value0, df=df)
  chi.attr=rep(chi.value,d1[2])

  #int.list=fun1.chi(data,class)
  int.list=int.list.start
  #the initial chi.stat
  #chi.stat=fun2.chi(int.list,mat.int)
  chi.stat=chi.stat.start

  flag=rep(TRUE,d1[2])
  #further
  calc=0
  while(TRUE)
  {
    if(all(!flag)) break
  for(irow in 1:d1[2])
  {
    if(flag[irow])
    {
      int.list1=int.list[[irow]]
      chi.stat1=chi.stat[[irow]]
      data1=data[,irow]
  while(TRUE)
  {
    #select the min element of chi.stat
    check=length(chi.stat[[irow]])
    if(check==0)
    {
      flag[irow]=FALSE
      break
    }
    icol=which.min(chi.stat[[irow]]) #which interval
    if(chi.stat[[irow]][icol]>chi.attr[irow])
    {
      break
    }
    vrem=int.list[[irow]]
    vrem[[icol]]$var=c(vrem[[icol]]$var,vrem[[icol+1]]$var)
    vrem[[icol]]$cl=c(vrem[[icol]]$cl,vrem[[icol+1]]$cl)

    if(icol!=(length(vrem)-1))
    {
      for( i in (icol+1):(length(vrem)-1))
      {
        vrem[[i]]=vrem[[i+1]]
        if(i==(length(vrem)-1)) break
        chi.stat[[irow]][i]=chi.stat[[irow]][i+1]
      }
    }

    vrem[[length(vrem)]]=NULL
    chi.stat[[irow]]=chi.stat[[irow]][-length(vrem)]

    #new chi values
    for(j in (icol-1):icol)
    {
      if((j>0)&&(j<length(vrem)))
      {
        for(ij in 1:2)
        {
          res=table(vrem[[j+(ij-1)]]$cl)
          mat.int[ij,]=rep(0,length(label))
          for(kj in 1:length(res))
          {
            mat.int[ij,names(res)[kj]]=res[kj]
          }
        }
        chi.stat[[irow]][j]=check.stat(mat.int)
      }
    }

    data[vrem[[icol]]$var,irow]=data[vrem[[icol]]$var[1],irow]

    int.list[[irow]]=vrem

    calc=calc+1
    #cat(paste("Iteration",calc,"\n"))

  }
  #check.incons
    incons=check.incons(data,vrem.nominal,class)

    if(incons<=threshold)
    {
      sig.attr[irow]=sig.attr[irow]-step


      chi.attr[irow]=qchisq(1-sig.attr[irow], df=df)
    }
    else
    {
      int.list[[irow]]=int.list1
      chi.stat[[irow]]=chi.stat1
      data[,irow]=data1
      flag[irow]=FALSE
    }
    }
  }
    if(shag==delta)
    {
      step=step*0.1
      delta=delta+9
    }
    shag=shag+1
  }
  rr=sapply(1:d1[2],function(z) length(unique(data[,z])))
  data.out=data[,which(rr>1),drop=FALSE]
  data.out=cbind(data.out,matrix[,attrs.nominal,drop=FALSE])
  
  return(list(data.out=data.out,subset=colnames(data.out)))
}


select.forward.Corr<- function(matrix,disc.method,attrs.nominal)
{
  evaluator <- function(subset) {
    #correaltion evaluate
    sum=0
    out=0
    len=length(subset)
    for(i in 1:len)
    {
      sum=sum+CalcGain(m3[,subset[i]],m3[,ncol(m3)],TRUE)
      if(len==1)
      {
        out=0
      }
      else
      {
      vrem=sapply(c(i:len),function(z) CalcGain(m3[,subset[z]],m3[,subset[i]],TRUE) )

      out=out+sum(vrem[-1])
      }
    }
    out=2*out #perhaps without len
    CFS=sum/sqrt(out+len) #according to the formula
    #cat(paste(subset,"\n"))
    return(CFS)
  }

  out=ProcessData(matrix,disc.method,attrs.nominal,FALSE)
  m3=out$m3
  sel.feature=out$sel.feature
  dd=dim(m3)
  if(dd[2]>1)
  {
  subset <- forward.search(names(m3)[-ncol(m3)], evaluator)
  }
  else
  {
    subset <- NULL
  }
}

select.forward.wrapper<- function(dattable)
{
  evaluator <- function(subset) {
    #k-fold cross validation
    results = sapply(1:k, function(i) {
      #test.idx <- (splits >= (i - 1) / k) & (splits < i / k)
      test.idx <- testind[,i]
      train.idx <- !test.idx
      test <- dattable[test.idx, , drop=FALSE]
      train <- dattable[train.idx, , drop=FALSE]
      tree <- rpart(as.simple.formula(subset, names(dattable)[ncol(dattable)]), train,method = "class")
      error.rate = sum(test[,ncol(dattable)] != predict(tree, test, type="c")) / nrow(test)
      return(1 - error.rate)
      #val=FUN(subset,dattable,train,test)
    })
    return(mean(results))
  }

  k <- 5
  splits <- runif(nrow(dattable))
  testind <- sapply(1:k, function(i) {(splits >= (i - 1) / k) & (splits < i / k)})
  subset <- forward.search(names(dattable)[-ncol(dattable)], evaluator)
}


CalcGain<-function(m1,m2,symm)
{
  dd=length(m1)
  fq1=table(m1)
  fq1=fq1/dd[1]

  entropyF1=-sapply(fq1, function(z) if(z==0) 0 else z*log(z))
  entropyF1=sum(entropyF1)

  fq2=table(m2)
  fq2=fq2/dd[1]

  entropyF2=-sapply(fq2, function(z) if(z==0) 0 else z*log(z))
  entropyF2=sum(entropyF2)

  fq=table(m1,m2)

  entropyF12=0

  for(i in 1:length(fq2))
  {
    fq0=fq[,i]/sum(fq[,i])
    vrem=-sapply(fq0,function(z) if(z==0) 0 else z*log(z))
    entropyF12=entropyF12+(fq2[i])*sum(vrem)
  }

  entropy=entropyF1-entropyF12
  if(symm)
  {
    if((entropyF1+entropyF2)==0)
    {
      entropy=0
    }
    else
    {
      entropy=2*entropy/(entropyF1+entropyF2)
    }
  }
  return(entropy)
}


ProcessData1<-function(matrix,disc.method,attrs.nominal)
{
  dd=dim(matrix)
  matrix=data.frame(matrix)

  matrix[,dd[2]]=as.factor(matrix[,dd[2]])
  #data=matrix

  if(disc.method=="MDL")
  {
    m3 <- Discretize(as.formula(paste(names(matrix)[dd[2]],"~.")), data = matrix)
    #m3<-mdlp(matrix)$Disc.data
  }

  if(disc.method=="equal frequency")
  {
    m3=matrix
    for(i in 1:(dd[2]-1))
    {
      if(!(i%in%attrs.nominal))
      {
        m3[,i] <- discretize(matrix[,i], method="frequency",categories=3)
      }
    }
  }

  if(disc.method=="equal interval width")
  {
    m3=matrix
    for(i in 1:(dd[2]-1))
    {
      if(!(i%in%attrs.nominal))
      {
        m3[,i] <- discretize(matrix[,i], categories=3)
      }
    }
  }
  #-------------------

  #extract the features with one interval

  sel.one=lapply(m3, function(z) (length(levels(z))==1)&&(levels(z)=="'All'"))

  sel.one=which(unlist(sel.one)==TRUE)

  #selected features
  sel.feature=1:dd[2]
  if(length(sel.one)>0)
  {
    sel.feature=sel.feature[-sel.one]

    matrix=matrix[,-sel.one,drop=FALSE]
    m3=m3[,-sel.one,drop=FALSE]
  }
  return (list(m3=m3,sel.feature=sel.feature))
}

ProcessData<-function(matrix,disc.method,attrs.nominal,flag=FALSE)
{
  dd=dim(matrix)
  matrix=data.frame(matrix)

  matrix[,dd[2]]=as.factor(matrix[,dd[2]])
  #data=matrix

  if(disc.method=="MDL")
  {
    m3 <- Discretize(as.formula(paste(names(matrix)[dd[2]],"~.")), data = matrix)
    #m3<-mdlp(matrix)$Disc.data
  }

  if(disc.method=="equal frequency")
  {
    m3=matrix
    for(i in 1:(dd[2]-1))
    {
      if(!(i%in%attrs.nominal))
      {
        m3[,i] <- discretize(matrix[,i], method="frequency",categories=3)
      }
    }
  }

  if(disc.method=="equal interval width")
  {
    m3=matrix
    for(i in 1:(dd[2]-1))
    {
      if(!(i%in%attrs.nominal))
      {
        m3[,i] <- discretize(matrix[,i], categories=3)
      }
    }
  }
  #-------------------
  sel.feature=1:dd[2]
  if(flag)
  {
  #extract the features with one interval

   sel.one=lapply(m3, function(z) (length(levels(z))==1)&&(levels(z)=="'All'"))

   sel.one=which(unlist(sel.one)==TRUE)

  #selected features

  if(length(sel.one)>0)
  {
    sel.feature=sel.feature[-sel.one]

    matrix=matrix[,-sel.one,drop=FALSE]
    m3=m3[,-sel.one,drop=FALSE]
  }
  }
  return (list(m3=m3,sel.feature=sel.feature))
}

select.inf.chi2<-function(matrix,disc.method,attrs.nominal)
{
  out=ProcessData(matrix,disc.method,attrs.nominal,FALSE)
  m3=out$m3
  sel.feature=out$sel.feature
  #algorithm
  dd=dim(m3)
  if(dd[2]>1)
  {
  #stat=sapply(1:(dd[2]-1), function(z) chisq.test(table(m3[,z],m3[,dd[2]]))$statistic)
  #names(stat)=colnames(m3[,1:(dd[2]-1)])
  #to compare
  weights <- chi.squared(as.formula(paste(names(m3)[dd[2]],"~.")), m3)

  #what features are selected
  res=sort(weights$attr_importance,decreasing = TRUE,index.return=TRUE)
  val=res$ix
  weights.sort=res$x
  num.feature=sel.feature[val] #val - sorting
  info=data.frame(names(m3)[val],weights.sort,num.feature)
  }
  else
  {
    info=data.frame(character(),numeric(),numeric())
  }
  names(info) <- c("Biomarker","ChiSquare","NumberFeature")
  return(info)
}

select.inf.symm<-function(matrix,disc.method,attrs.nominal)
{
  out=ProcessData(matrix,disc.method,attrs.nominal,FALSE)
  m3=out$m3
  sel.feature=out$sel.feature
  #algorithm
  dd=dim(m3)

  if(dd[2]>1)
  {
  #SU1=information.gain(names(matrix)[dd[2]]~., matrix) #package "FSelect"

  #entropy of feature
  entropy=c()
  class=m3[,dd[2]]

  for(j in 1:(dd[2]-1))
  {
    feature=m3[,j]

    #Function
    out=CalcGain(feature,class,TRUE)
    entropy=c(entropy,out)
    #--------
  }
  #what features are selected
  res=sort(entropy,decreasing = TRUE,index.return=TRUE)
  val=res$ix
  entropy.sort=res$x
  num.feature=sel.feature[val] #val - sorting
  info=data.frame(names(m3)[val],entropy.sort,num.feature)
  }
  else
  {
    info=data.frame(character(),numeric(),numeric())
  }
  names(info) <- c("Biomarker","SymmetricalUncertainty","NumberFeature")
  return(info)
}

select.inf.gain<-function(matrix,disc.method,attrs.nominal)
{
  out=ProcessData(matrix,disc.method,attrs.nominal,FALSE)
  m3=out$m3
  sel.feature=out$sel.feature
  #algorithm
  dd=dim(m3)
  if(dd[2]>1)
  {
  #SU1=information.gain(names(matrix)[dd[2]]~., matrix) #package "FSelect"

  #entropy of feature
  entropy=c()
  class=m3[,dd[2]]

  for(j in 1:(dd[2]-1))
  {
    feature=m3[,j]

    #Function
    out=CalcGain(feature,class,FALSE)
    entropy=c(entropy,out)
    #--------
  }
  #what features are selected
  res=sort(entropy,decreasing = TRUE,index.return=TRUE)
  val=res$ix
  entropy.sort=res$x
  num.feature=sel.feature[val] #val - sorting
  info=data.frame(names(m3)[val],entropy.sort,num.feature)
  }
  else
  {
    info=data.frame(character(),numeric(),numeric())
  }
  names(info) <- c("Biomarker","Information.Gain","NumberFeature")
  return(info)
}

select.fast.filter<-function(matrix,disc.method,threshold,attrs.nominal)
{

  #second package "RWeka"

  out=ProcessData(matrix,disc.method,attrs.nominal,FALSE)
  m3=out$m3
  sel.feature=out$sel.feature
  #algorithm
  dd=dim(m3)
  if(dd[2]>1)
  {
  #SU1=information.gain(names(matrix)[dd[2]]~., matrix) #package "FSelect"

  #entropy of feature
  entropy=c()
  class=m3[,dd[2]]

  for(j in 1:(dd[2]-1))
  {
  feature=m3[,j]

  #Function
  out=CalcGain(feature,class,FALSE)
  entropy=c(entropy,out)
  #--------
  }

  ind=sapply(entropy,function(z) z>=threshold)

  entropy=entropy[ind]
  m3=m3[,ind,drop=FALSE]

  index.F1=1

  res=sort(entropy,decreasing = TRUE,index.return=TRUE)
  val=res$ix
  entropy.sort=res$x

  while(index.F1<=length(val))
  {
    Fp=m3[,val[index.F1]]
    index.F2=index.F1+1
    while(index.F2<=length(val))
    {
    Fq=m3[,val[index.F2]]
    SUpq=CalcGain(Fp,Fq,FALSE)
    if(SUpq>=entropy.sort[index.F2])
    {
      val=val[-index.F2]
      entropy.sort=entropy.sort[-index.F2]
      index.F2=index.F2-1
    }
    index.F2=index.F2+1
    }
    index.F1=index.F1+1
  }

  #what features are selected, ind-features with SU(p,c)>threshold
  num.feature=sel.feature[ind]
  num.feature=num.feature[val] #val - sorting
  info=data.frame(names(m3)[val],entropy.sort,num.feature)
  }
  else
  {
    info=data.frame(character(),numeric(),numeric())
  }
  names(info) <- c("Biomarker","Information.Gain","NumberFeature")
  return(info)
}
