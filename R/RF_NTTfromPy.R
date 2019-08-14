#note to chanage the options to options("digits"=22) and options("scipen"=100)

library("ape")
library("phangorn")
library("readtext")
options("digits"=22)
options("scipen"=100)
#Beta function
Beta=function(m){
  if(m<=0){ans=1}
  else {
    ans=1
    for (i in 1:m) {
      ans=ans*(2*i+1)
    }
  }
  return(ans)
}

# the function for computing the number of internal edges
internaledges <- function(tree,ntip){
  intedges=array(0,c(1,ntip-1))
  edges=tree$edge
  for (i in (2*ntip-1):(ntip+1)) {
    children=which(edges[,1]==i)
    child1=edges[children[1],2]
    child2=edges[children[2],2]
    if(child1 <= ntip&child2 <= ntip){intedges[i-ntip]=0}
    else if(child1<= ntip & child2 > ntip){intedges[i-ntip]=intedges[child2-ntip]+1}
    else if(child2<= ntip & child1 > ntip){intedges[i-ntip]=intedges[child1-ntip]+1}
    else {intedges[i-ntip]=intedges[child2-ntip]+intedges[child1-ntip]+2}
  }
  return(intedges)
}

# the function for computing the number of internal children
internalchildren <- function(tree,v,ntip){
  edges=tree$edge
  children=which(edges[,1]==v)
  child1=edges[children[1],2]
  child2=edges[children[2],2]
  if(child1 > ntip & child2 > ntip){result=c(2,child1,child2)}
  else if(child1 > ntip & child2 <= ntip){result=c(1,child1)}
  else if(child2 > ntip & child1 <= ntip){result=c(1,child2)}
  else {result=0}
  return(result)
}

RF_Convolve=function(tree,n){
  tt=0
  t=(n-4)*(n-2)*3
  L= 2^ceiling(log(t)/log(2))
  #L=16384
  ntip=n-1
  N=tree$Nnode
  R=rep(list(matrix(0,(ntip-1),(ntip-1))),N)
  edges=internaledges(tree)
  B=c()
  for (k in 0:(n-2)) {
    B[k+1]=Beta(k)
  }
  for (v in N:1) {
    intchild=internalchildren(tree,v+ntip)
    intedges=edges[v]
    if(intchild[1]==0){
      R[[v]][1,1]=1
    }
    else if(intchild[1]==1){
      Rchild=R[[intchild[2]-ntip]]
      R[[v]][1,intedges+1]=1
      R[[v]][2:(ntip-1),1]=rowSums(t(t(Rchild[1:(ntip-2),])*B[1:(ntip-1)]))
      R[[v]][2:(ntip-1),2:(ntip-1)]=Rchild[2:(ntip-1),1:((ntip-2))]
    }
    else {
      Rchild1=R[[intchild[2]-ntip]]
      Rchild2=R[[intchild[3]-ntip]]
      R[[v]][1,intedges+1]=1
      R[[v]][3,1]=sum(t(t(Rchild1[1,])*B[1:(ntip-1)]))*sum(t(t(Rchild2[1,])*B[1:(ntip-1)]))
      for (s in 4:(ntip-1)) {
        R[[v]][s,1]=sum(rowSums(t(t(Rchild1[1:(s-2),])*B[1:(ntip-1)]))*rowSums(t(t(Rchild2[(s-2):1,])*B[1:(ntip-1)])))
      }
      sum1=matrix(0,(ntip-2),(ntip-2))
      sum1[1,1:(ntip-2)]=sum(t(t(Rchild1[1,])*B[1:(ntip-1)]))*Rchild2[1,1:(ntip-2)]
      for (s in 3:(ntip-1)) {
        temp=colSums(rowSums(t(t(Rchild1[1:(s-1),])*B[1:(ntip-1)]))*Rchild2[(s-1):1,1:(ntip-2)])
        sum1[s-1,1:(ntip-2)]=temp
      }
      sum2=matrix(0,(ntip-2),(ntip-2))
      sum2[1,1:(ntip-2)]=sum(t(t(Rchild2[1,])*B[1:(ntip-1)]))*Rchild1[1,1:(ntip-2)]
      for (s in 3:(ntip-1)) {
        temp=colSums(rowSums(t(t(Rchild2[1:(s-1),])*B[1:(ntip-1)]))*Rchild1[(s-1):1,1:(ntip-2)])
        sum2[s-1,1:(ntip-2)]=temp
      }

      R1=Rchild1[1:(ntip-1),1:(ntip-3)]
      #R1aug=numeric(nrow(R1)*(2*ncol(R1)-1))
      R1aug=numeric(L)
      t=1
      for(i in 1:ncol(R1)){
        R1aug[t:(t+nrow(R1)-1)]=R1[,i]
        t=t+3*nrow(R1)
      }
      R2=Rchild2[1:(ntip-1),1:(ntip-3)]
      #R2aug=numeric(nrow(R2)*(2*ncol(R2)-1))
      R2aug=numeric(L)
      t=1
      for(i in 1:ncol(R2)){
        R2aug[t:(t+nrow(R2)-1)]=R2[,i]
        t=t+3*nrow(R2)
      }

      #write(L,"testNTT.txt",ncolumns=L,append = TRUE)
      write(R1aug,"testNTT.txt",ncolumns=L,append = TRUE)
      write(R2aug,"testNTT.txt",ncolumns=L,append = TRUE)

      #read the output of Python

      system('python  /python_code/ntt_fromR.py')

      U=as.matrix(read.csv("outNTT.txt",header = FALSE, quote=""))

      Matc=ceiling(length(R1aug)/(3*nrow(R1)))
      sum3=matrix(c(U,numeric(Matc*3*nrow(R1)-length(R1aug))),nrow=3*nrow(R1))[1:nrow(R1),1:ncol(R1)]
      sum3=cbind(array(0, dim=c(nrow(R1)-1,1)),sum3[2:nrow(R1),])
      R[[v]][2:(ntip-1),2:(ntip-1)]=sum1+sum2+sum3
      file.remove("testNTT.txt")

    }
  }

  return(R)

}



#==========================================
RsT=function(R,n,s){
  B=c()
  for (k in 0:(n-2)) {
    B[k+1]=Beta(k)
  }
  rst =sum(t(t(R[[1]][s+1,1:(n-2-s)])*B[1:(n-2-s)]))
  return(rst)
}

#Compute the value of q_m(T)
qmT=function(R,n,m){
  qmt=0
  for (s in m:(n-3)) {
    rst=RsT(R,n,s)
    qmt=qmt+(factorial(s)/(factorial(m)*factorial(s-m)))*rst*(-1)^(s-m)
  }
  return(qmt)
}

#this function computes the RF distribution
ntt_polynomial=function(tree,n){
  Coef=numeric()
  R=RF_Convolve(tree,n)
  for (i in seq(0,2*(n-3),2)) {
    Coef=c(Coef,qmT(R,n,n-3-(i/2)))
  }
  return(Coef)
}


