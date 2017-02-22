#Using fft to compute the convolution
RConvolve=function(tree,n){
  ntip=n-1
  N=tree$Nnode
  R=rep(list(matrix(0,(ntip-1),(ntip-1))),N)
  edges=internaledges(tree,ntip)
  B=c()
  for (k in 0:(ntip-2)) {
    B[k+1]=Beta(k)        
  }
  for (v in N:1) {
    intchild=internalchildren(tree,v+ntip,ntip)
    intedges=edges[v]
    if(intchild[1]==0){
      R[[v]][1,1]=1
    }
    else if(intchild[1]==1){
      Rchild=R[[intchild[2]-ntip]]
      R[[v]][1,intedges+1]=Beta(intedges)
      R[[v]][2:(ntip-1),1]=rowSums(Rchild[1:(ntip-2),])
      R[[v]][2:(ntip-1),2:(ntip-1)]=t(t(Rchild[2:(ntip-1),1:((ntip-2))])*seq(3,(2*ntip-3),2))
      }
    else {
      Rchild1=R[[intchild[2]-ntip]]
      Rchild2=R[[intchild[3]-ntip]]
      R[[v]][1,intedges+1]=Beta(intedges)
      R[[v]][3,1]=sum(Rchild1[1,])*sum(Rchild2[1,])
      for (s in 4:(ntip-1)) {
        R[[v]][s,1]=sum(rowSums(Rchild1[1:(s-2),])*rowSums(Rchild2[(s-2):1,]))
        }
      sum1=matrix(0,(ntip-2),(ntip-2))
      sum1[1,1:(ntip-2)]=t(t(rowSums(t(Rchild1[1,]))*Rchild2[1,1:(ntip-2)])*seq(3,(2*ntip-3),2))
      for (s in 3:(ntip-1)) {
        temp=colSums(t(t(rowSums(Rchild1[1:(s-1),])*Rchild2[(s-1):1,1:(ntip-2)])*seq(3,(2*ntip-3),2)))
        sum1[s-1,1:(ntip-2)]=temp
       }
      sum2=matrix(0,(ntip-2),(ntip-2))
      sum2[1,1:(ntip-2)]=t(t(rowSums(t(Rchild2[1,]))*Rchild1[1,1:(ntip-2)])*seq(3,(2*ntip-3),2))
      for (s in 3:(ntip-1)) {
        temp=colSums(t(t(rowSums(Rchild2[1:(s-1),])*Rchild1[(s-1):1,1:(ntip-2)])*seq(3,(2*ntip-3),2)))
        sum2[s-1,1:(ntip-2)]=temp
     }
       
      R1=Rchild1[1:(ntip-1),1:(ntip-3)]
      R1=t(t(R1)/B[1:(ntip-3)])
      R2=Rchild2[1:(ntip-1),1:(ntip-3)]
      R2=t(t(R2)/B[1:(ntip-3)])
      R1aug=cbind(R1,matrix(0,nrow(R1),ncol(R1)))
      R2aug=cbind(R2,matrix(0,nrow(R2),ncol(R2)))
      U=round(convolve(t(R1aug), rev(t(R2aug)), type="open"),0)
      sum3=t(matrix(c(U,0),2*ncol(R1),2*nrow(R1))[1:ncol(R1),1:nrow(R1)])
      sum3=cbind(array(0, dim=c(nrow(R1)-1,1)),sum3[2:nrow(R1),])
      sum3=t(t(sum3)*B[2:(ntip-1)])
      R[[v]][2:(ntip-1),2:(ntip-1)]=sum1+sum2+sum3
     }
  }
  return(R)
}


#==========================================
RsT=function(R,n,s){
  rst =sum(R[[1]][s+1,1:(n-2-s)])
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

polynomial=function(tree,n){
  R=RConvolve(tree,n)
  for (i in seq(0,2*(n-3),2)) {
    print(qmT(R,n,n-3-(i/2)))
  }
}
