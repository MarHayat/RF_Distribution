#Beta function
Beta=function(m){
  if(m<0){return(1)}
  return(dfactorial(2*m+1))
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

#============================================
RFPaper=function(tree,n){
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
      sum3=matrix(0,(ntip-2),(ntip-2))
      for (s in 1:(ntip-2)){
        for (k in 2:(ntip-2)) {
          total3=0
          for (s1 in 0:(s)) {
            for (k1 in 0:(k-2)) {
              # if((s1+1) > 0 & (k1+1) > 0 &(s-s1+1) > 0 & (k-2-k1+1) > 0){
                total3=total3+Rchild1[s1+1,k1+1]*Rchild2[s-s1+1,k-2-k1+1]*Beta(k)/(Beta(k1)*Beta(k-2-k1))
              # }
            }
            sum3[s,k]=total3
          }
          
        }
       }
      
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
  R=RFPaper(tree,n)
  for (i in seq(0,2*(n-3),2)) {
    print(qmT(R,n,n-3-(i/2)))
  }
}
