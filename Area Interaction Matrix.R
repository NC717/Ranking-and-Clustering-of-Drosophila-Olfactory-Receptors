#Area Interaction Matrix

count =1
areamatrixS = matrix(,nrow =64 ,ncol=64)
areamatrixI = matrix(,nrow =64, ncol =64)
for(i in 1:64)
{  
  for(j in 1:64)
  { 
    areamatrixS[i,j] = sum(los[[count]])
    areamatrixI[i,j] = sum(loi[[count]])
    count = count +1
  }
}
