#Ranking The pairs of spatial fields based on spatial interactions

rank = matrix(, ncol =20,nrow =20)
for(i in 1:20)
{
  for(j in 1:20)
  {
    rank[i,j] = ((areamatrixI[i,j]/areamatrixS[i,j])*(min(dilationmatrix[i,j],erosionmatrix[i,j])/max(dilationmatrix[i,j],erosionmatrix[i,j])))
  }
}