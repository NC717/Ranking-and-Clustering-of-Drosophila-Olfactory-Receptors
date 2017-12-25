#Defining the pattern

library(stringr)
amino = c("A","B","C","D","E","F","G","H")
listofsequence= read.table("PPCA_E.txt", sep = "\t")
counter =1
listofsequence_clean = str_trim(listofsequence[[1]])
for(i in 1 :length(listofsequence_clean))
{
  listofsequence_clean[i] = str_replace_all(listofsequence_clean[i],"[DE]", "A") 
  listofsequence_clean[i] = str_replace_all(listofsequence_clean[i],"[RHK]", "B")
  listofsequence_clean[i] = str_replace_all(listofsequence_clean[i],"[YFW]", "C")
  listofsequence_clean[i] = str_replace_all(listofsequence_clean[i],"[ILVAG]", "D")
  listofsequence_clean[i] = str_replace_all(listofsequence_clean[i],"[P]", "E")
  listofsequence_clean[i] = str_replace_all(listofsequence_clean[i],"[MC]", "F")
  listofsequence_clean[i] = str_replace_all(listofsequence_clean[i],"[ST]", "G")
  listofsequence_clean[i] = str_replace_all(listofsequence_clean[i],"[QN]", "H")
}
count =1
pattern = c()
for(i in 1:length(amino))
{
  for(j in 1:length(amino))
  {
    pattern[count] = str_c(amino[i],amino[j])
    count =count +1
  }
}



#Finding the frequency

frequency =list()
LOF =list()
for(i in 1:length(listofsequence_clean)) {
  frequency[[i]]= str_count(listofsequence_clean[[i]],pattern)
}
for(i in 1:length(listofsequence_clean))
{
  LOF[[i]] = matrix(frequency[[i]], ncol =8, nrow=8)
}
for(i in 1:length(LOF))
{
  colnames(LOF[[i]]) = amino
  rownames(LOF[[i]]) = amino
}

#Findind Infima and Suprema


infima = matrix(,ncol =8 ,nrow=8)
suprema =matrix(,ncol =8, nrow=8)
loi= list()
los= list()
counter =1
count=1
count1=1
for(count in 1:length(LOF))
{
  for(count1 in 1: length(LOF)) 
  { for(i in 1:8)
  {
    for(j in 1:8)
    {
      if(LOF[[count]][i,j] > LOF[[count1]][i,j])
      {
        suprema[i,j] = LOF[[count]][i,j]
        infima[i,j] = LOF[[count1]][i,j]
      }
      else 
      {
        suprema[i,j] = LOF[[count1]][i,j]
        infima[i,j] = LOF[[count]][i,j]
      }
    }
  }
    los[[counter]] = suprema 
    loi[[counter]] = infima
    counter = counter+1
  }
}

#Area Interaction Matrix

count =1
areamatrixS = matrix(,nrow =length(LOF) ,ncol=length(LOF))
areamatrixI = matrix(,nrow =length(LOF), ncol =length(LOF))
for(i in 1:length(LOF))
{  
  for(j in 1:length(LOF))
  { 
    areamatrixS[i,j] = sum(los[[count]])
    areamatrixI[i,j] = sum(loi[[count]])
    count = count +1
  }
}


#Grayscale Morphological Dilation Distance

library(mmand)
counter = 1
B= matrix(c(0,1,0,1,1,1,0,1,0),nrow=3,ncol=3)
dilL = list()
ns =c()
xy = as.vector(areamatrixS) 
for(counter in 1:length(LOF))
{
  
  n=1
  s1= sum(dilate(loi[[counter]],B))
  dil = dilate(loi[[counter]],B)
  if(xy[counter]>=s1)
  {  
    while(!(xy[counter] <= s1))
    { 
      if(((sum(dilate(dil,B))==sum(dilate((dilate(dil,B)),B)))) & ((sum(dil))<=xy[counter])) 
      {
        n =n+1
        break
      }
      else
      {
        dil = dilate(dil, B)
        s1 = sum(dil)
        n=n+1
      }
    }     
    dilL[[counter]] = dil    
    ns[counter] =n
  }
  else
  {
    dilL[[counter]] =dil
    ns[counter]=n
  }
}
dilationmatrix =matrix(ns, nrow=length(LOF), ncol =length(LOF))

#Grayscale Morphological Erosion Distance

counter = 1
eroL = list()
ni =c()
xz = as.vector(areamatrixI) 
for(counter in 1:length(LOF))
{
  m=1
  s2= sum(erode(los[[counter]],B))
  ero = erode(los[[counter]],B)
  if(s2 >= xz[counter])
  { 
    while(!(s2 <= xz[counter]))
    { 
      if((sum(erode(ero,B)) == sum(erode(erode(ero,B),B))) & ((sum(ero))>=xz[counter]))
      {
        m= m+1
        break
      }
      else
      {
        ero = erode(ero, B)
        s2 = sum(ero)
        m=m+1
      }  
    }
    eroL[[counter]] = ero    
    ni[counter] = m
  }
  else
  {
    eroL[[counter]] =ero
    ni[counter] = m
  }
}
erosionmatrix = matrix(ni, nrow =length(LOF),ncol =length(LOF))

#Ranking The pairs of spatial fields based on spatial interactions

rank = matrix(, ncol =length(LOF),nrow =length(LOF))
for(i in 1:length(LOF))
{
  for(j in 1:length(LOF))
  {
    rank[i,j] = ((areamatrixI[i,j]/areamatrixS[i,j])*(min(dilationmatrix[i,j],erosionmatrix[i,j])/max(dilationmatrix[i,j],erosionmatrix[i,j])))
  }
}

#Phylogenetic Tree and Clustering

hcl = hclust(as.dist(rank), method ="single")
plot(hcl)