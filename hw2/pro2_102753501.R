######################################
# the reference code of program2 
######################################

######################################
# initial
######################################

library("Biostrings",verbose=F,quietly=T)



# read parameters
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("USAGE: Rscript pro2_<your student ID>.R --input test.fasta --score pam250.txt --aln global --gap_open -10 --gap_extend -2 --output result.fasta", call.=FALSE)
}

# parse parameters
i<-1 
while(i < length(args))
{
  if(args[i] == "--input"){
    i_f<-args[i+1]
    i<-i+1
  }else if(args[i] == "--score"){
    s_f<-args[i+1]
    i<-i+1
  }else if(args[i] == "--aln"){
    aln_mode <- args[i+1]
    i<-i+1
  }else if(args[i] == "--gap_open"){
    g_o<- as.integer(args[i+1]) 
    i<-i+1
  }else if(args[i] == "--gap_extend"){
    g_e<-as.integer(args[i+1])
    i<-i+1    
  }else if(args[i] == "--output"){
    o_f<-args[i+1]
    i<-i+1
  }else{
    stop(paste("Unknown flag", args[i]), call.=FALSE)
  }
  i<-i+1
}

print("PARAMETERS")
print(paste("input file         :", i_f))
print(paste("output file        :", o_f))
print(paste("score file         :", s_f))
print(paste("aln mode           :", aln_mode))
print(paste("gap open penalty   :", g_o))
print(paste("gap extend penalty :", g_e))

######################################
# main
######################################
# read fasta file
ff <- readAAStringSet(i_f)
seq_name = names(ff)
sequence = paste(ff)
sequence1 = paste(ff)

# aln length
aln_length <- nchar(sequence[1])

# read score file
s_m<-read.table(s_f)
s_m<-as.matrix(s_m)

aln_score<-0
#define an initial matrix
s <- array(0, dim=c(aln_length+1,aln_length+1))
#print(s)

#initial gap-open
for(i in 1:aln_length+1)
{
  if(i==1)s[i,i] = 0
  else s[i,1] = s[1,i] = s[i-1,1] + g_o
}

#print(s)

#best score
for(j in 2:aln_length+1)
{
  a<-substring(sequence[1], j-1, j-1)
  for(i in 2:aln_length+1)
  {
    
      
      #print(a)
      b<-substring(sequence[2], i-1, i-1)
      #print(b)
    
    
    #initial parameter
    s1 = s2 = s3 = 0
    
    # F(i,j) = F(i-1,j-1) + Mat(i,j)
    if((a != "-")&&(b != "-"))
    {
      s1 = s[i-1,j-1] + s_m[a,b]
      #print(s_m[a,b])
    }
    # if it's a gap-open
    else s1 = s[i-1,j-1] + g_o
    #print(s1)
    
    # F(i,j) = F(i,j-1) + gap_open
    s2 = s[i,j-1] + g_o
    #print(s2)
    
    # F(i,j) = F(i-1,j) + gap_open
    s3 = s[i-1,j] + g_o
    #print(s3)
    
    #global
    s[i,j] = s1
    if(s2 > s[i,j])s[i,j] = s2
    if(s3 > s[i,j])s[i,j] = s3
    
    #local
    if(aln_mode == "local")if(s[i,j]<0)s[i,j] = 0
    #print(s[i,j])
  }
}
#print(s)

# show score
#print(s[aln_length+1,aln_length+1])


print(sequence)

#traceback

i1 = aln_length+1
j1 = aln_length+1
a1<-substring(sequence[1],j1-1,j1-1)
substring(sequence1[1],j1-1,j1-1)<-a1
b1<-substring(sequence[2],i1-1,i1-1)
substring(sequence1[2],i1-1,i1-1)<-b1

while((i1>1)&&(j1>1))
{
  #go up and left if a1 or b1 is not -
   if((a1 != "-")&&(b1 != "-"))
   {
     if(s[i1,j1]==(s[i1-1,j1-1]+s_m[a1,b1]))
     {
       i1 = i1 - 1
       j1 = j1 - 1
       a1<-substring(sequence[1],j1-1,j1-1)
       substring(sequence1[1],j1-1,j1-1)<-a1
       b1<-substring(sequence[2],i1-1,i1-1)
       substring(sequence1[2],i1-1,i1-1)<-b1
       
     }
   }
   #go up and left if a1 or b1 is -
   # else if(s[i1,j1]==(s[i1-1,j1-1]+g_o))
   # {
   #   i1 = i1 - 1
   #   j1 = j1 - 1
   #   a1<-substring(sequence[1],j1-1,j1-1)
   #   substring(sequence1[1],j1-1,j1-1)<-a1
   #   b1<-substring(sequence[2],i1-1,i1-1)
   #   substring(sequence1[2],i1-1,i1-1)<-b1
   # }
  #go left
  else if(s[i1,j1]==(s[i1-1,j1]+g_o))
  {
    i1 = i1 - 1
    a1<-substring(sequence[1],j1-1,j1-1)
    substring(sequence1[1],j1-1,j1-1)<-a1
    b1 = "-"
    substring(sequence1[2],i1-1,i1-1)<-b1
  }
  # #go up
   else if(s[i1,j1]==(s[i1,j1-1]+g_o))
   {
     j1 = j1 - 1
     a1 = "-"
     substring(sequence1[1],j1-1,j1-1)<-a1
     b1<-substring(sequence[2],i1-1,i1-1)
     substring(sequence1[2],i1-1,i1-1)<-b1
   }
   #default
   else
   {
     i1 = i1 - 1
     j1 = j1 - 1
     a1 = "-"
     substring(sequence1[1],j1-1,j1-1)<-a1
     b1 = "-"
     substring(sequence1[2],i1-1,i1-1)<-b1
   }
}
print(sequence1)

    
    
    


# output
writeXStringSet(ff, o_f)
