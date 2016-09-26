# read PAM1 from data
pam1<-read.table("Z:/Users/Hu-2016/Desktop/BI16/hw1/pam1.txt")

# check PAM1 data
dim(pam1)
str(pam1)

pam1 <- data.matrix(pam1)
pam1 <- pam1 / 10000

#print(pam1)
# construct PAM250 from PAM1
pam250 <- pam1
for(i in 1:250){
  pam250 <- pam250 %*% pam1
}

# output PAM250 as a file
pam250 <- pam250 * 100
write.table(pam250,file = "Z:/Users/Hu-2016/Desktop/BI16/hw1/pam250.txt")

