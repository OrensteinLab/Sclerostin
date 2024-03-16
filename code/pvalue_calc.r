library(exactci)
df <- read.csv('all data by variants 1 mutations.csv')

lib1 <- df[["mutations1"]] 
lib2 <- df[["mutations2"]] 

Table<-matrix(NA,length(lib1),1)

for (i in 1:length(lib1)){
    Table[i,1]<-poisson.exact(c(lib1[i], lib2[i]), c(139568,141101),alternative ="two.sided",tsmethod="minlike",conf.level=0.95)$p.value
}

library(MASS)
write.matrix(Table,file="pvalue.csv")