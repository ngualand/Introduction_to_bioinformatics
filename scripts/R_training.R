# Exercise 1: Loading the Dataset
library(Biobase)
data(package = "Biobase")
data("geneData")
data("geneCovariate")

#Excercise 1: Transform data to data.frame and explore the dataframe, what they contains, what are rows and columns?
#Can you tell what geneData and GeneCovariates dataframes contain?
geneData = data.frame(geneData)
geneCovariate = data.frame(geneCovariate)

# Exercise 2: Basic Data Exploration, check what is the structure of data and count the number of samples (lines and rows) in each dataframe
str(geneData)
str(geneCovariate)
summary(geneData)
summary(geneCovariate)
dim(geneData)
dim(geneCovariate)

# Exercise 3: Retrieve expression data and covariates for the sample A, B, C, D and save them to 2 different objects
geneData.subset = geneData[, c("A", "B", "C", "D")]
geneCovariate.subset = geneCovariate[c("A", "B", "C", "D"), ]

# Exercise 4: Retrieve only probes (genes) with expression > 80 in sample A and check that the filter have worked (hint: for the test use the min funnction)
#How many probes pass this filter?
geneData.A.gt20 = geneData[geneData$A > 20,]
nrow(geneData.A.gt20)

#Excercise 4 bis: Create a variable that contains the sum of the expression of the probes AFFX-hum_alu_at in sample A and B
#What is the sum of expression of probes AFFX-hum_alu_at in sample A and B
exp.sum = geneData["AFFX-hum_alu_at", "A"] + geneData["AFFX-hum_alu_at", "B"]
#What is the higher expression for this probe
max.exp = max(geneData["AFFX-hum_alu_at", ])


# Exercise 5: Descriptive Statistics, find mean and median gene expression in sample A and B, which sample have the higher mean gene expression? 
mean.exp.A = mean(geneData$A)
mean.exp.B = mean(geneData$B)
mean.exp.A > mean.exp.B


# Exercise 6: Data Visualization, perform the correlation (scatterplot) of gene expression in sample A and B and sample A and L,
#what samples are similar (hint: perfomr a correlation test , cor.test) ? Can you hypothesize why? (Hint: look at the geneCovariates data)
plot(geneData$A,geneData$B, main="Gene expression A vs B", xlab="Gene expression in A", ylab="Gene expression in B", pch=19)
plot(geneData$A,geneData$C, main="Gene expression A vs L", xlab="Gene expression in A", ylab="Gene expression in L", pch=19)
cor.test(geneData$A,geneData$B)
cor.test(geneData$A,geneData$L)

# Exercise 7: test if the expression of probe (gene) AFFX-hum_alu_at, is significantly different between Case and Control, regardless of the Sex (hint: use a t.test)
t.test.result.group = t.test(as.numeric(geneData["AFFX-hum_alu_at",]) ~ geneCovariate$type)
print(t.test.result.group)
#and by Sex?
t.test.result.sex = t.test(as.numeric(geneData["AFFX-hum_alu_at",]) ~ geneCovariate$sex)
print(t.test.result.sex)

#Excercise 8: Plot the distribution of expression for gene AFFX-hum_alu_at in Case and controls and Male and Female 
#The boxplot are in accordance with the t,test performed before?
boxplot(as.numeric(geneData["AFFX-hum_alu_at",]) ~ geneCovariate$type)
boxplot(as.numeric(geneData["AFFX-hum_alu_at",]) ~ geneCovariate$sex)


# Exercise 10: Exporting Data, write both geneData and geneCovariates to 2 diffrent files and save them to disk for future use
write.csv(geneData, file = "geneData.csv", row.names = FALSE)
write.csv(geneCovariate, file = "geneCovariate.csv", row.names = FALSE)

