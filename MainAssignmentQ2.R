setwd("D:/TCD_MScDataScience/StatisticalModelling/LabFiles")
df1 <- read.csv("heart_data.csv")
df1$Class <- factor(df1$Class)
df1$Sex <- factor(df1$Sex)
df1$FastingBloodSugar <- factor(df1$FastingBloodSugar)
df1$ExerciseInduced <- factor(df1$ExerciseInduced)
df1$Slope <- factor(df1$Slope)
df1$MajorVessels <- factor(df1$MajorVessels)

df1$Age <- scale(df1$Age)
df1$RestBloodPressure <- scale(df1$RestBloodPressure)
df1$SerumCholestoral <- scale(df1$SerumCholestoral)
df1$MaxHeartRate <- scale(df1$MaxHeartRate)
head(df1)

apply(df1, 2, sd)

#Correlation scatterplots
pairs(subset(df1, select = c(Age, RestBloodPressure, SerumCholestoral, MaxHeartRate )))

#correlation matrix
cor(subset(df1, select = c(Age, RestBloodPressure, SerumCholestoral, MaxHeartRate )))

#boxplot b/w sex and response variable
boxplot(MaxHeartRate ~ Sex, data = df1)

##fitting the model
fit1 <- glm(Class ~ ., family = binomial(link = logit), data = df1)
summary(fit1)

summary(df1)
