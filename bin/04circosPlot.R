rm(list = ls())
# CIRCOS PLOT - TI-1 ------------------------------------------------------

library(circlize)
library(dplyr)

setwd("/home/clovis/Dropbox/Chumbo/git (cópia 1)/")

load("./Data/allTranscriptogramers80")

load("./Data/assocNoDup.RData")
load("./Data/clusters12.RData")

# Calculate gaps
temp1 <- clusters12[, -4]
temp1$na <- is.na(temp1$Clust1)

tmp <- 30
na <- c()
for (i in 1:nrow(temp1)) {
  if(temp1$na[i]) {
    na[i] <- tmp
  }
  else {
    na[i] <- temp1$Clust1[i]
    tmp <- tmp + 1
  }
}
temp1$clusters2 <- na

# Get association data
time1 <- transc[[1]]
association <- transc[[1]]@association
idx <- !duplicated(t(apply(association, 1, sort)))
assocNoDup <- association[idx,]

# Merge and organize dataframes
df1  <- merge(temp1, assocNoDup, by.x = "Protein", by.y = "p1")
colnames(df1)[1] <- "p1"
df2 <- unique(merge(temp1, df1, by.x = "Protein", by.y = "p2"))
df3 <- as.data.frame(table(df2$clusters2.x, df2$clusters2.y), stringsAsFactors = F)

# Calculate clusters size
size2 <- temp1 %>% 
  group_by(clusters2) %>% 
  summarise(size = n())

# Merge and prepare plot data
c_filter3 <- merge(df3, size2, by.x = "Var1", by.y = "clusters2")
c_filter4 <- c_filter3 %>% 
  group_by(Var1) %>% 
  mutate(k = Freq / size) %>% 
  select(Var1, Var2, k)

c_filter5 <- c_filter4 %>%
  ungroup() %>% 
  mutate(Var1 = as.numeric(Var1), Var2 = as.numeric(Var2)) %>% 
  filter(Var1 < 30 & Var2 < 30) %>% 
  arrange(Var1)

col_grid <- c("#ff2222ff", "#009900ff", "#ffff00ff", "#ff9e30ff", "#1818fcff", "#b2a77fff", 
              "#000055ff", "#00aaaaff", "#00ffffff", "#8000aaff", "#fd11c2ff")

# Plot
circos.clear()
chordDiagram(c_filter5, grid.col = col_grid)





# CIRCOS PLOT - TI-2 ------------------------------------------------------

library(circlize)
library(dplyr)


# Calculate gaps
temp1 <- clusters12[, -3]
temp1$na <- is.na(temp1$Clust2)

tmp <- 30
na <- c()
for (i in 1:nrow(temp1)) {
  if(temp1$na[i]) {
    na[i] <- tmp
  }
  else {
    na[i] <- temp1$Clust2[i]
    tmp <- tmp + 1
  }
}
temp1$clusters2 <- na

# Get association data
time1 <- transc[[2]]
association <- transc[[2]]@association
idx <- !duplicated(t(apply(association, 1, sort)))
assocNoDup <- association[idx,]

# Merge dataframes
df1  <- merge(temp1, assocNoDup, by.x = "Protein", by.y = "p1")
colnames(df1)[1] <- "p1"
df2 <- unique(merge(temp1, df1, by.x = "Protein", by.y = "p2"))
df3 <- as.data.frame(table(df2$clusters2.x, df2$clusters2.y), stringsAsFactors = F)

# Calculate clusters size
size2 <- temp1 %>% 
  group_by(clusters2) %>% 
  summarise(size = n())

# Merge and prepare plot data
c_filter3 <- merge(df3, size2, by.x = "Var1", by.y = "clusters2")
c_filter4 <- c_filter3 %>% 
  group_by(Var1) %>% 
  mutate(k = Freq / size) %>% 
  select(Var1, Var2, k)

c_filter5 <- c_filter4 %>%
  ungroup() %>% 
  mutate(Var1 = as.numeric(Var1), Var2 = as.numeric(Var2)) %>% 
  filter(Var1 < 30 & Var2 < 30) %>% 
  arrange(Var1)

col_grid <- c("#ff2222ff", "#ff9e30ff", "#01ff40ff", "#0000ffff", "#ada27aff", "#000055ff", 
              "#00aaaaff", "#8000aaff", "#ff05c0ff", "#6666ffff", "#3333aaff")

# Plot
circos.clear()
chordDiagram(c_filter5, grid.col = col_grid)






