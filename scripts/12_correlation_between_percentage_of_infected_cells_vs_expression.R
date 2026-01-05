##Supplementary figure S12
###running 2_cell_type_cluster_NV_DAV_tropism_with_infection_percentage_in_cell_types.R script first
# df is your original dataframe
df_wide <- df %>%
  pivot_wider(names_from = virus, values_from = percent)

# View resulting dataframe
head(df_wide)
colnames(df_wide) <- c("cluster.ident", "condition", "NV" ,"DAV")

df_wide <- df_wide %>% 
  filter(condition == "Virgin")
 

meta<- fatbody_v.harmony@meta.data
# ---- Feature extraction ----
# Fetch expression values for specific genes/TEs/viral features from the Seurat object.
DAV<-FetchData(fatbody_v.harmony, vars = "KP969946.1-Drosophila-A-virus-isolate-LJ35") 
NV<-FetchData(fatbody_v.harmony, vars = "JX220408.1-Nora-virus-isolate-FR1") 
meta$DAV<- DAV$`KP969946.1-Drosophila-A-virus-isolate-LJ35`
meta$NV<- NV$`JX220408.1-Nora-virus-isolate-FR1`
meta <- meta[ , c("Identity", "type","cluster.ident","DAV","NV")]

meta <- meta %>% 
  filter(type == "Virgin", NV > 0)

meta <- meta %>% 
  group_by(cluster.ident) %>% 
  summarise(mean_NV = mean(NV, na.rm = TRUE))

combined_df <- merge(meta, df_wide,
                     by = "cluster.ident")

ggscatter(combined_df, 
          x = "NV", 
          y = "mean_NV",
          color = "cluster.ident",
          add = "reg.line",       
          add.params = list(color = "black", linetype = "solid"), # global line in black dashed style
          conf.int = F,        
          cor.coef = T,        
          cor.method = "spearman",
          xlab = "% of NV infected cells", 
          ylab = "NV Normalized reads (viral RNA)", 
          legend = "right")



# df is your original dataframe
df_wide <- df %>%
  pivot_wider(names_from = virus, values_from = percent)

# View resulting dataframe
head(df_wide)
colnames(df_wide) <- c("cluster.ident", "condition", "NV" ,"DAV")

df_wide <- df_wide %>% 
  filter(condition == "Mated")


meta<- fatbody_v.harmony@meta.data
DAV<-FetchData(fatbody_v.harmony, vars = "KP969946.1-Drosophila-A-virus-isolate-LJ35") 
NV<-FetchData(fatbody_v.harmony, vars = "JX220408.1-Nora-virus-isolate-FR1") 
meta$DAV<- DAV$`KP969946.1-Drosophila-A-virus-isolate-LJ35`
meta$NV<- NV$`JX220408.1-Nora-virus-isolate-FR1`
meta <- meta[ , c("Identity", "type","cluster.ident","DAV","NV")]

meta <- meta %>% 
  filter(type == "Mated", NV > 0, DAV >0)

meta <- meta %>%
  group_by(cluster.ident) %>%
  summarise(mean_NV = mean(NV, na.rm = TRUE),
            mean_DAV = mean(DAV, na.rm = TRUE))


combined_df <- merge(meta, df_wide,
                     by = "cluster.ident")

ggscatter(combined_df, 
          x = "NV", 
          y = "mean_NV",
          color = "cluster.ident",
          add = "reg.line",       
          add.params = list(color = "black", linetype = "solid"), # global line in black dashed style
          conf.int = F,        
          cor.coef = T,        
          cor.method = "spearman",
          xlab = "% of NV infected cells", 
          ylab = "NV Normalized reads (viral RNA)", 
          legend = "right")

ggscatter(combined_df, 
          x = "DAV", 
          y = "mean_DAV",
          color = "cluster.ident",
          add = "reg.line",       
          add.params = list(color = "black", linetype = "solid"), # global line in black dashed style
          conf.int = F,        
          cor.coef = T,        
          cor.method = "spearman",
          xlab = "% of DAV infected cells", 
          ylab = "DAV Normalized reads (viral RNA)", 
          legend = "right")
