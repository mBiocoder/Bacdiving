#For all 12 Meta IBS datasets extract the taxtable from phyloseq object and save it as tsv file

######################################## AGP ##################################
#phyloseq object
agp_phyloseq_object <- readRDS("./input_data/phyloseq_objects/physeq_agp.rds")
agp_taxtab <- as.data.frame(agp_phyloseq_object@tax_table@.Data)
nrow(agp_taxtab)
head(agp_taxtab)

#Write to tsv file
write.table(agp_taxtab, file = "./input_data/taxtable_from_phyloseq/agp_taxtab.tsv", quote=FALSE, row.names = TRUE, col.names=NA, sep='\t',)
d <- read.table("./input_data/taxtable_from_phyloseq/agp_taxtab.tsv")
nrow(d)
head(d)


######################################## FUKUI ##################################
#phyloseq object
fukui_phyloseq_object <- readRDS("./input_data/phyloseq_objects/physeq_fukui.rds")
fukui_taxtab <- as.data.frame(fukui_phyloseq_object@tax_table@.Data)
nrow(fukui_taxtab)
head(fukui_taxtab)

#Write to tsv file
write.table(fukui_taxtab, file = "./input_data/taxtable_from_phyloseq/fukui_taxtab.tsv", quote=FALSE, row.names = TRUE, col.names=NA, sep='\t',)
d <- read.table("./input_data/taxtable_from_phyloseq/fukui_taxtab.tsv")
nrow(d)
head(d)



######################################## HUGERTH ##################################
#phyloseq object
hugerth_phyloseq_object <- readRDS("./input_data/phyloseq_objects/physeq_hugerth.rds")
hugerth_taxtab <- as.data.frame(hugerth_phyloseq_object@tax_table@.Data)
nrow(hugerth_taxtab)
head(hugerth_taxtab)

#Write to tsv file
write.table(hugerth_taxtab, file = "./input_data/taxtable_from_phyloseq/hugerth_taxtab.tsv", quote=FALSE, row.names = TRUE, col.names=NA, sep='\t',)
d <- read.table("./input_data/taxtable_from_phyloseq/hugerth_taxtab.tsv")
nrow(d)
head(d)



######################################## LABUS ##################################
#phyloseq object
labus_phyloseq_object <- readRDS("./input_data/phyloseq_objects/physeq_labus.rds")
labus_taxtab <- as.data.frame(labus_phyloseq_object@tax_table@.Data)
nrow(labus_taxtab)
head(labus_taxtab)

#Write to tsv file
write.table(labus_taxtab, file = "./input_data/taxtable_from_phyloseq/labus_taxtab.tsv", quote=FALSE, row.names = TRUE, col.names=NA, sep='\t',)
d <- read.table("./input_data/taxtable_from_phyloseq/labus_taxtab.tsv")
nrow(d)
head(d)


######################################## LIU ##################################
#phyloseq object
liu_phyloseq_object <- readRDS("./input_data/phyloseq_objects/physeq_liu.rds")
liu_taxtab <- as.data.frame(liu_phyloseq_object@tax_table@.Data)
nrow(liu_taxtab)
head(liu_taxtab)

#Write to tsv file
write.table(liu_taxtab, file = "./input_data/taxtable_from_phyloseq/liu_taxtab.tsv", quote=FALSE, row.names = TRUE, col.names=NA, sep='\t',)
d <- read.table("./input_data/taxtable_from_phyloseq/liu_taxtab.tsv")
nrow(d)
head(d)


######################################## LOPRESTI ##################################
#phyloseq object
lopresti_phyloseq_object <- readRDS("./input_data/phyloseq_objects/physeq_lopresti.rds")
lopresti_taxtab <- as.data.frame(lopresti_phyloseq_object@tax_table@.Data)
nrow(lopresti_taxtab)
head(lopresti_taxtab)

#Write to tsv file
write.table(lopresti_taxtab, file = "./input_data/taxtable_from_phyloseq/lopresti_taxtab.tsv", quote=FALSE, row.names = TRUE, col.names=NA, sep='\t',)
d <- read.table("./input_data/taxtable_from_phyloseq/lopresti_taxtab.tsv")
nrow(d)
head(d)

######################################## MARS ##################################
#phyloseq object
mars_phyloseq_object <- readRDS("./input_data/phyloseq_objects/physeq_mars.rds")
mars_taxtab <- as.data.frame(mars_phyloseq_object@tax_table@.Data)
nrow(mars_taxtab)
head(mars_taxtab)

#Write to tsv file
write.table(mars_taxtab, file = "./input_data/taxtable_from_phyloseq/mars_taxtab.tsv", quote=FALSE, row.names = TRUE, col.names=NA, sep='\t',)
d <- read.table("./input_data/taxtable_from_phyloseq/mars_taxtab.tsv")
nrow(d)
head(d)


######################################## NAGEL ##################################
#phyloseq object
nagel_phyloseq_object <- readRDS("./input_data/phyloseq_objects/physeq_nagel.rds")
nagel_taxtab <- as.data.frame(nagel_phyloseq_object@tax_table@.Data)
nrow(nagel_taxtab)
head(nagel_taxtab)

#Write to tsv file
write.table(nagel_taxtab, file = "./input_data/taxtable_from_phyloseq/nagel_taxtab.tsv", quote=FALSE, row.names = TRUE, col.names=NA, sep='\t',)
d <- read.table("./input_data/taxtable_from_phyloseq/nagel_taxtab.tsv")
nrow(d)
head(d)


######################################## POZUELO ##################################
#phyloseq object
pozuelo_phyloseq_object <- readRDS("./input_data/phyloseq_objects/physeq_pozuelo.rds")
pozuelo_taxtab <- as.data.frame(pozuelo_phyloseq_object@tax_table@.Data)
nrow(pozuelo_taxtab)
head(pozuelo_taxtab)

#Write to tsv file
write.table(pozuelo_taxtab, file = "./input_data/taxtable_from_phyloseq/pozuelo_taxtab.tsv", quote=FALSE, row.names = TRUE, col.names=NA, sep='\t',)
d <- read.table("./input_data/taxtable_from_phyloseq/pozuelo_taxtab.tsv")
nrow(d)
head(d)


######################################## RINGEL ##################################
#phyloseq object
ringel_phyloseq_object <- readRDS("./input_data/phyloseq_objects/physeq_ringel.rds")
ringel_taxtab <- as.data.frame(ringel_phyloseq_object@tax_table@.Data)
nrow(ringel_taxtab)
head(ringel_taxtab)

#Write to tsv file
write.table(ringel_taxtab, file = "./input_data/taxtable_from_phyloseq/ringel_taxtab.tsv", quote=FALSE, row.names = TRUE, col.names=NA, sep='\t',)
d <- read.table("./input_data/taxtable_from_phyloseq/ringel_taxtab.tsv")
nrow(d)
head(d)


######################################## ZEBER ##################################
#phyloseq object
zeber_phyloseq_object <- readRDS("./input_data/phyloseq_objects/physeq_zeber.rds")
zeber_taxtab <- as.data.frame(zeber_phyloseq_object@tax_table@.Data)
nrow(zeber_taxtab)
head(zeber_taxtab)

#Write to tsv file
write.table(zeber_taxtab, file = "./input_data/taxtable_from_phyloseq/zeber_taxtab.tsv", quote=FALSE, row.names = TRUE, col.names=NA, sep='\t',)
d <- read.table("./input_data/taxtable_from_phyloseq/zeber_taxtab.tsv")
nrow(d)
head(d)


######################################## ZHU ##################################
#phyloseq object
zhu_phyloseq_object <- readRDS("./input_data/phyloseq_objects/physeq_zhu.rds")
zhu_taxtab <- as.data.frame(zhu_phyloseq_object@tax_table@.Data)
nrow(zhu_taxtab)
head(zhu_taxtab)

#Write to tsv file
write.table(zhu_taxtab, file = "./input_data/taxtable_from_phyloseq/zhu_taxtab.tsv", quote=FALSE, row.names = TRUE, col.names=NA, sep='\t',)
d <- read.table("./input_data/taxtable_from_phyloseq/zhu_taxtab.tsv")
nrow(d)
head(d)


######################################## ZHUANG ##################################
#phyloseq object
zhuang_phyloseq_object <- readRDS("./input_data/phyloseq_objects/physeq_zhuang.rds")
zhuang_taxtab <- as.data.frame(zhuang_phyloseq_object@tax_table@.Data)
nrow(zhuang_taxtab)
head(zhuang_taxtab)

#Write to tsv file
write.table(zhuang_taxtab, file = "./input_data/taxtable_from_phyloseq/zhuang_taxtab.tsv", quote=FALSE, row.names = TRUE, col.names=NA, sep='\t',)
d <- read.table("./input_data/taxtable_from_phyloseq/zhuang_taxtab.tsv")
nrow(d)
head(d)


