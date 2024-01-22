library(ggplot2)
library(ggrepel)

fnmCSV = '/zDsk/Manuscripts/Baydoun/UK-BioBank/2024-01-05—12—ProteomicADPRS/Volcano/ADPRS_PROTEOME_VOLCANOPLOTDATA.csv'

datObj = read.csv(fnmCSV, as.is=T)
nrow(datObj)
summary(datObj)
head(datObj)

pltDat = within(datObj, {
	lblOutcome	= sub('^ z', '', zOutcome, perl=T)
	})

summary(pltDat)
head(pltDat)

options(ggrepel.max.overlaps = Inf)	# avoids error message about colliding labels

Volcano = function(volDat) {
	ggObj =	ggplot(data=volDat, aes(x=estimate, y=mlog10p, label=lblOutcome)) +
		geom_point(shape=20, size=1) +
                geom_point(shape=20, size=3, col='red', data=subset(volDat, selected==1)) +
		geom_point(shape=20, size=3, col='orange', data=subset(volDat, selected==1 & signif==1 & estimate>0.05)) +
                geom_point(shape=20, size=3, col='blue', data=subset(volDat, selected==1 & signif==1 & estimate< -0.05)) +
		xlim(-.15, .15) + ylim(0, 125) +
		geom_text_repel(size=3, col='red',  data=subset(volDat, selected==1 & estimate >  0.03), aes(label=lblOutcome)) +
		geom_text_repel(size=3, col='blue', data=subset(volDat, selected==1 & estimate < -0.03), aes(label=lblOutcome))
                labs(x='Effect size') +
		theme_minimal() +
		geom_vline(xintercept=-.05, col="gray", linetype=2) +
		geom_vline(xintercept=.05,  col="gray", linetype=2) +
		geom_hline(yintercept=-log10(0.05/1463), col="gray", linetype=2)

	ggObj
	}

(pltObj = Volcano(pltDat))

dirName = '/zDsk/Manuscripts/Baydoun/UK-BioBank/2024-01-05—12—ProteomicADPRS/Volcano/'
imgName = '2024-01-14-UKB12-Volcano'
pdfFile = paste0(dirName, imgName, '.pdf')
pngFile = paste0(dirName, imgName, '.png')

pdf(pdfFile)
pltObj
jnk = dev.off()

png(pngFile)
pltObj
jnk = dev.off()
