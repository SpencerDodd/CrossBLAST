summary_data <- read.csv('/Users/spencerdodd/Documents/Research/Khrapko Lab/Mitochondrial Subspecies/Summaries/Summary_15_12_07_(23h_1m_3s).csv')

levels = c('Subspecies', 'Species', 'Genus', 'Subfamily', 'Family', 'Other')

other <- subset(summary_data$Percent_difference, summary_data$Level == 'Other')
family <- subset(summary_data$Percent_difference, summary_data$Level == 'Family')
subfamily <- subset(summary_data$Percent_difference, summary_data$Level == 'Subfamily')
genus <- subset(summary_data$Percent_difference, summary_data$Level == 'Genus')
species <- subset(summary_data$Percent_difference, summary_data$Level == 'Species')
subspecies <- subset(summary_data$Percent_difference, summary_data$Level == 'Subspecies')

par(mfrow=c(3, 2))
hist(other, main = 'Other', xlab = 'Distance to common ancestor (% divergence / 2)')
hist(family, main = 'Family', xlab = 'Distance to common ancestor (% divergence / 2)')
hist(subfamily, main = 'Subfamily', xlab = 'Distance to common ancestor (% divergence / 2)')
hist(genus, main = 'Genus', xlab = 'Distance to common ancestor (% divergence / 2)')
hist(species, main = 'Species', xlab = 'Distance to common ancestor (% divergence / 2)')
hist(subspecies, main = 'Subspecies', xlab = 'Distance to common ancestor (% divergence / 2)')