#nuevo proyecto
#Script methane data management

library(data.table)
library(dplyr)
library(plyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(GGally)

# Set the directory
setwd("C:/Users/Ester/Documents/curso metano")
################################## READ FILES ############################################
# Read the sniffer output
bd=read.table("output_methane.txt",sep=",",header=T) #500 rows & 20 cols
bd_id=unique(bd[,1]) #63 individuals
summary(bd)
# Read test day file
test=read.table("control.csv",sep=";",header=T)#107 rows & 21 cols
test_id=unique(test[,1])#64 #july and sept
#obtain date from test day records
test$test_date1=dmy(test$test_date)

# Obtain kgm of protein and fat 

test$kgmfat=test$fat*test$milk/100                  
test$kgmprotein=test$protein*test$milk/100 
################################### Join both tables by ID and closest date ##################
# Obtain date from output
bd<-tidyr::separate(data = bd, col ="date",into =  c("day", "month","date","time","z","year"), sep = " ", extra = "merge")
bd$sniffer_date=paste(bd$date, bd$month, bd$year, sep = "/")
bd$sniffer_date1<-dmy(bd$sniffer_date)

#colSums(is.na(test))
# Join the two datatables-
bd_full <- lapply(intersect(bd$cow,test$cow),function(id) {
  d1 <- subset(bd,cow==id)
  d2 <- subset(test,cow==id)
  
  d1$indices <- sapply(as.Date(d1$sniffer_date1),function(d) which.min(abs(as.Date(d2$test_date1) - d)))
  d2$indices <- 1:nrow(d2)
  
  merge(d1,d2,by=c('cow','indices'))
}) ###

bd_full <- do.call(rbind,bd_full) #342
bd_full$indices <- NULL
#Check animals not joining 
bd_antijoin=anti_join(bd,bd_full,by="cow") #158, 3 individuals
bd_full$dif_days=abs(bd_full$sniffer_date1-bd_full$test_date1)
#check the difference of days ###choose less than 40 days
table(bd_full$dif_days)

############################### CALVING AND MILKING ###################################

# Days in milking
bd_full=tidyr::separate(data=bd_full,col ="calving_date", into=c('calving_date', 'calving_time'),sep=' ')
bd_full$calving_date1=dmy(bd_full$calving_date)
bd_full$daysinmilking=bd_full$sniffer_date1-bd_full$calving_date1
bd_full$daysinmilking=gsub("[a-z]","",bd_full$daysinmilking)#retain only numbers
bd_full$daysinmilking=as.numeric(bd_full$daysinmilking)
# Check difference of days
table(bd_full$daysinmilking)
# Choose The threshold the difference of days
bd_full1=bd_full %>% filter(daysinmilking <= 365)#336
# Obtain week of lactation 
bd_full1$week_lactation=floor(bd_full1$daysinmilking / 7)
# Codify days in milking as levels 
bd_full1<-mutate(bd_full1,state_lactation=case_when(
  daysinmilking < 91 ~ "1",
  daysinmilking > 90 & daysinmilking < 151 ~ "2",
  daysinmilking >150 ~ "3"))
table(bd_full1$state_lactation)

# Codify number of calvings 
hist(bd_full1$numpar)
table(bd_full1$numpar)
# If there is data of more than 3 calving
bd_full1<-mutate(bd_full1,num_calving=case_when(
  numpar <= 1 ~ "1",
  numpar > 1 & numpar < 3  ~ "2",
  numpar > 2 ~ "3"))
table(bd_full1$num_calving)
########################### Ratio and grams per day ####################################
# Obtain ratio of mean CH4 and mean CO2
bd_full1$ratioCH4CO2=bd_full1$meanCH4/bd_full1$meanCO2
mean(bd_full1$ratioCH4CO2)
# Obtain grams per day 
# Madsen equation
bd_full1$ECM=bd_full1$milk*(0.25+0.122*bd_full1$fat+0.077*bd_full1$protein)
bd_full1$days_inpregnancy=0
bd_full1$gd_madsen<-(0.714*bd_full1$ratioCH4CO2)*180*24*0.001*(5.6*bd_full1$weight**0.75+22*bd_full1$ECM+
                                                                 1.6*0.00001*bd_full1$days_inpregnancy**3)
mean(bd_full1$gd_madsen,na.rm=T)
  
# Tier 2 equation
Cf<-0.386*1.2 #coefficient for feeding situation, lactating
Ca<-0 #Coefficient for activity
Cp<-0.10 #coefficient for pregnancy
DE<-80 #Digestible energy per gross energy in cows fed low quality forage and concentrate

NEm<-Cf*bd_full1$weight**0.75
NEa<-NEm*Ca
NEl<-bd_full1$milk*(1.47+0.4*bd_full1$fat)
NEp<-ifelse(bd_full1$daysinmilking>100, Cp*NEm,0)
REM<-1.123-(0.004092)*DE+0.00001126*(DE)**2-25.4/DE
GE<- 1000*((NEm+NEa+NEl+NEp)/REM)/(DE/100) #g/d
bd_full1$CH4_tier2<-GE*0.065/55.65

mean(bd_full1$CH4_tier2,na.rm=T)
######################################### Distribution of data #######################
summary(bd_full1)
# Plot the phenotypes distribution

# By trait
ggplot(bd_full1, aes(x = meanCH4)) +
  geom_histogram( fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(title = "distribution of phenotypes", x = "meanCH4(ppm)", y = "Frequency")

# All traits 
names(bd_full1)
df_plot <- bd_full1 %>%
  pivot_longer(cols = c(9:23,49,52,53), names_to = "Variable", values_to = "Valor")

ggplot(df_plot, aes(x = Valor)) +
  geom_histogram( fill = "skyblue", color = "black") +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal() +
  labs(title = "distribution of phenotypes", x = "value", y = "Frequency")

##################################### Outlier detection ##################################
ggplot(df_plot, aes(x = Valor)) +
  geom_boxplot(fill = "skyblue", color = "black", outlier.color = "red", outlier.shape = 16) +
  facet_wrap(~ Variable, scales = "free")+
  theme_minimal() +
  labs(title = "Box and whisker plot", x = "Variable", y = "Values")

####################################### Pattern by hour ##################################

# Plot the distribution by hour 
bd_full1$hour=substr(bd_full1$time, 1, 2) #susbtract the hour to plot
bd_full1$hour=as.numeric(bd_full1$hour)
table(bd_full1$hour)

# By individual
animal_9765=subset(bd_full1,bd_full1$cow=="9765")
ggplot(data=animal_9765, aes(x= hour, y=gd_madsen)) +
  geom_line()+geom_point()+
  geom_smooth(method = "loess",  span = 0.2,se=TRUE)+
  facet_wrap(facets = vars(date))+
  ggtitle("9765")

  # All individuals
ggplot(data=bd, aes(x= hour, y=gd_madsen)) +
  # geom_point() +  # Puntos individuales
  geom_line()+    #línea de tendencia
  geom_smooth(method = "loess",  span = 0.3,se=TRUE)+ 
  facet_wrap(facets = vars(cow))+
  ggtitle("Pattern by hour")


################################# DATA FILTERING ########################################
# Threshold by CO2 concentration
bd_full2 <- bd_full1 %>%
  dplyr::mutate(across(c(meanCH4, meanCH4_5s,meanCO2,meanRatioCH4_CO2,AUC_CH4,AUC_Ratio,
                         Sum_of_PeaksCH4,Sum_of_PeaksCH4_5s,Sum_of_PeaksCO2,Sum_of_PeaksRatio,
                         Mean_of_PeaksCH4,Mean_of_PeaksRatio,Sum_MaxPeak, ratioCH4CO2, gd_madsen ,
                         CH4_tier2), ~ ifelse(meanCO2 < 2500, NA, .)))

#View(bd_full1)
colSums(is.na(bd_full1))
colSums(is.na(bd_full2))


############################################ DATA CORRECTION BY SD ##########################
# Function for correction of data ±3 SD
function_outliers <- function(x) {
  mean_x <- mean(x, na.rm = TRUE)
  sd_x <- sd(x, na.rm = TRUE)
  
  # Establecer el límite superior e inferior (± 3 SD)
  lower_lim <- mean_x - 3 * sd_x
  upper_lim <- mean_x + 3 * sd_x
  
  # Reemplazar valores fuera de los límites por NA 
  x <- ifelse(x < lower_lim | x > upper_lim, NA, x)
  
  return(x)
}

data_corrected <- bd_full2 %>%
  dplyr::mutate(across(
    c(meanCH4, meanCH4_5s,meanCO2,meanRatioCH4_CO2,AUC_CH4,AUC_Ratio,
      Sum_of_PeaksCH4,Sum_of_PeaksCH4_5s,Sum_of_PeaksCO2,Sum_of_PeaksRatio,
      Mean_of_PeaksCH4,Mean_of_PeaksRatio,Sum_MaxPeak, ratioCH4CO2, gd_madsen ,
      CH4_tier2), #Apply the function to these variables  
    function_outliers    # Function
  ))
colSums(is.na(data_corrected))

##################################### Check the data after correction ##################### 
summary(data_corrected)
# Plot again to compare
# Plot the distribution again, but without outliers
names(data_corrected)
df_plot2 <- data_corrected %>%
  pivot_longer(cols = c(9:23,49,52,53), names_to = "Variable", values_to = "Valor")

ggplot(df_plot2, aes(x = Valor)) +
  geom_histogram( fill = "skyblue", color = "black") +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal() +
  labs(title = "distribution of phenotypes", x = "value", y = "Frequency")

############################# Check the pattern by hour after correction ####################### 
# By individual
animal_7436=subset(data_corrected,data_corrected$cow=="7436")
ggplot(data=animal_7436, aes(x= hour, y=meanCH4)) +
  geom_line()+
  geom_smooth(method = "loess",  span = 0.3,se=TRUE)+
  ggtitle("pattern by hour")

 # All cows 
ggplot(data=data_corrected, aes(x= hour, y=meanCH4)) +
  geom_line()+
  geom_smooth(method = "loess",  span = 0.3,se=TRUE)+
  facet_wrap(facets = vars(cow))+
  ggtitle("pattern by hour")

############################# Distribution by week of lactation ############################ 
table(data_corrected$week_lactation)

ggplot(data=data_corrected, aes(x= week_lactation, y=meanCH4)) +
  geom_smooth(method = "loess",  span = 0.6,se=TRUE)+
  #facet_wrap(facets = vars(cow))+
  ggtitle("week of lactation")

################################# Plot Correlations between phenotypes ################

names(data_corrected)
correlaciones=data_corrected[,c(1,9,16,21,49,52,53,31,41,42)]
GGally::ggpairs(correlaciones, columns = 2:10,ggplot2::aes(colour="farm"))

#################################### MEAN PER WEEK ##################################
# Obtain week with lubridate
data_corrected$epiweek=epiweek(data_corrected$sniffer_date1)
data_corrected$epiyear=epiyear(data_corrected$sniffer_date1)
data_corrected$robot="1"
data_corrected$farm="1"
mean_week <- data_corrected %>%
  group_by(cow,epiweek,epiyear) %>% #some farms could be studied in two or more periods
  dplyr::summarise(across((where(is.numeric)), list(mean = ~mean(.,na.rm=T)), .names = "{.col}_mean_week"),
                   count = n()) %>%
  dplyr::left_join(data_corrected, by = c("cow","epiweek","epiyear"))#keep the other variables
# Filter by count per week, keep count > 3/5 records per week
table(mean_week$count)
mean_week1=mean_week[mean_week$count > 3, ]

# Choose the variables to keep 
names(mean_week1)
mean_week2=mean_week1[,c(1:3,5:38,86,87,94,95)]
mean_week3=mean_week2[!duplicated(mean_week2),]
#sort the columns
names(mean_week3)
mean_week4=mean_week3[,c(1:3,41,40,39,38,4:30,33,36,37)] #38 rows

########################## Check the number of records by group and animal #####################
# Check the number of animals by comparison group. i.e =herd-season-year or herd-robot-week-year
mean_week5=mean_week5%>%
  group_by(farm,robot,epiweek,epiyear)%>%
  dplyr::mutate(N_group=n())
table(mean_week5$N_group)

  # Check the number of records per animal (to be considered in the model, repeated measurements)
mean_week5=mean_week5%>%
  group_by(cow)%>%
  dplyr::mutate(N_group=n())
table(mean_week5$N_anim) 


##################################### EXTRA #########################################
 # Join columns to create mixed effects
mean_week5=unite(mean_week5,col = "random_effect",c("farm","robot","epiweek","epiyear"),sep="-",remove = F)  # create a new column "random_effect",combining 4 columns
#keeping the original column in the dataset
