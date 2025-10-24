library(tidyverse)
library(arrow)
library(sads)
library(INLA)

all_fileINLAall_files <- list.files('./model_input_data',full.names=T, all.files=T,recursive=T)

hyper2.rw = list(prec = list(prior='pc.prec', param=c(0.3, 0.01))) # medium

# ds.ls <- lapply(all_files,function(X) read_csv(X) )
# 
# ds <- bind_rows(ds.ls)
# write_parquet(ds,'./data/raw.parquet')


ds <- read_parquet('./data/raw.parquet') %>%
  mutate(date = as.Date(paste(year, month,'01',sep='-') ),
         offset1 = pop_total/100000
  )
  
new_rows <- ds %>%
  distinct(fcode) %>%
  mutate(dummy = 1) %>%
  crossing(date = new_months) %>%
  select(-dummy)

# Get the latest offset1 value per fcode
latest_offsets <- ds %>%
  group_by(fcode) %>%
  filter(date == max(date)) %>%
  dplyr::select(fcode, offset1)


ds2 <- new_rows %>%
  left_join(latest_offsets, by = "fcode") %>%
  bind_rows(ds, .) %>%
  arrange(fcode, date) %>%
  dplyr::select(fcode, date, obs_dengue_cases,offset1)


## SEASONAL DECOMPOSITION IN INLA TO GET THE BASELINE
ts_decomposition_inla <- function(forecast_year, fcode.select, dsname=ds2){
  set.seed(8123)  # Global R seed for reproducibility
  
  new_months <- seq.Date(as.Date("2025-10-01"), as.Date("2025-12-01"), by = "month")
  
 
  c1 <- dsname %>%
    filter( date>='2004-09-01')%>%
    arrange(fcode, date) %>%
    mutate( t = lubridate::interval(min(date), date) %/% months(1) + 1) %>%
    group_by(fcode) %>%
    mutate(fcode2=fcode,
        
           month=as.factor(month(date)),
           monthN=month(date),
           #log_offset=log(pop_total/100000)
    ) %>%
    ungroup() %>%
    mutate(
      fcodeID2 = fcode,
      fcodeID3 = fcode,
      fcodeID4 = fcode,
      t = t - min(t, na.rm = TRUE) + 1, #make sure timeID starts at 1
      
      time_id1= t , 
      time_id2=t,
      time_id3= t,
      time_id4= t,
      year=year(date),
      yearN= as.numeric(as.factor(year))) %>%
    arrange(date,fcode) %>% #SORT FOR SPACE_TIME
    mutate(fcodeIDpad=str_pad(fcode, 3, pad = "0", side='left'),
           timeIDpad=str_pad(time_id1, 5, pad = "0", side='left')
           
    )
  
  c2 <- c1 %>%
    filter(fcode==fcode.select & year <=forecast_year) %>%
    mutate(obs_dengue_cases_fit = if_else(year>=forecast_year, NA_real_, obs_dengue_cases))
  
  offset1 <- c2$offset1
  
  #single fcode version
  form2 <- as.formula( 'obs_dengue_cases_fit ~ 1+
            f(time_id1, model="ar1",constr=TRUE) +
            f(yearN, model="rw2", constr=TRUE)+
            f(monthN, model="rw1", hyper=hyper2.rw, cyclic=TRUE, scale.model=TRUE, constr=TRUE)'
  )
  
  n_times <- length(unique(c2$time_id1))
  n_years <- length(unique(c2$year))
  
  pred.time_id1 <- c2$time_id1[(n_times-11):n_times]
  pred.yearN <- c2$yearN[(n_times-11):n_times]
  pred.monthN <- c2$monthN[(n_times-11):n_times]
  
  time.mat <- model.matrix(~ -1 + as.factor(time_id1) , data=c2)
  year.mat <- model.matrix(~ -1 + as.factor(yearN) , data=c2)
  month.mat <- model.matrix(~ -1 + as.factor(monthN) , data=c2)
  
  #confirm linear combs gives same results as summary.linear.predictor
  # lc.lin.pred = inla.make.lincombs("(Intercept)"=rep(1,12),
  #                              'time_id1' = time.mat[(n_times-11):n_times,], 
  #                              'yearN' = year.mat[(n_times-11):n_times,],
  #                              'monthN'=month.mat[(n_times-11):n_times,])
  
  lc.lin.pred.no.ar1 = inla.make.lincombs("(Intercept)"=rep(1,nrow(c2)),
                                          'time_id1' = rep(0,nrow(c2)), 
                                          'yearN' = year.mat,
                                          'monthN'=month.mat)
  
  mod1 <- inla(form2, data = c2,  family = "poisson",E=offset1,
               lincomb = lc.lin.pred.no.ar1,
               control.compute = list(dic = FALSE, 
                                      waic = FALSE, 
                                      config = T,
                                      return.marginals=F
               ),
               # save predicted values on response scale
               control.predictor = list(compute=TRUE, link=1),
               control.fixed = list(mean.intercept=0, 
                                    prec.intercept=1e-4, # precision 1
                                    mean=0, 
                                    prec=1), # weakly regularising on fixed effects (sd of 1)
               inla.mode = "experimental", # new version of INLA algorithm (requires R 4.1 and INLA testing version)
               num.threads=8
  )    
  mod.family <- mod1$.args$family
  
 # samps = inla.posterior.sample(n=10000,mod1)
  
  #need the linear predictor minus the AR(1) using lincomb; ar(1) random effect should be uncorrelated from rest of model;
  lin.pred.no.ar1 <- mod1$summary.lincomb.derived[,c('mean','sd')] %>%
    dplyr::rename(mean.no.ar1_inc=mean, sd.no.ar1_inc=sd) %>%
    mutate(
          mean_log_baseline = mean.no.ar1_inc + log(offset1) , #lambda, mean of cases
           sd_log_baseline = sd.no.ar1_inc # SD unchaned because just adding a constant to mean
           )
  
  ##dataset with the baseline expected values!
  baseline <- cbind.data.frame('date'=c2$date ,'obs_dengue_cases'=c2$obs_dengue_cases,c2$obs_dengue_cases, lin.pred.no.ar1) %>%
  
    mutate(year=lubridate::year(date),
           fcode=fcode.select) %>%
    filter(year==forecast_year) %>%
    saveRDS( paste0('./output/baseline_', forecast_year,'_',fcode.select,'.rds'))
  
  # lin.pred <- mod1$summary.linear.predictor[,c('mean','sd')]%>%
  #   rename(mean.lin.pred=mean, sd.lin.pred=sd) 
  # 
  # comb.pred <- cbind.data.frame('date'=c2$date ,lin.pred.no.ar1, lin.pred) %>%
  #   mutate(t=row_number())
  # 
  # ggplot(comb.pred, aes(x=date, y=mean.no.ar1))+
  #   geom_line(color='red')+
  #   geom_line(aes(x=date, y=mean.lin.pred))
}



##calculate baseline for 2025

for(i in 2025:max(ds$year)){
  for(j in unique(ds$fcode)){
    print(i)
    print(j)
    ts_decomposition_inla(forecast_year=i, fcode.select=j)
  }
}


##Check baseline
ds.out <- readRDS('./output/baseline_2025_ED_AN_GIANG_AN_PHU_DISTRICT.rds') %>%
  mutate(base_cases = exp(mean_log_baseline),
          base_case_2sd =  exp(mean_log_baseline+2*sd_log_baseline), #ucl based on uncertainty intevral
         
         base_cases_ucl_v2 = base_case_2sd + 2*sqrt(base_case_2sd) #uncertainty based on Poisson 95% CI
          )

p1 <- ggplot(ds.out)+
  geom_line(aes(x=date,y=base_cases))+
  geom_line(aes(x=date,y=base_case_2sd), lty=2 ,color='gray')+
  geom_line(aes(x=date,y=base_cases_ucl_v2), lty=2)+
  
        geom_line(aes(x=date, y=obs_dengue_cases), color='red') + 
  theme_classic()


plotly::ggplotly(p1)
p1 +ggtitle('ED_AN_GIANG_AN_PHU_DISTRICT')
#base_case_2sd =  exp(mean_log_baseline)*exp(2*sd_log_baseline)*offset,
            