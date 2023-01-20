## Title: Multiscale spatially varying coefficient modelling using a Geographical Gaussian Process GAM
# Alexis Comber, Paul Harris and Chris Brunsdon
# If you have any problems with data / code / versions etc please contact Lex Comber a.comber@leeds.ac.uk

## 1. load packages
library(tidyverse)
library(cowplot)
library(cols4all)
library(mgcv)
library(GWmodel)
library(parlitools)
library(broom)
library(spdep)
library(stringr)

## 2. Create simulated data in the manner of Fotheringham et al 2017 Part I
# 25 * 25 lattice 
# v = vertical, coordinate u = horizontal coordinate both with increment of 1
# b0 = 3		# zero heterogeneity
# b1 = 1 + ( (1/12) * (u + v) )		# low heterogeneity
# b2 = 1 + ( (1/324) * (36-( 6- (u/2))^2) * (36-( 6- (v/2))^2) )		# high heterogeneity
gr = expand.grid(u = 1:25, v = 25:1)
gr$b0 = 3
gr$b1 = 1 + ( (1/12) * (gr$u + gr$v) )	
# gr$b12 = 1 + ( (1/12) * (gr$u + rev(gr$v)) )
gr$b2 = 1 + ( (1/324) * (36-( 6- (gr$v/2))^2) * (36-( 6- (gr$u/2))^2) )	
# head(gr)

## 3. Create Figure 1
tit =expression(paste(""*beta[0]*""))
p1 = 
  ggplot(gr, aes(x = u, y = v, fill = b0)) +
  geom_raster()+
  coord_fixed() +
  scale_fill_continuous_c4a_seq(palette="scico.lajolla",  limits = c(0, 5.4)) + 
  theme_bw() + xlab("") + ylab("") +  
  ggtitle(tit) +
  theme(legend.position = "none", axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), 
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) 

tit =expression(paste(""*beta[1]*""))
p2 = 
  ggplot(gr, aes(x = u, y = v, fill = b1)) +
  geom_raster()+
  coord_fixed() +
  scale_fill_continuous_c4a_seq(palette="scico.lajolla", limits = c(0, 5.4)) + 
  theme_bw() + xlab("") + ylab("") +
  ggtitle(tit) +
  theme(legend.position = "none",axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), 
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) 

tit =expression(paste(""*beta[2]*""))
p3 = 
  ggplot(gr, aes(x = u, y = v, fill = b2)) +
  geom_raster()+
  coord_fixed() +
  scale_fill_continuous_c4a_seq(palette="scico.lajolla", name = "", limits = c(0, 5.4)) + 
  theme_bw() + xlab("") + ylab("") +
  ggtitle(tit) +
  theme(legend.position = "bottom", legend.key.width = unit(2.1, "cm"),axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), 
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())

legend_betas <- get_legend(p3)
p3 = p3 + theme(legend.position = "none")

# open plot window and plot 
if (.Platform$GUI == "AQUA") {
		quartz(w=10,h=5.5) } else  {
			x11(w=10,h=5.5) } 
plot_grid(p1, p2, p3, NULL, legend_betas, NULL, nrow = 2, ncol = 3, rel_heights = c(1, 0.2))

## 4. Create simulated data in the manner of Fotheringham et al 2017 Part II
# x1 and x2 were generated randomly from a normal distribution N(0,1) # error term was generated from a normal distributed E~ N(0, 0.5^2) # variable y was computed# done 100 times

MakeSim = F
if(MakeSim) {
  ## Create Simulations
  make_x = function(){
  	x = rnorm(nrow(gr)) 
  	x = (x - min(x))/(max(x) - min(x))
  	x
  }	
  y_sim = x1_sim = x2_sim = e_sim = vector()
  for(i in 1:100){
  	x1 = make_x()
  	x2 = make_x()
  	e = make_x()*.25
  	y.i = gr$b0 + (gr$b1*x1) + (gr$b2*x2) + e
  	y_sim = cbind(y_sim, y.i)
  	x1_sim= cbind(x1_sim, x1)
  	x2_sim= cbind(x2_sim, x2)
  	e_sim= cbind(e_sim, e)
  }
  save(list = c("y_sim", "x1_sim", "x2_sim", "e_sim"), file = "sim1_fullpaper.RData")
} else {
  load("sim1_fullpaper.RData")
}

## 5. Evaluate the simulated data with GGP-GAM and MGWR
# some fit functions
rsq <- function (x, y) cor(x, y) ^ 2
rmse <- function(x, y) sqrt(mean((x - y)^2))
mae <- function(x, y){
	n = length(x) 
	sum = 0
	for (i in 1:n){
		sum = abs(x[i] - y[i]) + sum
	}
	sum/n
}
# Big loop - this takes a few hours to run because of the MGWR
EvalSim = F
if(EvalSim) { 
  # define some outputs
  gam_res_tab = mgwr_res_tab = vector()
  gb0.all = gb1.all = gb2.all = vector()
  mb0.all = mb1.all = mb2.all = vector()
  
  for (i in 1:100) {
    ## 1. set up the data
    df.i = data.frame(y = y_sim[,i],
                      x0 = 3,
                      x1 = x1_sim[,i],
                      x2 = x2_sim[,i],
                      #e = e_sim[,i],
                      Intercept = 1, 
                      X = gr$u, Y = gr$v)
    ## 2. GGP-GAM: make the model
    gam.i <- gam(y~ 0 + 
                Intercept + s(X,Y,bs='gp',by=Intercept) + 
                x1 + s(X,Y,bs='gp',by=x1) +
                x2 + s(X,Y,bs='gp',by=x2), data = df.i)             
    ## 3. GGP-GAM SVCs: setting each beta to 1 others to 0
    df.ii = df.i
    b0 <- df.i %>% mutate(x1 = 0, x2 = 0, Intercept = 1)
    df.ii <- df.ii %>% mutate(se0 = predict(gam.i, se = TRUE, newdata = b0)$se.fit,
                              b0=predict(gam.i,newdata=b0))
    b1 <- df.i %>% mutate(x1 = 1, x2 = 0, Intercept = 0)
    df.ii <- df.ii %>% mutate(se1 = predict(gam.i, se = TRUE, newdata = b1)$se.fit,
                              b1=predict(gam.i,newdata=b1))
    b2 <- df.i %>% mutate(x1 = 0, x2 = 1, Intercept = 0)
    df.ii <- df.ii %>% mutate(se2 = predict(gam.i, se = TRUE, newdata = b2)$se.fit,
                              b2=predict(gam.i,newdata=b2))
    ## 4. GGP-GAM: evaluate betas prediction accuracy
    rsq.i = c(gam.i$aic, rsq(gr$b1, df.ii$b1), rsq(gr$b2, df.ii$b2))
    rmse.i = c(rmse(gr$b0, df.ii$b0), rmse(gr$b1, df.ii$b1), rmse(gr$b2, df.ii$b2))
    mae.i = c(mae(gr$b0, df.ii$b0), mae(gr$b1, df.ii$b1), mae(gr$b2, df.ii$b2))
    gam_res = c(rsq.i, rmse.i, as.vector(mae.i))
    gam_res_tab = rbind(gam_res_tab, gam_res)
    # save coefficients
    gb0.all = cbind(gb0.all, df.ii$b0)
    gb1.all = cbind(gb1.all, df.ii$b1)
    gb2.all = cbind(gb2.all, df.ii$b2)
    ## 5. MGWR
    # make some sp data
    df.i.sp = SpatialPointsDataFrame(coords = df.i[, c("X", "Y")], 
    		data = data.frame(df.i[, c("y", "x1", "x2")]))
    # determine MGWR bandwidths
    bws.i = gwr.multiscale(y~x1+x2,  
  	                        data = df.i.sp,
  	                        adaptive = T, max.iterations = 1000,
  	                        criterion="CVR",
  	                        kernel = "bisquare",
  	                        bws0=c(10,10,10),
  	                        verbose = F, predictor.centered=rep(T, 3))
  	# run with bandwidths to confirm
    msgwr_bws.i = as.vector(bws.i$GW.arguments$bws)
    bws.i = gwr.multiscale(y~x1+x2,  
  	                        data = df.i.sp,
  	                        adaptive = T, max.iterations = 1000,
  	                        criterion="CVR",
  	                        kernel = "bisquare",
                            bws0 = msgwr_bws.i,
                            bw.seled=rep(T, 3),
  	                        verbose = F, predictor.centered=rep(T, 3))
    ## 6. MGWR: evaluate betas prediction accuracy
    df.ii = bws.i$SDF@data[,1:3]
  	names(df.ii) = c("b0", "b1", "b2")
    rsq.i = c(bws.i$GW.diagnostic$AIC, rsq(gr$b1, df.ii$b1), rsq(gr$b2, df.ii$b2))
    rmse.i = c(rmse(gr$b0, df.ii$b0), rmse(gr$b1, df.ii$b1), rmse(gr$b2, df.ii$b2))
    mae.i = c(mae(gr$b0, df.ii$b0), mae(gr$b1, df.ii$b1), mae(gr$b2, df.ii$b2))
  	mgwr_res = c(rsq.i, rmse.i, as.vector(mae.i))
    mgwr_res_tab = rbind(mgwr_res_tab, mgwr_res)
    # save coefficients
    mb0.all = cbind(mb0.all, df.ii$b0)
    mb1.all = cbind(mb1.all, df.ii$b1)
    mb2.all = cbind(mb2.all, df.ii$b2)
    #save GGP-GAM SPs and MGWR BWS
    bws_sp = rbind(bws_sp, c(as.vector(gam.i$sp), bws.i$GW.arguments$bws) )
    
    if(i %% 10 == 0) cat(i, "\t")
    
  }
  save(list = c("gam_res_tab", "gb0.all", "gb1.all", "gb2.all",
                "mgwr_res_tab", "mb0.all", "mb1.all", "mb2.all",
                "bws_sp"), file = "simulation_result_fullpaper.RData")
} else {
  load("simulation_result_fullpaper.RData")
}

## 6. Create Figure 2
# name the results 
colnames(gam_res_tab) = c("AIC", "Rsq_B1", "Rsq_B2", "RMSE_B0", "RMSE_B1", "RMSE_B2", "MAE_B0", "MAE_B1", "MAE_B2")
colnames(mgwr_res_tab) = c("AIC", "Rsq_B1", "Rsq_B2", "RMSE_B0", "RMSE_B1", "RMSE_B2", "MAE_B0", "MAE_B1", "MAE_B2")
# combine to a dingle data frame
gam_res_tab = data.frame(gam_res_tab, res = "GAM")
mgwr_res_tab = data.frame(mgwr_res_tab, res = "MGWR")
df = rbind(gam_res_tab, mgwr_res_tab)
rownames(df) = 1:nrow(df)
# lengthen
df %>% 
  pivot_longer(-res) -> df
# rename
df$name = recode(df$name,
  "MAE_B0" = '"MAE "*beta[0]',
  "MAE_B1" = '"MAE "*beta[1]',
  "MAE_B2" = '"MAE "*beta[2]',
  "RMSE_B0" = '"RMSE "*beta[0]',
  "RMSE_B1" = '"RMSE "*beta[1]',
  "RMSE_B2" = '"RMSE "*beta[2]',
  "Rsq_B1"= 'R^2*" "*beta[1]',
  "Rsq_B2"= 'R^2*" "*beta[2]')
p4 = 
  df %>%   
  ggplot(aes(y = value,  fill = res)) +
  geom_boxplot() +
  scale_fill_discrete_c4a_cat(palette="tableau.classic_color_blind", name = "SVC model") + 
  facet_wrap(~name, scales = "free", labeller=label_parsed, dir = "h", as.table = T) +theme_bw()+
  theme(#legend.position = c(0.85, 0.15),
       legend.position = "bottom", 
        legend.key.size = unit(1.4, "cm"),
        legend.text=element_text(size=10))
# open plot window and plot 
if (.Platform$GUI == "AQUA") {
		quartz(w=10,h=8) } else  {
			x11(w=10,h=8) } 
p4

## 7. Create Figure 3
# using the 4th simulation run
gr = cbind(gr, gam_b0 = gb0.all[,4])
gr = cbind(gr, gam_b1 = gb1.all[,4])
gr = cbind(gr, gam_b2 = gb2.all[,4])
gr = cbind(gr, mgw_b0 = mb0.all[,4])
gr = cbind(gr, mgw_b1 = mb1.all[,4])
gr = cbind(gr, mgw_b2 = mb2.all[,4])
# head(gr)
tit =expression(paste("GAM "*beta[0]*""))
p5 = 
  ggplot(gr, aes(x = u, y = v, fill = gam_b0)) +
  geom_raster()+
  coord_fixed() +
  scale_fill_continuous_c4a_seq(palette="scico.lajolla", limits = c(0, 5.4)) + 
  theme_bw() + xlab("") + ylab("") +  
  ggtitle(tit) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), 
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
tit =expression(paste("GAM "*beta[1]*""))
p6 = 
  ggplot(gr, aes(x = u, y = v, fill = gam_b1)) +
  geom_raster()+
  coord_fixed() +
  scale_fill_continuous_c4a_seq(palette="scico.lajolla", limits = c(0, 5.4)) + 
  theme_bw() + xlab("") + ylab("") +  
  ggtitle(tit) +
  theme(legend.position = "none",axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), 
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) 
tit =expression(paste("GAM "*beta[2]*""))
p7 = 
  ggplot(gr, aes(x = u, y = v, fill = gam_b2)) +
  geom_raster()+
  coord_fixed() +
  scale_fill_continuous_c4a_seq(palette="scico.lajolla", limits = c(0, 5.4)) + 
  theme_bw() + xlab("") + ylab("") +  
  ggtitle(tit) +
  theme(legend.position = "none",axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), 
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) 

tit =expression(paste("MGWR "*beta[0]*""))
p8 = 
  ggplot(gr, aes(x = u, y = v, fill = mgw_b0)) +
  geom_raster()+
  coord_fixed() +
  scale_fill_continuous_c4a_seq(palette="scico.lajolla", limits = c(0, 5.4)) + 
  theme_bw() + xlab("") + ylab("") +  
  ggtitle(tit) +
  theme(legend.position = "none",axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), 
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) 
tit =expression(paste("MGWR "*beta[1]*""))
p9 = 
  ggplot(gr, aes(x = u, y = v, fill = mgw_b1)) +
  geom_raster()+
  coord_fixed() +
  scale_fill_continuous_c4a_seq(palette="scico.lajolla", limits = c(0, 5.4)) + 
  theme_bw() + xlab("") + ylab("") +  
  ggtitle(tit) +
  theme(legend.position = "none",axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), 
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) 
tit =expression(paste("MGWR "*beta[2]*""))
p10 = 
  ggplot(gr, aes(x = u, y = v, fill = mgw_b2)) +
  geom_raster()+
  coord_fixed() +
  scale_fill_continuous_c4a_seq(palette="scico.lajolla", limits = c(0, 5.4), name = "") + 
  theme_bw() + xlab("") + ylab("") +  
  ggtitle(tit) +
  theme(legend.position = "bottom", legend.key.width = unit(2.1, "cm"),axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), 
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
legend_betas2 <- get_legend(p10)
p10 = p10 + theme(legend.position = "none")
if (.Platform$GUI == "AQUA") {
		quartz(w=10,h=8) } else  {
			x11(w=10,h=8) } 
plot_grid(p5, p6, p7, p8, p9, p10, NULL, 
	legend_betas2, NULL, nrow = 3, ncol = 3, rel_heights = c(1, 1, 0.2))

## 8. Create Table 1 
bws_sp = data.frame(bws_sp)
colnames(bws_sp) = c("gam_b0", "gam_b1", "gam_b2", "mgw_b0", "mgw_b1", "mgw_b2")
bws_sp %>% 
  mutate(ID = 1:n()) %>%
  pivot_longer(-ID) -> df

# extract the GAM SP and MGWR BW summaries
tab = t(apply(bws_sp, 2, summary))
tab = data.frame(tab)
rownames(tab) = c("GAM SP $x$~0~", "GAM SP $x$~1~", "GAM SP $x$~2~", "MGWR BW $x$~0~", "MGWR BW $x$~1~", "MGWR BW $x$~2~")
# format the GAM SPs to scientific notation
tab[1:3,] = format(tab[1:3,], digits = 2, scientific = T, big.mark = ",")
# for the mean
tab = tab[,-4]
colnames(tab) = c("Min.", "Q1", "Median", "Q3", "Max.")
knitr::kable(tab, digits = 0, caption = "\\label{tab:tab1}Sumamries of MGWR bandwidths (BW) and GAM GP spline smoothing paramters (SP).", linesep = "") 


## 8. Brexit case study  
# Merge the data for analysis
# Focus on Great Britain,  not NI
NI <- str_detect(west_hex_map$gss_code,'^N')
n <- sum(!NI) ## Cheeky recast of logicals to numeric 0/1
west_hex_gb <- west_hex_map[!NI,] ## Subset of map
st_crs(west_hex_gb) = st_crs(west_hex_gb)
# Join the other data tables
yabda <- west_hex_gb %>% rename(ons_const_id=gss_code) %>% 
  left_join(bes_2015,by = 'ons_const_id') %>%
  left_join(leave_votes_west, by = 'ons_const_id') %>% 
  left_join(census_11, by = 'ons_const_id') 
## Select variables for analysis from Beecham et al 
# Beecham, R., Williams, N. and Comber, A., 2020. Regionally-structured explanations behind area-level populism: An update to recent ecological analyses. Plos one, 15(3), p.e0229974.
# % in leisure & hospitality jobs net 
# % change in manufacturing jobs 
# % in transport, trade & utilities jobs
# % with Bachelors level or higher median household income
# % not in good health | living in poverty
# % not UK-born | not US-born 
# % age 65+ years
yabda %>% 
  mutate(over65 = age_65_to_74+ age_75_to_84	+ age_85_to_89 + age_90_plus, 
         ttu = nssecsemi_routine + nssecroutine + nsseclower_supervisor, 
         lhosp = industry_accommodation,
         manu = industry_manufacturing,
         degree = qual_level_4,
         badhealth = health_bad + health_very_bad,
         bornuk = born_uk,
         leave = figure_to_use *100, 
         to15 = turnout_15) %>% 
  select(ons_const_id, leave, to15, over65,  ttu, lhosp, manu, degree, badhealth, bornuk, region_name) -> hex.sp  
hex.sp$region_name = gsub(" Euro Region", "", hex.sp$region_name) 
hex.sp$region_name = gsub(" and the Humber", "", hex.sp$region_name) 
hex.sp$region_name = gsub("Eastern", "East", hex.sp$region_name) 
# create region layer for context / interpretation
hex.sp %>% 
    group_by(region_name) %>%
    summarise(geometry = sf::st_union(geometry)) %>%
    ungroup() -> regions.sp
hex.sp %>% select(-region_name) -> hex.sp

## 9. Create Table 2
df = data.frame(Variable = names(hex.sp)[-c(1,11)],
                Description = c("Leave vote (%)",
                                "Turn out in the 2015 general election (%)",
                                "Over 65s (%)",
                                "Transport trade and utilities employment (%)",
                                "Leisure & hospitality employment (%)",
                                "Manufacturing employment (%)",
                                "Degree (level 4) qualification (%)",
                                "Bad health (%)",
                                "UK born (%)"))
  
knitr::kable(df, caption = "\\label{tab:tab2}The variables used to construct the GAM with SP SVC regression model of Leave voting rates, for each parliamentary consituency.", linesep = "",
      booktabs = T)

## 10. Create Figure 4
p1 = 
  ggplot() +
  geom_sf(data = hex.sp,aes(fill=(leave-50)), col = NA) + 
  scale_fill_continuous_c4a_div(palette="tol.sunset",name='Leave (%) - 50%' ) + 
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=8))
p1a = 
  ggplot() +
  geom_sf(data = regions.sp, aes(fill = region_name), col = NA) +
  scale_fill_discrete_c4a_cat(palette="brewer.paired",name='Regions') + 
  theme_bw() +
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=8))
if (.Platform$GUI == "AQUA") {
		quartz(w=10,h=6) } else  {
			x11(w=10,h=6) } 
plot_grid(p1, p1a, nrow = 1)

## 11. Create Figure 5
plot_var_func = function(var.name, tit) {
  ggplot(hex.sp, aes_string(fill=var.name)) + 
  geom_sf(col = NA) + 
  scale_fill_continuous_c4a_seq(palette="scico.tokyo", name= tit ) + 
  theme_bw() +
  theme(legend.position = "bottom", 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text=element_text(size=8))   
}
p2 = plot_var_func("to15", "to15")
p3 = plot_var_func("over65", "over65")
p4 = plot_var_func("ttu", "ttu")
p5 = plot_var_func("lhosp", "lhosp")
p6 = plot_var_func("manu", "manu")
p7 = plot_var_func("degree", "degree")
p8 = plot_var_func("badhealth", "badhealth")
p9 = plot_var_func("bornuk", "bornuk")
if (.Platform$GUI == "AQUA") {
		quartz(w=12,h=8) } else  {
			x11(w=12,h=8) } 
plot_grid(p2, p3, p4, p5,p6,p7,p8,p9, ncol = 4)

## 12. Create Table 3
m = lm(leave~., data = hex.sp %>% st_drop_geometry() %>% select(-ons_const_id))
# summary(m)
m.pred = predict(m, newdata = hex.sp %>% st_drop_geometry() %>% select(-ons_const_id))    
r2 = round(summary(m)$r.squared, 3)
aic = round(AIC(m),1)
rm.se = round(rmse(hex.sp$leave, m.pred),3)
knitr::kable(tidy(m), booktabs = T, digits = 3,
      caption = paste0("\\label{tab:tab3}The results of standard OLS regression of the Brexit leave vote (RMSE = ", rm.se, "; $R^{2}$ = ",r2, "; AIC = ", aic, ")."))

## 13. Create GGP-GAM
hex.gb = st_transform(hex.sp, 27700)
coords = st_coordinates(st_centroid(hex.gb))
df.gam = 
  hex.sp %>% 
  mutate(Intercept = 1, 
         X = coords[,"X"]/1000, 
         Y = coords[,"Y"]/1000) %>%
  st_drop_geometry() %>%
  as_tibble()
gam.m = gam(leave~ 0 + 
              Intercept + s(X,Y,bs='gp',by=Intercept) + 
              to15 + s(X,Y,bs='gp',by=to15) +
              over65 + s(X,Y,bs='gp',by=over65) +
              ttu + s(X,Y,bs='gp',by=ttu) +
              lhosp + s(X,Y,bs='gp',by=lhosp) +
              manu + s(X,Y,bs='gp',by=manu) +
              degree + s(X,Y,bs='gp',by=degree) +
              badhealth + s(X,Y,bs='gp',by=badhealth) + 
              bornuk + s(X,Y,bs='gp',by=bornuk), 
            data = df.gam)            
# GGP-GAM SVCs: setting each beta to 1 others to 0
df.ii = df.gam %>% select(-ons_const_id)
b0 <- df.gam %>% mutate(to15 = 0, over65 = 0, ttu = 0, lhosp = 0,
	manu = 0, degree = 0, badhealth = 0, bornuk = 0, Intercept = 1)
df.ii <- df.ii %>% mutate(se0 = predict(gam.m, se = TRUE, newdata = b0)$se.fit,
	b0=predict(gam.m,newdata=b0))
bto15 <- df.gam %>% mutate(to15 = 1, over65 = 0, ttu = 0, lhosp = 0, 
	manu = 0, degree = 0, badhealth = 0, bornuk = 0, Intercept = 0)
df.ii <- df.ii %>% mutate(seto15 = predict(gam.m, se = TRUE, newdata = bto15)$se.fit,
	bto15=predict(gam.m,newdata=bto15))
bover65 <- df.gam %>% mutate(to15 = 0, over65 = 1, ttu = 0, lhosp = 0, 
	manu = 0, degree = 0, badhealth = 0, bornuk = 0, Intercept = 0)
df.ii <- df.ii %>% mutate(seover65 = predict(gam.m, se = TRUE, newdata = bover65)$se.fit,
	bover65=predict(gam.m,newdata=bover65))
bttu <- df.gam %>% mutate(to15 = 0, over65 = 0, ttu = 1, lhosp = 0, 
	manu = 0, degree = 0, badhealth = 0, bornuk = 0, Intercept = 0)
df.ii <- df.ii %>% mutate(settu = predict(gam.m, se = TRUE, newdata = bttu)$se.fit,
	bttu=predict(gam.m,newdata=bttu))
blhosp <- df.gam %>% mutate(to15 = 0, over65 = 0, ttu = 0, lhosp = 1, 
	manu = 0, degree = 0, badhealth = 0, bornuk = 0, Intercept = 0)
df.ii <- df.ii %>% mutate(selhosp= predict(gam.m, se = TRUE, newdata = blhosp)$se.fit,
	blhosp=predict(gam.m,newdata=blhosp))
bmanu <- df.gam %>% mutate(to15 = 0, over65 = 0, ttu = 0, lhosp = 0, 
	manu = 1, degree = 0, badhealth = 0, bornuk = 0, Intercept = 0)
df.ii <- df.ii %>% mutate(semanu= predict(gam.m, se = TRUE, newdata = bmanu)$se.fit,
	bmanu=predict(gam.m,newdata=bmanu))
bdegree <- df.gam %>% mutate(to15 = 0, over65 = 0, ttu = 0, lhosp = 0, 
	manu = 0, degree = 1, badhealth = 0, bornuk = 0, Intercept = 0)
df.ii <- df.ii %>% mutate(sedegree= predict(gam.m, se = TRUE, newdata = bdegree)$se.fit,
	bdegree=predict(gam.m,newdata=bdegree))
bbadhealth <- df.gam %>% mutate(to15 = 0, over65 = 0, ttu = 0, lhosp = 0,
	manu = 0, degree = 0, badhealth = 1, bornuk = 0, Intercept = 0)
df.ii <- df.ii %>% mutate(sebadhealth= predict(gam.m, se = TRUE, newdata = bbadhealth)$se.fit,
	bbadhealth=predict(gam.m,newdata=bbadhealth))
bbornuk <- df.gam %>% mutate(to15 = 0, over65 = 0, ttu = 0, lhosp = 0, 
	manu = 0, degree = 0, badhealth = 0, bornuk = 1, Intercept = 0)
df.ii <- df.ii %>% mutate(sebornuk= predict(gam.m, se = TRUE, newdata = bbornuk)$se.fit,
	bbornuk=predict(gam.m,newdata=bbornuk))
# create some fit measures
gam_hex.sp = df.ii
gam.pred = predict(gam.m, newdata = gam_hex.sp)    
r2 = round(rsq(hex.sp$leave, gam.pred),3)
rm.se = round(rmse(hex.sp$leave, gam.pred),3)
#mae(hex.sp$leave, gam.pred)
aic = round(AIC(gam.m), 1)
st_geometry(gam_hex.sp) <- st_geometry(hex.sp)

## 14. Create Tables 4 and 5
tab = tidy(gam.m, digits = 3) 
tab2 <- gam_hex.sp %>% st_drop_geometry() %>% 
  select(b0, bto15, bover65, bttu, blhosp, bmanu, bdegree, bbadhealth, bbornuk)
tab2 = round(apply(tab2,2,summary),3)
colnames(tab2) = c("Intercept", names(hex.sp)[3:10])
tab2 = t(round(tab2, 2))
knitr::kable(tab, booktabs = T, digits = 3, row.names = F, linesep = "",
      caption = paste0("\\label{tab:tab4} The smooth terms of the GGP-GAM model."))
knitr::kable(tab2, booktabs = T, digits = 3, row.names = T, linesep = "",
      caption = paste0("\\label{tab:tab5}The distributions of the GGP-GAM spatially varying coefficient estimates (RMSE = ", rm.se, "; $R^{2}$ = ",r2, "; AIC = ", aic, ")."))

## 15. Create Figure 6 
plot_vgam_coef_func = function(var.name = "b0", tit) {
  var  = gam_hex.sp %>% 
    st_drop_geometry() %>% 
    dplyr::select(all_of(var.name)) %>%
    unlist() %>% 
    as.vector()
  if (sign(max(var)) * sign(min(var)) == 1) flip = F
  if (sign(max(var)) * sign(min(var)) == -1) flip = T

  if(flip) {
    ggplot(gam_hex.sp, aes_string(fill=var.name)) + 
    geom_sf(col = NA) + 
    scale_fill_continuous_c4a_div(palette="scico.vik", name= tit, mid = 0, reverse = T ) + 
    theme_bw() +
    theme(legend.position = "bottom", 
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), 
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          text=element_text(size=8))   
  } else {
    ggplot(gam_hex.sp, aes_string(fill=var.name)) + 
    geom_sf(col = NA) + 
    scale_fill_continuous_c4a_seq(palette="scico.nuuk", name= tit) + 
    theme_bw() +
    theme(legend.position = "bottom", 
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), 
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          text=element_text(size=8))
  }
}

p1 = plot_vgam_coef_func("b0", tit = "Intercept")
p2 = plot_vgam_coef_func("bto15", tit = "to15")
p3 = plot_vgam_coef_func("bover65", "over65")
p4 = plot_vgam_coef_func("bttu", "ttu")
p5 = plot_vgam_coef_func("blhosp", "lhosp")
p6 = plot_vgam_coef_func("bmanu", "manu")
p7 = plot_vgam_coef_func("bdegree", "degree")
p8 = plot_vgam_coef_func("bbadhealth", "badhealth")
p9 = plot_vgam_coef_func("bbornuk", "bornuk")
if (.Platform$GUI == "AQUA") {
		quartz(w=9,h=11) } else  {
			x11(w=9,h=11) } 
plot_grid(p1, p2, p3, p4, p5,p6,p7,p8,p9, ncol = 3)

### END ###
