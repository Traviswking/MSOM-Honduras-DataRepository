library(camtrapR)
library(purrr)
library(DT)
library(knitr)
library(ggplot2)
library(sf)
library(rlist)
library(rjags)
library(coda)
library(gridExtra)

mod.jagsT18_2 <- readRDS("mod.jagsT18_BPV.rds")

############ import in a saved text file of the model structure from the tmp model text file on main computer 

modelFile = 'MSOM_T18_ModelText_2.txt'

params = c("mu.psi", "mu.p", mod.jagsT18_2@params)




n.iter = 60000
n.burnin = 30000
thin = 50
chains = 5



########

mod <- rjags::jags.model(file =  modelFile, 
                         data =  mod.jagsT18_2@data, 
                         inits =  mod.jagsT18_2@inits_fun(),
                         n.chain=chains, 
                         n.adapt=0,
                         quiet = TRUE)

out <- rjags::coda.samples(model = mod,
                           variable.names =  params, 
                           n.iter	= n.iter, 
                           thin = thin)




out_mcmclist <- coda::mcmc.list(out)

fit.jagsT18_2 <- window(out_mcmclist, 
                        start=n.burnin+1, 
                        end = n.iter)




print(fit.jagsT18_2, dig = 3)
list.save(fit.jagsT18_2, "fit.jagsT18_2.rds")


print(mod.jagsT18_2, dig = 3)
list.save(mod.jagsT18_2, "mod.jagsT18_2.rds")






fit_summary <- summary(fit.jagsT18_2)


MSOM_T18_2_Statistics <- DT::datatable(round(fit_summary$statistics, 3))

write.csv(MSOM_T18_2_Statistics$x$data, file = "MSOM_T18_2_Statistics.csv")





MSOM_T18_2_Quantiles <-DT::datatable(round(fit_summary$quantiles, 3))

write.csv(MSOM_T18_2_Quantiles$x$data, file = "MSOM_T18_2_Quantiles.csv")


MSOM_T18_2_OccEffectPlot <- plot_effects(mod.jagsT18_2,
                                         fit.jagsT18_2,
                                         submodel = "state")

ggsave(
  filename = "MSOM_T18_2_OccEffectPlots.pdf", 
  plot = marrangeGrob(MSOM_T18_2_OccEffectPlot, nrow=1, ncol=1), 
  width = 15, height = 9
)


MSOM_T18_2_DetEffectPlot <- plot_effects(mod.jagsT18_2,
                                         fit.jagsT18_2,
                                         submodel = "det")

ggsave(
  filename = "MSOM_T18_2_DetEffectPlots.pdf", 
  plot = marrangeGrob(MSOM_T18_2_DetEffectPlot, nrow=1, ncol=1), 
  width = 15, height = 9
)


MSOM_T18_2_OccCoefPlot_Combined <- plot_coef(mod.jagsT18_2,
                                             fit.jagsT18_2,
                                             submodel = "state",
                                             combine = T)


ggsave(
  filename = "MSOM_T18_2_OccCoefPlot_Combined.pdf", 
  plot = MSOM_T18_2_OccCoefPlot_Combined, 
  width = 15, height = 9
)



MSOM_T18_2_OccCoefPlots <- plot_coef(mod.jagsT18_2,
                                     fit.jagsT18_2,
                                     submodel = "state",
                                     ordered = TRUE)

ggsave(
  filename = "MSOM_T18_2_OccCoefPlots.pdf", 
  plot = marrangeGrob(MSOM_T18_2_OccCoefPlots, nrow=1, ncol=1), 
  width = 15, height = 9
)


MSOM_T18_2_DetCoefPlots <- plot_coef(mod.jagsT18_2,
                                     fit.jagsT18_2,
                                     submodel = "det")

ggsave(
  filename = "MSOM_T18_2_DetCoefPlots.pdf", 
  plot = marrangeGrob(MSOM_T18_2_DetCoefPlots, nrow=1, ncol=1), 
  width = 15, height = 9
)



gd <- coda::gelman.diag(fit.jagsT18_2,  multivariate = FALSE)
gd


MSOM_T18_2_GelmanD <- DT::datatable(round(gd$psrf, 3))

write.csv(MSOM_T18_2_GelmanD$x$data, file = "MSOM_T18_2_GelmanD.csv")


gp <- coda::gelman.plot(fit.jagsT18_2,  multivariate = FALSE)
gp


print(gp, dig = 3)
list.save(gp, "MSOM_T18_2_gp.rds")


mcmc.df <- as.data.frame(as.matrix(fit.jagsT18_2))

write.csv(mcmc.df, file = "FullMCMC_T18_2.csv")