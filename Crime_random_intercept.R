library(rstan)
library(dplyr)
library(Ecdat)
library(ggplot2)

#load Crime Data from Ecdat
data(Crime)
Crime

#data preprocessing 
#transform all variables to log rates, delete counties 185 and 115 to get rid of faulty log odds of arrest data > 1 
#delete prob of conviction all around since numerous values > 1 

Crime_transformed = Crime %>%
  mutate(log_crime_rate = log(crmrte))%>%
  mutate(log_odd_arrest = log((prbarr/(1-prbarr))))%>%
  mutate(log_odd_prison = log(prbpris/(1-prbpris)))%>%
  mutate(log_police_ratio = log(polpc/mean(polpc)))%>%
  mutate(log_pop_density_ratio = log(density/mean(density)))%>% 
  filter(county != 185 & county != 115) %>% 
  select(county, year, log_crime_rate, log_odd_arrest,log_odd_prison,log_police_ratio, log_pop_density_ratio, region)
Crime_transformed

#Make first 6 years 81-86 training data 
#grouping first 6 years by county averages 
Crime.grouped.train = Crime_transformed %>% 
  filter(year!= 87) %>% 
  group_by(county, region) %>%  
  summarize(avg_crime = mean(log_crime_rate), avg_arrest = mean(log_odd_arrest),
            avg_prison = mean(log_odd_prison),avg_police = mean(log_police_ratio), avg_pop_dens = mean(log_pop_density_ratio))
Crime.grouped.train

#year 87 grouped by county test data 
Crime.grouped.test = Crime_transformed %>% 
  filter(year== 87) %>% 
  group_by(county, region) %>%  
  summarize(avg_crime = mean(log_crime_rate), avg_arrest = mean(log_odd_arrest),
            avg_prison = mean(log_odd_prison),avg_police = mean(log_police_ratio), avg_pop_dens = mean(log_pop_density_ratio))  
Crime.grouped.test


#make first 6 yearss x1, x2, x3, x4 variables
x1 = Crime.grouped.train$avg_arrest
x2 = Crime.grouped.train$avg_prison
x3 = Crime.grouped.train$avg_police
x4 = Crime.grouped.train$avg_pop_dens

y = Crime.grouped.train$avg_crime
N = length(y)

y_test = Crime.grouped.test$avg_crime
x1_test = Crime.grouped.test$avg_arrest
x2_test = Crime.grouped.test$avg_prison
x3_test = Crime.grouped.test$avg_police
x4_test = Crime.grouped.test$avg_pop_dens

#all priors/hyper priors
sd_alpha = 1
mu_beta1 = -.4
sd_beta1 = .3
mu_beta2 = .5
sd_beta2 = .3
mu_beta3 = .3
sd_beta3 = .2
mu_beta4 = .3
sd_beta4 = .2

nu = 7
A = sd(y)
county = seq(1, 88, by = 1)
J = length(county)
n_test <- length(y_test)
x_grid = seq(-2.5, 1, by = 0.04)
n_grid = length(x_grid)
length(x_grid)

mlr.random = list(N = N, J = J, y = y, group_id = county, group_id_test = county, x1 = x1, x2 = x2, x3 = x3, x4 = x4,
                  y_test = y_test, x1_test = x1_test, x2_test = x2_test,
                  x3_test = x3_test, x4_test = x4_test, mu_alpha = mu_alpha, sd_alpha = sd_alpha,
                  mu_beta1 = mu_beta1, sd_beta1 = sd_beta1,
                  mu_beta2 = mu_beta2, sd_beta2 = sd_beta2,
                  mu_beta3 = mu_beta3, sd_beta3 = sd_beta3,
                  mu_beta4 = mu_beta4, sd_beta4 = sd_beta4,
                  nu = nu, A = A, n_test = n_test, n_grid = n_grid, x_grid = x_grid)

multiple.model <- stan_model(file = "~/Crime_random_intercept.stan")
mlr.fit <- sampling(object = multiple.model,
                    data = mlr.random,
                    pars = c("alpha", "beta1", "beta2", "beta3", "beta4", "ystar", "post_line", "prob_grid"))

mlr.fit

ystar_samples <- extract(object = mlr.fit, pars = "ystar")[["ystar"]]
ystar_mean <- apply(ystar_samples, MARGIN = 2, FUN = mean)
ystar_l95 <- apply(ystar_samples, MARGIN = 2, FUN = quantile, probs = 0.025)
ystar_u95 <- apply(ystar_samples, MARGIN = 2, FUN = quantile, probs = 0.975)

#Model Diagnostics 
mean((y_test - ystar_mean)^2 )/mean((y_test - mean(y))^2)
test_in_interval <- ((y_test >= ystar_l95) & (y_test <= ystar_u95))  
mean(test_in_interval) 

#PLOTTING Posterior Intervals 
sorted_test <- sort(y_test, decreasing = FALSE, index.return = TRUE)
sorted_index <- sorted_test$ix 

png("MLR_OOS_hierarchical.png", width = 9, height = 4, units = "in", res= 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
plot(1, type = "n", xlim = c(0, 90), ylim = range(c(ystar_samples, y_test)),
     main = "Hierarchical Model Posterior Predictive Uncertainty", ylab = "Log Crime Rate", xlab = "County", xaxt = "n")
for(i in 1:n_test){
  lines(x = c(i,i), y = c(ystar_l95[sorted_index[i]], ystar_u95[sorted_index[i]]), col = 'gray')
  points(x = i, y = ystar_mean[sorted_index[i]], col = 'black', pch = 16)
  points(x = i, y = y_test[sorted_index[i]], pch = 17,
         col = ifelse(test_in_interval[sorted_index[i]],'green', 'red'))
}
legend("bottomright", legend = c("Prediction", "Actual", "Uncertainty interval"),
       pch = c(16, 17, NA), lty = c(NA, NA, 1))
dev.off()

#coefficient distributions together
library(bayesplot)
color_scheme_set("blue")

plot_title = ggtitle("Posterior Distributions with Random Intercepts")
mcmc_areas(mlr.fit, pars = c("beta1", "beta2", "beta3", "beta4")) +
  scale_y_discrete(expand = c(0, 0)) + plot_title
mcmc_dens(mlr.fit, pars = c("beta3", "beta4", "beta1", "beta2")) +
  scale_y_discrete(expand = c(0, 0)) + plot_title
mcmc_areas(mlr.fit, pars = c("beta1", "beta2", "beta3", "beta4"), prob = 0.89) +
  scale_y_discrete(expand = c(0, 0)) + plot_title

#boxplot
std_alpha_samples <- extract(mlr.fit, pars = "alpha")[["alpha"]]
colnames(std_alpha_samples) <- county

png("county_baseline_probs.png", width = 9, height = 4, units = "in", res = 400)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
boxplot(std_alpha_samples, pch = 16, cex = 0.5, medlwd = 0.5,
        main = "Hierarchical County-Specific Log Crime Rate Intercepts", xlab = "County #", ylab = "Log Crime Rate")
dev.off()

#Individual coefficient distributions (NOT USED)
alpha_samples <- extract(object = mlr.fit, pars = "alpha")[["alpha"]]
beta1_samples <- extract(object = mlr.fit, pars = "beta1")[["beta1"]]
beta2_samples <- extract(object = mlr.fit, pars = "beta2")[["beta2"]]
beta3_samples <- extract(object = mlr.fit, pars = "beta3")[["beta3"]]
beta4_samples <- extract(object = mlr.fit, pars = "beta4")[["beta4"]]

png("beta1_params_post.png", width = 9, height = 4, units = "in", res = 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(1,2))
hist(alpha_samples, breaks = 100, main = "Posterior draws of alpha", xlab = "alpha")
hist(beta1_samples, breaks = 100, main = "Posterior draws of beta1", xlab = "beta1")
dev.off()
png("beta2_params_post.png", width = 9, height = 4, units = "in", res = 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(1,2))
hist(alpha_samples, breaks = 100, main = "Posterior draws of alpha", xlab = "alpha")
hist(beta2_samples, breaks = 100, main = "Posterior draws of beta2", xlab = "beta2")
dev.off()
png("beta3_params_post.png", width = 9, height = 4, units = "in", res = 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(1,2))
hist(alpha_samples, breaks = 100, main = "Posterior draws of alpha", xlab = "alpha")
hist(beta3_samples, breaks = 100, main = "Posterior draws of beta3", xlab = "beta3")
dev.off()
png("beta4_params_post.png", width = 9, height = 4, units = "in", res = 300)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(1,2))
hist(alpha_samples, breaks = 100, main = "Posterior draws of alpha", xlab = "alpha")
hist(beta4_samples, breaks = 100, main = "Posterior draws of beta4", xlab = "beta4")
dev.off()

#Not used
post_pred_samples <- extract(mlr.fit, pars = "prob_grid")[["prob_grid"]]
dimnames(post_pred_samples)[[2]] <- c(county, "new_county") 
new_county_grid_samples <- post_pred_samples[,"new_county",]

new_county_mean <- apply(new_county_grid_samples, MAR = 2, FUN = mean)
new_county_l95 <- apply(new_county_grid_samples, MAR = 2, FUN = quantile, probs = 0.025)
new_county_u95 <- apply(new_county_grid_samples, MAR = 2, FUN = quantile, probs = 0.975)

plot(1, type= "n", xlim = c(-3, 1), ylim = c(-5,-2),
     main = "Posterior Predictive Log Crime Rate", xlab = "County", ylab = "Log Crime Rate")
polygon(x = c(x_grid, rev(x_grid)), 
        y = c(new_county_l95, rev(new_county_u95)),
        col = rgb(0, 0, 1, 1/4),
        border = NA)
lines(x = x_grid, y = new_county_mean, col = 'blue', lwd = 2)
legend("bottomleft", legend = c("New county"),
       col = c("blue"), lty = 1)
dev.off()


