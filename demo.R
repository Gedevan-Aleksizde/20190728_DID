pacman::p_load(CausalImpact, tidyverse, ggplot2, ggthemes, patchwork)
xts2tbl <- function(data) data %>% fortify.zoo %>% as_tibble
tbl2xts <- function(data) xts(data %>% dplyr::select(-Index) %>% as.matrix)

plot_raw <- function(data){
  ggplot(data %>% xts2tbl %>% gather(key="variable", value="value", -Index),
         aes(x=Index, y=value, group=variable, color=variable, linetype=variable)) +
    geom_line() +
    facet_wrap(~variable, ncol=1, scales="free_y", strip.position='right') +
    geom_vline(xintercept=70, linetype=2) +
    theme_bw() + theme(legend.position="none") + scale_color_pander() +
    labs(x="time", y="", title='Original Series')
}

plot_CI_valid <- function(result, y0, y1){
  temp <- result$series %>% xts2tbl
  pre_periods <- result$model$pre.period[1]:result$model$pre.period[2]
  post_periods <- result$model$post.period[1]:result$model$post.period[2]
  d <- bind_rows(
    temp %>% dplyr::select(starts_with("point.pred"), Index) %>% rename_all(function(x) sub("point\\.pred", "value", x)) %>% mutate(series="difference"),
    temp %>% dplyr::select(starts_with("response"), Index) %>% rename(value=response) %>% mutate(series="difference", type="0_real"),
    temp %>% dplyr::select(starts_with("point.effect"), Index) %>% rename_all(function(x) sub("point\\.effect", "value", x)) %>%
      mutate(series="pointwise effect"),
    temp %>% dplyr::select(starts_with("cum.effect"), Index) %>% rename_all(function(x) sub("cum\\.effect", "value", x)) %>%
      mutate(series="cumulative effect"),
    temp %>% dplyr::select(starts_with("point.pred"), Index) %>% rename_all(function(x) sub("point\\.pred", "value", x)) %>% mutate(series="y0 fitting"),
    tibble(Index=post_periods, value=y0[post_periods], series="y0 fitting", type="2_counterfactual")
  ) %>% mutate(type=if_else(is.na(type), "1_pred", type)) %>% 
    mutate(series=factor(series, levels=c("y0 fitting", "difference", "pointwise effect", "cumulative effect")))
  
  ggplot(d, aes(x=Index, y=value, group=type, linetype=type)) + geom_line() +
    geom_ribbon(aes(ymin=value.lower, ymax=value.upper), alpha=.1, fill=4) +
    facet_wrap(~series, scales="free_y", ncol=1, strip.position="right") +
    geom_hline(data=filter(d, series=="pointwise effect"), aes(yintercept=0), color="grey") +
    geom_hline(data=filter(d, series=="cumulative effect"), aes(yintercept=0), color="grey") +
    geom_vline(xintercept=result$model$pre.period[2], linetype=2) +
    theme_bw() + scale_color_pander(guide=F) + scale_fill_pander(guide=F) +
    theme(axis.title.y=element_blank(), legend.position="none") + labs(x="time")
}

effect <- 10
pre.intervention <- c(from=1, to=70)
post.intervention <- c(from=71, to=100)
set.seed(42)
x1 <- 100 + arima.sim(model = list(ar = 0.999), n = 100) + .5 * 1:100
y0 <- 1.2 * x1 + rnorm(100)
data <- cbind(y1=y0, x1)

###### CASE1: Tutorial#####
data1 <- data
data1[71:100, "y1"] <- data1[71:100, "y1"] + effect + rnorm(length(71:100))

plot_raw(data1)
impact1 <- CausalImpact(data1, pre.intervention, post.intervention)
summary(impact1)
plot_CI_valid(impact1, y0, y1) + labs(title="CASE 1: Tutorial")

##### CASE2: covariate shift #####
set.seed(42)
x1 <- data[, "x1"]
x1[71:100] <- arima.sim(model=list(ar=c(.9)), n=30) + .5 * 1:30
x1[71:100] <- x1[71:100] + (x1[70] - x1[71])
y0 <- 1.2 * x1 + rnorm(100)
y1 <- y0
y1[71:100] <- y1[71:100] + effect
data2 <- cbind(y1, x1)

plot_raw(data2)
impact2 <- CausalImpact(data2, pre.intervention, post.intervention)
plot_CI_valid(impact2, y0, y1) + labs(title="CASE 2: Covariate Shift")

##### CASE3: nonlinearity misspecification #####
set.seed(42)
x1 <- data[, "x1"]
y0 <- 0.01 * x1^2 + rnorm(100)
y1 <- y0
y1[71:100] <- y1[71:100] + effect
data3 <- cbind(y1, x1)

plot_raw(data3)
impact3 <- CausalImpact(data3, pre.intervention, post.intervention)
plot_CI_valid(impact3, y0, y1) + labs(title="")

# 結果は省略するが, covariate shift を追加するとさらに差が開く

##### CASE4: multiplicative effect #####
set.seed(42)
y <- 1.2 * data[, "x1"] + rnorm(100)

data3 <- cbind(y, data[, "x1"])

plot_raw(data3)

impact3 <- CausalImpact(data1, pre.intervention, post.intervention)
plot(impact3)

######
data_ <- data1
y0 <- data_[, "y1"] + 5 * sin((1:100) * pi / 10)
y1 <- y0
y1[71:100] <- y1[71:100] + effect
data_[, "y1"] <- y1
data_[71:100, "y1"] <- NA
data_ <- zoo(data_)

plot_raw(data_)

ss <- list() %>% bsts::AddLocalLevel(data_$y1) %>% AddSeasonal(data_$y1, nseasons=20)
m <- bsts(y1 ~ x1, ss, niter=1000, data = data_)
imp <- CausalImpact(bsts.model=m, post.period.response=y1[71:100])

plot_CI_valid(imp, y0, y1)

require(bsts)
data(iclaims)
ss <- list() %>% AddLocalLevel(., initial.claims$iclaimsNSA) %>% AddSeasonal(initial.claims$iclaimsNSA, nseasons=52)
