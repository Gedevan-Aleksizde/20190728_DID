pacman::p_load(CausalImpact, tidyverse, tsibble, ggthemes, lubridate)
pacman::p_load_gh("thomasp85/patchwork")

##### draw figure 1; DID framework ####
tribble(
  ~time, ~type, ~KPI, ~varname,
  0, "処置群 (平行トレンド仮定)", 3, "hat(y)[b](1)",
  0, "処置群", 3, "y[b](1)",
  0, "対照群", 1, "y[b](0)",
  1, "処置群", 5, "y[a](1)",
  1, "処置群 (平行トレンド仮定)", 4, "hat(y)[a](1)",
  1, "対照群", 2, "y[a](0)"
) %>% ggplot(aes(x=time, y=KPI, group=type, linetype=type, color=type, label=varname)) +
  geom_line(size=1) + geom_point() + geom_label(parse=T, nudge_y = .2, size=5) +
  labs(x="t") + theme_classic() +
  theme(legend.position="bottom", legend.title=element_blank(), axis.text.y=element_blank(),
        axis.title.y=element_text(size=20, angle=0), axis.title.x=element_text(size=20),
        axis.text.x=element_text(size=15),
        legend.text=element_text(size=15),
        panel.grid.major.y=element_line(color="grey")
        ) +
  scale_x_continuous(breaks = c(0, 1), labels=c("before", "after")) +
  coord_cartesian(xlim=c(-.5, 1.5)) + scale_y_continuous(breaks=0:5)

##### draw figure 2 ####
set.seed(42)
df_fig2 <- tibble(time=1:50) %>% mutate(
  y0=rnorm(n=n(), sd=.05) + .01 * time,
  y1=if_else(time <35, y0, y0 + rnorm(mean=.2, sd=.05, n=n()))
  )
p1 <- ggplot(gather(df_fig2, key=series, value=y, -time), aes(x=time, y=y, group=series, color=series)) + geom_line() +
  theme_clean() + scale_color_pander() + theme(legend.position="bottom", legend.title=element_blank())
p2 <- ggplot(df_fig2, aes(x=time, ymin=y0, ymax=y1, fill="diff")) + geom_ribbon() +
  theme_clean() + theme(legend.position="bottom", legend.title=element_blank())  + scale_fill_pander()
p3 <- ggplot(mutate(df_fig2, pointwise=y1-y0, cumulative=cumsum(pointwise)) %>% dplyr::select(-y0, -y1) %>% gather(key=series, value=effect, -time) %>% mutate(series=factor(series, levels=c("pointwise", "cumulative"))), aes(x=time, y=effect, group=series, color=series)) + geom_line() + facet_wrap(~series, ncol=1, scales="free_y") +
  theme_clean() + theme(legend.position="bottom", legend.title=element_blank())  + scale_color_pander()
p1 / p2 / p3

##### define some utility functions #####
tsibble2zoo <- function(x){
  stopifnot(inherits(x, "tbl_ts"))
  zoo::read.zoo(dplyr::select(x, index_var(x), everything()))
}

zoo2tsibble <- function(x){
  stopifnot(inherits(x, "zoo"))
  fortify.zoo(x) %>% as_tsibble(index=Index)
}

print_bsts_ss_components <- function(x){
  if(inherits(x, "CausalImpact"))
    ss <- x$model$bsts.model$state.specification 
  else if(inherits(x, "bsts"))
    ss <- x$state.specification
  map(ss, function(x) attr(x, "class"))
}

plot_raw <- function(data, intervention=NULL){
  if(!inherits(data, "tbl_ts")){
    data <- as_tsibble(data, index=Index)
  }
  data <- data %>% gather(key="series", value="value", -Index) %>% mutate(xy=substr(series, 1, 1))
  plot_x <- filter(data, xy=="x") %>% ggplot(aes(x=Index, y=value, group=series, color=series, linetype=series)) +
    geom_line(size=1) + theme_clean() + scale_color_pander() +
    theme(title=element_text(size=15), text=element_text(size=15), axis.title=element_text(size=15),
          legend.position="bottom", legend.title=element_blank(), axis.title.x=element_blank()) + labs(y="covariates")
  plot_y <- filter(data, xy=="y") %>% ggplot(aes(x=Index, y=value, group=series, color=series, linetype=series)) +
    geom_line(size=1) + theme_clean() + scale_color_pander() +
    theme(title=element_text(size=15), text=element_text(size=15), axis.title=element_text(size=15),
          legend.position="bottom", legend.title=element_blank()) + labs(x="time", y="outcomes")
  if(!is.null(intervention)){
    plot_x <- plot_x + geom_vline(xintercept=intervention, linetype=2)
    plot_y <- plot_y + geom_vline(xintercept=intervention, linetype=2)
  }
  if("POSIXct" %in% class(data$Index)){
    plot_x <- plot_x + scale_x_datetime(date_labels="%Y-%m-%d")
    plot_y <- plot_y + scale_x_datetime(date_labels="%Y-%m-%d")
  }
  else if("Date" %in% class(data$Index)){
    plot_x <- plot_x + scale_x_date(date_labels="%Y-%m-%d")
    plot_y <- plot_y + scale_x_date(date_labels="%Y-%m-%d")
  }
  plot_x / plot_y
}

plot_CI_valid <- function(result, data){
  # result は CI の結果, data y0 を含む CI への入力データ
  # y0 <- as.numeric(y0)
  # temp <- result$series %>% zoo2df
  d <- left_join(
    zoo2tsibble(result$series) %>% rename_at(vars(-(1:2)), .funs=function(x) paste0(x, "_estimated")) %>% rename(y1_true=response),
    dplyr::select(data, y0) %>% rename(y0_true=y0)
    ) %>%
    mutate(point.effect_true=y1_true - y0_true,
           cum.effect_true=cumsum(point.effect_true)) %>% gather() %>% as_tibble %>% mutate(vlu=factor(case_when(
             str_detect(key, "\\.lower_") ~ "lower",
             str_detect(key, "\\.upper_") ~ "upper",
             T ~ "value"
           )),
           type=factor(str_replace(key, ".+_(.+)$", "\\1"), levels=c("true", "estimated")),
           key=str_replace(key, "_.+$", "") %>% str_replace("(\\.upper|\\.lower)$", "")
           ) %>% filter(!key %in% c("cum.response", "cum.pred")) %>%
    mutate(
      key=case_when(
        key=="point.effect" ~ "pointwise effect",
        key=="cum.effect" ~ "cumulative effect",
        key=="point.pred" ~ "y0",
        key=="response" ~ "y1",
        T ~ key
      )
    ) %>% mutate(
      y=if_else(key %in% c("y0", "y1"), key, ""),
      key=if_else(key %in% c("y0", "y1"), "y", key) %>% factor(levels=c("y", "pointwise effect", "cumulative effect"))
    ) %>%
    filter(!(type=="estimated" & ! key=="y" & Index %within% pre_span)) %>%
    spread(key=vlu, value=value)
  
  g <- ggplot(d, aes(x=Index, y=value, group=paste(y, type), linetype=type)) + geom_line(size=1, aes(color=type)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=.1, fill=4) +
    facet_wrap(~key, scales="free_y", ncol=1, strip.position="top") +
    geom_hline(data=filter(d, key=="pointwise effect"), aes(yintercept=0), color="grey") +
    geom_hline(data=filter(d, key=="cumulative effect"), aes(yintercept=0), color="grey") +
    geom_vline(xintercept=result$model$pre.period[2], linetype=2) +
    scale_color_pander() + scale_fill_pander(guide=F) +
    scale_linetype_discrete() + theme_clean() +
    theme(title=element_text(size=15), text=element_text(size=15),
          axis.title.x=element_blank(), axis.title.y=element_blank(),
          legend.position="bottom", legend.title=element_blank()) + labs(x="time")
  if("POSIXct" %in% class(d$Index)){
    g <- g + scale_x_datetime(date_labels="%Y-%m-%d")
  }
  else if("Date" %in% class(d$Index)){
    g <- g + scale_x_date(date_labels="%Y-%m-%d")
  }
  g
}


##### make dataset ####
effect <- 10
set.seed(42)
n <- 100
x1 <- arima.sim(model = list(ar = 0.999), n = n) + .5 * 1:n
y0 <- 1.2 * x1 + rnorm(n)
data <- tibble(y0, x1) %>% mutate(Index=seq(ymd("2019-01-01"), by="1 day", length.out=n())) %>% as_tsibble(index=Index)

pre_intervention <- c(start=data$Index[1], end=data$Index[70])
post_intervention <- c(start=data$Index[71], end=data$Index[n])
pre_span <- do.call(interval, as.list(pre_intervention))
post_span <- do.call(interval, as.list(post_intervention))


###### CASE1: Tutorial#####
set.seed(42)
data1 <- mutate(data,
  y=if_else(Index %within% pre_span, y0, y0 + effect + rnorm(n=n()))
  )

plot_raw(data1, pre_intervention[2]) + plot_annotation(title="CASE 1: Tutorial")
impact1 <- CausalImpact(dplyr::select(data1, y, everything(), -y0) %>% tsibble2zoo,
                        pre_intervention, post_intervention)
summary(impact1)
plot(impact1)
plot_CI_valid(impact1, data1) + labs(title="CASE 1: Tutorial")


##### CASE2: covariate shift #####
set.seed(42)
data2 <- data %>% mutate(
  x1=c(filter(., Index %within% pre_span)$x1,
       arima.sim(model=list(ar=c(.9)), n=NROW(filter(., Index %within% post_span)))+ last(filter(., Index  %within% pre_span)$x1) +
         .5 * 1:NROW(filter(., Index %within% post_span))
       ),
  y0= 1.2 * x1 + rnorm(n()),
  y= y0 + if_else(Index %within% pre_span, 0, effect + rnorm(n=n()))
  )

plot_raw(data2, pre_intervention[2]) + plot_annotation(title="CASE2 Raw Data")
impact2 <- CausalImpact(dplyr::select(data2, y, everything(), -y0) %>% tsibble2zoo,
                        pre_intervention, post_intervention)
summary(impact2)
plot(impact2)
plot_CI_valid(impact2, data2) + labs(title="CASE 2: Covariate Shift")


##### CASE3-1: nonlinearity misspecification #####
set.seed(42)
data3 <- data %>% mutate(
  y0=0.01 * x1^2 + rnorm(n()),
  y=if_else(Index %within% pre_span, y0, y0 + effect + rnorm(n=n()))
  )

plot_raw(data3, pre_intervention[2]) + plot_annotation(title="CASE 3-1 Raw Data")
impact3 <- CausalImpact(dplyr::select(data3, y, everything(), -y0) %>% tsibble2zoo,
                        pre_intervention, post_intervention)
summary(impact3)
plot(impact3)
plot_CI_valid(impact3, data3) + labs(title="CASE 3-1 Nonlinearity Misspecification")

# covariate shift を追加するとさらに差が開く
set.seed(42)
data3_2 <- data %>% mutate(
  x1=c(filter(., Index %within% pre_span)$x1,
       arima.sim(model=list(ar=c(.9)), n=NROW(filter(., Index %within% post_span)))+ last(filter(., Index  %within% pre_span)$x1) +
         .5 * 1:NROW(filter(., Index %within% post_span))
  ),
  y0= .01 * x1^2 + rnorm(n()),
  y=if_else(Index %within% pre_span, y0, y0 + effect + rnorm(n=n()))
)

plot_raw(data3_2, pre_intervention[2]) + plot_annotation(title="CASE 3-2 Raw Data")
impact3_2 <- CausalImpact(dplyr::select(data3_2, y, everything(), -y0) %>% tsibble2zoo,
                        pre_intervention, post_intervention)
summary(impact3_2)
plot_CI_valid(impact3_2, data3_2) + labs(title="CASE 3-2 Nonlinearity Misspecification + Covariate Shift")

set.seed(42)
data3_3 <- data %>% mutate(
  seasonal = 5 * sin(row_number()*pi/5),
  y0=1.2 * x1 + as.numeric(arima.sim(model=list(ar=c(.5, .3)), n=n())) + seasonal,
  y=if_else(Index %within% pre_span, y0, y0+effect + rnorm(n=n()))
)
plot_raw(data3_3, pre_intervention[2]) + plot_annotation(title="CASE 3-3 Raw Data")
impact3_3 <- CausalImpact(dplyr::select(data3_3, y, everything(), -y0, -seasonal) %>% tsibble2zoo,
                          pre_intervention, post_intervention)
summary(impact3_3)
plot_CI_valid(impact3_3, data3_3) + labs(title="CASE 3-3 Autoregressive + Seasonal")



impact3_3$series %>% zoo2tsibble() %>% mutate(residual=response - point.pred) %>%
  ggplot(aes(x=Index, y=residual)) + geom_line() + labs(title="Residuals") + theme_clean() +
  theme(legend.title=element_blank(), axis.text.y=element_blank(),
        axis.title.y=element_text(size=20), axis.title.x=element_text(size=20),
        axis.text.x=element_text(size=15),
        legend.text=element_text(size=15),
        panel.grid.major.y=element_line(color="grey")
  )
impact3_3$series %>% zoo2tsibble() %>% mutate(resid=response - point.pred) %>% dplyr::select(resid) %>% acf

# model.args で周期成分を追加
impact3_3 <- CausalImpact(dplyr::select(data3_3, y, everything(), -y0, -seasonal) %>% tsibble2zoo,
                          pre_intervention, post_intervention, model.args=list(nseasons=10))
summary(impact3_3)
plot_CI_valid(impact3_3, data3_3) + labs(title="CASE 3-3 Autoregressive + Seasonal (Fitted with Seasonal Local Level Model")

# bsts を与える場合少しめんどくさい
data3_3_input <- data3_3 %>% mutate(y=if_else(Index %within% pre_span, y, NA_real_))
ss3_3 <- list() %>% AddAr(y=data3_3_input$y, lags=2) %>%
  AddSeasonal(y=data3_3_input$y, nseasons=10)
model3_3 <- bsts(formula=y~x1, state.specification=ss3_3, timestamps=data3_3$Index,
                 data=data3_3_input %>% as_tibble, niter=1000, seed=42)
plot(model3_3)
summary(model3_3)
plot(model3_3, "comp")

impact3_3_mod <- CausalImpact(bsts.model=model3_3, post.period.response=filter(data3_3, Index %within% post_span)$y)
summary(impact3_3_mod)
plot_CI_valid(impact3_3_mod, data3_3) + labs(title="CASE 3-3 Autoregressive + Seasonal (Fitted with AR(2) + Seasonal)")


##### CASE4: structural change #####
set.seed(42)
data4 <- data %>% mutate(
  y0=if_else(Index %within% pre_span, y0, (filter(data, Index==pre_intervention[2])$y0 - filter(data, Index==pre_intervention[2])$x1) + 0.9 * x1 + rnorm(n=n())),
  y=if_else(Index %within% pre_span, y0, y0 + effect + rnorm(n=n()))
)

plot_raw(data4, pre_intervention[2]) + plot_annotation(title="CASE 4 Raw Data")
impact4 <- CausalImpact(dplyr::select(data4, y, everything(), -y0) %>% tsibble2zoo,
                        pre_intervention, post_intervention)
summary(impact4)
plot(impact4)
plot_CI_valid(impact4, data4) + labs(title="CASE 4: Structural Change")

##### CASE 5: nonstationarity #####
set.seed(42)
data5 <- data %>% mutate(
  y0=1.2 + x1 + cumsum(rnorm(n=n(), sd=2)),
  y=if_else(Index %within% pre_span, y0, y0 + effect + rnorm(n=n())))
plot_raw(data5, pre_intervention[2]) + plot_annotation(title="CASE 5 Raw Data")
impact5 <- CausalImpact(dplyr::select(data5, y, everything(), -y0) %>% tsibble2zoo,
                        pre_intervention, post_intervention)
summary(impact5)
plot(impact5)
plot_CI_valid(impact5, data5) + labs(title="CASE 5: Random Walk Noise")