## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup,echo=FALSE,warning=FALSE,message=FALSE-----------------------------
library(tsarma)
library(data.tree)
library(data.table)
library(flextable)
library(xts)

## ---- highlight=T,message=FALSE,warning=FALSE---------------------------------
spy <- tsdatasets::spy
spyr <- diff(log(spy$Close))[-1]
spec <- arma_modelspec(spyr[1:2000], order = c(2,1), distribution = "jsu")
mod <- estimate(spec)
summary(mod, vcov_type = "QMLE")

## -----------------------------------------------------------------------------
model_summary <- summary(mod, vcov_type = "QMLE")
out <- print(model_summary, format = "flextable", table.caption = "SPY ARMA(2,1)~JSU")
out

## ---- fig.height=8, fig.width=7-----------------------------------------------
plot(mod)

## ----highlight=T--------------------------------------------------------------
p <- predict(mod, h = 25, bootstrap = TRUE, nsim = 5000)
head(cbind(p$mean, p$sigma))

## ---- fig.width=7,fig.height = 4----------------------------------------------
par(mar = c(3,3,3,3))
plot(p, n_original = 100, main = "SPY ARMA(2,1)~JSU 25-period prediction", ylab = "", xlab = "", gradient_color = "azure3",interval_color = "orange")

## -----------------------------------------------------------------------------
parmatrix <- copy(mod$parmatrix)
spec$parmatrix <- parmatrix

## -----------------------------------------------------------------------------
filter_spec <- tsfilter(spec, y = spyr[2001:2010])
filter_model <- tsfilter(mod, y = spyr[2001:2010])
all.equal(fitted(filter_spec), fitted(filter_model))

## -----------------------------------------------------------------------------
filter_spec <- tsfilter(spec)
all.equal(fitted(mod), fitted(filter_spec))

## -----------------------------------------------------------------------------
sim_spec <- spec
sim_spec$parmatrix <- copy(mod$parmatrix)
filter_spec <- tsfilter(sim_spec)
maxpq <- max(spec$model$order)
series_init <- tail(sim_spec$target$y_orig, maxpq)
innov_init <- as.numeric(tail(residuals(filter_spec), maxpq))
simulate_spec <- simulate(sim_spec, nsim = 2000, seed = 100, h = 10, series_init = series_init, innov_init = innov_init)
simulate_model <- simulate(mod, nsim = 2000, seed = 100, h = 10, series_init = series_init, innov_init = innov_init)
set.seed(100)
p <- predict(mod, h = 10, nsim = 2000)
# equivalence between estimation and specification simulated distributions
all.equal(simulate_spec$simulated, simulate_model$simulated)
# equivalence between simulation method and prediction method distributions
all.equal(simulate_spec$simulated, p$distribution)

## ----echo=FALSE,warning=FALSE,message=FALSE-----------------------------------
x <- read.csv('function_map.csv')
colnames(x)[1] <- "(...)"
tsarma_map <- as.Node(x)
map_df <- ToDataFrameTree(tsarma_map, "type","input","output")
map_df[,1] <- gsub(" ", "\t", map_df[,1], fixed = TRUE)
names(map_df)[1] <- c("(...)")
map_df[1, 2:4] <- " "
out <- flextable(map_df) |> theme_alafoli() |> color(i = 1, j = 1, color = "purple")
out <- out |> color(i = 2, j = 1, color = "steelblue")
out <- out |> color(i = 5, j = 1, color = "steelblue")
out <- out |> color(i = 34, j = 1, color = "steelblue")
out <- out |> color(i = 32, j = 1, color = "steelblue")
out <- out |> color(i = 6, j = 1, color = "cadetblue")
out <- out |> color(i = 22, j = 1, color = "cadetblue")
out <- out |> set_caption(caption = "Table: tsarma function map")
out

