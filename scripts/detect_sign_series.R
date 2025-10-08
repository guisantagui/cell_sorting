if(!require("tseries", quietly = T)){
        install.packages("tseries")
}
library(tseries)
if(!require("urca", quietly = T)){
        install.packages("urca")
}
if(!require("signal", quietly = T)){
        install.packages("signal")
}
library(signal)
library(tseries)
library(urca)
library(Kendall)
library(plotUtils)

outDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/neurodeg_aging_project/results/detect_sign_series/"
create_dir_if_not(outDir)


# Combines Ljung-Box test p-value and Mann-Kendall p-value, keeping the 
# minimum of them, as Ljung Box tends to give a non-significant result when
# changes happen at the end of the series, and Mann-Kendall cannot detect
# cyclic patterns
test_series <- function(series){
        box_pval <- Box.test(series,
                             lag = floor(sqrt(length(series))),
                             type = "Ljung-Box")$p.value
        mk_pval <- as.numeric(MannKendall(series)$sl)
        print(sprintf("%s p-value used",
                      c("Ljung-Box", "Mann-Kendall")[which.min(c(box_pval, mk_pval))]))
        pval <- c(box_pval, mk_pval)[which.min(c(box_pval, mk_pval))]
        return(pval)
}

plot_prog <- function(x, y, get_pval = F){
        df <- data.frame(x = x,
                         y = y)
        plt <- ggplot(data = df, mapping = aes(x = x, y = y)) +
                geom_line()
        if (get_pval){
                box_pval <- Box.test(y,
                                     lag = floor(sqrt(length(y))),
                                     type = "Ljung-Box")$p.value
                mk_pval <- as.numeric(MannKendall(y)$sl)
                plt <- plt +
                        labs(caption = sprintf("Ljung-Box p-value = %s; Mann-Kendall p-value = %s",
                                               box_pval,
                                               mk_pval))
                out <- list(plot = plt, ljung_box_pval = box_pval,
                            mann_kendall_pval = mk_pval)
        }else{
                out <- list(plot = plt)
        }
        return(out)
}
#funcs <- 

ages <- 20:100       

func1 <- function(ages, baseline = 0){
        return(log2(ages) + baseline)
}


func2 <- function(ages, x){
        return(ages ^ x)
}

func3 <- function(ages, baseline = 0, periodicity = 5){
        return(sin(ages/periodicity) + baseline)
}

func4 <- function(ages, baseline = 0){
        return(ages ^ 3 + ages ^ 2  + ages + baseline)
}

noise <- function(ages, sd = .01, heterosk = 0){
        if (heterosk == 0){
                sdev <- rep(sd, length(ages))
        }else if (heterosk > 0){
                sdev <- ages * heterosk * sd
        }else{
                sdev <- ages[length(ages):1] * abs(heterosk) * sd
        }
        noise <- sapply(sdev, function(x) rnorm(1, mean = 0, sd = x))
        return(noise)
}

ac_res_noise <- acf(noise(ages, sd = 10, heterosk = 2))

func5 <- function(ages, x = 2, age1 = 50, slope = 1){
        return(ifelse(ages < age1,
                      ages ^ x,
                      ((ages - age1) * slope) + age1 ^ x))
}
func5(ages)

#ac_res$acf

plot(x = ages, y = func3(ages, baseline = 10, periodicity = 3) + func1(ages) + noise(ages, sd = .1, heterosk = .1))

plot(x = ages, y = func5(ages, age1 = 40, slope = 1, x = .2) + func3(ages) * 5 + noise(ages, sd = 10, heterosk = .1))


set.seed(123)
prog1 <- func5(ages, age1 = 40, slope = 1, x = .2) + func3(ages) * 5 + noise(ages, sd = 0.1, heterosk = 0)
set.seed(234)
prog2 <- func5(ages, age1 = 40, slope = 1, x = .2) + func3(ages) * 5 + noise(ages, sd = 1, heterosk = 0)
set.seed(345)
prog3 <- func5(ages, age1 = 40, slope = 1, x = .2) + func3(ages) * 5 + noise(ages, sd = 10, heterosk = 0)
set.seed(456)
prog4 <- func5(ages, age1 = 40, slope = 1, x = .2) + func3(ages) * 5 + noise(ages, sd = 100, heterosk = 0)
set.seed(567)
prog5 <- func5(ages, age1 = 40, slope = 1, x = .2) + func3(ages) * 5 + noise(ages, sd = 20, heterosk = 0)
set.seed(678)
prog6 <- func5(ages, age1 = 40, slope = 1, x = .2) + func3(ages) * 5 + noise(ages, sd = 30, heterosk = 0)
set.seed(789)
prog7 <- func5(ages, age1 = 40, slope = 1, x = .2) + func3(ages) * 5 + noise(ages, sd = 40, heterosk = 0)
set.seed(890)
prog8 <- func5(ages, age1 = 40, slope = 1, x = .2) + func3(ages) * 5 + noise(ages, sd = 50, heterosk = 0)
set.seed(990)
prog9 <- func5(ages, age1 = 40, slope = 1, x = .2) + func3(ages) * 5 + noise(ages, sd = 60, heterosk = 0)
set.seed(891)
prog10 <- func5(ages, age1 = 40, slope = 1, x = .2) + func3(ages) * 5 + noise(ages, sd = 70, heterosk = 0)
set.seed(892)
prog11 <- func5(ages, age1 = 40, slope = 1, x = .2) + func3(ages) * 5 + noise(ages, sd = 80, heterosk = 0)
set.seed(893)
prog12 <- func5(ages, age1 = 40, slope = 1, x = .2) + func3(ages) * 5 + noise(ages, sd = 90, heterosk = 0)
set.seed(894)
prog13 <- func5(ages, age1 = 40, slope = 1, x = .2) + func3(ages) * 5 + noise(ages, sd = 100, heterosk = 0)
set.seed(896)
prog14 <- func5(ages, age1 = 40, slope = 1, x = .2) + func3(ages) * 5 + noise(ages, sd = 110, heterosk = 0)
set.seed(2)
prog15 <- func3(ages) * 5 + noise(ages, sd = .1, heterosk = .1)
set.seed(1)
prog16 <- func3(ages) * 5 + noise(ages, sd = 1., heterosk = .1)
set.seed(3)
prog17 <- func3(ages) * 5 + noise(ages, sd = 2., heterosk = .1)
set.seed(4)
prog18 <- func3(ages) * 5 + noise(ages, sd = 3., heterosk = .1)

proglist <- list(prog1 = prog1,
                 prog2 = prog2,
                 prog3 = prog3,
                 prog4 = prog4,
                 prog5 = prog5,
                 prog6 = prog6,
                 prog7 = prog7,
                 prog8 = prog8,
                 prog9 = prog9,
                 prog10 = prog10,
                 prog11 = prog11,
                 prog12 = prog12,
                 prog13 = prog13,
                 prog14 = prog14,
                 prog15 = prog15,
                 prog16 = prog16,
                 prog17 = prog17,
                 prog18 = prog18)

plot_prog(20:100, prog18, get_pval = T)
MannKendall(func3(ages))

pltList <- lapply(proglist, plot_prog, x = ages, get_pval = T)
mapply(function(p, n) ggsave(filename = sprintf("%s%s.pdf", outDir, n),
                          plot = p,
                          height = 5, width = 5),
       lapply(pltList, function(x) x$plot),
       names(pltList),
       SIMPLIFY = F)

plot_prog(ages, prog1, get_pval = T)
plot_prog(ages, prog2, get_pval = T)
plot_prog(ages, prog3, get_pval = T)
plot_prog(ages, prog4, get_pval = T)
plot_prog(ages, prog5, get_pval = T)
plot_prog(ages, prog6, get_pval = T)
plot_prog(ages, prog7, get_pval = T)
plot_prog(ages, prog8, get_pval = T)
plot_prog(ages, prog9, get_pval = T)
plot_prog(ages, prog10, get_pval = T)
plot_prog(ages, prog11, get_pval = T)
plot_prog(ages, prog12, get_pval = T)
plot_prog(ages, prog13, get_pval = T)
plot_prog(ages, prog14, get_pval = T)
plot_prog(ages, prog15, get_pval = T)

plot_prog(20:100, noise(20:100, sd = 1, heterosk = 2), get_pval = T)
res_tst <- MannKendall(noise(20:100, sd = 1, heterosk = 1))
summ <- summary(res_tst)
summary(ur.df(prog8, type = "drift", selectlags = "AIC"))
names(res_tst)

Box.test(func5(ages, age1 = 40, slope = 1, x = .2) + func3(ages) * 5 + noise(ages, sd = 10, heterosk = .06),
         lag = 10, type = "Ljung-Box")
adf.test(func5(ages, age1 = 40, slope = 1, x = .2) + func3(ages) * 5 + noise(ages, sd = 10, heterosk = .1))
adf.test(noise(ages, sd = 10, heterosk = .1))


t <- 1:100
d <- ts(10 + 0.1 * t + 5 * sin(2 * pi * t / 12) + rnorm(100), frequency = 12)

t <- 20:100
d <- ts(prog10, frequency = 1)
plot(d)
stl_result <- stl(d, s.window = NULL)
plot(stl_result)

Box.test((stl_result$time.series[, "remainder"]), lag = 10, type = "Ljung-Box")
MannKendall(stl_result$time.series[, "trend"])
plot(stl_result$time.series[, "remainder"])
