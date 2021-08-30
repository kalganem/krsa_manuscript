library(tidyverse)

data <- tibble(
  Method = rep("original",20),
  values = c(34,23,45,56,57,87,92,77,55,64,22,34,45,89,76,34,67,87,54,41)
  
)


jit_fun <- function(a, x, fct = 1) {
  l <- length(x)
  
  tibble(
    Method = rep(as.character(a), l),
    values = jitter(x = x, amount = a, factor = fct)
  )
  
}


jit_fun(c(34,23,45,56,57,87,92,77,55,64,22,34,45,89,76,34,67,87,54,41), a = 0)

map_df(c(0.5,0.9), jit_fun, x = data$values) %>% rbind(data) %>% 
  ggplot(aes(values, fill = Method)) + geom_density( alpha = 1/5)


jitter(c(34,23,45,56,57,87,92,77,55,64,22,34,45,89,76,34,67,87,54,41))

max(data$values) - min(data$values)


a_values <- seq(2, 50, by = 1)

map_df(.x = a_values, jit_fun, x = data$values) %>% 
  rbind(data) %>% 
  ggplot(aes(values, fill = Method)) + geom_density( alpha = 1/5)

map(a_values, jit_fun, x = c(34,23,45,56,57,87,92,77,55,64,22,34,45,89,76,34,67,87,54,41))


jitter(c(34,23,45,56,57,87,92,77,55,64,22,34,45,89,76,34,67,87,54,41), a = 50)


# new method
library(tidyverse)
library(KRSA)
# signal_to_noise_ratio = 15
# data = c(34,23,45,56,57,87,92,77,55,64,22,34,45,89,76,34,67,87,54,41)
# noise = rnorm(data) # generate standard normal errors
# k <- sqrt(var(data)/(signal_to_noise_ratio*var(noise)))
# data_wNoise = data + k*noise 

LFC_values <- read_delim("DLPFC_LFC.txt", delim = "\t")


chipCov <- KRSA_coverage_STK_PamChip_87102_v1
KRSA_file <- KRSA_Mapping_STK_PamChip_87102_v1

library(furrr)
plan(multisession)


add_noise <- function(x) {
  input_data <- LFC_values$LFC
  signal_to_noise_ratio = x
  set.seed(123)
  noise <- rnorm(input_data) # generate standard normal errors
  k <- sqrt(var(input_data)/(signal_to_noise_ratio*var(noise)))
  LFC_values$LFC_with_noise <- input_data + k*noise 
  LFC_values$SNR <- rep(x, length(input_data))
  LFC_values
}

krsa2 <- function (peptides, itr = 2000, return_count = F, 
          map_file = KRSA_file, cov_file = chipCov) {
  message("Running KRSA ...")
  temp <- purrr::map_df(1:itr, krsa_sampling, cov_file, map_file, 
                        length(peptides))
  temp2 <- temp %>% dplyr::group_by(Kin) %>% dplyr::summarise(SamplingAvg = mean(counts), 
                                                              SD = stats::sd(counts))
  temp3 <- cov_file %>% dplyr::group_by(Kin) %>% dplyr::summarise(Observed = sum(Substrates %in% 
                                                                                   peptides))
  fin <- dplyr::left_join(temp2, temp3) %>% dplyr::mutate(Z = (Observed - 
                                                                 SamplingAvg)/SD) %>% dplyr::arrange(dplyr::desc(abs(Z))) %>% 
    dplyr::filter(!Kin %in% c("BARK1", "VRK2")) %>% dplyr::select(Kin, 
                                                                  Observed, SamplingAvg, SD, Z) %>% dplyr::rename(Kinase = Kin)
  if (return_count == T) {
    return(list(count_mtx = temp, KRSA_Table = fin))
  }
  else {
    fin
  }
}

run_krsa_with_noise <- function(x) {
  lfc_tbl_with_noise <- add_noise(x)
  
  krsa_get_diff(lfc_tbl_with_noise,LFC_with_noise ,c(0.2,0.3,0.4)) %>% list("meanLFC" = .) -> sigPeps
  sigPeps_total <- list(sigPeps) %>% unlist(recursive = F) %>%  unlist(recursive = F)
  
  print(map_dbl(sigPeps_total, length))
  
  future_map(sigPeps_total, krsa2) -> mutiple_krsa_outputs
  
  df <- data.frame(matrix(unlist(mutiple_krsa_outputs), ncol = max(lengths(mutiple_krsa_outputs)), byrow = TRUE))
  df <- setNames(do.call(rbind.data.frame, mutiple_krsa_outputs), names(mutiple_krsa_outputs$meanLFC.0.2))
  
  df %>% rownames_to_column("method") %>% select(Kinase, Z, method) %>% 
    mutate(method = str_extract(method, "\\w+\\.\\w+\\.\\w+")) %>% 
    mutate(method = gsub("(^\\w+)[\\.]", "\\1>", method)) %>% 
    mutate_if(is.numeric, round, 2) %>% 
    filter(grepl("mean", method)) %>% 
    select(Kinase, Z, method) %>% group_by(Kinase) %>% mutate(AvgZ = mean(Z)) %>% 
    ungroup() %>% 
    select(Kinase, AvgZ) %>% distinct() %>% 
    mutate(SNR = x) -> AvgZTable
  
  AvgZTable
  
}


map_df(seq(1, to = 3, by = 1), run_krsa_with_noise) -> final_res


org_krsa_table <- read_delim("acrossChip_KRSA_Table_multipleComp.txt", delim = "\t") %>% 
  mutate_if(is.numeric, round, 2) %>% 
  filter(grepl("mean", method)) %>% 
  select(Kinase, Z, method) %>% group_by(Kinase) %>% mutate(AvgZ = mean(Z)) %>% 
  ungroup() %>% 
  select(Kinase, AvgZ) %>% distinct() %>% 
  mutate(SNR = "Raw_Data") -> org_krsa_table

rbind(final_res, org_krsa_table) -> final_res2


final_res2 %>% pivot_wider(names_from = SNR, values_from = AvgZ) %>% 
  column_to_rownames("Kinase") -> final_res3


M <-cor(final_res3)
library(corrplot)
library(RColorBrewer)
corrplot(M, type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu")) 

M %>% as.data.frame() %>% 
  rownames_to_column("SNR") %>% 
  select(SNR, everything()) %>% 
  pivot_longer(2:ncol(.), names_to = "Method", values_to = "Cor") %>% 
  filter(Method == "Raw_Data") %>% 
  ggplot(aes(SNR, Cor)) + geom_point()

final_res %>% pivot_wider(names_from = SNR, values_from = AvgZ) %>%
  ggplot(aes(`0.05`, `0.5`)) + geom_point()

