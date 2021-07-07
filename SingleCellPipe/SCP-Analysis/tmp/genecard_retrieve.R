setwd("/home/zhanghao/scData/")
## run the command in the console
## java -jar ~/Program/SystemTools/selenium-server-standalone-3.141.59.jar
library(rvest)
library(httr)
library(magrittr)
library(RSelenium)
library(stringr)
remDr <- remoteDriver(
  remoteServerAddr = "127.0.0.1",
  port = 4444,
  browserName = "firefox"
)
remDr$open()

url_main <- "https://discovery.lifemapsc.com/in-vivo-development/cellular"
remDr$navigate(url_main)
pageid <- remDr$getPageSource()[[1]][1] %>% read_html() %>% html_nodes("div.page-sizer label") %>% html_attr("for")
element <- remDr$findElement(using = "id",value = pageid)
element$sendKeysToElement(sendKeys = list("All",key = "enter"))
element$executeScript(paste0("document.getElementById('",pageid,"').click()"))
df_cell <- element$getPageSource()[[1]][1] %>% read_html() %>% html_nodes("div.data-container")  %>% .[[1]] %>% html_table() %>% as.data.frame()
colnames(df_cell) <- df_cell[1,]
df_cell <- df_cell[-1,]
lines <- element$getPageSource()[[1]][1] %>% read_html() %>% html_nodes("div.data-container") %>% .[[1]] %>% html_nodes("tbody tr")
list_href <- lapply(lines,function(x) x %>% html_children() %>% html_node("a") %>% html_attr("href"))
df_cell_href <- as.data.frame(do.call(rbind, list_href))
saveRDS(df_cell,file = "df_cell.rds")

cell_gene_list <- list()
compartment_gene_list <- list()
organ_gene_list <- list()
for (i in 367:nrow(df_cell_href)) {
  cat(i,"/",nrow(df_cell_href),"\r")
  
  cell <- df_cell[i,1]
  if (!cell %in% names(cell_gene_list)) {
    url_cell <- df_cell_href[i, 1]
    remDr$navigate(url_cell)
    panel <- remDr$getPageSource()[[1]][1] %>% read_html() %>% html_nodes("div.content div.grid-panel")
    id <- panel %>% html_attr("id")
    id_use <- which(id == "GeneExpression")
    pageid <- panel[[id_use]] %>% html_nodes("div.page-sizer label") %>% html_attr("for")
    element <- remDr$findElement(using = "id",value = pageid)
    element$sendKeysToElement(sendKeys = list("All",key = "enter"))
    element$executeScript(paste0("document.getElementById('",pageid,"').click()"))
    panel <- element$getPageSource()[[1]][1] %>% read_html() %>% html_nodes("div.content div.grid-panel")
    container <-  panel %>% .[[id_use]] %>% html_node("div.data-container")
    if (html_text(container) != "No Data Found") {
      df_gene <- container %>% html_table() %>% as.data.frame()
      lines <- container %>% html_nodes("tbody tr")
      nrow <- length(lines)
      ncol <- lines[[1]] %>% html_children() %>% length()
      info <- as.data.frame(matrix(nrow = nrow,ncol = ncol))
      for (j in 1:ncol) {
        info[,j] <- sapply(lines,function(x) html_children(x)[[j]] %>% html_nodes("span") %>% html_attr("title") %>% paste0(collapse = ","))
      }
      list_href <- lapply(lines,function(x) x %>% html_children() %>% html_node("a") %>% html_attr("href"))
      info_id <- as.data.frame(do.call(rbind, list_href))
      info[,3] <- str_extract(info_id[,3],pattern = "(?<=symbols\\=).*")
      info[,6] <- str_extract(info_id[,6],pattern = "(?<=gene\\=).*") 
      info[,7] <- str_extract(info_id[,7],pattern = "(?<=gene/).*") 
      index <- which(apply(info,2,function(x)any(x != "",na.rm = TRUE)))
      df_gene[,index] <- info[,index]
      df_gene[,"cell"] <- cell
      colnames(df_gene)[c(1,6,7)] <- c("Species","GeneCards-symbol","ENTREZID")
      cell_gene_list[[cell]] <- df_gene
    }
  }
  
  compartment <- df_cell[i,2]
  if (!compartment %in% names(compartment_gene_list)) {
    url_compartment <- df_cell_href[i, 2]
    remDr$navigate(url_compartment)
    panel <- remDr$getPageSource()[[1]][1] %>% read_html() %>% html_nodes("div.content div.grid-panel")
    id <- panel %>% html_attr("id")
    id_use <- which(id == "GeneExpression")
    pageid <- panel[[id_use]] %>% html_nodes("div.page-sizer label") %>% html_attr("for")
    element <- remDr$findElement(using = "id",value = pageid)
    element$sendKeysToElement(sendKeys = list("All",key = "enter"))
    element$executeScript(paste0("document.getElementById('",pageid,"').click()"))
    panel <- element$getPageSource()[[1]][1] %>% read_html() %>% html_nodes("div.content div.grid-panel")
    container <-  panel %>% .[[id_use]] %>% html_node("div.data-container")
    if (html_text(container) != "No Data Found") {
      df_gene <- container %>% html_table() %>% as.data.frame()
      lines <- container %>% html_nodes("tbody tr")
      colnames(df_gene) <- df_gene[1,]
      df_gene <- df_gene[-1,]
      nrow <- length(lines)
      ncol <- lines[[1]] %>% html_children() %>% length()
      info <- as.data.frame(matrix(nrow = nrow,ncol = ncol))
      for (j in 1:ncol) {
        info[,j] <- sapply(lines,function(x) html_children(x)[[j]] %>% html_nodes("span") %>% html_attr("title") %>% paste0(collapse = ","))
      }
      list_href <- lapply(lines,function(x) x %>% html_children() %>% html_node("a") %>% html_attr("href"))
      info_id <- as.data.frame(do.call(rbind, list_href))
      info[,3] <- str_extract(info_id[,3],pattern = "(?<=symbols\\=).*")
      info[,5] <- str_extract(info_id[,5],pattern = "(?<=gene\\=).*")
      info[,6] <- str_extract(info_id[,6],pattern = "(?<=gene/).*")
      index <- which(apply(info,2,function(x)any(x != "")))
      df_gene[,index] <- info[,index]
      df_gene[,"compartment"] <- compartment
      colnames(df_gene)[c(1,5,6)] <- c("Species","GeneCards-symbol","ENTREZID")
      compartment_gene_list[[compartment]] <- df_gene
    }
  }
  
  organ <- df_cell[i,3]
  if (!organ %in% names(organ_gene_list)) {
    url_organ <- df_cell_href[i, 3]
    remDr$navigate(url_organ)
    panel <- remDr$getPageSource()[[1]][1] %>% read_html() %>% html_nodes("div.content div.grid-panel")
    id <- panel %>% html_attr("id")
    id_use <- which(id == "GeneExpression")
    pageid <- panel[[id_use]] %>% html_nodes("div.page-sizer label") %>% html_attr("for")
    element <- remDr$findElement(using = "id",value = pageid)
    element$sendKeysToElement(sendKeys = list("All",key = "enter"))
    element$executeScript(paste0("document.getElementById('",pageid,"').click()"))
    panel <- element$getPageSource()[[1]][1] %>% read_html() %>% html_nodes("div.content div.grid-panel")
    container <-  panel %>% .[[id_use]] %>% html_node("div.data-container")
    if (html_text(container) != "No Data Found") {
      df_gene <- container %>% html_table() %>% as.data.frame()
      lines <- container %>% html_nodes("tbody tr")
      nrow <- length(lines)
      ncol <- lines[[1]] %>% html_children() %>% length()
      info <- as.data.frame(matrix(nrow = nrow,ncol = ncol))
      for (j in 1:ncol) {
        info[,j] <- sapply(lines,function(x) html_children(x)[[j]] %>% html_nodes("span") %>% html_attr("title") %>% paste0(collapse = ","))
      }
      list_href <- lapply(lines,function(x) x %>% html_children() %>% html_node("a") %>% html_attr("href"))
      info_id <- as.data.frame(do.call(rbind, list_href))
      info[,3] <- str_extract(info_id[,3],pattern = "(?<=symbols\\=).*")
      info[,6] <- str_extract(info_id[,6],pattern = "(?<=gene\\=).*")
      info[,7] <- str_extract(info_id[,7],pattern = "(?<=gene/).*")
      index <- which(apply(info,2,function(x)any(x != "")))
      df_gene[,index] <- info[,index]
      df_gene[,"organ"] <- organ
      colnames(df_gene)[c(1,6,7)] <- c("Species","GeneCards-symbol","ENTREZID")
      organ_gene_list[[organ]] <- df_gene
    }
  }
}
saveRDS(cell_gene_list,"cell_gene_list.rds")
saveRDS(compartment_gene_list,"compartment_gene_list.rds")
saveRDS(organ_gene_list,"organ_gene_list.rds")

















##### PTM info for every protein
##### 
library(rvest)
library(iptmnetr)

library(httr)
library(data.table)
## Roles: enz/AND-OR/sub
role="enzORsub" 
## All PTMs: ac,gn,go,gc,gs,me,my,p,su,ub,sno
mod="ac,gn,go,gc,gs,me,my,p,su,ub,sno" 
## Species: human-9606,mouse-10090
species=10090
int_url=paste0("https://research.bioinformatics.udel.edu/iptmnet/browse/role-",role,"&mod-",mod,"&taxon-",species)
web<-read_html(x = int_url,encoding="UTF-8")

protein_num <- web %>% html_nodes("div#search-title strong") %>%extract(3) %>% html_text() %>%as.numeric()
attrs<-  web %>% html_nodes("li a")  %>% html_attrs()
pages_num <- attrs[which(grepl(pattern = "\\d",x = attrs,perl=T))] %>%
  gsub(pattern = "\\?page=",replacement = "") %>% as.numeric() %>% max()

start.time <- Sys.time()
uniprot_id <- c()
for (page in 1:pages_num) {
  cat("\r pages: ",page," out of ",pages_num)
  url=paste0(int_url,"?page=",page)
  web=read_html(httr::RETRY("GET", url,times=1000,timeout(1000)))
  while(TRUE){
    web <- try(read_html(RETRY("GET", url,times=100,timeout(100))), silent=TRUE)
    if(!is(web, 'try-error')) break
  }
  uniprot_id <- c(uniprot_id,
                  web %>% html_nodes("div#iptm-search-report table tbody tr input") %>% html_attr("value"))
}
names(uniprot_id) <- uniprot_id
end.time <- Sys.time()
(end.time-start.time)


start.time <- Sys.time()
dat<-list()
.to_dataframe <- function(data) {
  con <- textConnection(data)
  dataframe <- utils::read.csv(con)
  return <- dataframe
}
for (i in 1:length(uniprot_id)){
  cat("\r proteins: ",i," out of ",length(uniprot_id))
  id <- uniprot_id[i]
  url <- sprintf("%s/%s/substrate",get_host_url(),id)
  while(TRUE){
    result <- try(RETRY("GET", url,httr::add_headers("Accept"="text/plain"),
                        times=100,timeout(100)), silent=TRUE)
    if(!is(result, 'try-error')) break
  }
  dat[[id]]  <- .to_dataframe(httr::content(result,"text",encoding = "UTF-8"))
}
end.time <- Sys.time()
(end.time-start.time)
