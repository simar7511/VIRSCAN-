library(shiny)
library(shinyjs)
library(shinyBS)
library(readxl)
library(DT)
library(dplyr)
library(ggplot2)
library(reshape)
library(shinydashboard)
library(tidyr)
library(ggforce)
library(scales)
library(stringr)



workingDir <- "/Users/Simar/VIRSCAN/data/"
setwd(workingDir)
load("qcInput.RData")


qcInput <- qcInput %>% 
  mutate(Assay = case_when(Plate %in% c(1, 2) ~ 1,
                           Plate %in% c(3, 4) ~ 2))
qcInput$Sample.Annotation[qcInput$Sample.Annotation == "Emperical"] <- "Empirical"

# Define a function to calculate the quality call based on criteria
calculate_quality_call_replicate <- function(originalQCreport) {
  quality_call <- originalQCreport %>% 
    # select(SampleID, reads_mapped, readDepth, pctDetected, CPMcor, ZscoreCor, sumAllSpecies) %>% 
    mutate(alignedReadsQC = case_when(reads_mapped >= 200000 ~ "Good",
                                      reads_mapped < 200000 & reads_mapped >= 160000 ~ "PossibleGood",
                                      reads_mapped < 160000 & reads_mapped >= 100001 ~ "Questionable",
                                      reads_mapped <= 100000 ~ "Fail"),
           readDepthQC = case_when(readDepth >= 4 ~ "Good",
                                   readDepth < 4 & readDepth >= 3.5 ~ "PossibleGood",
                                   readDepth < 3.5 & readDepth >= 2.5 ~ "Questionable",
                                   readDepth < 2.5 ~ "Fail"),
           # need to consider between reps for Questionable and Fail!!! ----
           pctDetectedQC = case_when(pctDetected >= 40 ~ "Good",
                                     pctDetected < 40 & pctDetected >= 30 ~ "Questionable",
                                     pctDetected < 30 ~ "Fail"),
           CPMcorQC = case_when(CPMcor >= .25 ~ "Good",
                                CPMcor < .25 ~ "Questionable",
                                is.na(CPMcor) ~ "Fail"),
           ZscoreCorQC = case_when(ZscoreCor >= .25 ~ "Good",
                                   ZscoreCor < .25 ~ "Fail",
                                   is.na(ZscoreCor) ~ "Fail"),
           sumAllSpeciesQC = case_when(sumAllSpecies >= 100 ~ "Good",
                                       sumAllSpecies < 100 & sumAllSpecies >= 50 ~ "PossibleGood",
                                       sumAllSpecies < 50 ~ "Fail")) %>% 
    group_by(SampleID) %>% 
    mutate(lowestQCvalTmp = paste(alignedReadsQC, readDepthQC, pctDetectedQC, CPMcorQC, ZscoreCorQC, sumAllSpeciesQC, collapse = " ")) %>% 
    data.frame() %>% 
    mutate(lowestQCval = case_when(grepl("Fail", lowestQCvalTmp) ~ "Fail",
                                   grepl("Questionable", lowestQCvalTmp) ~ "Questionable",
                                   grepl("PossibleGood", lowestQCvalTmp) ~ "PossibleGood",
                                   grepl(" Good|Good ", lowestQCvalTmp) ~ "Good")) %>% 
    dplyr::rename(SampleAnnotation = Sample.Annotation) %>% 
    # select(SampleID, alignedReadsQC, readDepthQC, pctDetectedQC, CPMcorQC, ZscoreCorQC, sumAllSpeciesQC, lowestQCval)
    select(SampleID, 
           alignedReadsQC, readDepthQC, pctDetectedQC, CPMcorQC, ZscoreCorQC, sumAllSpeciesQC, lowestQCval,
           raw_total_sequences,
           reads_mapped,
           pctReadsAligned,
           error_rate,
           average_quality,
           peptidesGTE_15,
           peptidesDetected,
           pctDetected,
           alignments,
           readDepth,
           CPMcor,
           ZscoreCor,
           Human.herpesvirus.1, Human.herpesvirus.2, Human.herpesvirus.3, Human.herpesvirus.4, Human.herpesvirus.5,
           Human.respiratory.syncytial.virus, Influenza.A.virus, Influenza.B.virus, Rhinovirus.A, Rhinovirus.B, Streptococcus.pneumoniae,
           Measles.virus, Simian.foamy.virus, Aravan.virus,
           sumAllSpecies,
           SampleAnnotation,
           PoolID)
  return(quality_call)
}

# Apply the replicate quality assessment function
qualityCallReplicate <- calculate_quality_call_replicate(qcInput)


# Sample level QC scores
# function to calculate quality scores for samples
calculate_quality_call_sample <- function(originalQCreport) {
  quality_call <- originalQCreport %>% 
    select(PoolID, Sample.Annotation, VirScanID, SampleID, reads_mapped, readDepth, pctDetected, CPMcor, ZscoreCor, sumAllSpecies) %>% 
    group_by(VirScanID) %>% 
    mutate(minMappedReads = min(reads_mapped),
           maxMappedReads = max(reads_mapped),
           minReadDepth = min(readDepth),
           maxReadDepth = max(readDepth),
           corCPM = min(CPMcor),
           corZ = min(ZscoreCor),
           minPctDetected = min(pctDetected),
           maxPctDetected = max(pctDetected),
           SumAllSpecies = min(sumAllSpecies),
           # mappedReadsFC = round(maxMappedReads/minMappedReads, 1),
           pctDetectedFC = round(maxPctDetected/minPctDetected, 1)) %>% 
    data.frame() %>% 
    select(VirScanID,
           PoolID,
           Sample.Annotation,
           minMappedReads,
           maxMappedReads,
           # mappedReadsFC,
           minReadDepth,
           maxReadDepth,
           minPctDetected,
           maxPctDetected,
           pctDetectedFC,
           corCPM, corZ,
           SumAllSpecies) %>% 
    unique() %>% 
    dplyr:::rename(CPMcor = corCPM,
                   ZscoreCor = corZ,
                   sumAllSpecies = SumAllSpecies,
                   SampleAnnotation = Sample.Annotation) %>%
    mutate(FinalQualityCall = case_when(minMappedReads >= 200000 & maxMappedReads >= 200000 & 
                                          minReadDepth >= 4 & maxReadDepth >= 4 & 
                                          minPctDetected >= 40 & maxPctDetected >= 40 &
                                          pctDetectedFC <= 2 &
                                          CPMcor >= .25 & 
                                          sumAllSpecies >= 100 ~ "Good",
                                        minMappedReads >= 200000 & maxMappedReads >= 200000 & 
                                          minReadDepth >= 4 & maxReadDepth >= 4 & 
                                          minPctDetected < 40 &
                                          pctDetectedFC <= 2 &
                                          CPMcor >= .25 & 
                                          sumAllSpecies >= 100 ~ "Good",
                                        sumAllSpecies < 50 ~ "Fail",
                                        minMappedReads < 100000 ~ "Fail",
                                        sumAllSpecies < 50 & CPMcor < .25 ~ "Fail",
                                        minReadDepth >= 2.5 & minReadDepth < 3.5 & CPMcor < .6 ~ "Fail",
                                        minReadDepth >= 2.5 & minReadDepth < 3.5 & maxReadDepth >= 4 & CPMcor >= .6 & pctDetectedFC <= 2 ~ "Questionable",
                                        minReadDepth >= 2.5 & minReadDepth < 3.5 & maxReadDepth >= 4 & CPMcor >= .6 & pctDetectedFC > 2 ~ "Fail",
                                        minReadDepth < 2.5 ~ "Fail",
                                        minPctDetected < 30 & maxPctDetected < 30 & pctDetectedFC <= 2 ~ "Good",
                                        minPctDetected < 30 & maxPctDetected < 30 & pctDetectedFC > 2 ~ "Fail",
                                        minReadDepth >= 3.5 & minReadDepth < 4 & pctDetectedFC <= 2 & CPMcor >= .6 ~ "Good",
                                        minMappedReads >= 160000 & minMappedReads < 200000 & maxMappedReads >= 200000 & CPMcor >= .6 & pctDetectedFC <= 2 ~ "Good",
                                        minMappedReads >= 160000 & minMappedReads < 200000 & maxMappedReads >= 200000 & pctDetectedFC > 2 ~ "Fail",
                                        minMappedReads >= 160000 & minMappedReads < 200000 & maxMappedReads >= 200000 & CPMcor < .6 ~ "Fail",
                                        CPMcor < .25 ~ "Questionable",
                                        minMappedReads >= 200000 & maxMappedReads >= 200000 & minReadDepth >= 4 & maxReadDepth >= 4 & CPMcor >= .25 & pctDetectedFC > 2 ~ "Questionable",
                                        sumAllSpecies < 100 & sumAllSpecies >= 50 ~ "Questionable",
                                        (minMappedReads < 200000 & minMappedReads >= 160000 & maxMappedReads < 200000 & maxMappedReads >= 160000) | 
                                          (minReadDepth < 3.5 & minReadDepth >= 2.5 & maxReadDepth < 3.5 & maxReadDepth >= 2.5) ~ "Questionable")) %>% 
    mutate(FinalQualityCallNotes = case_when(minMappedReads >= 200000 & maxMappedReads >= 200000 & 
                                               minReadDepth >= 4 & maxReadDepth >= 4 & 
                                               minPctDetected >= 40 & maxPctDetected >= 40 &
                                               pctDetectedFC <= 2 &
                                               CPMcor >= .25 & 
                                               sumAllSpecies >= 100 ~ "minMappedReads >= 200000; maxMappedReads >= 200000; minReadDepth >= 4; maxReadDepth >= 4; minPctDetected >= 40 & maxPctDetected >= 40; pctDetectedFC <= 2; CPMcor >= .25",
                                             minMappedReads >= 200000 & maxMappedReads >= 200000 & 
                                               minReadDepth >= 4 & maxReadDepth >= 4 & 
                                               minPctDetected < 40 &
                                               pctDetectedFC <= 2 &
                                               CPMcor >= .25 & 
                                               sumAllSpecies >= 100 ~ "minMappedReads >= 200000; maxMappedReads >= 200000; minReadDepth >= 4; maxReadDepth >= 4; minPctDetected < 40; pctDetectedFC <= 2; CPMcor >= .25",
                                             sumAllSpecies < 50 ~ "sumAllSpecies < 50",
                                             minMappedReads < 100000 ~ "minMappedReads < 100000",
                                             sumAllSpecies < 50 & CPMcor < .25 ~ "sumAllSpecies < 50; CPMcor < .25",
                                             minReadDepth >= 2.5 & minReadDepth < 3.5 & CPMcor < .6 ~ "minReadDepth >= 2.5 & < 3.5; CPMcor < .6",
                                             minReadDepth >= 2.5 & minReadDepth < 3.5 & maxReadDepth >= 4 & CPMcor >= .6 & pctDetectedFC <= 2 ~ "minReadDepth >= 2.5 & < 3.5; maxReadDepth >=4; CPMcor < .6; pctDetectedFC <= 2",
                                             minReadDepth >= 2.5 & minReadDepth < 3.5 & maxReadDepth >= 4 & CPMcor >= .6 & pctDetectedFC > 2 ~ "minReadDepth >= 2.5 & < 3.5; maxReadDepth >=4; CPMcor < .6; pctDetectedFC > 2",
                                             minReadDepth < 2.5 ~ "minReadDepth < 2.5",
                                             minPctDetected < 30 & maxPctDetected < 30 & pctDetectedFC <= 2 ~ "minPctDetected < 30; maxPctDetected < 30; pctDetectedFC <= 2",
                                             minPctDetected < 30 & maxPctDetected < 30 & pctDetectedFC > 2 ~ "minPctDetected < 30; maxPctDetected < 30; pctDetectedFC > 2",
                                             minReadDepth >= 3.5 & minReadDepth < 4 & pctDetectedFC <= 2 & CPMcor >= .6 ~ "minReadDepth >= 3.5 & < 4; pctDetectedFC <= 2; CPMcor >= .6",
                                             minMappedReads >= 160000 & minMappedReads < 200000 & maxMappedReads >= 200000 & CPMcor >= .6 & pctDetectedFC <= 2 ~ "minMappedReads >= 160000 & < 200000; maxMappedReads >= 200000; CPMcor >= .6; pctDetectedFC <= 2",
                                             minMappedReads >= 160000 & minMappedReads < 200000 & maxMappedReads >= 200000 & pctDetectedFC > 2 ~ "minMappedReads >= 160000 & < 200000; maxMappedReads >= 200000; pctDetectedFC > 2",
                                             minMappedReads >= 160000 & minMappedReads < 200000 & maxMappedReads >= 200000 & CPMcor < .6 ~ "minMappedReads >= 160000 & < 200000; maxMappedReads >= 200000; CPMcor < .6",
                                             CPMcor < .25 ~ "CPMcor < .25",
                                             minMappedReads >= 200000 & maxMappedReads >= 200000 & minReadDepth >= 4 & maxReadDepth >= 4 & CPMcor >= .25 & pctDetectedFC > 2 ~ " minMappedReads >= 200000; maxMappedReads >= 200000; minReadDepth >= 4; maxReadDepth >= 4; CPMcor >= .25; pctDetectedFC > 2",
                                             sumAllSpecies < 100 & sumAllSpecies >= 50 ~ "sumAllSpecies < 100 & >= 50",
                                             (minMappedReads < 200000 & minMappedReads >= 160000 & maxMappedReads < 200000 & maxMappedReads >= 160000) | 
                                               (minReadDepth < 3.5 & minReadDepth >= 2.5 & maxReadDepth < 3.5 & maxReadDepth >= 2.5) ~ "(minMappedReads < 200000 & >= 160000; maxMappedReads < 200000 & >= 160000) | (minReadDepth < 3.5 & >= 2.5; maxReadDepth < 3.5 & >= 2.5)")) %>% 
    select(VirScanID, 
           minMappedReads,
           maxMappedReads,
           minReadDepth,
           maxReadDepth,
           minPctDetected,
           maxPctDetected,
           pctDetectedFC,
           CPMcor,
           ZscoreCor,
           sumAllSpecies,
           FinalQualityCall,
           FinalQualityCallNotes,
           SampleAnnotation,
           PoolID)
}

qualityCallSample <- calculate_quality_call_sample(qcInput)


# Apply the quality assessment function
# qualityCalls <- calculate_quality_call(qcInput)

# Accessing the results
# overallQuality <- qualityCalls$overall
# sampleQuality <- qualityCalls$sample

# Control Samples Analysis
control_analysis <- qcInput %>% 
  filter(grepl("Control", Sample.Annotation)) %>% 
  group_by(PoolID) %>% 
  mutate(
    mean_SumAllSpecies = mean(sumAllSpecies),
    mean_ViralScore = mean(Human.herpesvirus.1),  # Choose a relevant QC virus
    mean_mapped_reads = mean(reads_mapped),
    mean_percent_reads_aligned = mean(pctReadsAligned),
    mean_corr_counts = mean(CPMcor),
    mean_corr_pct_detected = mean(pctDetected)
  ) %>% 
  select(Sample.Annotation, starts_with("mean_")) %>% 
  unique()

# Samples Analysis
samples_analysis <- qcInput %>% 
  filter(grepl("Emperical", Sample.Annotation)) %>% 
  group_by(PoolID) %>% 
  mutate(
    mean_mapped_reads = mean(reads_mapped),
    median_pct_reads_aligned = median(pctReadsAligned),
    mean_corr_counts = mean(CPMcor),
    mean_corr_pct_detected = mean(pctDetected)
  ) %>% 
  select(Sample.Annotation, starts_with("mean_")) %>% 
  unique()

# Mock IP Analysis
mock_analysis <- qcInput %>% 
  filter(grepl("Mock-IP", Sample.Annotation)) %>% 
  group_by(PoolID) %>% 
  mutate(
    mean_mapped_reads = mean(reads_mapped),
    median_pct_reads_aligned = median(pctReadsAligned),
    mean_corr_counts = mean(CPMcor),
    mean_corr_pct_detected = mean(pctDetected)
  ) %>% 
  select(Sample.Annotation, starts_with("mean_")) %>% 
  unique()

# Batch Level Comparison
batch_comparison <- qcInput %>% 
  group_by(PoolID) %>% 
  summarize(
    num_controls = sum(grepl("Control", Sample.Annotation)),
    num_samples = sum(grepl("Empirical", Sample.Annotation)),
    num_mock = sum(grepl("Mock-IP", Sample.Annotation)),
    mean_SumAllSpecies = mean(sumAllSpecies),
    mean_ViralScore = mean(Human.herpesvirus.1),  
    mean_mapped_reads = mean(reads_mapped),
    median_pct_reads_aligned = median(pctReadsAligned),
    mean_CPMcor = mean(CPMcor),
    mean_pct_detected = mean(pctDetected)
  ) %>% 
  ungroup() %>% 
  mutate(mean_SumAllSpecies = round(mean_SumAllSpecies, 1),
         mean_ViralScore = round(mean_ViralScore, 1),
         mean_mapped_reads = round(mean_mapped_reads, 1),
         mean_CPMcor = round(mean_CPMcor, 1),
         mean_pct_detected = round(mean_pct_detected, 1))
batchComparison <- melt(data.frame(batch_comparison))


# Assay level comparison ----
assay_comparison <- qcInput %>% 
  mutate(
    Assay = case_when(
      Plate %in% c(1, 2) ~ 1,
      Plate %in% c(3, 4) ~ 2)) %>% 
  group_by(PoolID, Assay) %>% 
  summarize(
    num_controls = sum(grepl("Control", Sample.Annotation)),
    num_samples = sum(grepl("Empirical", Sample.Annotation)),
    num_mock = sum(grepl("Mock-IP", Sample.Annotation)),
    mean_SumAllSpecies = mean(sumAllSpecies),
    mean_ViralScore = mean(Human.herpesvirus.1),  # Choose a relevant QC virus ----
    mean_mapped_reads = mean(reads_mapped),
    median_pct_reads_aligned = median(pctReadsAligned),
    mean_CPMcor = mean(CPMcor),
    mean_pct_detected = mean(pctDetected)
  ) %>% 
  ungroup() %>% 
  mutate(mean_SumAllSpecies = round(mean_SumAllSpecies, 1),
         mean_ViralScore = round(mean_ViralScore, 1),
         mean_mapped_reads = round(mean_mapped_reads, 1),
         mean_CPMcor = round(mean_CPMcor, 1),
         mean_pct_detected = round(mean_pct_detected, 1))

assayComparison <- assay_comparison %>%
  mutate(AssayID = paste(PoolID, Assay, sep = "_")) %>%
  select(AssayID, everything()) %>%
  as.data.frame() %>%
  melt() %>%
  mutate(Assay = as.numeric(sub("^.+_", "", AssayID))) %>%
  select(AssayID, PoolID, Assay, variable, value)


# Plate level comparison ----
plate_comparison <- qcInput %>% 
  group_by(PoolID, Plate) %>% 
  summarize(
    num_controls = sum(grepl("Control", Sample.Annotation)),
    num_samples = sum(grepl("Empirical", Sample.Annotation)),
    num_mock = sum(grepl("Mock-IP", Sample.Annotation)),
    mean_SumAllSpecies = mean(sumAllSpecies),
    mean_ViralScore = mean(Human.herpesvirus.1),  # Choose a relevant QC virus
    mean_mapped_reads = mean(reads_mapped),
    median_pct_reads_aligned = median(pctReadsAligned),
    mean_CPMcor= mean(CPMcor),
    mean_pct_detected = mean(pctDetected)
  ) %>%
  ungroup() %>% 
  mutate(mean_SumAllSpecies = round(mean_SumAllSpecies, 1),
         mean_ViralScore = round(mean_ViralScore, 1),
         mean_mapped_reads = round(mean_mapped_reads, 1),
         mean_CPMcor = round(mean_CPMcor, 1),
         mean_pct_detected = round(mean_pct_detected, 1))

plateComparison <- plate_comparison %>% 
  mutate(PlateID = paste(PoolID, Plate, sep = "_")) %>% 
  select(PlateID, everything()) %>% 
  data.frame() %>% 
  melt() %>% 
  mutate(Plate = sub("^.+_", "",PlateID)) %>% 
  select(PlateID, PoolID, Plate, variable, value)


# Shiny UI ----
ui <- dashboardPage(
  dashboardHeader(),
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("Main Tab", tabName = "main"),
      menuItem(
        checkboxGroupInput(
          inputId = "selectedPoolID",
          label = "Select Batch(es)",
          choices = unique(qcInput$PoolID)
        )
      ),
      menuItem("Batch", tabName = "Batch"),
      menuItem("Assay", tabName = "Assay"),
      menuItem("Plate", tabName = "Plate"),
      menuItem("Plate Map", tabName = "PlateMap"),
      menuItem("Sample", tabName = "Sample"),
      menuItem("Replicate", tabName = "Replicate")
    )
  ),
  dashboardBody(
    # tags$head(tags$style(HTML('
    #   .main-header .logo {
    #     font-family: "Georgia", Times, "Times New Roman", serif;
    #     font-weight: bold;
    #     font-size: 24px;
    #   }
    # '))),
    useShinyjs(),
    tabItems(
      tabItem(
        tabName = "main",
        "Main Tab Content"
      ),
      tabItem(
        tabName = "Batch",
        body <- dashboardBody(
          fluidRow(
            box(width = 3,
                radioButtons(inputId = "metricBatch",
                             label = "Select Metric",
                             choices = c("num_controls",
                                         "num_samples",
                                         "num_mock",
                                         "mean_SumAllSpecies",
                                         "mean_ViralScore",
                                         "mean_mapped_reads",
                                         "median_pct_reads_aligned",
                                         "mean_CPMcor",
                                         "mean_pct_detected"))),
            box(plotOutput("batchPlot"), width = 9)
          ),
          fluidRow(
            dataTableOutput("batchQC")
          )
        )
      ),
      tabItem(
        tabName = "Assay",
        body <- dashboardBody(
          fluidRow(
            box(width = 3,
                checkboxGroupInput(inputId = "selectedAssay",
                                   label = "Selected Assay(s)",
                                   choices = unique(qcInput$Assay)),
                radioButtons(inputId = "metricAssay",
                             label = "Select Metric",
                             choices = c("num_controls",
                                         "num_samples",
                                         "num_mock",
                                         "mean_SumAllSpecies",
                                         "mean_ViralScore",
                                         "mean_mapped_reads",
                                         "median_pct_reads_aligned",
                                         "mean_CPMcor",
                                         "mean_pct_detected"))),
            box(plotOutput("AssayPlot"), width = 9)
          ),
          fluidRow(
            dataTableOutput("assayQC")
          )
        )
      ),
      tabItem(
        tabName = "Plate",
        body <- dashboardBody(
          fluidRow(
            box(width = 3,
                checkboxGroupInput(inputId = "selectedPlate",
                                   label = "Selected Plate(s)",
                                   choices = unique(qcInput$Plate)),
                radioButtons(inputId = "metricPlate",
                             label = "Select Metric",
                             choices = c("num_controls",
                                         "num_samples",
                                         "num_mock",
                                         "mean_SumAllSpecies",
                                         "mean_ViralScore",
                                         "mean_mapped_reads",
                                         "median_pct_reads_aligned",
                                         "mean_CPMcor",
                                         "mean_pct_detected"))),
            box(plotOutput("platePlot"), width = 9)
          ),
          fluidRow(
            dataTableOutput("plateQC")
          )
        )
      ),
      tabItem(
        tabName = "PlateMap",
        sidebarPanel(width = 3,
                     radioButtons(inputId = "selectedPoolIDplateMap",
                                  label = "Select Pool",
                                  choices = unique(qcInput$PoolID)),
                     radioButtons(inputId = "metricPlateMap",
                                  label = "Select Metric",
                                  choices = names(qcInput[6:36]))
        ),
        mainPanel(width = 9,
                  fluidRow(
                    splitLayout(cellWidths = c("50%", "50%"), plotOutput("plate96_1"), plotOutput("plate96_2"))
                  ),
                  fluidRow(
                    splitLayout(cellWidths = c("50%", "50%"), plotOutput("plate96_3"), plotOutput("plate96_4"))
                  )
        )
      ),
      tabItem(
        tabName = "Sample",
        body <- dashboardBody(
          fluidRow(
            box(width = 3,
                checkboxGroupInput(inputId = "sampleAnnotationType",
                                   label = "Sample Type",
                                   choices = unique(qcInput$Sample.Annotation)),
                radioButtons(inputId = "metricSample",
                             label = "Select Metric",
                             choices = c("minMappedReads",
                                         "maxMappedReads",
                                         "minReadDepth",
                                         "maxReadDepth",
                                         "minPctDetected",
                                         "maxPctDetected",
                                         "pctDetectedFC",
                                         "CPMcor",
                                         "ZscoreCor",
                                         "sumAllSpecies"))),
            box(plotOutput("qualityCallSamplePlot"), width = 9)
          ),
          fluidRow(
            dataTableOutput("qualityCallSampleTable")
          )
        )
      ),
      tabItem(
        tabName = "Replicate",
        body <- dashboardBody(
          fluidRow(
            box(width = 3,
                checkboxGroupInput(inputId = "replicateAnnotationType",
                                   label = "Sample Type",
                                   choices = unique(qcInput$Sample.Annotation)),
                radioButtons(inputId = "metricReplicate",
                             label = "Select Metric",
                             choices = c("raw_total_sequences",
                                         "reads_mapped",
                                         "pctReadsAligned",
                                         "error_rate",
                                         "average_quality",
                                         "peptidesGTE_15",
                                         "peptidesDetected",
                                         "pctDetected",
                                         "alignments",
                                         "readDepth",
                                         "CPMcor",
                                         "ZscoreCor"))),
            box(plotOutput("qualityCallReplicatePlot"), width = 9)
          ),
          fluidRow(
            dataTableOutput("qualityCallReplicateTable")
          )
        )
      ) 
    )   
  )
)

# Shiny server ----
# Function to filter data based on selected PoolID
# Helper function to filter data based on selected PoolID
filterDataByPoolID <- function(data, selectedPoolID) {
  data[data$PoolID %in% selectedPoolID, ]
}

# Helper function to render assay level plots
renderAssayPlots <- function(assayData, variable, x, y) {
  ggplot(assayData, aes_string(x, y)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    ylab(variable)
}

renderBatchPlots <- function(batchData, variable, x, y) {
  ggplot(batchData, aes_string(x, y)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    ylab(variable)
}


# Helper function to render plate level plots
renderPlatePlots <- function(plateData, variable, x, y) {
  ggplot(plateData, aes_string(x, y)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    ylab(variable)
}


# Helper function to render 96 well plate plots
renderPlate96Plots <- function(selectTableMap, title, subtitle, metricPlateMap) {
  ggplot(data = selectTableMap()) +
    labs(title = paste(title, ": ", subtitle, sep = " "), subtitle = metricPlateMap, x = "", y = "") +
    geom_circle(aes(x0 = col, y0 = row, r = 0.45, fill = value)) +
    coord_equal() +
    scale_x_continuous(breaks = 1:12, expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(breaks = 1:8, labels = LETTERS[1:8], expand = expansion(mult = c(0.01, 0.01)), trans = reverse_trans()) +
    scale_fill_gradient(low = "white", high = "purple") +
    geom_text(aes(x = col, y = row, label = paste0(value)), size = 3) +
    labs(fill = "")
}

# Helper function to render sample level plots
renderSamplePlots <- function(sampleData, x, y) {
  ggplot(sampleData, aes_string(x, y)) +
    geom_violin() +
    theme_bw() +
    ylab(y)
}

# Generic function to render plots
renderPlots <- function(plotData, x, y) {
  ggplot(plotData, aes_string(x, y)) +
    geom_violin() +
    theme_bw() +
    ylab(y)
}


server <- function(input, output, session) {
  
  # Assay level table ----
  assayQCselected <- reactive({
    assay_comparison[assay_comparison$PoolID %in% input$selectedPoolID,] %>% filter(Assay %in% input$selectedAssay)
  })
  
  output$assayQC <- renderDataTable({
    req(input$selectedPoolID)
    DT::datatable(assayQCselected(),
                  extensions = c('FixedColumns',"FixedHeader"),
                  rownames = FALSE,
                  class = 'cell-border stripe',
                  filter = 'top',
                  options = list(autoWidth = TRUE,
                                 scrollX = TRUE,
                                 fixedHeader = TRUE,
                                 fixedColumns = list(leftColumns = 1, rightColumns = 0)))
  })
  
  
  
  # Assay level plots ----
  output$AssayPlot <- renderPlot({
    renderAssayPlots(
      filterDataByPoolID(assayComparison, input$selectedPoolID) %>% filter(variable == input$metricAssay,
                                                                           Assay %in% input$selectedAssay),
      input$metricAssay,
      "AssayID",
      "value"
    )
  })
  
  
  # Batch level table ----
  batchQCselected <- reactive({
    batch_comparison[batch_comparison$PoolID %in% input$selectedPoolID,]
  })
  
  output$batchQC <- renderDataTable({
    req(input$selectedPoolID)
    DT::datatable(batchQCselected(),
                  extensions = c('FixedColumns',"FixedHeader"),
                  rownames = FALSE,
                  class = 'cell-border stripe',
                  filter = 'top',
                  options = list(autoWidth = TRUE,
                                 scrollX = TRUE,
                                 fixedHeader = TRUE,
                                 fixedColumns = list(leftColumns = 1, rightColumns = 0)))
  })
  
  
  # Batch level plots ----
  output$batchPlot <- renderPlot({
    renderBatchPlots(
      filterDataByPoolID(batchComparison, input$selectedPoolID) %>% filter(variable == input$metricBatch),
      input$metricBatch,
      "PoolID",
      "value"
    )
  })
  
  
  
  # Plate level table ----
  plateQCselected <- reactive({
    plate_comparison[plate_comparison$PoolID %in% input$selectedPoolID,] %>% filter(Plate %in% input$selectedPlate)
  })
  
  output$plateQC <- renderDataTable({
    req(input$selectedPoolID)
    DT::datatable(plateQCselected(),
                  extensions = c('FixedColumns',"FixedHeader"),
                  rownames = FALSE,
                  class = 'cell-border stripe',
                  filter = 'top',
                  options = list(autoWidth = TRUE,
                                 scrollX = TRUE,
                                 fixedHeader = TRUE,
                                 fixedColumns = list(leftColumns = 1, rightColumns = 0)))
  })
  
  # Plate level plots ----
  output$platePlot <- renderPlot({
    renderPlatePlots(
      filterDataByPoolID(plateComparison, input$selectedPoolID) %>% filter(variable == input$metricPlate,
                                                                           Plate %in% input$selectedPlate),
      input$metricPlate,
      "PlateID",
      "value"
    )
  })
  
  
  # 96 well plate plots ----
  selectTableMap_1 <- reactive({
    # p1 <- qcInput[qcInput$PoolID %in% input$selectedPoolID,] %>% 
    p1 <- qcInput[qcInput$PoolID %in% input$selectedPoolIDplateMap,] %>%
      # filter(Plate == input$selectedPlate) %>%
      filter(Plate == 1) %>% 
      mutate(sampleIndex = as.numeric(sub("^.+_", "", VirScanID))) %>%
      arrange(rep, sampleIndex) %>% 
      select(input$metricPlateMap)
    p1.matrix <- matrix(p1[,1], nrow = 8, ncol = 12, byrow = TRUE)
    p1.df <- as.data.frame(p1.matrix)
    p1.df %>%
      mutate(row = 1:8) %>%
      pivot_longer(-row, names_to = "col", values_to = "value") %>%
      mutate(col = as.integer(str_remove(col, "V"))) %>%
      mutate(well = paste0(LETTERS[row], col)) %>% 
      data.frame()
  })
  
  selectTableMap_2 <- reactive({
    # p1 <- qcInput[qcInput$PoolID %in% input$selectedPoolID,] %>% 
    p1 <- qcInput[qcInput$PoolID %in% input$selectedPoolIDplateMap,] %>%
      # filter(Plate == input$selectedPlate) %>%
      filter(Plate == 2) %>% 
      mutate(sampleIndex = as.numeric(sub("^.+_", "", VirScanID))) %>%
      arrange(rep, sampleIndex) %>% 
      select(input$metricPlateMap)
    p1.matrix <- matrix(p1[,1], nrow = 8, ncol = 12, byrow = TRUE)
    p1.df <- as.data.frame(p1.matrix)
    p1.df %>%
      mutate(row = 1:8) %>%
      pivot_longer(-row, names_to = "col", values_to = "value") %>%
      mutate(col = as.integer(str_remove(col, "V"))) %>%
      mutate(well = paste0(LETTERS[row], col)) %>% 
      data.frame()
  })
  
  selectTableMap_3 <- reactive({
    # p1 <- qcInput[qcInput$PoolID %in% input$selectedPoolID,] %>% 
    p1 <- qcInput[qcInput$PoolID %in% input$selectedPoolIDplateMap,] %>%
      # filter(Plate == input$selectedPlate) %>%
      filter(Plate == 3) %>% 
      mutate(sampleIndex = as.numeric(sub("^.+_", "", VirScanID))) %>%
      arrange(rep, sampleIndex) %>% 
      select(input$metricPlateMap)
    p1.matrix <- matrix(p1[,1], nrow = 8, ncol = 12, byrow = TRUE)
    p1.df <- as.data.frame(p1.matrix)
    p1.df %>%
      mutate(row = 1:8) %>%
      pivot_longer(-row, names_to = "col", values_to = "value") %>%
      mutate(col = as.integer(str_remove(col, "V"))) %>%
      mutate(well = paste0(LETTERS[row], col)) %>% 
      data.frame()
  })
  
  selectTableMap_4 <- reactive({
    # p1 <- qcInput[qcInput$PoolID %in% input$selectedPoolID,] %>% 
    p1 <- qcInput[qcInput$PoolID %in% input$selectedPoolIDplateMap,] %>%
      # filter(Plate == input$selectedPlate) %>%
      filter(Plate == 4) %>% 
      mutate(sampleIndex = as.numeric(sub("^.+_", "", VirScanID))) %>%
      arrange(rep, sampleIndex) %>% 
      select(input$metricPlateMap)
    p1.matrix <- matrix(p1[,1], nrow = 8, ncol = 12, byrow = TRUE)
    p1.df <- as.data.frame(p1.matrix)
    p1.df %>%
      mutate(row = 1:8) %>%
      pivot_longer(-row, names_to = "col", values_to = "value") %>%
      mutate(col = as.integer(str_remove(col, "V"))) %>%
      mutate(well = paste0(LETTERS[row], col)) %>% 
      data.frame()
  })
  
  
  
  output$plate96_1 <- renderPlot({
    renderPlate96Plots(selectTableMap_1, input$selectedPoolIDplateMap, "Plate 1", input$metricPlateMap)
  })
  
  output$plate96_2 <- renderPlot({
    renderPlate96Plots(selectTableMap_2, input$selectedPoolIDplateMap, "Plate 2", input$metricPlateMap)
  })
  
  output$plate96_3 <- renderPlot({
    renderPlate96Plots(selectTableMap_3, input$selectedPoolIDplateMap, "Plate 3", input$metricPlateMap)
  })
  
  output$plate96_4 <- renderPlot({
    renderPlate96Plots(selectTableMap_4, input$selectedPoolIDplateMap, "Plate 4", input$metricPlateMap)
  })
  
  # Sample level plots ----
  output$qualityCallSamplePlot <- renderPlot({
    renderSamplePlots(
      filterDataByPoolID(qualityCallSample, input$selectedPoolID) %>% filter(SampleAnnotation %in% input$sampleAnnotationType),
      "PoolID",
      input$metricSample
    )
  })
  
  
  # Sample level tables ----
  sampleQCselected <- reactive({
    x <- qualityCallSample[qualityCallSample$PoolID %in% input$selectedPoolID,]
    x %>% filter(SampleAnnotation %in% input$sampleAnnotationType)
  })
  
  output$qualityCallSampleTable <- renderDataTable({
    req(input$selectedPoolID)
    DT::datatable(sampleQCselected(),
                  extensions = c('FixedColumns',"FixedHeader"),
                  rownames = FALSE,
                  class = 'cell-border stripe',
                  filter = 'top',
                  options = list(autoWidth = TRUE,
                                 scrollX = TRUE,
                                 fixedHeader = TRUE,
                                 fixedColumns = list(leftColumns = 1, rightColumns = 0)))
  })
  
  
  # Replicate level plots ----
  output$qualityCallReplicatePlot <- renderPlot({
    renderPlots(
      filterDataByPoolID(qualityCallReplicate, input$selectedPoolID) %>% filter(SampleAnnotation %in% input$replicateAnnotationType),
      "PoolID",
      input$metricReplicate
    )
  })  
  
  # Replicate level tables ----
  replicateQCselected <- reactive({
    x <- qualityCallReplicate[qualityCallReplicate$PoolID %in% input$selectedPoolID,]
    x %>% filter(SampleAnnotation %in% input$replicateAnnotationType)
  })
  
  output$qualityCallReplicateTable <- renderDataTable({
    req(input$selectedPoolID)
    DT::datatable(replicateQCselected(),
                  extensions = c('FixedColumns',"FixedHeader"),
                  rownames = FALSE,
                  class = 'cell-border stripe',
                  filter = 'top',
                  options = list(autoWidth = TRUE,
                                 scrollX = TRUE,
                                 fixedHeader = TRUE,
                                 fixedColumns = list(leftColumns = 1, rightColumns = 0)))
  })
  
  
}

shinyApp(ui, server)