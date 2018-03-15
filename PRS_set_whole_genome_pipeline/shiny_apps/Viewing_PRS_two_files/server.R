#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(scales)
library(data.table)
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(shiny)
library(stringr)
library(grid)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

    My_data <- reactive({ req(input$file1)
      
      ## Read in data
      Full_data <- fread(input$file1$datapath)
      
      ## Create new columns parsing the identifiers in the Full_data score column and input to the shiny app
      Full_data[, Genesets := gsub(pattern = ".*_SCORE_(.*)_.*", replacement = "\\1", x = Full_data$score,perl = T)]
      Full_data[, Gene_regions := gsub(pattern = "^(.*)_geneset_SCORE_.*", replacement = "\\1", x = Full_data$score,perl = T)]
      Full_data[, Significance_thresholds := gsub(pattern = ".*_(.*$)", replacement = "\\1", x = Full_data$score,perl = T)]
      
      ## Create arguments to shiny app
      Gene.sets.input <- unique(Full_data$Genesets)
      significance_threshold.input <- unique(Full_data$Significance_thresholds)
      DSM.input <- unique(Full_data$samples.i.)
      
      ## Identify which rows in the data table contain whole genome information
      whole_genome_genic_positions_Full_data <- grep(x = Full_data$Gene_regions, pattern = "genic.genome_SCORE_whole_genome",perl = T)
      whole_genome_all_genome_positions_Full_data <- grep(x = Full_data$Gene_regions, pattern = "All.genome_SCORE_whole_genome",perl = T)
      
      whole_genome_plot_all_positions <- c(whole_genome_genic_positions_Full_data,whole_genome_all_genome_positions_Full_data)
      
      ## Change the Gene regions identifier to enable comparison of gene-set PRS to whole genome PRS
      Full_data[whole_genome_genic_positions_Full_data, Gene_regions := gsub(pattern = "^(.*).genic.genome_SCORE_.*", replacement = "\\1", x = Full_data$Gene_regions[whole_genome_genic_positions_Full_data],perl = T)]
      Full_data[whole_genome_all_genome_positions_Full_data, Gene_regions := gsub(pattern = ".*", replacement = "Full", x = Full_data$Gene_regions[whole_genome_all_genome_positions_Full_data],perl = T)]
      
      Full_data[whole_genome_plot_all_positions, Type := "Whole_genome"]
      Full_data[!whole_genome_plot_all_positions, Type:= "Pathway"]
      Full_data$logp <- -log10(Full_data$p)
      Full_data$SE_higher <- Full_data$estimate + Full_data$SE
      Full_data$SE_lower <- Full_data$estimate - Full_data$SE
      Full_data$r2_dir <- 100 * (as.numeric(Full_data$r.squared) *
                                   (sign(as.numeric(Full_data$estimate))))
      Full_data$p_value_text <- paste("p =", scientific(Full_data$p, digits = 2), sep = " ")
      
      
      ## Add alterations column to create "human readable" formats of the data
      alterations <- Full_data$score
      alterations[-whole_genome_plot_all_positions] <- str_replace(string = alterations[-whole_genome_plot_all_positions], pattern = ".*SCORE_(.*)_.*",replacement = "\\1")
      alterations[-whole_genome_plot_all_positions] <- str_replace(string = alterations[-whole_genome_plot_all_positions], pattern = "long_term_potentiation",replacement = "LTP")
      alterations[-whole_genome_plot_all_positions] <- str_replace(string = alterations[-whole_genome_plot_all_positions], pattern = "action_potential",replacement = "AP")
      alterations[-whole_genome_plot_all_positions] <- str_replace(string = alterations[-whole_genome_plot_all_positions], pattern = "depolarization",replacement = "DP")
      alterations[whole_genome_genic_positions_Full_data] <- "Whole Genome PRS GENIC"
      alterations[whole_genome_all_genome_positions_Full_data] <- "Whole Genome PRS ALL"
      
      Full_data[, alterations := alterations]
      
      
      output$Significance_threshold <- renderUI({
        checkboxGroupInput("Significance_threshold", label = "PRS P Value Threshold:",
                           choices = significance_threshold.input, selected = "0.05")
      })
      
      output$DSM <- renderUI({
        selectInput("DSM", label = "DSM type:",
                    choices = DSM.input)
      })
      
      output$geneset <- renderUI({
        checkboxGroupInput("geneset", label = "Geneset PRS to include:",
                           choices = Gene.sets.input)
      })
      
      Full_data
      
    })

    
  output$PvalPlot <- renderPlot({
      
      My_data()
      
      if (is.null(input$Significance_threshold)) {
        return(NULL)
      }    
      if (is.null(input$geneset)) {
        return(NULL)
      }    
      if (is.null(input$Gene_regions)) {
        return(NULL)
      }    
      
      ## Limit data table to input arguments and pipe to limiting columns and ordering based on significance
      sample_analysis <- My_data() %>%
        filter(samples.i. == input$DSM,
               Gene_regions %in% input$Gene_regions,
               Significance_thresholds %in% input$Significance_threshold,
               Genesets %in% input$geneset
        )
      
      ## Format DF to DT and apply fixes to the number of decimal points, format "g" = change to nn.dde-dd only if required
      
      ## Okay, so faceting was not meant to have differing axis lables, but in order to place the axis in the right order, I need to specify just one threshold and repeat across all facets
      ## I've used a short-cut here, the line 108 sorts the alterations column by the score and type and then only selects the unique labels for these columns so that the structure is "repeated" across all thresholds
      ## despite not knowing how many thresholds are in the analysis...i've saved a few lines of code and thought here.
      
      Sample_analysis_2 <- as.data.table(sample_analysis)
      Sample_analysis_2$score <- factor(Sample_analysis_2$score, levels = Sample_analysis_2$score[order(Sample_analysis_2$score, Sample_analysis_2$Type)])
      Sample_analysis_2$alterations <- factor(Sample_analysis_2$alterations, levels = unique(Sample_analysis_2$alterations[order(Sample_analysis_2$score, Sample_analysis_2$Type)]))
      Sample_analysis_2
      
      # Plot the resulting table for comparisons
      p <- ggplot(Sample_analysis_2, aes(x=score, y=logp, fill = Type, group=Significance_thresholds))
      
      p <- p +
        geom_point(aes(colour = Type))
      
      p <- p + scale_x_discrete(labels= levels(Sample_analysis_2$alterations))
      
      p <- p + facet_grid(. ~ Significance_thresholds,scales = "free_x", space = "free_x") +
        theme(strip.text.x = element_text(size = 10))
      
      p <- p + geom_hline(aes(yintercept=1.30103), colour = "red", linetype= "solid", alpha = 0.25)
      p <- p + scale_fill_brewer(palette = "Paired")
      p <- p + theme(axis.text.x = element_text(angle = 90,size = 10,hjust = 1,vjust = 0.5))
      p <- p + ggtitle(Sample_analysis_2$.id[1])
      #p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5,vjust = 0.5,angle = 90))
      p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
      p <- p + ylab(label = expression(-'log'[10]*'(p)'))
      p <- p + xlab(label = "Polygenic risk score")
      p
      #ggplotly(p) %>% 
      #  layout(height = input$plotHeight, autosize=TRUE)
      
      # Possible improvements:
      # Implement in switch from whole genome to gene-sets
      # Implement data-table of the raw results
      # Implement output file of the plots
      # Colour rows for significant values
      # Incorporate into its own app
      # 
      
    })
  output$Beta_plot <- renderPlot({
      My_data()
      
      if (is.null(input$Significance_threshold)) {
        return(NULL)
      }    
      if (is.null(input$geneset)) {
        return(NULL)
      }    
      if (is.null(input$Gene_regions)) {
        return(NULL)
      }    
      
      ## Limit data table to input arguments and pipe to limiting columns and ordering based on significance
      sample_analysis <- My_data() %>%
        filter(samples.i. == input$DSM,
               Gene_regions %in% input$Gene_regions,
               Significance_thresholds %in% input$Significance_threshold,
               Genesets %in% input$geneset
        )
      
      ## Format DF to DT and apply fixes to the number of decimal points, format "g" = change to nn.dde-dd only if required
      
      ## Okay, so faceting was not meant to have differing axis lables, but in order to place the axis in the right order, I need to specify just one threshold and repeat across all facets
      ## I've used a short-cut here, the line 108 sorts the alterations column by the score and type and then only selects the unique labels for these columns so that the structure is "repeated" across all thresholds
      ## despite not knowing how many thresholds are in the analysis...i've saved a few lines of code and thought here.
      
      Sample_analysis_2 <- as.data.table(sample_analysis)
      Sample_analysis_2$score <- factor(Sample_analysis_2$score, levels = Sample_analysis_2$score[order(Sample_analysis_2$score, Sample_analysis_2$Type)])
      Sample_analysis_2$alterations <- factor(Sample_analysis_2$alterations, levels = unique(Sample_analysis_2$alterations[order(Sample_analysis_2$score, Sample_analysis_2$Type)]))
      Sample_analysis_2
      
      # Put in the code below above, removing all of the excess alterations work to create the pdf plots...
      
      p <- ggplot(Sample_analysis_2, aes(x=score, y=estimate, fill = Type, group=Significance_thresholds))
      
      p <- p +
        geom_errorbar(aes(ymin = upper, ymax = lower), position = "dodge", width = 0.25) +
        geom_point(aes(colour = Type))
      
      p <- p + scale_x_discrete(labels= levels(Sample_analysis_2$alterations))
      p <- p + facet_grid(. ~ Significance_thresholds, scales = "free_x", space = "free_x") +
        theme(strip.text.x = element_text(size = 10))
      p <- p + geom_hline(aes(yintercept=0), colour = "red", linetype= "solid", alpha = 0.25)
      p <- p + scale_fill_brewer(palette = "Paired")
      p <- p + theme(axis.text.x = element_text(angle = 90,size = 10,hjust = 1,vjust = 0.5))
      p <- p + ggtitle(Sample_analysis_2$.id[1])
      p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
      p <- p + ylab(label = "BETA")
      p <- p + xlab(label = "Polygenic risk score")
      p
      # ggplotly(p) %>% 
      # layout(height = input$plotHeight, autosize=TRUE)
      
      # Possible improvements:
      # Implement in switch from whole genome to gene-sets
      # Implement data-table of the raw results
      # Implement output file of the plots
      # Colour rows for significant values
      # Incorporate into its own app
      # 
      
    })
  output$R2_plot <- renderPlot({
      My_data()
      
      if (is.null(input$Significance_threshold)) {
        return(NULL)
      }    
      if (is.null(input$geneset)) {
        return(NULL)
      }    
      if (is.null(input$Gene_regions)) {
        return(NULL)
      }    
      
      ## Limit data table to input arguments and pipe to limiting columns and ordering based on significance
      sample_analysis <- My_data() %>%
        filter(samples.i. == input$DSM,
               Gene_regions %in% input$Gene_regions,
               Significance_thresholds %in% input$Significance_threshold,
               Genesets %in% input$geneset
        )
      
      ## Format DF to DT and apply fixes to the number of decimal points, format "g" = change to nn.dde-dd only if required
      
      ## Okay, so faceting was not meant to have differing axis lables, but in order to place the axis in the right order, I need to specify just one threshold and repeat across all facets
      ## I've used a short-cut here, the line 108 sorts the alterations column by the score and type and then only selects the unique labels for these columns so that the structure is "repeated" across all thresholds
      ## despite not knowing how many thresholds are in the analysis...i've saved a few lines of code and thought here.
      
      Sample_analysis_2 <- as.data.table(sample_analysis)
      Sample_analysis_2$score <- factor(Sample_analysis_2$score, levels = Sample_analysis_2$score[order(Sample_analysis_2$score, Sample_analysis_2$Type)])
      Sample_analysis_2$alterations <- factor(Sample_analysis_2$alterations, levels = unique(Sample_analysis_2$alterations[order(Sample_analysis_2$score, Sample_analysis_2$Type)]))
      Sample_analysis_2
      
      p <- ggplot(Sample_analysis_2, aes(x=score, y=r2_dir, fill = Type, group=Significance_thresholds))
      p <- p +
        geom_bar(stat = "identity", aes(colour = Type), position = "dodge") +
        geom_text(data=subset(Sample_analysis_2, p < 0.05),
                  aes(x=score,y=r2_dir,label=p_value_text, hjust=ifelse(sign(r2_dir)>0, 0, 0)), angle = 90, position = position_dodge(width = 1), size = 2.9)
      
      #Problem with labels with a workaround
      # I use the score column in the format of factors and reference each relevant dataset for ggplot.
      # However this relies on having 0.05 and 0.5 in the value name.
      # scale_x_discrete accepts functions, but I also need to convert SCORE_0.05 and Score_0.5 into a "Whole_genome_PRS" which is almost impossible to write"
      # However as the labels function accepts key:value pairs, I wrote a vector in R that maps the original names of the pathways to "human readable" format using names function in R
      # This should work for most instances
      
      p <- p + scale_x_discrete(labels= levels(Sample_analysis_2$alterations))
      p <- p + scale_y_continuous(expand = expand_scale(mult = c(0.2,.6)))
      p <- p + facet_grid(. ~ Significance_thresholds, scales = "free_x", space = "free_x") +
        theme(strip.text.x = element_text(size = 10))
      p <- p + theme(axis.text.x = element_text(angle = 90, size = 10, hjust = 1,vjust = 0.5))
      p <- p + ggtitle(Sample_analysis_2$.id[1])
      p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
      p <- p + ylab(label = "R2_dir (%)")
      p <- p + xlab(label = "Polygenic risk score")
      p
      
      #ggplotly(p) %>% 
      #  layout(height = input$plotHeight, autosize=TRUE)
      
      # Possible improvements:
      # Implement in switch from whole genome to gene-sets
      # Implement data-table of the raw results
      # Implement output file of the plots
      # Colour rows for significant values
      # Incorporate into its own app
      # 
      
    })



  output$PvalPlot_step <- renderPlot({
    
    My_data()
    
    if (is.null(input$Significance_threshold)) {
      return(NULL)
    }    
    if (is.null(input$geneset)) {
      return(NULL)
    }    
    if (is.null(input$Gene_regions)) {
      return(NULL)
    }    
    
    ## Limit data table to input arguments and pipe to limiting columns and ordering based on significance
    sample_analysis <- My_data() %>%
      filter(samples.i. == input$DSM,
             Gene_regions %in% input$Gene_regions,
             Significance_thresholds %in% input$Significance_threshold,
             Genesets %in% input$geneset
      )
    
    ## Format DF to DT and apply fixes to the number of decimal points, format "g" = change to nn.dde-dd only if required
    
    ## Okay, so faceting was not meant to have differing axis lables, but in order to place the axis in the right order, I need to specify just one threshold and repeat across all facets
    ## I've used a short-cut here, the line 108 sorts the alterations column by the score and type and then only selects the unique labels for these columns so that the structure is "repeated" across all thresholds
    ## despite not knowing how many thresholds are in the analysis...i've saved a few lines of code and thought here.
    
    Sample_analysis_2 <- as.data.table(sample_analysis)
    Sample_analysis_2$score <- factor(Sample_analysis_2$score, levels = Sample_analysis_2$score[order(Sample_analysis_2$score, Sample_analysis_2$Type)])
    Sample_analysis_2$alterations <- factor(Sample_analysis_2$alterations, levels = unique(Sample_analysis_2$alterations[order(Sample_analysis_2$score, Sample_analysis_2$Type)]))
    Sample_analysis_2
    
    # Plot the resulting table for comparisons
    p <- ggplot(Sample_analysis_2, aes(x=score, y=logp, fill = Type, group=Significance_thresholds))
    
    p <- p +
      geom_point(aes(colour = Type))
    
    p <- p + scale_x_discrete(labels= levels(Sample_analysis_2$alterations))
    
    p <- p + facet_grid(. ~ Significance_thresholds,scales = "free_x", space = "free_x") +
      theme(strip.text.x = element_text(size = 10))
    
    p <- p + geom_hline(aes(yintercept=1.30103), colour = "red", linetype= "solid", alpha = 0.25)
    p <- p + scale_fill_brewer(palette = "Paired")
    p <- p + theme(axis.text.x = element_text(angle = 90,size = 10,hjust = 1,vjust = 0.5))
    p <- p + ggtitle(Sample_analysis_2$.id[1])
    #p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5,vjust = 0.5,angle = 90))
    p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
    p <- p + ylab(label = expression(-'log'[10]*'(p)'))
    p <- p + xlab(label = "Polygenic risk score")
    p
    #ggplotly(p) %>% 
    #  layout(height = input$plotHeight, autosize=TRUE)
    
    # Possible improvements:
    # Implement in switch from whole genome to gene-sets
    # Implement data-table of the raw results
    # Implement output file of the plots
    # Colour rows for significant values
    # Incorporate into its own app
    # 
    
  })
  output$Beta_plot_step <- renderPlot({
    My_data()
    
    if (is.null(input$Significance_threshold)) {
      return(NULL)
    }    
    if (is.null(input$geneset)) {
      return(NULL)
    }    
    if (is.null(input$Gene_regions)) {
      return(NULL)
    }    
    
    ## Limit data table to input arguments and pipe to limiting columns and ordering based on significance
    sample_analysis <- My_data() %>%
      filter(samples.i. == input$DSM,
             Gene_regions %in% input$Gene_regions,
             Significance_thresholds %in% input$Significance_threshold,
             Genesets %in% input$geneset
      )
    
    
    ## Format DF to DT and apply fixes to the number of decimal points, format "g" = change to nn.dde-dd only if required
    
    ## Okay, so faceting was not meant to have differing axis lables, but in order to place the axis in the right order, I need to specify just one threshold and repeat across all facets
    ## I've used a short-cut here, the line 108 sorts the alterations column by the score and type and then only selects the unique labels for these columns so that the structure is "repeated" across all thresholds
    ## despite not knowing how many thresholds are in the analysis...i've saved a few lines of code and thought here.
    
    Sample_analysis_2 <- as.data.table(sample_analysis)
    Sample_analysis_2$score <- factor(Sample_analysis_2$score, levels = Sample_analysis_2$score[order(Sample_analysis_2$score, Sample_analysis_2$Type)])
    Sample_analysis_2$alterations <- factor(Sample_analysis_2$alterations, levels = unique(Sample_analysis_2$alterations[order(Sample_analysis_2$score, Sample_analysis_2$Type)]))
    Sample_analysis_2
    
    if(any(grep(input$DSM, Sample_analysis_2$predictors_retained))){
      rows_to_be_highlighted <- grep(input$DSM, Sample_analysis_2$predictors_retained)
      Sample_analysis_2[, colour_change := "black"]
      Sample_analysis_2[rows_to_be_highlighted,colour_change := "green"]
      testing_sample_analysis <- Sample_analysis_2[order(Sample_analysis_2$Significance_thresholds,Sample_analysis_2$alterations)]
      d1 <- split(testing_sample_analysis$colour_change, ceiling(seq_along(testing_sample_analysis$colour_change)/length(levels(testing_sample_analysis$alterations))))
    }
    
    # Put in the code below above, removing all of the excess alterations work to create the pdf plots...
    
    p <- ggplot(Sample_analysis_2, aes(x=score, y=estimate, fill = Type, group=Significance_thresholds))
    
    p <- p +
      geom_errorbar(aes(ymin = upper, ymax = lower), position = "dodge", width = 0.25) +
      geom_point(aes(colour = Type))
    
    p <- p + scale_x_discrete(labels= levels(Sample_analysis_2$alterations))
    p <- p + facet_grid(. ~ Significance_thresholds, scales = "free_x", space = "free_x") +
      theme(strip.text.x = element_text(size = 10))
    p <- p + geom_hline(aes(yintercept=0), colour = "red", linetype= "solid", alpha = 0.25)
    p <- p + scale_fill_brewer(palette = "Paired")
    p <- p + theme(axis.text.x = element_text(angle = 90,size = 10,hjust = 1,vjust = 0.5))
    p <- p + ggtitle(Sample_analysis_2$.id[1])
    p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
    p <- p + ylab(label = "BETA")
    p <- p + xlab(label = "Polygenic risk score")
    pt <- ggplotGrob(p)
    
    integer <- length(input$Significance_threshold)*2 +1
    for (i in 1:length(input$Significance_threshold)){
    pt$grobs[[integer+i]]$children[[2]]$grobs[[2]]$children[[1]]$gp$col <- d1[[i]]
    }
    
    grid::grid.draw(pt)
    
    # ggplotly(p) %>% 
    # layout(height = input$plotHeight, autosize=TRUE)
    
    # Possible improvements:
    # Implement in switch from whole genome to gene-sets
    # Implement data-table of the raw results
    # Implement output file of the plots
    # Colour rows for significant values
    # Incorporate into its own app
    # 
    
  })
  output$R2_plot_step <- renderPlot({
    My_data()
    
    if (is.null(input$Significance_threshold)) {
      return(NULL)
    }    
    if (is.null(input$geneset)) {
      return(NULL)
    }    
    if (is.null(input$Gene_regions)) {
      return(NULL)
    }    
    
    ## Limit data table to input arguments and pipe to limiting columns and ordering based on significance
    sample_analysis <- My_data() %>%
      filter(samples.i. == input$DSM,
             Gene_regions %in% input$Gene_regions,
             Significance_thresholds %in% input$Significance_threshold,
             Genesets %in% input$geneset
      )
    
    ## Format DF to DT and apply fixes to the number of decimal points, format "g" = change to nn.dde-dd only if required
    
    ## Okay, so faceting was not meant to have differing axis lables, but in order to place the axis in the right order, I need to specify just one threshold and repeat across all facets
    ## I've used a short-cut here, the line 108 sorts the alterations column by the score and type and then only selects the unique labels for these columns so that the structure is "repeated" across all thresholds
    ## despite not knowing how many thresholds are in the analysis...i've saved a few lines of code and thought here.
    
    Sample_analysis_2 <- as.data.table(sample_analysis)
    Sample_analysis_2$score <- factor(Sample_analysis_2$score, levels = Sample_analysis_2$score[order(Sample_analysis_2$score, Sample_analysis_2$Type)])
    Sample_analysis_2$alterations <- factor(Sample_analysis_2$alterations, levels = unique(Sample_analysis_2$alterations[order(Sample_analysis_2$score, Sample_analysis_2$Type)]))
    Sample_analysis_2
    
    p <- ggplot(Sample_analysis_2, aes(x=score, y=r2_dir, fill = Type, group=Significance_thresholds))
    p <- p +
      geom_bar(stat = "identity", aes(colour = Type), position = "dodge") +
      geom_text(data=subset(Sample_analysis_2, p < 0.05),
                aes(x=score,y=r2_dir,label=p_value_text, hjust=ifelse(sign(r2_dir)>0, 0, 0)), angle = 90, position = position_dodge(width = 1), size = 2.9)
    
    #Problem with labels with a workaround
    # I use the score column in the format of factors and reference each relevant dataset for ggplot.
    # However this relies on having 0.05 and 0.5 in the value name.
    # scale_x_discrete accepts functions, but I also need to convert SCORE_0.05 and Score_0.5 into a "Whole_genome_PRS" which is almost impossible to write"
    # However as the labels function accepts key:value pairs, I wrote a vector in R that maps the original names of the pathways to "human readable" format using names function in R
    # This should work for most instances
    
    p <- p + scale_x_discrete(labels= levels(Sample_analysis_2$alterations))
    p <- p + scale_y_continuous(expand = expand_scale(mult = c(0.2,.6)))
    p <- p + facet_grid(. ~ Significance_thresholds, scales = "free_x", space = "free_x") +
      theme(strip.text.x = element_text(size = 10))
    p <- p + theme(axis.text.x = element_text(angle = 90, size = 10, hjust = 1,vjust = 0.5))
    p <- p + ggtitle(Sample_analysis_2$.id[1])
    p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
    p <- p + ylab(label = "R2_dir (%)")
    p <- p + xlab(label = "Polygenic risk score")
    p
    
    #ggplotly(p) %>% 
    #  layout(height = input$plotHeight, autosize=TRUE)
    
    # Possible improvements:
    # Implement in switch from whole genome to gene-sets
    # Implement data-table of the raw results
    # Implement output file of the plots
    # Colour rows for significant values
    # Incorporate into its own app
    # 
    
  })

  output$summary_table <- renderDataTable({
      
      My_data()
      
      # These are required in case no tick boxes are selected
      if (is.null(input$Significance_threshold)) {
        return(NULL)
      }    
      if (is.null(input$geneset)) {
        return(NULL)
      }    
      if (is.null(input$Gene_regions)) {
        return(NULL)
      }    
      
      # Select columns you wish to output
      cols <- c("estimate", "SE","r.squared","p")
      
      ## Limit data table to input arguments and pipe to limiting columns and ordering based on significance
      sample_analysis <- My_data() %>%
        filter(samples.i. == input$DSM,
               Gene_regions %in% input$Gene_regions,
               Significance_thresholds %in% input$Significance_threshold,
               Genesets %in% input$geneset
        )  %>%
        select(c(score:SE,p,r.squared)) %>%
        arrange(p)
      
      ## Format DF to DT and apply fixes to the number of decimal points, format "g" = change to nn.dde-dd only if required
      Sample_analysis_2 <- as.data.table(sample_analysis)
      Sample_analysis_2[, (cols) := lapply(.SD, formatC, digits = 3, format = "g"), .SDcols = cols]
      
      ## leave datatable function in for "prettyfying" the final result    
      datatable(data = Sample_analysis_2,
                options = list(pageLength = 10),
                rownames = F)
      
      # Possible improvements:
      # colour rows for significant values
      # incorporate into its own app
      # 
    })
  })