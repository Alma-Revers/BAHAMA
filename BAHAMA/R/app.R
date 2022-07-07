
plotBAHAMA <- function(tidy_rr_pt, tidy_rr_hlt, tidy_rr_hlgt, tidy_rr_soc){
    tidy_rr_pt <- as.data.frame(tidy_rr_pt)
    tidy_rr_hlt <- as.data.frame(tidy_rr_hlt)
    tidy_rr_hlgt <- as.data.frame(tidy_rr_hlgt)
    tidy_rr_soc <- as.data.frame(tidy_rr_soc)

# Define UI for application
ui <- shiny::fluidPage(
    sidebarLayout(
        sidebarPanel(
            selectInput(inputId = "soc",
                        "SOC", choices = c("All",levels(as.factor(tidy_rr_soc$label)))),
            selectInput(inputId = "hlgt",
                        "HLGT", choices = c("All",levels(as.factor(tidy_rr_hlgt$label)))),
            selectInput(inputId = "hlt",
                        "HLT", choices = c("All",levels(as.factor(tidy_rr_hlt$label)))),
            selectInput(inputId = "pt",
                        "PT", choices = c("All",levels(as.factor(tidy_rr_pt$label)))),
            selectInput(inputId = "table_input",
                        "Table to show", choices = c("SOC", "HLGT", "HLT", "PT")),
        ),
        mainPanel(plotlyOutput(outputId = "plotOutput_pt"),
                  plotlyOutput(outputId = "plotOutput_hlt"),
                  plotlyOutput(outputId = "plotOutput_hlgt"),
                  plotlyOutput(outputId = "plotOutput_soc"),
                  dataTableOutput(outputId = "table"))
    )
)

# Define server logic required
server <- function(input, output) {
    output$plotOutput_pt <- renderPlotly({
        data_plot_pt <- data.frame(RR_PT_x = tidy_rr_pt$estimate,
                                RR_PT_y = abs(tidy_rr_pt$estimate/tidy_rr_pt$std.error),
                                PT_label = tidy_rr_pt$label,
                                soc = tidy_rr_pt$soc,
                                hlgt = tidy_rr_pt$hlgt,
                                hlt = tidy_rr_pt$hlt,
                                size = tidy_rr_pt$abs_diff)
        p_pt <- ggplot(data_plot_pt, aes(x = RR_PT_x, y = RR_PT_y, colour=soc, label=PT_label, size=size))+
            theme_bw() +
            labs(x = "log(RR)", y = "log(RR)/SD", tag = "PT-level")+
            geom_vline(xintercept = 0, linetype = "dashed", alpha=0.4) +
            scale_colour_hue(h=c(25, 275), limits = levels(tidy_rr_pt$soc))

        if(input$soc == "All" & input$hlgt == "All" & input$hlt == "All" & input$pt == "All"){
            p_pt <- p_pt + geom_point()
            plotly::ggplotly(p_pt, tooltip = c("label"))
        }
        else if( (!input$soc == "All") & input$hlgt == "All" & input$hlt == "All" & input$pt == "All"){
            alpha_pt <- rep(0.5, length(as.factor(tidy_rr_pt$soc)))
            alpha_pt[data_plot_pt$soc == as.character(input$soc)] <- 1
            data_plot_pt$alpha  <- alpha_pt
            p_pt <- p_pt + geom_point(aes(alpha = alpha), data = data_plot_pt)
            plotly::ggplotly(p_pt, tooltip = c("label"))
        }
        else if( (!input$hlgt == "All") & input$hlt == "All" & input$pt == "All" ){
            alpha_pt <- rep(0.5, length(as.factor(tidy_rr_pt$hlgt)))
            alpha_pt[as.factor(data_plot_pt$hlgt) == input$hlgt] <- 1
            data_plot_pt$alpha  <- alpha_pt
            p_pt <- p_pt + geom_point(aes(alpha = alpha), data = data_plot_pt)
            plotly::ggplotly(p_pt, tooltip = c("label"))
        }
        else if( (!input$hlt == "All") & input$pt == "All" ){
            alpha <- rep(0.5, length(as.factor(tidy_rr_pt$hlt)))
            alpha[as.factor(data_plot_pt$hlgt) == input$hlt] <- 1
            data_plot_pt$alpha  <- alpha
            p_pt <- p_pt + geom_point(aes(alpha = alpha), data = data_plot_pt)
            plotly::ggplotly(p_pt, tooltip = c("label"))
        }
        else {
            alpha <- rep(0.5, length(as.factor(tidy_rr_pt$pt)))
            alpha[as.factor(data_plot_pt$hlgt) == input$pt] <- 1
            data_plot_pt$alpha  <- alpha
            p_pt <- p_pt + geom_point(aes(alpha = alpha), data = data_plot_pt)
            plotly::ggplotly(p_pt, tooltip = c("label"))
        }

    })
    output$plotOutput_hlt <- renderPlotly({
        data_plot_hlt <- data.frame(RR_HLT_x = tidy_rr_hlt$estimate,
                                    RR_HLT_y = abs(tidy_rr_hlt$estimate/tidy_rr_hlt$std.error),
                                    HLT_label = tidy_rr_hlt$label,
                                    soc = tidy_rr_hlt$soc,
                                    hlgt = as.factor(tidy_rr_hlt$hlgt),
                                    size = tidy_rr_hlt$abs_diff)
        p_hlt <- ggplot(data_plot_hlt, aes(x = RR_HLT_x, y = RR_HLT_y, color=soc, label=HLT_label, size=size))+
            theme_bw() +
            labs(x = "log(RR)", y = "log(RR)/SD", tag = "HLT-level")+
            geom_vline(xintercept = 0, linetype = "dashed", alpha=0.4) +
            guides(alpha = FALSE, size=FALSE, colour = guide_legend(override.aes = list(size=10)))+
            scale_colour_hue(h=c(25, 275), limits = levels(tidy_rr_hlt$soc))
        if(input$soc == "All" & input$hlgt == "All" & input$hlt == "All" & input$pt == "All"){
            p_hlt <- p_hlt + geom_point()
            plotly::ggplotly(p_hlt, tooltip = c("label"))
        }
        else if( (!input$soc == "All") & input$hlgt == "All" & input$hlt == "All" & input$pt == "All"){
            alpha_hlt <- rep(0.5, length(as.factor(tidy_rr_hlt$soc)))
            alpha_hlt[as.factor(data_plot_hlt$soc) == input$soc] <- 1
            data_plot_hlt$alpha  <- alpha_hlt
            p_hlt <- p_hlt + geom_point(aes(alpha = alpha), data = data_plot_hlt)
            plotly::ggplotly(p_hlt, tooltip = c("label"))
        }
        else if( (!input$hlgt == "All") & input$hlt == "All" & input$pt == "All" ){
            alpha_hlt <- rep(0.5, length(as.factor(tidy_rr_hlt$hlgt)))
            alpha_hlt[as.factor(data_plot_hlt$hlgt) == input$hlgt] <- 1
            data_plot_hlt$alpha  <- alpha_hlt
            p_hlt <- p_hlt + geom_point(aes(alpha = alpha), data = data_plot_hlt)
            plotly::ggplotly(p_hlt, tooltip = c("label"))
        }
        else if( (!input$hlt == "All") & input$pt == "All" ){
            alpha_hlt <- rep(0.5, length(as.factor(tidy_rr_hlt$label)))
            alpha_hlt[as.factor(data_plot_hlt$hlgt) == input$hlt] <- 1
            data_plot_hlt$alpha  <- alpha_hlt
            p_hlt <- p_hlt + geom_point(aes(alpha = alpha), data = data_plot_hlt)
            plotly::ggplotly(p_hlt, tooltip = c("label"))
        }
        else {
            alpha_hlt <- rep(0.5, length(as.factor(tidy_rr_hlt$label)))
            pt_hlt <- unique(tidy_rr_pt$hlt[as.factor(tidy_rr_pt$pt) == input$pt])
            alpha_hlt[as.factor(data_plot_hlt$label) %in% pt_hlt] <- 1
            data_plot_hlt$alpha  <- alpha_hlt
            p_hlt <- p_hlt + geom_point(aes(alpha = alpha), data = data_plot_hlt)
            plotly::ggplotly(p_hlt, tooltip = c("label"))
        }

    })

    output$plotOutput_hlgt <- renderPlotly({
        data_plot_hlgt <- data.frame(RR_HLGT_x = tidy_rr_hlgt$estimate,
                                     RR_HLGT_y = abs(tidy_rr_hlgt$estimate/tidy_rr_hlgt$std.error),
                                     HLGT_label = tidy_rr_hlgt$label,
                                     soc = tidy_rr_hlgt$soc,
                                     size = tidy_rr_hlgt$abs_diff)
        p_hlgt <- ggplot(data_plot_hlgt, aes(x = RR_HLGT_x, y = RR_HLGT_y, color=soc, label=HLGT_label, size=size))+
            theme_bw() +
            labs(x = "log(RR)", y = "log(RR)/SD", tag = "HLGT-level")+
            geom_vline(xintercept = 0, linetype = "dashed", alpha=0.4) +
            guides(alpha = FALSE, size=FALSE, colour = guide_legend(override.aes = list(size=10)))+
            scale_colour_hue(h=c(25, 275), limits = levels(tidy_rr_hlgt$soc))

        if(input$soc == "All" & input$hlgt == "All" & input$hlt == "All" & input$pt == "All"){
            p_hlgt <- p_hlgt + geom_point()
            plotly::ggplotly(p_hlgt, tooltip = c("label"))
        }
        else if( (!input$soc == "All") & input$hlgt == "All" & input$hlt == "All" & input$pt == "All"){
            alpha_hlgt <- rep(0.5, length(as.factor(tidy_rr_hlgt$soc)))
            alpha_hlgt[as.factor(data_plot_hlgt$soc) == input$soc] <- 1
            data_plot_hlgt$alpha  <- alpha_hlgt
            p_hlgt <- p_hlgt + geom_point(aes(alpha = alpha), data = data_plot_hlgt)
            plotly::ggplotly(p_hlgt, tooltip = c("label"))
        }
        else if( (!input$hlgt == "All") & input$hlt == "All" & input$pt == "All" ){
            alpha_hlgt <- rep(0.5, length(as.factor(tidy_rr_hlgt$label)))
            alpha_hlgt[as.factor(data_plot_hlgt$label) == input$hlgt] <- 1
            data_plot_hlgt$alpha  <- alpha_hlgt
            p_hlgt <- p_hlgt + geom_point(aes(alpha = alpha), data = data_plot_hlgt)
            plotly::ggplotly(p_hlgt, tooltip = c("label"))
        }
        else if( (!input$hlt == "All") & input$pt == "All" ){
            alpha_hlgt <- rep(0.5, length(as.factor(tidy_rr_hlgt$label)))
            HLT_hlgt <- unique(tidy_rr_hlt$hlgt[as.factor(tidy_rr_hlt$hlt) == input$hlt])
            alpha_hlgt[data_plot_hlgt$label %in% HLT_hlgt] <- 1
            data_plot_hlgt$alpha  <- alpha_hlgt
            p_hlgt <- p_hlgt + geom_point(aes(alpha = alpha), data = data_plot_hlgt)
            plotly::ggplotly(p_hlgt, tooltip = c("label"))
        }
        else {
            alpha_hlgt <- rep(0.5, length(tidy_rr_hlgt$label))
            pt_hlgt <- unique(tidy_rr_pt$hlgt[as.factor(tidy_rr_pt$pt) == input$pt])
            alpha_hlgt[data_plot_hlgt$label %in% pt_hlgt] <- 1
            data_plot_hlgt$alpha  <- alpha_hlgt
            p_hlgt <- p_hlgt + geom_point(aes(alpha = alpha), data = data_plot_hlgt)
            plotly::ggplotly(p_hlgt, tooltip = c("label"))
        }

    })

    output$plotOutput_soc <- renderPlotly({
        data_plot_soc <- data.frame(RR_SOC_x = tidy_rr_soc$estimate,
                                    RR_SOC_y = abs(tidy_rr_soc$estimate/tidy_rr_soc$std.error),
                                    SOC_label = tidy_rr_soc$label,
                                    size = tidy_rr_soc$abs_diff)
        p_soc <- ggplot(data_plot_soc, aes(x = RR_SOC_x, y = RR_SOC_y, color=SOC_label, label=SOC_label, size=size))+
            theme_bw() +
            labs(x = "log(RR)", y = "log(RR)/SD", tag = "SOC-level")+
            geom_vline(xintercept = 0, linetype = "dashed", alpha=0.4) +
            guides(alpha = FALSE, size=FALSE, colour = guide_legend(override.aes = list(size=10)))+
            scale_colour_hue(h=c(25, 275), limits = levels(tidy_rr_soc$label))

        if(input$soc == "All" & input$hlgt == "All" & input$hlt == "All" & input$pt == "All"){
            p_soc <- p_soc + geom_point()
            plotly::ggplotly(p_soc, tooltip = c("label"))
        }
        else if( (!input$soc == "All") & input$hlgt == "All" & input$hlt == "All" & input$pt == "All"){
            alpha_soc <- rep(0.5, length(as.factor(tidy_rr_soc$label)))
            alpha_soc[as.factor(data_plot_soc$SOC_label) == input$soc] <- 1
            data_plot_soc$alpha  <- alpha_soc
            p_soc <- p_soc + geom_point(aes(alpha = alpha), data = data_plot_soc)
            plotly::ggplotly(p_soc, tooltip = c("label"))
        }
        else if( (!input$hlgt == "All") & input$hlt == "All" & input$pt == "All" ){
            alpha_soc <- rep(0.5, length(as.factor(tidy_rr_soc$label)))
            HLGT_SOC <- unique(tidy_rr_hlgt$soc[as.factor(tidy_rr_hlgt$label) == input$hlgt])
            alpha_soc[as.factor(data_plot_soc$SOC_label) %in% HLGT_SOC] <- 1
            data_plot_soc$alpha  <- alpha_soc
            p_soc <- p_soc + geom_point(aes(alpha = alpha), data = data_plot_soc)
            plotly::ggplotly(p_soc, tooltip = c("label"))
        }
        else if( (!input$hlt == "All") & input$pt == "All" ){
            alpha_soc <- rep(0.5, length(as.factor(tidy_rr_soc$label)))
            HLT_SOC <- unique(tidy_rr_hlt$soc[as.factor(tidy_rr_hlt$label) == input$hlt])
            alpha_soc[as.factor(data_plot_soc$SOC_label) %in% HLT_SOC] <- 1
            data_plot_soc$alpha  <- alpha_soc
            p_soc <- p_soc + geom_point(aes(alpha = alpha), data = data_plot_soc)
            plotly::ggplotly(p_soc, tooltip = c("label"))
        }
        else {
            alpha_soc <- rep(0.5, length(tidy_rr_soc$label))
            pt_soc <- unique(tidy_rr_pt$soc[as.factor(tidy_rr_pt$label) == input$pt])
            alpha_soc[data_plot_soc$SOC_label %in% pt_soc] <- 1
            data_plot_soc$alpha  <- alpha_soc
            p_soc <- p_soc + geom_point(aes(alpha = alpha), data = data_plot_soc)
            plotly::ggplotly(p_soc, tooltip = c("label"))
        }

    })

    output$table <- renderDataTable({
        if(input$table_input == "PT"){
            dataset_toPrint <- data.frame("PT" =  tidy_rr_pt$label,
                                          "Y control" = tidy_rr_pt$y_pt_c,
                                          "Y treatment" = tidy_rr_pt$y_pt_t,
                                          "log(RR)" = round(tidy_rr_pt$estimate, digits = 3),
                                          "SD log(RR)" = round(tidy_rr_pt$std.error, digits = 3),
                                          "rhat" = round(tidy_rr_pt$rhat, digits = 2),
                                          "ess" = tidy_rr_pt$ess)
            if(input$soc == "All" & input$hlgt == "All" & input$hlt == "All" & input$pt == "All"){
                dataset_toPrint
            }
            else if(input$soc != "All" & input$hlgt == "All" & input$hlt == "All" & input$pt == "All"){
                dataset_toPrint[,as.factor(tidy_rr_pt$soc) == input$soc]
                dataset_toPrint
            }
            else if(input$hlgt != "All" & input$hlt == "All" & input$pt == "All"){
                dataset_toPrint[,as.factor(tidy_rr_pt$hlgt) == input$hlgt]
                dataset_toPrint
            }
            else if(input$hlt != "All" & input$pt == "All"){
                dataset_toPrint[,as.factor(tidy_rr_pt$hlt) == input$hlt]
                dataset_toPrint
            }
            else{
                dataset_toPrint[,as.factor(tidy_rr_pt$pt) == input$pt]
                dataset_toPrint
            }
        }

        if(input$table_input == "HLT"){
            dataset_toPrint <- data.frame("HLT" =  tidy_rr_hlt$label,
                                          "Y control" = tidy_rr_hlt$y_hlt_c,
                                          "Y treatment" = tidy_rr_hlt$y_hlt_t,
                                          "log(RR)" = round(tidy_rr_hlt$estimate, digits = 3),
                                          "SD log(RR)" = round(tidy_rr_hlt$std.error, digits = 3),
                                          "rhat" = round(tidy_rr_hlt$rhat, digits = 2),
                                          "ess" = tidy_rr_hlt$ess)

            if(input$soc == "All" & input$hlgt == "All" & input$hlt == "All" & input$pt == "All"){
                dataset_toPrint
            }
            else if(input$soc != "All" & input$hlgt == "All" & input$hlt == "All" & input$pt == "All"){
                dataset_toPrint <- dataset_toPrint[tidy_rr_hlt$soc == input$soc,]
                dataset_toPrint
            }
            else if(input$hlgt != "All" & input$hlt == "All" & input$pt == "All"){
                dataset_toPrint <- dataset_toPrint[tidy_rr_hlt$hlgt == input$hlgt,]
                dataset_toPrint
            }
            else if(input$hlt != "All" & input$pt == "All"){
                dataset_toPrint <- dataset_toPrint[tidy_rr_hlt$label == input$hlt,]
                dataset_toPrint
            }
            else{
                pt_hlt <- unique(tidy_rr_pt$hlt[tidy_rr_pt$pt == input$pt])
                dataset_toPrint <- dataset_toPrint[tidy_rr_hlt$label %in% pt_hlt,]
                dataset_toPrint
            }
        }
        if(input$table_input == "HLGT"){
            dataset_toPrint <- data.frame("HLGT" =  tidy_rr_hlgt$label,
                                          "Y control" = tidy_rr_hlgt$y_hlgt_c,
                                          "Y treatment" = tidy_rr_hlgt$y_hlgt_t,
                                          "log(RR)" = round(tidy_rr_hlgt$estimate, digits = 3),
                                          "SD log(RR)" = round(tidy_rr_hlgt$std.error, digits = 3),
                                          "rhat" = round(tidy_rr_hlgt$rhat, digits = 2),
                                          "ess" = tidy_rr_hlgt$ess)

            if(input$soc == "All" & input$hlgt == "All" & input$hlt == "All" & input$pt == "All"){
                dataset_toPrint
            }
            else if(input$soc != "All" & input$hlgt == "All" & input$hlt == "All" & input$pt == "All"){
                dataset_toPrint <- dataset_toPrint[tidy_rr_hlgt$soc == input$soc,]
                dataset_toPrint
            }
            else if(input$hlgt != "All" & input$hlt == "All" & input$pt == "All"){
                dataset_toPrint <- dataset_toPrint[tidy_rr_hlgt$label == input$hlgt,]
                dataset_toPrint
            }
            else if(input$hlt != "All" & input$pt == "All"){
                hlt_hlgt <- unique(tidy_rr_hlt$hlgt[tidy_rr_hlt$hlt == input$hlt])
                dataset_toPrint <- dataset_toPrint[tidy_rr_hlgt$label %in% hlt_hlgt,]
                dataset_toPrint
            }
            else{
                pt_hlgt <- unique(tidy_rr_pt$hlgt[tidy_rr_pt$pt == input$pt])
                dataset_toPrint <- dataset_toPrint[tidy_rr_hlgt$label %in% pt_hlgt,]
                dataset_toPrint
            }
        }

        if(input$table_input == "SOC"){
            dataset_toPrint <- data.frame("SOC" =  tidy_rr_soc$label,
                                          "Y control" = tidy_rr_soc$y_soc_c,
                                          "Y treatment" = tidy_rr_soc$y_soc_t,
                                          "log(RR)" = round(tidy_rr_soc$estimate, digits = 3),
                                          "SD log(RR)" = round(tidy_rr_soc$std.error, digits = 3),
                                          "rhat" = round(tidy_rr_soc$rhat, digits = 2),
                                          "ess" = tidy_rr_soc$ess)

            if(input$soc == "All" & input$hlgt == "All" & input$hlt == "All" & input$pt == "All"){
                dataset_toPrint
            }
            else if(input$soc != "All" & input$hlgt == "All" & input$hlt == "All" & input$pt == "All"){
                dataset_toPrint <- dataset_toPrint[as.factor(tidy_rr_soc$label) == input$soc,]
                dataset_toPrint
            }
            else if(input$hlgt != "All" & input$hlt == "All" & input$pt == "All"){
                hlgt_soc <- unique(tidy_rr_hlgt$soc[as.factor(tidy_rr_hlgt$label) == input$hlgt])
                dataset_toPrint <- dataset_toPrint[tidy_rr_soc$label %in% hlgt_soc,]
                dataset_toPrint
            }
            else if(input$hlt != "All" & input$pt == "All"){
                hlt_soc <- unique(tidy_rr_hlt$soc[as.factor(tidy_rr_hlt$label) == input$hlt])
                dataset_toPrint <- dataset_toPrint[tidy_rr_soc$label %in% hlt_soc,]
                dataset_toPrint
            }
            else{
                pt_soc <- unique(tidy_rr_pt$soc[as.factor(tidy_rr_pt$label) == input$pt])
                dataset_toPrint <- dataset_toPrint[tidy_rr_soc$label %in% pt_soc,]
                dataset_toPrint
            }
        }
    })
}

# Run the application
shiny::shinyApp(ui = ui, server = server)

}
