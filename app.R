library(shiny)
library(shinydashboard)
library(rgdal)
library(DT)
library(leaflet)
library(leaflet.extras)
library(RColorBrewer)
library(devtools)
library(spdep)
library(dashboardthemes)
library(spdep)
library(rgeos)
library(mboost)
library(spatialreg)
library(httr)

#install_github("nik01010/dashboardthemes")
#setwd("D:/Skripsweet/Shiny/Pemilihan Matrix/data")
#IGI17<-readOGR(".","jateng_IDI17")

#setwd("~/Skripsweet/Shiny/Pemilihan___Matrix")
IGI17 <- readOGR("data/jateng_IDI17.shp",
                 layer = "jateng_IDI17",
                 GDAL1_integer64_policy = TRUE)
IGI17$id <- 1:nrow(IGI17)

colnames(IGI17@data)[1] <- "KODE"
colnames(IGI17@data)[2] <- "PROV"
colnames(IGI17@data)[3] <- "KABKOT"

dataIGI <- IGI17@data
rownames(dataIGI) <- dataIGI$KABKOT
vy <- dataIGI$IGI
vx <- dataIGI[, 5:9]
coords <- cbind(dataIGI$x_long, dataIGI$y_lat)

dataYBox <- readxl::read_xlsx("www/dataYBox.xlsx")
dataBox <- readxl::read_xlsx("www/dataBox.xlsx")
skemaBox <- readxl::read_xlsx("www/skemaBox.xlsx")

migrasi <- read.csv("data/migrasi_selisih.csv")
rownames(migrasi) <- colnames(migrasi)
w_migrasi <- as.matrix(migrasi)
tpt <- read.csv("data/Matriks_TPT_inv.csv")
rownames(tpt) <- colnames(tpt)
w_tpt <- as.matrix(tpt)
pdrb <- read.csv("data/Matriks_PDRB_Berlaku_inverse.csv")
rownames(pdrb) <- colnames(pdrb)
w_pdrb <- as.matrix(pdrb)

kk_migrasi <- read.csv("data/kk_migrasi_selisih.csv")
rownames(kk_migrasi) <- colnames(kk_migrasi)

kk_tpt <- read.csv("data/kk_Matriks_TPT_inv.csv")
rownames(kk_tpt) <- colnames(kk_tpt)

kk_pdrb <- read.csv("data/kk_Matriks_PDRB_Berlaku_inverse.csv")
rownames(kk_pdrb) <- colnames(kk_pdrb)


w1 <- get(load("data/IGI_jateng.w1.Rdata"))
w2 <- get(load("data/IGI_jateng.w2.Rdata"))
w3 <- get(load("data/IGI_jateng.w3.Rdata"))
w4 <- get(load("data/IGI_jateng.w4.Rdata"))
w5 <- get(load("data/IGI_jateng.w5.Rdata"))
w6 <- get(load("data/IGI_jateng.w6.Rdata"))
w7 <- get(load("data/IGI_jateng.w7.Rdata"))
w8 <- get(load("data/IGI_jateng.w8.Rdata"))
w9 <- get(load("data/IGI_jateng.w9.Rdata"))
w10 <- get(load("data/IGI_jateng.w10.Rdata"))
w11 <- get(load("data/IGI_jateng.w11.Rdata"))
w12 <- get(load("data/IGI_jateng.w12.Rdata"))

w1 <- as.matrix(w1)
w2 <- as.matrix(w2)
w3 <- as.matrix(w3)
w4 <- as.matrix(w4)
w5 <- as.matrix(w5)
w6 <- as.matrix(w6)
w7 <- as.matrix(w7)
w8 <- as.matrix(w8)
w9 <- as.matrix(w9)
w10 <- as.matrix(w10)
w11 <- as.matrix(w11)
w12 <- as.matrix(w12)

skema123 <- function(n, w, w_custom) {
    #Neighborhood
    knn <- knearneigh(coords, k = n)
    knn_nb <- knn2nb(knn)
    dlist <- nbdists(knn_nb, coords, longlat = TRUE)
    indis <- function(x)
        1 / (x ^ w)
    dlist <- lapply(dlist, indis)
    bobot1 <-
        nb2listw(knn_nb,
                 glist = dlist,
                 style = "B",
                 zero.policy = T)
    bobot1_mat <- listw2mat(bobot1)
    bobot_custom <- w_custom * bobot1_mat
    b_cust_tot <- rowSums(bobot_custom, na.rm = TRUE)
    bobot_custom <- bobot_custom / b_cust_tot
    w_final <- mat2listw(bobot_custom, style = "W")
    w_final_nb <- w_final[["neighbours"]]
    #Moran's I
    moran <- moran.test(dataIGI$IGI, w_final, randomisation = TRUE)
    #Spatial Lag Model
    r.ols <-
        lm(
            dataIGI$IGI ~ dataIGI$INFLASI + dataIGI$PMTB + dataIGI$UMK + dataIGI$PPS +
                dataIGI$PP
        )
    r.lag <- lagsarlm(r.ols, dataIGI, w_final, zero.policy = TRUE)
    slagm <- summary(r.lag, Nagelkerke = TRUE)
    #Impact Effect
    W <- as(w_final, "CsparseMatrix")
    trMatc <- trW(W, type = "mult")
    set.seed(1)
    impactEffect <-
        summary(impacts(r.lag, tr = trMatc, R = 99),
                zstats = TRUE,
                short = TRUE)
    return(
        list(
            "Listw" = w_final,
            "Neighborhood" = w_final_nb,
            "Moran's I" = moran,
            "Spatial Lag Model" = slagm,
            "Direct & Indirect Effect" = impactEffect,
            "Final W" = bobot_custom
        )
    )
}

skema456 <- function(w, w_custom) {
    #Neighborhood
    k <- dim(dataIGI)[1] - 1
    knn <- knearneigh(coords, k)
    knn_nb <- knn2nb(knn)
    dlist <- nbdists(knn_nb, coords, longlat = TRUE)
    indis <- function(x)
        1 / (x ^ w)
    dlist <- lapply(dlist, indis)
    bobot1 <-
        nb2listw(knn_nb,
                 glist = dlist,
                 style = "B",
                 zero.policy = T)
    bobot1_mat <- listw2mat(bobot1)
    bobot_custom <- w_custom * bobot1_mat
    b_cust_tot <- rowSums(bobot_custom, na.rm = TRUE)
    bobot_custom <- bobot_custom / b_cust_tot
    w_final <- mat2listw(bobot_custom, style = "W")
    w_final_nb <- w_final[["neighbours"]]
    #Moran's I
    moran <- moran.test(dataIGI$IGI, w_final, randomisation = TRUE)
    #Spatial Lag Model
    r.ols <-
        lm(
            dataIGI$IGI ~ dataIGI$INFLASI + dataIGI$PMTB + dataIGI$UMK + dataIGI$PPS +
                dataIGI$PP
        )
    r.lag <- lagsarlm(r.ols, dataIGI, w_final, zero.policy = TRUE)
    slagm <- summary(r.lag, Nagelkerke = TRUE)
    #Impact Effect
    W <- as(w_final, "CsparseMatrix")
    trMatc <- trW(W, type = "mult")
    set.seed(1)
    impactEffect <-
        summary(impacts(r.lag, tr = trMatc, R = 99),
                zstats = TRUE,
                short = TRUE)
    return(
        list(
            "Listw" = w_final,
            "Neighborhood" = w_final_nb,
            "Moran's I" = moran,
            "Spatial Lag Model" = slagm,
            "Direct & Indirect Effect" = impactEffect,
            "Final W" = bobot_custom
        )
    )
}

skema789 <- function(o, w_custom) {
    #Neighborhood
    if (o == 1) {
        nb <- poly2nb(IGI17, queen = TRUE)
        w_queen <- nb2listw(nb, style = "B")
    }
    else {
        nb <- poly2nb(IGI17, queen = TRUE)
        nb_1 <- nblag(nb, o)
        nb_cum <- nblag_cumul(nb_1)
        w_queen <- nb2listw(nb_cum, style = "B")
    }
    w_queen <- listw2mat(w_queen)
    m_cust <- w_custom * w_queen
    m_cust_tot <- rowSums(m_cust, na.rm = TRUE)
    m_custom <- m_cust / m_cust_tot
    w_final <- mat2listw(m_custom, style = "W")
    w_final_nb <- w_final[["neighbours"]]
    #Moran's I
    moran <- moran.test(dataIGI$IGI, w_final, randomisation = TRUE)
    #Spatial Lag Model
    r.ols <-
        lm(
            dataIGI$IGI ~ dataIGI$INFLASI + dataIGI$PMTB + dataIGI$UMK + dataIGI$PPS +
                dataIGI$PP
        )
    r.lag <- lagsarlm(r.ols, dataIGI, w_final, zero.policy = TRUE)
    slagm <- summary(r.lag, Nagelkerke = TRUE)
    #Impact Effect
    W <- as(w_final, "CsparseMatrix")
    trMatc <- trW(W, type = "mult")
    set.seed(1)
    impactEffect <-
        summary(impacts(r.lag, tr = trMatc, R = 99),
                zstats = TRUE,
                short = TRUE)
    return(
        list(
            "Listw" = w_final,
            "Neighborhood" = w_final_nb,
            "Moran's I" = moran,
            "Spatial Lag Model" = slagm,
            "Direct & Indirect Effect" = impactEffect,
            "Final W" = m_custom
        )
    )
}

skema101112 <- function(o, w_custom) {
    #Neighborhood
    if (o == 1) {
        nb <- poly2nb(IGI17, queen = FALSE)
        w_queen <- nb2listw(nb, style = "B")
    }
    else {
        nb <- poly2nb(IGI17, queen = FALSE)
        nb_1 <- nblag(nb, o)
        nb_cum <- nblag_cumul(nb_1)
        w_queen <- nb2listw(nb_cum, style = "B")
    }
    w_queen <- listw2mat(w_queen)
    m_cust <- w_custom * w_queen
    m_cust_tot <- rowSums(m_cust, na.rm = TRUE)
    m_custom <- m_cust / m_cust_tot
    w_final <- mat2listw(m_custom, style = "W")
    w_final_nb <- w_final[["neighbours"]]
    #Moran's I
    moran <- moran.test(dataIGI$IGI, w_final, randomisation = TRUE)
    #Spatial Lag Model
    r.ols <-
        lm(
            dataIGI$IGI ~ dataIGI$INFLASI + dataIGI$PMTB + dataIGI$UMK + dataIGI$PPS +
                dataIGI$PP
        )
    r.lag <- lagsarlm(r.ols, dataIGI, w_final, zero.policy = TRUE)
    slagm <- summary(r.lag, Nagelkerke = TRUE)
    #Impact Effect
    W <- as(w_final, "CsparseMatrix")
    trMatc <- trW(W, type = "mult")
    set.seed(1)
    impactEffect <-
        summary(impacts(r.lag, tr = trMatc, R = 99),
                zstats = TRUE,
                short = TRUE)
    return(
        list(
            "Listw" = w_final,
            "Neighborhood" = w_final_nb,
            "Moran's I" = moran,
            "Spatial Lag Model" = slagm,
            "Direct & Indirect Effect" = impactEffect,
            "Final W" = m_custom
        )
    )
}

ui <- dashboardPage(
    dashboardHeader(title = "Analisis Spasial dengan Optimasi W"),
    dashboardSidebar(
        sidebarMenu(
            menuItem("Tentang",
                     tabName = "tentang",
                     icon = icon("book")),
            menuItem(
                "Data",
                tabName = "data",
                icon = icon("table"),
                startExpanded = TRUE,
                menuSubItem("Tabel",
                            tabName = "tabelTab"),
                menuSubItem("Peta",
                            tabName = "petaTab")
            ),
            menuItem(
                "Matriks",
                tabName = "matriks",
                icon = icon("align-justify")
            ),
            menuItem(
                "Analisis",
                tabName = "analisis",
                icon = icon("sign-in-alt")
            ),
            menuItem("Hasil",
                     tabName = "hasil",
                     icon = icon("sign-out-alt"))
            
        )
    ),
    dashboardBody(
        shinyDashboardThemes(theme = "poor_mans_flatly"),
        tabItems(
            tabItem(
                tabName = "tentang",
                fluidRow(
                    box(
                        title = "Selamat Datang",
                        status = "success",
                        collapsible = FALSE,
                        width = 12,
                        includeText("www/tentangBox.txt")
                    )
                ),
                fluidRow(
                    box(
                        title = "Data",
                        status = "success",
                        collapsible = FALSE,
                        width = 12,
                        tableOutput(outputId = "dataYBox"),
                        tableOutput(outputId = "dataBox")
                    )
                ),
                fluidRow(
                    box(
                        title = "Metode",
                        status = "success",
                        collapsible = FALSE,
                        width = 12,
                        includeText("www/metodeBox.txt")
                    )
                ),
                fluidRow(
                    box(
                        title = "Keterangan Skema",
                        status = "success",
                        collapsible = FALSE,
                        width = 12,
                        tableOutput(outputId = "skemaBox")
                    )
                )
            ),
            tabItem(tabName = "tabelTab",
                    h3("Tabel Data"),
                    fluidRow(
                        box(
                            status = "success",
                            solidHeader = FALSE,
                            DT::dataTableOutput(outputId = "tabel1"),
                            style = "height:600px;overflow-y: scroll;overflow-x: scroll;",
                            width = 600
                        )
                    )),
            tabItem(
                tabName = "petaTab",
                h3("Peta Choropleth"),
                leafletOutput(outputId = "peta1",
                              height = 550),
                absolutePanel(
                    id = "control1",
                    class = "panel panel-default",
                    fixed = TRUE,
                    draggable = FALSE,
                    top = 180,
                    left = 256.5,
                    right = "auto",
                    bottom = "auto",
                    width = 200,
                    height = "auto",
                    selectInput(
                        inputId = "var",
                        label = "Pilih variabel",
                        choices = c("IGI", "Inflasi", "PMTB", "UMK", "PPS", "PP"),
                        selected = "IGI"
                    )
                )
            ),
            tabItem(tabName = "analisis",
                    tabsetPanel(
                        tabPanel("Seleksi Matriks",
                                 sidebarLayout(
                                     sidebarPanel(
                                         width = 3,
                                         h4("Pemilihan Skema Kustomisasi Ketetanggaan"),
                                         radioButtons(
                                             inputId = "sk",
                                             label = "Pilih skema ketetanggaan",
                                             choices = list("KNN",
                                                            "Inverse Jarak",
                                                            "Queen",
                                                            "Rook"),
                                             selected = "KNN"
                                         ),
                                         radioButtons(
                                             inputId = "mc",
                                             label = "Pilih variabel sosioekonomi",
                                             choices = list("Migrasi",
                                                            "TPT",
                                                            "PDRB"),
                                             selected = "Migrasi"
                                         ),
                                         actionButton(inputId = "submit", label = "OK"),
                                     ),
                                     mainPanel(
                                         h3("Matriks Terpilih"),
                                         verbatimTextOutput("sko"),
                                         verbatimTextOutput("mco"),
                                         verbatimTextOutput("rangkuman"),
                                         verbatimTextOutput("aic"),
                                         verbatimTextOutput("gMDL"),
                                         verbatimTextOutput("chosen"),
                                         tags$hr(),
                                         downloadButton("dlMat", label = "Unduh Matriks Terpilih"),
                                         tags$hr(),
                                         fluidRow(
                                             box(
                                                 title = "Moran Plot",
                                                 width = 12,
                                                 collapsible = TRUE,
                                                 status = "success",
                                                 plotOutput("moranplot")
                                             )
                                         ),
                                         verbatimTextOutput("moran")
                                     )
                                 )),
                        tabPanel("Peta Ketetanggaan",
                                 sidebarLayout(
                                     sidebarPanel(
                                         width = 3,
                                         h4("Skema Ketetanggaan"),
                                         h6("Klik pada poligon untuk mengetahui wilayah tetangganya.")
                                     ),
                                     mainPanel(leafletOutput(
                                         outputId = "peta2",
                                         width = 760,
                                         height = 590
                                     ))
                                 ))
                    ),),
            tabItem(
                tabName = "hasil",
                h3("Hasil Pemodelan"),
                fluidRow(
                    box(
                        title = "Keterangan Skema",
                        status = "success",
                        collapsible = TRUE,
                        width = 10,
                        verbatimTextOutput("sko1"),
                        verbatimTextOutput("mco1")
                    )
                ),
                fluidRow(
                    box(
                        title = "Spatial Lag Model",
                        status = "success",
                        collapsible = TRUE,
                        width = 10,
                        verbatimTextOutput(outputId = "slm")
                    ),
                ),
                fluidRow(
                    box(
                        title = "Direct & Indirect Effect",
                        status = "success",
                        collapsible = TRUE,
                        width = 10,
                        verbatimTextOutput(outputId = "impact")
                    )
                )
            ),
            tabItem(
                tabName = "matriks",
                h3("Matriks Variabel Sosioekonomi"),
                sidebarLayout(
                    sidebarPanel(
                        width = 3,
                        radioButtons(
                            inputId = "matSel",
                            label = "Pilih matriks untuk ditampilkan",
                            choices = c("Migrasi", "TPT", "PDRB"),
                            selected = "Migrasi"
                        ),
                        tags$hr(),
                        h5("Keterangan"),
                        textOutput(outputId = "ketMat")
                    ),
                    mainPanel(
                        fluidRow(column(
                            width = 12,
                            box(
                                width = NULL,
                                status = "primary",
                                div(DT::dataTableOutput(outputId = "matVis", width = "100%"))
                            )
                        ))
                    )
                )
            )
        )
    )
)

choro_map <- function(var1, palet1, legend1) {
    shades <- brewer.pal(5, palet1)
    bin_pal <- colorBin(shades, var1, 5, pretty = TRUE)
    labels1 <- sprintf(
        "<div style = 'overflow-wrap: anywhere;'>
                  %s
                  <br/>IGI: %s</div>
                  <br/>Inflasi: %s&#37 </div>
                  <br/>PMTB: %s&#37 </div>
                  <br/>UMK: %s Ratus ribu Rp </div>
                  <br/>PP: %s&#37 </div>
                  <br/>PPS: %s&#37 </div>",
        IGI17@data$KABKOT,
        IGI17@data$IGI,
        IGI17@data$INFLASI,
        IGI17@data$PMTB,
        IGI17@data$UMK,
        IGI17@data$PP,
        IGI17@data$PPS
    ) %>% lapply(htmltools::HTML)
    
    leaflet(options = leafletOptions(zoomControl = FALSE)) %>%
        htmlwidgets::onRender("function(el,x) {
                L.control.zoom({ position: 'topright' }).addTo(this)
                          }") %>%
        addProviderTiles(providers$CartoDB.DarkMatter,
                         options = providerTileOptions((opacity = 1))) %>%
        setView(
            lng = mean(IGI17@data$x_long),
            lat = mean(IGI17@data$y_lat),
            zoom = 7.5
        ) %>%
        addPolygons(
            data = IGI17,
            stroke = TRUE,
            color = "white",
            fillColor = ~ bin_pal(var1),
            weight = 1.3,
            fillOpacity = 0.7,
            highlightOptions = highlightOptions(weight = 6,
                                                color = "yellow"),
            label = labels1
        ) %>%
        addLegend(
            data = IGI17,
            pal = bin_pal,
            title = legend1,
            values = ~ (var1),
            position = 'bottomright'
        ) %>%
        addResetMapButton()
}



server <- function(input, output) {
    matDl <- reactive({
        if (input$sk == "KNN" & input$mc == "Migrasi") {
            n <- substr(terpilih(), 5, 5)
            n <- as.numeric(n)
            w <- substr(terpilih(), 7, nchar(terpilih()))
            w <- as.numeric(w)
            finalw <-
                skema123(n = n,
                         w = w,
                         w_custom = w_migrasi)[["Final W"]]
            
        }
        else if (input$sk == "KNN" & input$mc == "TPT") {
            n <- substr(terpilih(), 5, 5)
            n <- as.numeric(n)
            w <- substr(terpilih(), 7, nchar(terpilih()))
            w <- as.numeric(w)
            finalw <-
                skema123(n = n,
                         w = w,
                         w_custom = w_tpt)[["Final W"]]
            
        }
        else if (input$sk == "KNN" & input$mc == "PDRB") {
            n <- substr(terpilih(), 5, 5)
            n <- as.numeric(n)
            w <- substr(terpilih(), 7, nchar(terpilih()))
            w <- as.numeric(w)
            finalw <-
                skema123(n = n,
                         w = w,
                         w_custom = w_pdrb)[["Final W"]]
            
        }
        else if (input$sk == "Inverse Jarak" &
                 input$mc == "Migrasi") {
            w <- substr(terpilih(), 5, 5)
            w <- as.numeric(w)
            finalw <-
                skema456(w = w, w_custom = w_migrasi)[["Final W"]]
            
        }
        else if (input$sk == "Inverse Jarak" &
                 input$mc == "TPT") {
            w <- substr(terpilih(), 5, 5)
            w <- as.numeric(w)
            finalw <-
                skema456(w = w, w_custom = w_tpt)[["Final W"]]
            
        }
        else if (input$sk == "Inverse Jarak" &
                 input$mc == "PDRB") {
            w <- substr(terpilih(), 5, 5)
            w <- as.numeric(w)
            finalw <-
                skema456(w = w, w_custom = w_pdrb)[["Final W"]]
            
        }
        else if (input$sk == "Queen" & input$mc == "Migrasi") {
            o <- substr(terpilih(), 5, nchar(terpilih()))
            o <- as.numeric(o)
            finalw <-
                skema789(o = o, w_custom = w_migrasi)[["Final W"]]
            
        }
        else if (input$sk == "Queen" & input$mc == "TPT") {
            o <- substr(terpilih(), 5, nchar(terpilih()))
            o <- as.numeric(o)
            finalw <-
                skema789(o = o, w_custom = w_tpt)[["Final W"]]
            
        }
        else if (input$sk == "Queen" & input$mc == "PDRB") {
            o <- substr(terpilih(), 5, nchar(terpilih()))
            o <- as.numeric(o)
            finalw <-
                skema789(o = o, w_custom = w_pdrb)[["Final W"]]
            
        }
        else if (input$sk == "Rook" & input$mc == "Migrasi") {
            o <- substr(terpilih(), 6, nchar(terpilih()))
            o <- as.numeric(o)
            finalw <-
                skema101112(o = o, w_custom = w_migrasi)[["Final W"]]
            
        }
        else if (input$sk == "Rook" & input$mc == "TPT") {
            o <- substr(terpilih(), 6, nchar(terpilih()))
            o <- as.numeric(o)
            finalw <-
                skema101112(o = o, w_custom = w_tpt)[["Final W"]]
        }
        else if (input$sk == "Rook" & input$mc == "PDRB") {
            o <- substr(terpilih(), 6, nchar(terpilih()))
            o <- as.numeric(o)
            finalw <-
                skema101112(o = o, w_custom = w_pdrb)[["Final W"]]
        }
        return(finalw)
    })
    
    output$dlMat <- downloadHandler(
        filename = function() {
            paste("finalw", ".csv", sep = "")
        },
        content = function(file) {
            write.csv(matDl(), file, row.names = FALSE)
        }
    )
    
    matSel1 <- reactive({
        switch(
            input$matSel,
            "Migrasi" = kk_migrasi,
            "TPT" = kk_tpt,
            "PDRB" = kk_pdrb
        )
    })
    
    matSel2 <- reactive({
        switch(
            input$matSel,
            "Migrasi" = paste("Matriks dengan elemen berupa selisih arus migrasi risen"),
            "TPT" = paste("Matriks dengan elemen berupa inverse selisih nilai TPT"),
            "PDRB" = paste("Matriks dengan elemen berupa inverse selisih nilai PDRB")
        )
    })
    
    output$matVis <- renderDataTable({
        DT::datatable(
            matSel1(),
            rownames = TRUE,
            extensions = c('FixedColumns', 'FixedHeader'),
            options = list(
                dom = 't',
                ordering = FALSE,
                paging = FALSE,
                fixedHeader = TRUE,
                scrollX = TRUE,
                scrollY = 600,
                fixedColumns = list(leftColumns = 1, rightColumns = 0)
            ),
            fillContainer = TRUE,
            class = "display"
        )
    })
    
    output$ketMat <- renderText(matSel2())
    
    output$sko <- renderText({
        input$submit
        isolate(paste("Skema ketetanggaan: ", input$sk))
    })
    
    output$mco <- renderText({
        input$submit
        isolate(paste("Matriks variabel sosioekonomi: ", input$mc))
    })
    
    output$rangkuman <- renderText({
        input$submit
        isolate({
            if (input$sk == "KNN" & input$mc == "Migrasi") {
                paste(
                    "Keterangan: \n- KNN dengan range disesuaikan\n- Inverse distance dengan power [0,4:4] dengan selisih 0,1\n- Matriks selisih migrasi"
                )
            }
            else if (input$sk == "KNN" & input$mc == "TPT") {
                paste(
                    "Keterangan: \n- KNN dengan range disesuaikan\n- Inverse distance dengan power [0,4:4] dengan selisih 0,1\n- Matriks inverse selisih TPT"
                )
            }
            else if (input$sk == "KNN" & input$mc == "PDRB") {
                paste(
                    "Keterangan: \n- KNN dengan range disesuaikan\n- Inverse distance dengan power [0,4:4] dengan selisih 0,1\n- Matriks inverse selisih PDRB"
                )
            }
            else if (input$sk == "Inverse Jarak" &
                     input$mc == "Migrasi") {
                paste(
                    "Keterangan: \n- Inverse distance dengan power [0,4:4] dengan selisih 0,1\n- Matriks selisih migrasi"
                )
            }
            else if (input$sk == "Inverse Jarak" &
                     input$mc == "TPT") {
                paste(
                    "Keterangan: \n- Inverse distance dengan power [0,4:4] dengan selisih 0,1\n- Matriks inverse selisih TPT"
                )
            }
            else if (input$sk == "Inverse Jarak" &
                     input$mc == "PDRB") {
                paste(
                    "Keterangan: \n- Inverse distance dengan power [0,4:4] dengan selisih 0,1\n- Matriks inverse selisih PDRB"
                )
            }
            else if (input$sk == "Queen" & input$mc == "Migrasi") {
                paste("Keterangan: \n- Queen ordo 1 s.d 3\n- Matriks selisih migrasi")
            }
            else if (input$sk == "Queen" & input$mc == "TPT") {
                paste("Keterangan: \n- Queen ordo 1 s.d 3\n- Matriks inverse selisih TPT")
            }
            else if (input$sk == "Queen" & input$mc == "PDRB") {
                paste("Keterangan: \n- Queen ordo 1 s.d 3\n- Matriks inverse selisih PDRB")
            }
            else if (input$sk == "Rook" & input$mc == "Migrasi") {
                paste("Keterangan: \n- Rook ordo 1 s.d 3\n- Matriks selisih migrasi")
            }
            else if (input$sk == "Rook" & input$mc == "TPT") {
                paste("Keterangan: \n- Rook ordo 1 s.d 3\n- Matriks inverse selisih TPT")
            }
            else if (input$sk == "Rook" & input$mc == "PDRB") {
                paste("Keterangan: \n- Rook ordo 1 s.d 3\n- Matriks inverse selisih PDRB")
            }
        })
    })
    
    output$sko1 <- renderText({
        input$submit
        isolate(paste("Skema ketetanggaan: ", input$sk))
        
    })
    
    output$mco1 <- renderText({
        input$submit
        isolate(paste("Matriks variabel sosioekonomi: ", input$mc))
    })
    
    w_final_nb <- reactive({
        if (input$sk == "KNN" & input$mc == "Migrasi") {
            n <- substr(terpilih(), 5, 5)
            n <- as.numeric(n)
            w <- substr(terpilih(), 7, nchar(terpilih()))
            w <- as.numeric(w)
            wfinal_nb <-
                skema123(n = n,
                         w = w,
                         w_custom = w_migrasi)[["Neighborhood"]]
        }
        else if (input$sk == "KNN" & input$mc == "TPT") {
            n <- substr(terpilih(), 5, 5)
            n <- as.numeric(n)
            w <- substr(terpilih(), 7, nchar(terpilih()))
            w <- as.numeric(w)
            wfinal_nb <-
                skema123(n = n,
                         w = w,
                         w_custom = w_tpt)[["Neighborhood"]]
        }
        else if (input$sk == "KNN" & input$mc == "PDRB") {
            n <- substr(terpilih(), 5, 5)
            n <- as.numeric(n)
            w <- substr(terpilih(), 7, nchar(terpilih()))
            w <- as.numeric(w)
            wfinal_nb <-
                skema123(n = n,
                         w = w,
                         w_custom = w_pdrb)[["Neighborhood"]]
        }
        else if (input$sk == "Inverse Jarak" &
                 input$mc == "Migrasi") {
            w <- substr(terpilih(), 5, 5)
            w <- as.numeric(w)
            wfinal_nb <-
                skema456(w = w, w_custom = w_migrasi)[["Neighborhood"]]
        }
        else if (input$sk == "Inverse Jarak" & input$mc == "TPT") {
            w <- substr(terpilih(), 5, 5)
            w <- as.numeric(w)
            wfinal_nb <-
                skema456(w = w, w_custom = w_tpt)[["Neighborhood"]]
        }
        else if (input$sk == "Inverse Jarak" & input$mc == "PDRB") {
            w <- substr(terpilih(), 5, 5)
            w <- as.numeric(w)
            wfinal_nb <-
                skema456(w = w, w_custom = w_pdrb)[["Neighborhood"]]
        }
        else if (input$sk == "Queen" & input$mc == "Migrasi") {
            o <- substr(terpilih(), 5, nchar(terpilih()))
            o <- as.numeric(o)
            wfinal_nb <-
                skema789(o = o, w_custom = w_migrasi)[["Neighborhood"]]
        }
        else if (input$sk == "Queen" & input$mc == "TPT") {
            o <- substr(terpilih(), 5, nchar(terpilih()))
            o <- as.numeric(o)
            wfinal_nb <-
                skema789(o = o, w_custom = w_tpt)[["Neighborhood"]]
        }
        else if (input$sk == "Queen" & input$mc == "PDRB") {
            o <- substr(terpilih(), 5, nchar(terpilih()))
            o <- as.numeric(o)
            wfinal_nb <-
                skema789(o = o, w_custom = w_pdrb)[["Neighborhood"]]
        }
        else if (input$sk == "Rook" & input$mc == "Migrasi") {
            o <- substr(terpilih(), 6, nchar(terpilih()))
            o <- as.numeric(o)
            wfinal_nb <-
                skema101112(o = o, w_custom = w_migrasi)[["Neighborhood"]]
        }
        else if (input$sk == "Rook" & input$mc == "TPT") {
            o <- substr(terpilih(), 6, nchar(terpilih()))
            o <- as.numeric(o)
            wfinal_nb <-
                skema101112(o = o, w_custom = w_tpt)[["Neighborhood"]]
        }
        else if (input$sk == "Rook" & input$mc == "PDRB") {
            o <- substr(terpilih(), 6, nchar(terpilih()))
            o <- as.numeric(o)
            wfinal_nb <-
                skema101112(o = o, w_custom = w_pdrb)[["Neighborhood"]]
        }
        return(wfinal_nb)
    })
    
    m1_boost <- reactive({
        if (input$sk == "KNN" & input$mc == "Migrasi") {
            m1 <-
                glmboost(w1,
                         vy,
                         control = boost_control(mstop = 1000, nu = 0.1),
                         center = FALSE)
        }
        else if (input$sk == "KNN" & input$mc == "TPT") {
            m1 <-
                glmboost(w2,
                         vy,
                         control = boost_control(mstop = 1000, nu = 0.1),
                         center = FALSE)
        }
        else if (input$sk == "KNN" & input$mc == "PDRB") {
            m1 <-
                glmboost(w3,
                         vy,
                         control = boost_control(mstop = 1000, nu = 0.1),
                         center = FALSE)
        }
        else if (input$sk == "Inverse Jarak" &
                 input$mc == "Migrasi") {
            m1 <-
                glmboost(w4,
                         vy,
                         control = boost_control(mstop = 1000, nu = 0.1),
                         center = FALSE)
        }
        else if (input$sk == "Inverse Jarak" & input$mc == "TPT") {
            m1 <-
                glmboost(w5,
                         vy,
                         control = boost_control(mstop = 1000, nu = 0.1),
                         center = FALSE)
        }
        else if (input$sk == "Inverse Jarak" & input$mc == "PDRB") {
            m1 <-
                glmboost(w6,
                         vy,
                         control = boost_control(mstop = 1000, nu = 0.1),
                         center = FALSE)
        }
        else if (input$sk == "Queen" & input$mc == "Migrasi") {
            m1 <-
                glmboost(w7,
                         vy,
                         control = boost_control(mstop = 1000, nu = 0.1),
                         center = FALSE)
        }
        else if (input$sk == "Queen" & input$mc == "TPT") {
            m1 <-
                glmboost(w8,
                         vy,
                         control = boost_control(mstop = 1000, nu = 0.1),
                         center = FALSE)
        }
        else if (input$sk == "Queen" & input$mc == "PDRB") {
            m1 <-
                glmboost(w9,
                         vy,
                         control = boost_control(mstop = 1000, nu = 0.1),
                         center = FALSE)
        }
        else if (input$sk == "Rook" & input$mc == "Migrasi") {
            m1 <-
                glmboost(w10,
                         vy,
                         control = boost_control(mstop = 1000, nu = 0.1),
                         center = FALSE)
        }
        else if (input$sk == "Rook" & input$mc == "TPT") {
            m1 <-
                glmboost(w11,
                         vy,
                         control = boost_control(mstop = 1000, nu = 0.1),
                         center = FALSE)
        }
        else if (input$sk == "Rook" & input$mc == "PDRB") {
            m1 <-
                glmboost(w12,
                         vy,
                         control = boost_control(mstop = 1000, nu = 0.1),
                         center = FALSE)
        }
        return(m1)
    })
    
    aicOutput <- reactive({
        m1 <- m1_boost()
        aic1 <- AIC(m1, method = "corrected")
        mbest1aic <- m1[mstop(aic1)]
        hasil1 <- names(coef(mbest1aic)[abs(coef(mbest1aic)) > 0])
        return(hasil1)
    })
    
    output$aic <- renderText({
        input$submit
        isolate(paste("By AIC: ", aicOutput()))
    })
    
    gMDLOutput <- reactive({
        m1 <- m1_boost()
        gMDL1 <- AIC(m1, method = "gMDL")
        mbest1gMDL = m1[mstop(gMDL1)]
        hasil2 <- names(coef(mbest1gMDL)[abs(coef(mbest1gMDL)) > 0])
        return(hasil2)
    })
    
    output$gMDL <- renderText({
        input$submit
        isolate(paste("By gmDL: ", gMDLOutput()))
    })
    
    terpilih <- reactive({
        if (aicOutput() == gMDLOutput()) {
            hasil <- aicOutput()
        }
        else {
            hasil <- "Not Found"
        }
        return(hasil)
    })
    
    notifstat <- function(){
        withProgress(message = 'Processing...', value = 0, {
            incProgress(0.1)
            aicOutput()
            incProgress(0.1)
            gMDLOutput()
            incProgress(0.1)
            terpilih()
            incProgress(0.1)
        })
    }
    
    observeEvent(input$submit, {
        notifstat()
    })
    
    output$chosen <- renderText({
        input$submit
        isolate(paste("Selected W: ", terpilih()))
    })
    
    output$moranplot <- renderPlot({
        input$submit
        isolate({
            if (input$sk == "KNN" & input$mc == "Migrasi") {
                n <- substr(terpilih(), 5, 5)
                n <- as.numeric(n)
                w <- substr(terpilih(), 7, nchar(terpilih()))
                w <- as.numeric(w)
                moran.plot(
                    x = dataIGI$IGI,
                    listw = skema123(
                        n = n,
                        w = w,
                        w_custom = w_migrasi
                    )[["Listw"]],
                    labels = as.character(dataIGI$KABKOT),
                    pch = 19
                )
            }
            else if (input$sk == "KNN" & input$mc == "TPT") {
                n <- substr(terpilih(), 5, 5)
                n <- as.numeric(n)
                w <- substr(terpilih(), 7, nchar(terpilih()))
                w <- as.numeric(w)
                moran.plot(
                    x = dataIGI$IGI,
                    listw = skema123(
                        n = n,
                        w = w,
                        w_custom = w_tpt
                    )[["Listw"]],
                    labels = as.character(dataIGI$KABKOT),
                    pch = 19
                )
            }
            else if (input$sk == "KNN" & input$mc == "PDRB") {
                n <- substr(terpilih(), 5, 5)
                n <- as.numeric(n)
                w <- substr(terpilih(), 7, nchar(terpilih()))
                w <- as.numeric(w)
                moran.plot(
                    x = dataIGI$IGI,
                    listw = skema123(
                        n = n,
                        w = w,
                        w_custom = w_pdrb
                    )[["Listw"]],
                    labels = as.character(dataIGI$KABKOT),
                    pch = 19
                )
            }
            else if (input$sk == "Inverse Jarak" &
                     input$mc == "Migrasi") {
                w <- substr(terpilih(), 5, 5)
                w <- as.numeric(w)
                moran.plot(
                    x = dataIGI$IGI,
                    listw = skema456(w = w, w_custom = w_migrasi)[["Listw"]],
                    labels = as.character(dataIGI$KABKOT),
                    pch = 19
                )
            }
            else if (input$sk == "Inverse Jarak" &
                     input$mc == "TPT") {
                w <- substr(terpilih(), 5, 5)
                w <- as.numeric(w)
                moran.plot(
                    x = dataIGI$IGI,
                    listw = skema456(w = w, w_custom = w_tpt)[["Listw"]],
                    labels = as.character(dataIGI$KABKOT),
                    pch = 19
                )
            }
            else if (input$sk == "Inverse Jarak" &
                     input$mc == "PDRB") {
                w <- substr(terpilih(), 5, 5)
                w <- as.numeric(w)
                moran.plot(
                    x = dataIGI$IGI,
                    listw = skema456(w = w, w_custom = w_pdrb)[["Listw"]],
                    labels = as.character(dataIGI$KABKOT),
                    pch = 19
                )
            }
            else if (input$sk == "Queen" & input$mc == "Migrasi") {
                o <- substr(terpilih(), 5, nchar(terpilih()))
                o <- as.numeric(o)
                moran.plot(
                    x = dataIGI$IGI,
                    listw = skema789(o = o, w_custom = w_migrasi)[["Listw"]],
                    labels = as.character(dataIGI$KABKOT),
                    pch = 19
                )
            }
            else if (input$sk == "Queen" & input$mc == "TPT") {
                o <- substr(terpilih(), 5, nchar(terpilih()))
                o <- as.numeric(o)
                moran.plot(
                    x = dataIGI$IGI,
                    listw = skema789(o = o, w_custom = w_tpt)[["Listw"]],
                    labels = as.character(dataIGI$KABKOT),
                    pch = 19
                )
            }
            else if (input$sk == "Queen" & input$mc == "PDRB") {
                o <- substr(terpilih(), 5, nchar(terpilih()))
                o <- as.numeric(o)
                moran.plot(
                    x = dataIGI$IGI,
                    listw = skema789(o = o, w_custom = w_pdrb)[["Listw"]],
                    labels = as.character(dataIGI$KABKOT),
                    pch = 19
                )
            }
            else if (input$sk == "Rook" & input$mc == "Migrasi") {
                o <- substr(terpilih(), 6, nchar(terpilih()))
                o <- as.numeric(o)
                moran.plot(
                    x = dataIGI$IGI,
                    listw = skema101112(o = o, w_custom = w_migrasi)[["Listw"]],
                    labels = as.character(dataIGI$KABKOT),
                    pch = 19
                )
            }
            else if (input$sk == "Rook" & input$mc == "TPT") {
                o <- substr(terpilih(), 6, nchar(terpilih()))
                o <- as.numeric(o)
                moran.plot(
                    x = dataIGI$IGI,
                    listw = skema101112(o = o, w_custom = w_tpt)[["Listw"]],
                    labels = as.character(dataIGI$KABKOT),
                    pch = 19
                )
            }
            else if (input$sk == "Rook" & input$mc == "PDRB") {
                o <- substr(terpilih(), 6, nchar(terpilih()))
                o <- as.numeric(o)
                moran.plot(
                    x = dataIGI$IGI,
                    listw = skema101112(o = o, w_custom = w_pdrb)[["Listw"]],
                    labels = as.character(dataIGI$KABKOT),
                    pch = 19
                )
            }
        })
    })
    
    output$moran <- renderPrint({
        input$submit
        isolate({
            if (input$sk == "KNN" & input$mc == "Migrasi") {
                n <- substr(terpilih(), 5, 5)
                n <- as.numeric(n)
                w <- substr(terpilih(), 7, nchar(terpilih()))
                w <- as.numeric(w)
                skema123(n = n,
                         w = w,
                         w_custom = w_migrasi)[["Moran's I"]]
            }
            else if (input$sk == "KNN" & input$mc == "TPT") {
                n <- substr(terpilih(), 5, 5)
                n <- as.numeric(n)
                w <- substr(terpilih(), 7, nchar(terpilih()))
                w <- as.numeric(w)
                skema123(n = n,
                         w = w,
                         w_custom = w_tpt)[["Moran's I"]]
            }
            else if (input$sk == "KNN" & input$mc == "PDRB") {
                n <- substr(terpilih(), 5, 5)
                n <- as.numeric(n)
                w <- substr(terpilih(), 7, nchar(terpilih()))
                w <- as.numeric(w)
                skema123(n = n,
                         w = w,
                         w_custom = w_pdrb)[["Moran's I"]]
            }
            else if (input$sk == "Inverse Jarak" &
                     input$mc == "Migrasi") {
                w <- substr(terpilih(), 5, 5)
                w <- as.numeric(w)
                skema456(w = w, w_custom = w_migrasi)[["Moran's I"]]
            }
            else if (input$sk == "Inverse Jarak" &
                     input$mc == "TPT") {
                w <- substr(terpilih(), 5, 5)
                w <- as.numeric(w)
                skema456(w = w, w_custom = w_tpt)[["Moran's I"]]
            }
            else if (input$sk == "Inverse Jarak" &
                     input$mc == "PDRB") {
                w <- substr(terpilih(), 5, 5)
                w <- as.numeric(w)
                skema456(w = w, w_custom = w_pdrb)[["Moran's I"]]
            }
            else if (input$sk == "Queen" & input$mc == "Migrasi") {
                o <- substr(terpilih(), 5, nchar(terpilih()))
                o <- as.numeric(o)
                skema789(o = o, w_custom = w_migrasi)[["Moran's I"]]
            }
            else if (input$sk == "Queen" & input$mc == "TPT") {
                o <- substr(terpilih(), 5, nchar(terpilih()))
                o <- as.numeric(o)
                skema789(o = o, w_custom = w_tpt)[["Moran's I"]]
            }
            else if (input$sk == "Queen" & input$mc == "PDRB") {
                o <- substr(terpilih(), 5, nchar(terpilih()))
                o <- as.numeric(o)
                skema789(o = o, w_custom = w_pdrb)[["Moran's I"]]
            }
            else if (input$sk == "Rook" & input$mc == "Migrasi") {
                o <- substr(terpilih(), 6, nchar(terpilih()))
                o <- as.numeric(o)
                skema101112(o = o, w_custom = w_migrasi)[["Moran's I"]]
            }
            else if (input$sk == "Rook" & input$mc == "TPT") {
                o <- substr(terpilih(), 6, nchar(terpilih()))
                o <- as.numeric(o)
                skema101112(o = o, w_custom = w_tpt)[["Moran's I"]]
            }
            else if (input$sk == "Rook" & input$mc == "PDRB") {
                o <- substr(terpilih(), 6, nchar(terpilih()))
                o <- as.numeric(o)
                skema101112(o = o, w_custom = w_pdrb)[["Moran's I"]]
            }
        })
    })
    
    output$slm <- renderPrint({
        input$submit
        isolate({
            if (input$sk == "KNN" & input$mc == "Migrasi") {
                n <- substr(terpilih(), 5, 5)
                n <- as.numeric(n)
                w <- substr(terpilih(), 7, nchar(terpilih()))
                w <- as.numeric(w)
                skema123(n = n,
                         w = w,
                         w_custom = w_migrasi)[["Spatial Lag Model"]]
            }
            else if (input$sk == "KNN" & input$mc == "TPT") {
                n <- substr(terpilih(), 5, 5)
                n <- as.numeric(n)
                w <- substr(terpilih(), 7, nchar(terpilih()))
                w <- as.numeric(w)
                skema123(n = n,
                         w = w,
                         w_custom = w_tpt)[["Spatial Lag Model"]]
            }
            else if (input$sk == "KNN" & input$mc == "PDRB") {
                n <- substr(terpilih(), 5, 5)
                n <- as.numeric(n)
                w <- substr(terpilih(), 7, nchar(terpilih()))
                w <- as.numeric(w)
                skema123(n = n,
                         w = w,
                         w_custom = w_pdrb)[["Spatial Lag Model"]]
            }
            else if (input$sk == "Inverse Jarak" &
                     input$mc == "Migrasi") {
                w <- substr(terpilih(), 5, 5)
                w <- as.numeric(w)
                skema456(w = w, w_custom = w_migrasi)[["Spatial Lag Model"]]
            }
            else if (input$sk == "Inverse Jarak" &
                     input$mc == "TPT") {
                w <- substr(terpilih(), 5, 5)
                w <- as.numeric(w)
                skema456(w = w, w_custom = w_tpt)[["Spatial Lag Model"]]
            }
            else if (input$sk == "Inverse Jarak" &
                     input$mc == "PDRB") {
                w <- substr(terpilih(), 5, 5)
                w <- as.numeric(w)
                skema456(w = w, w_custom = w_pdrb)[["Spatial Lag Model"]]
            }
            else if (input$sk == "Queen" & input$mc == "Migrasi") {
                o <- substr(terpilih(), 5, nchar(terpilih()))
                o <- as.numeric(o)
                skema789(o = o, w_custom = w_migrasi)[["Spatial Lag Model"]]
            }
            else if (input$sk == "Queen" & input$mc == "TPT") {
                o <- substr(terpilih(), 5, nchar(terpilih()))
                o <- as.numeric(o)
                skema789(o = o, w_custom = w_tpt)[["Spatial Lag Model"]]
            }
            else if (input$sk == "Queen" & input$mc == "PDRB") {
                o <- substr(terpilih(), 5, nchar(terpilih()))
                o <- as.numeric(o)
                skema789(o = o, w_custom = w_pdrb)[["Spatial Lag Model"]]
            }
            else if (input$sk == "Rook" & input$mc == "Migrasi") {
                o <- substr(terpilih(), 6, nchar(terpilih()))
                o <- as.numeric(o)
                skema101112(o = o, w_custom = w_migrasi)[["Spatial Lag Model"]]
            }
            else if (input$sk == "Rook" & input$mc == "TPT") {
                o <- substr(terpilih(), 6, nchar(terpilih()))
                o <- as.numeric(o)
                skema101112(o = o, w_custom = w_tpt)[["Spatial Lag Model"]]
            }
            else if (input$sk == "Rook" & input$mc == "PDRB") {
                o <- substr(terpilih(), 6, nchar(terpilih()))
                o <- as.numeric(o)
                skema101112(o = o, w_custom = w_pdrb)[["Spatial Lag Model"]]
            }
        })
    })
    
    output$impact <- renderPrint({
        input$submit
        isolate({
            if (input$sk == "KNN" & input$mc == "Migrasi") {
                n <- substr(terpilih(), 5, 5)
                n <- as.numeric(n)
                w <- substr(terpilih(), 7, nchar(terpilih()))
                w <- as.numeric(w)
                skema123(n = n,
                         w = w,
                         w_custom = w_migrasi)[["Direct & Indirect Effect"]]
            }
            else if (input$sk == "KNN" & input$mc == "TPT") {
                n <- substr(terpilih(), 5, 5)
                n <- as.numeric(n)
                w <- substr(terpilih(), 7, nchar(terpilih()))
                w <- as.numeric(w)
                skema123(n = n,
                         w = w,
                         w_custom = w_tpt)[["Direct & Indirect Effect"]]
            }
            else if (input$sk == "KNN" & input$mc == "PDRB") {
                n <- substr(terpilih(), 5, 5)
                n <- as.numeric(n)
                w <- substr(terpilih(), 7, nchar(terpilih()))
                w <- as.numeric(w)
                skema123(n = n,
                         w = w,
                         w_custom = w_pdrb)[["Direct & Indirect Effect"]]
            }
            else if (input$sk == "Inverse Jarak" &
                     input$mc == "Migrasi") {
                w <- substr(terpilih(), 5, 5)
                w <- as.numeric(w)
                skema456(w = w, w_custom = w_migrasi)[["Direct & Indirect Effect"]]
            }
            else if (input$sk == "Inverse Jarak" &
                     input$mc == "TPT") {
                w <- substr(terpilih(), 5, 5)
                w <- as.numeric(w)
                skema456(w = w, w_custom = w_tpt)[["Direct & Indirect Effect"]]
            }
            else if (input$sk == "Inverse Jarak" &
                     input$mc == "PDRB") {
                w <- substr(terpilih(), 5, 5)
                w <- as.numeric(w)
                skema456(w = w, w_custom = w_pdrb)[["Direct & Indirect Effect"]]
            }
            else if (input$sk == "Queen" & input$mc == "Migrasi") {
                o <- substr(terpilih(), 5, nchar(terpilih()))
                o <- as.numeric(o)
                skema789(o = o, w_custom = w_migrasi)[["Direct & Indirect Effect"]]
            }
            else if (input$sk == "Queen" & input$mc == "TPT") {
                o <- substr(terpilih(), 5, nchar(terpilih()))
                o <- as.numeric(o)
                skema789(o = o, w_custom = w_tpt)[["Direct & Indirect Effect"]]
            }
            else if (input$sk == "Queen" & input$mc == "PDRB") {
                o <- substr(terpilih(), 5, nchar(terpilih()))
                o <- as.numeric(o)
                skema789(o = o, w_custom = w_pdrb)[["Direct & Indirect Effect"]]
            }
            else if (input$sk == "Rook" & input$mc == "Migrasi") {
                o <- substr(terpilih(), 6, nchar(terpilih()))
                o <- as.numeric(o)
                skema101112(o = o, w_custom = w_migrasi)[["Direct & Indirect Effect"]]
            }
            else if (input$sk == "Rook" & input$mc == "TPT") {
                o <- substr(terpilih(), 6, nchar(terpilih()))
                o <- as.numeric(o)
                skema101112(o = o, w_custom = w_tpt)[["Direct & Indirect Effect"]]
            }
            else if (input$sk == "Rook" & input$mc == "PDRB") {
                o <- substr(terpilih(), 6, nchar(terpilih()))
                o <- as.numeric(o)
                skema101112(o = o, w_custom = w_pdrb)[["Direct & Indirect Effect"]]
            }
        })
    })
    
    output$dataYBox <- renderTable(dataYBox)
    
    output$dataBox <- renderTable(dataBox)
    
    output$skemaBox <- renderTable(skemaBox)
    
    output$tabel1 <- DT::renderDataTable({
        datatable(IGI17@data,
                  options = list(
                      scrollX = TRUE,
                      scrollCollapse = TRUE,
                      columnDefs = list(list(
                          target = c(2),
                          visible = TRUE,
                          width = '90'
                      ))
                  ))
    })
    
    output$peta1 <- renderLeaflet({
        var1 <- switch(
            input$var,
            "IGI" = IGI17@data$IGI,
            "Inflasi" = IGI17@data$INFLASI,
            "PMTB" = IGI17@data$PMTB,
            "UMK" = IGI17@data$UMK,
            "PPS" = IGI17@data$PPS,
            "PP" = IGI17@data$PP
        )
        
        palet1 <- switch(
            input$var,
            "IGI" = "YlOrRd",
            "Inflasi" = "YlGnBu",
            "PMTB" = "YlGn",
            "UMK" = "RdPu",
            "PPS" = "Purples",
            "PP" = "Reds"
        )
        
        legend1 <- switch(
            input$var,
            "IGI" = "IGI",
            "Inflasi" = "Inflasi (%)",
            "PMTB" = "PMTB (%)",
            "UMK" = "UMK (Ratus Ribu Rp)",
            "PPS" = "PPS (%)",
            "PP" = "PP (%)"
        )
        
        choro_map(var1, palet1, legend1)
    })
    
    output$peta2 <- renderLeaflet({
        labels2 <-
            sprintf("<div style = 'overflow-wrap: anywhere;'>%s</div>",
                    IGI17@data$KABKOT) %>% lapply(htmltools::HTML)
        
        peta2 <- leaflet(IGI17) %>%
            addProviderTiles(providers$CartoDB.Positron,
                             options = providerTileOptions((opacity = 1))) %>%
            addPolygons(
                layerId = ~ id,
                color = "blue",
                fillColor = "transparent",
                weight = 0.5,
                smoothFactor = 0.1,
                label = labels2
            ) %>%
            setView(
                lng = mean(IGI17@data$x_long),
                lat = mean(IGI17@data$y_lat),
                zoom = 7.5
            ) %>%
            addResetMapButton()
        
        peta2
    })
    
    click_tract <- eventReactive(input$peta2_shape_click, {
        return(input$peta2_shape_click$id)
    })
    
    focal_tract <- reactive({
        req(click_tract())
        return(IGI17[IGI17$id == click_tract(),])
    })
    
    neighbors <- reactive({
        req(click_tract())
        if ((input$sk == "KNN") & (input$mc == "Migrasi")) {
            return(IGI17[IGI17$id %in% w_final_nb()[[click_tract()]],])
        }
        else if ((input$sk == "KNN") & (input$mc == "TPT")) {
            return(IGI17[IGI17$id %in% w_final_nb()[[click_tract()]],])
        }
        else if ((input$sk == "KNN") & (input$mc == "PDRB")) {
            return(IGI17[IGI17$id %in% w_final_nb()[[click_tract()]],])
        }
        else if ((input$sk == "Inverse Jarak") &
                 (input$mc == "Migrasi")) {
            return(IGI17[IGI17$id %in% w_final_nb()[[click_tract()]],])
        }
        else if ((input$sk == "Inverse Jarak") &
                 (input$mc == "TPT")) {
            return(IGI17[IGI17$id %in% w_final_nb()[[click_tract()]],])
        }
        else if ((input$sk == "Inverse Jarak") &
                 (input$mc == "PDRB")) {
            return(IGI17[IGI17$id %in% w_final_nb()[[click_tract()]],])
        }
        else if ((input$sk == "Queen") & (input$mc == "Migrasi")) {
            return(IGI17[IGI17$id %in% w_final_nb()[[click_tract()]],])
        }
        else if ((input$sk == "Queen") & (input$mc == "TPT")) {
            return(IGI17[IGI17$id %in% w_final_nb()[[click_tract()]],])
        }
        else if ((input$sk == "Queen") & (input$mc == "PDRB")) {
            return(IGI17[IGI17$id %in% w_final_nb()[[click_tract()]],])
        }
        else if ((input$sk == "Rook") & (input$mc == "Migrasi")) {
            return(IGI17[IGI17$id %in% w_final_nb()[[click_tract()]],])
        }
        else if ((input$sk == "Rook") & (input$mc == "TPT")) {
            return(IGI17[IGI17$id %in% w_final_nb()[[click_tract()]],])
        }
        else if ((input$sk == "Rook") & (input$mc == "PDRB")) {
            return(IGI17[IGI17$id %in% w_final_nb()[[click_tract()]],])
        }
    })
    
    observe({
        req(click_tract())
        proxy <- leafletProxy('peta2')
        if (!is.null(neighbors())) {
            proxy %>%
                removeShape('focal') %>%
                clearGroup('neighbors') %>%
                addPolygons(
                    data = neighbors(),
                    fill = FALSE,
                    color = '#FFFF00',
                    group = 'neighbors',
                    opacity = 1
                ) %>%
                addPolygons(
                    data = focal_tract(),
                    color = '#00FFFF',
                    opacity = 1,
                    layerId = 'focal',
                    fillColor = 'transparent'
                )
        } else {
            proxy %>%
                removeShape('focal') %>%
                clearGroup('neighbors') %>%
                addPolygons(
                    data = focal_tract(),
                    color = '#00FFFF',
                    opacity = 1,
                    layerId = 'focal',
                    fillColor = 'transparent'
                )
        }
    })
}

shinyApp(ui, server)
