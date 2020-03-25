require(shiny)
require(spatstat)
require(maps)
require(maptools)
#require(ggmap)
require(ggplot2)
require(reshape)

# load in data------------------------------------------------
results = read.csv("twitt_data.csv",header =TRUE,comment.char = "",check.names=FALSE)
#-------------------------------------------------------------------------
# exclude non US twiit user------------------------------------------------
nonUS<-which(is.na(map.where(database = "usa", results$longtitude, results$latitude)))
reduce_results<-results[-nonUS,]
n<-dim(reduce_results)[1]
p<-dim(reduce_results)[2]-3

# Create an owin object from the USA outline(s)
usmap <- map('usa', fill=TRUE, col="transparent", plot=FALSE)

IDs <- sapply(strsplit(usmap$names, ":"), function(x) x[1])
usa <- map2SpatialPolygons(usmap, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))
spatstat.options(checkpolygons=FALSE)
usowin <- as.owin.SpatialPolygons(usa)
spatstat.options(checkpolygons=TRUE)

# Create a spatstat ppp object using USA window
location<-data.frame(x=reduce_results$longtitude,y=reduce_results$latitude)
pts <- ppp(location$x,location$y,window=usowin)
raw_density<-density(pts) # raw data density function
raw_density$v<-raw_density$v/sum(raw_density$v,na.rm = TRUE) # raw data normalize density function

data_raw<-data.frame(density=melt(raw_density$v)$value,Longtitude=sort(rep(raw_density$xcol,128)),Latitude=rep(raw_density$yrow,128))

choice_word<-seq(1,p)
names(choice_word)<-names(reduce_results)[seq(1,p)]

print(length(choice_word))

ui<-fluidPage(
  
  titlePanel("Word map Explorer"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("selectword","Input a word:",choices = choice_word),
      selectInput("selectword_pred","Input words for location prediction:",choices = choice_word,multiple = TRUE,selected =c(32,193,229,686)),
      actionButton("button", "Show Prediction"),
      h3( verbatimTextOutput("pred_word")),
      h3( verbatimTextOutput("hellinger")),
      h3( verbatimTextOutput("KL")),
      h3( verbatimTextOutput("wordpred_location"))
    )
    ,
    mainPanel(
      plotOutput("word_density"),
      plotOutput("pred_word_map")
    )
  ))

server <-function(input,output){
  
  colmnum<-reactive({
    as.numeric(input$selectword)
  })
  
  word_density_func<-function(c){
    inds <- reduce_results[,c] >0
    xypp_inds <- ppp(location$x[inds],location$y[inds],window = usowin)
    word_density<-density(xypp_inds)
    word_density$v<-density(xypp_inds)$v/sum(density(xypp_inds)$v,na.rm = TRUE)
    word_density_hellinger<-sum((word_density$v-raw_density$v)^2,na.rm = TRUE) # specific word data density function
    word_density_KL<-sum(word_density$v*log(word_density$v/raw_density$v),na.rm = TRUE) # specific word data density function
    return(list(word_density=word_density,inds=inds,word_density_hellinger=word_density_hellinger,word_density_KL=word_density_KL))
  }
  #----------------------------------------------------------
  # Naive Bayes Classifier
  #----------------------------------------------------------
  observeEvent(input$button,{
    cols <- as.numeric(input$selectword_pred)
    word_cond_prob<-lapply(cols,function(x) word_density_func(x)$word_density)
    post_prob<-raw_density$v
    
    for(i in 1:length(cols)) {
      post_prob<-post_prob*word_cond_prob[[i]]$v
    }
    
    data_pred<-na.omit(data.frame(post_prob=melt(post_prob)$value,Longtitude=sort(rep(raw_density$xcol,128)),Latitude=rep(raw_density$yrow,128)))
    pred_id<-which(data_pred[,"post_prob"]==max(data_pred[,"post_prob"]))
    
    pred_location<-data.frame(pred_long=data_pred$Longtitude[pred_id],pred_lat=data_pred$Latitude[pred_id])
   
    output$wordpred_location <- renderText({
    paste("Predicted Longtitude=", round(pred_location$pred_long,3),"\n","Predicted Latitude=", round(pred_location$pred_lat,3))
  })
    output$pred_word_map <- renderPlot({
      m<-map_data('state')
      ggplot(m,aes(x=long,y=lat))+
        scale_color_gradientn(colours = rev(rainbow(4)))+
        geom_point(aes_string(x="Longtitude",y="Latitude",col="post_prob"),size=2, data = data_pred,alpha=1)+
        geom_point(aes_string(x="pred_long",y="pred_lat"),size=1, data = pred_location,alpha=1)+
        theme(plot.title = element_text(hjust = 0.5))+
        labs(x ="Longtitude", y = "Latitude")+coord_fixed(ratio=1.5)+
        ggtitle("Posterior probability")
  })
  })   
  
  
  #----------------------------------------------------------
  # output plots
  #----------------------------------------------------------
  
  #output$raw_density <- renderPlot({
  #  m<-map_data('state')
  #  ggplot(m,aes(x=long,y=lat))+
  #    scale_color_gradientn(colours = rev(rainbow(4)))+
  #    geom_point(aes_string(x="Longtitude",y="Latitude",col="density"),size=2, data = data_raw,alpha=1)+
  #    geom_point(aes_string(x="longtitude",y="latitude"),size=1, data = reduce_results,alpha=1)+
  #    theme(plot.title = element_text(hjust = 0.5))+
  #    labs(x ="Longtitude", y = "Latitude")+coord_fixed(ratio=1.5)+
  #    ggtitle("overall data density")
  #})
  
  output$word_density <- renderPlot({
    require(reshape)
    word_density<-word_density_func(c=colmnum())$word_density
    inds<-word_density_func(c=colmnum())$inds
    data_word<-na.omit(data.frame(density=melt(word_density$v)$value,Longtitude=sort(rep(raw_density$xcol,128)),Latitude=rep(raw_density$yrow,128)))
    
    m<-map_data('state')
    ggplot(m,aes(x=long,y=lat))+
      scale_color_gradientn(colours = rev(rainbow(4)))+
      geom_point(aes_string(x="Longtitude",y="Latitude",col="density"),size=2, data = data_word,alpha=1)+
      #geom_point(aes_string(x="x",y="y"),size=1, data = location[inds,],alpha=1)+
      theme(plot.title = element_text(hjust = 0.5))+
      labs(x ="Longtitude", y = "Latitude")+coord_fixed(ratio=1.5)+
      ggtitle(paste(names(reduce_results)[colmnum()],"kernel density"))
  })
  
  output$hellinger = renderText({  
    word_density_hellinger<-word_density_func(c=colmnum())$word_density_hellinger
    paste("The hellinger distance:",round(word_density_hellinger,4)) 
  })
  
  output$KL = renderText({  
    word_density_KL<-word_density_func(c=colmnum())$word_density_KL
    paste("The K-L distance:",round(word_density_KL,4))    
  })
  
}

shinyApp(ui, server)