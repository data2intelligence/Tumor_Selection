countSelect <- 3
ratioExclude <- 0.1
ratioCutoff <- 0.5
ratioDiffCutoff <- 0.2

xlsx <- rio::import_list("Sample_data.xlsx")
tumorSize <- xlsx[[1]]


v1 <- tumorSize[tumorSize[,2]==1,5]
v2 <- tumorSize[tumorSize[,2]==2,5]

names(v1) <- paste0("ID",tumorSize[tumorSize[,2]==1,1])
names(v2) <- paste0("ID",tumorSize[tumorSize[,2]==2,1])

v1_sorted <- sort(v1)
v2_sorted <- sort(v2)

# Exclude top and bottom tumors
countExclude <- ceiling(length(v1) * 0.1)
v1_selected_pool <- v1_sorted[ (1+countExclude):(length(v1)-countExclude) ]
v2_selected_pool <- v2_sorted[ (1+countExclude):(length(v2)-countExclude) ]


set.seed(123)
sampleTimes <- 0
for(rpt in 1:2000)
{
  sampleTimes <- sampleTimes+1

  v1_selected <- sample(v1_selected_pool, countSelect)
  v2_selected <- sample(v2_selected_pool, countSelect)

  v1_rest <- v1[setdiff(names(v1),names(v1_selected))]
  v2_rest <- v2[setdiff(names(v2),names(v2_selected))]


  # condition 1: with each tumor size within all-tumor mean within 50%
  satisfyFlag1 <- TRUE
  for(i in 1:countSelect)
  {
    singleRatio1 <- v1_selected[i]/mean(v1)
    singleRatio2 <- v2_selected[i]/mean(v1)

    if(singleRatio1 > 1+ratioCutoff |
       singleRatio1 < 1-ratioCutoff |
       singleRatio2 > 1+ratioCutoff |
       singleRatio2 < 1-ratioCutoff
       )
    {
      satisfyFlag1 <- FALSE
      break
    }
  }

  if(satisfyFlag1)
  {
    print("condition 1 success")
  }else{
    next
  }

  # condition 2: After removing these 3 tumors,
  # the all-tumor mean ratio between two groups should not change beyond 20%
  satisfyFlag2 <- TRUE

  ratio_before_select <- mean(v1)/mean(v2)
  ratio_after_select <- mean(v1_rest)/mean(v2_rest)

  ratioDiff <- ratio_after_select/ratio_before_select

  if(ratioDiff > 1+ratioDiffCutoff |
     ratioDiff < 1-ratioDiffCutoff
  )
  {
    next
  }else{
    print("condition 2 success")
  }

  # condition 3: For the 3 sample tumors,
  # the 3-tumor mean ratio be similar to all-tumor mean ratio, within +- 20%
  satisfyFlag3 <- TRUE

  ratio_before_select <- mean(v1)/mean(v2)
  ratio_selected <- mean(v1_selected)/mean(v2_selected)

  ratioDiff <- ratio_selected/ratio_before_select

  if(ratioDiff > 1+ratioDiffCutoff |
     ratioDiff < 1-ratioDiffCutoff
  )
  {
    next
  }else{
    print("condition 3 success")
  }

  if(satisfyFlag1 & satisfyFlag2 & satisfyFlag3)
  {
    print("***********************")
    print("Group 1")
    #print(v1_selected)
    names(v1_sorted)[names(v1_sorted)%in%names(v1_selected)] <- paste0("<",names(v1_selected),">")
    print(v1_sorted)

    print("***********************")
    print("Group 2")
    #print(v2_selected)
    names(v2_sorted)[names(v2_sorted)%in%names(v2_selected)] <- paste0("<",names(v2_selected),">")
    print(v2_sorted)
    break
  }

}

