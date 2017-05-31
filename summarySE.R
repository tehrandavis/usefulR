summary_se = function(data, measurevar, idvar, betweenvars) {
  require('dplyr')
  require('lazyeval')
  
  summarized = data %>%
    group_by_(.dots=c(idvar, betweenvars)) %>%
    dplyr::summarise_(.dots = list(
      MeasureVar = interp( ~mean( MEASUREVAR, na.rm=TRUE), MEASUREVAR = as.name(measurevar)  )
    )) %>%
    ungroup() %>%
    group_by_(.dots=betweenvars) %>%
    dplyr::summarise(N    = n(),
                     Mean = mean(MeasureVar, na.rm=TRUE),
                     SD   = sd(MeasureVar, na.rm=TRUE),
                     SE   = SD / sqrt(N) 
    ) %>%
    ungroup()
  
  colnames(summarized)[colnames(summarized)=='Mean'] = measurevar
  
  return( summarized )
  
}

summary_se_within = function(data, measurevar, idvar, withinvars= NULL, betweenvars= NULL) {
  require('dplyr')
  require('lazyeval')
  
  # Really just between subjects?
  if (is.null(withinvars)) {
    return( summary_se(data, measurevar, idvar, betweenvars) )
  }
  
  # Warn about impossibility of errorbars in mixed designs:
  if (!is.null(betweenvars)) {
    warning('Error bars cannot be accurately represented in mixed designs. ',
            'Treating each level of any between-subjects factors as separate experiment.')
  }
  
  # Collapse Multiple Observations in the Same Participant x Design Cell, Get Normed MeasureVar:
  normed = collapse_and_norm(data, measurevar, idvar, withinvars, betweenvars)
  
  # Get Correction factor:
  num_within_groups = prod( apply(data[,withinvars,drop=FALSE], MARGIN = 2, FUN = function(x) length(unique(x))) )
  correction_factor = sqrt( num_within_groups / (num_within_groups-1) )
  
  # Get Means, SDs, Etc:
  summarized = normed %>%
    group_by_(.dots= c(betweenvars, withinvars) ) %>%
    dplyr::summarise(N    = n(),
                     Mean = mean(MeasureVar, na.rm= TRUE),
                     SD   = sd(MeasureVarNormed, na.rm= TRUE),
                     SE   = SD / sqrt( N ),
                     CI   = SE * qt(.975, df = N-1) 
    ) %>%
    dplyr::mutate(SD = SD*correction_factor,
                  SE = SE*correction_factor,
                  CI = CI*correction_factor)
  
  colnames(summarized)[colnames(summarized)=='Mean'] = measurevar
  
  summarized
  
}

collapse_and_norm = function(data, measurevar, idvar, withinvars, betweenvars= NULL) {
  require('dplyr')
  require('lazyeval')
  
  # Collapse Multiple Observations in the Same Participant x Design Cell : 
  collapsed = data %>%
    group_by_(.dots= c(idvar, betweenvars, withinvars)) %>%
    dplyr::summarise_(.dots = list(MeasureVar = interp( ~mean( MEASUREVAR,na.rm=TRUE), MEASUREVAR = as.name(measurevar) )
    )) %>%
    ungroup() %>%
    group_by_(.dots= c(idvar, betweenvars) ) %>%
    dplyr::mutate(SubMean = mean(MeasureVar, na.rm=TRUE),
                  Diff    = MeasureVar - SubMean
    ) %>%
    ungroup() 
  
  # Get Normed Data:
  normed = collapsed %>%
    dplyr::mutate(MeasureVarNormed = (MeasureVar-SubMean) + mean(MeasureVar, na.rm=TRUE) ) 
}