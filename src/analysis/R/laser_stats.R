library(dplyr)
library(lmerTest)
library(rlang)


sem = function(x) sqrt( var(x, na.rm=TRUE) / length(x))

create.animal.summary = function(df, var) {
  var = enquo(var)
  group_by(df, animal, exp, laserOn, stage_desc, channelLocation) %>% 
    dplyr::summarise(var.mean=mean(!!var, na.rm=TRUE),
                     var.sem=sem(!!var))
}

plot.laser.change = function(df, var) {
  var = enquo(var)
  summary.df = create.animal.summary(filter(df, stage_desc != 'after_stim'), !!var)
  xlabels = c('off', 'on')
  summary.df %>%
    ggplot() +
    geom_point(aes(x=laserOn, y=var.mean, color=animal)) +
    geom_line(aes(x=laserOn, y=var.mean, group=animal)) +
    geom_ribbon(aes(x=laserOn, 
                    ymin=var.mean - var.sem, 
                    ymax=var.mean + var.sem,
                    group=animal), 
                alpha=0.1) +
    facet_grid(. ~ exp + channelLocation, scales = 'free') +
    xlab('Stimulation') +
    scale_x_discrete(labels=xlabels) +
    gtheme
}

calc.laser.fixed.effects = function(tested.df, var) {
  var = enquo(var)
  var.name = quo_name(var)
  form = glue::glue('{varname} ~ laserOn + exp + laserOn * exp + (1 + laserOn | exp_animal)', varname=var.name)
  m.full = lmerTest::lmer(form, 
                          data=tested.df,
                          REML = TRUE) 
  print(summary(m.full))
  #print(coef(m.full)$animal)
  print(anova(m.full, refit=FALSE, ddf='Satterthwaite'))
  plot.model.diagnostics(m.full, tested.df$animal, tested.df$laserOn)
}

plot.model.diagnostics = function(m, animals, laser_states) {
  df = data.frame(res=unname(residuals(m)), 
                  fitted=unname(fitted(m)), 
                  animal=animals, 
                  laser=laser_states) 
  
  ggplot(df) +
    geom_point(aes(x=res, y=fitted, colour=animal, shape=laser)) +
    theme(legend.position = 'none') -> g1
  
  X=qqnorm(residuals(m), plot.it=FALSE)
  data.frame(x=X$x, y=X$y, animal=animals, laser=laser_states) %>%
    ggplot() +
    geom_point(aes(x=x, y=y, colour=animal, shape=laser)) +
    stat_qq_line(mapping=aes(sample=res), data=df) +
    xlab('Theoretical quantiles') +
    ylab('Sample quantiles') +
    theme(legend.position = 'top') -> g2
  
  plot_grid(g1,g2)
}