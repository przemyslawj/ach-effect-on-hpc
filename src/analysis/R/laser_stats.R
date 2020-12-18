library(dplyr)
library(lmerTest)
library(rlang)


sem = function(x) sqrt( var(x, na.rm=TRUE) / length(x))

create.animal.summary = function(df, var, group_vars=vars(animal, exp, exp_animal, correct, laserOn, stage_desc, channelLocation)) {
  var = enquo(var)
  if (!'correct' %in% colnames(df)) {
    df$correct = 1
  }
  group_by(df, !!!group_vars) %>% 
    dplyr::summarise(var.mean=mean(!!var, na.rm=TRUE),
                     var.sem=sem(!!var),
                     n=n())
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
    facet_grid(. ~ channelLocation + exp, scales = 'free') +
    xlab('Stimulation') +
    scale_x_discrete(labels=xlabels) +
    gtheme
}

cbPalette = c("#777777", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
expPalette = c("#777777", "#CC79A7")

plot.samples = function(df, var, xvar=laserOn, xlabels = c('off', 'on'),
                        group_vars=vars(animal, exp, exp_animal, correct, laserOn, stage_desc, channelLocation)) {
  var = enquo(var)
  xvar = enquo(xvar)
  summary.df = create.animal.summary(filter(df, stage_desc != 'after_stim'), !!var, group_vars)
  ggplot(df, aes(x=!!xvar, y=!!var)) +
    geom_jitter(aes(color=exp), shape=1, width=0.25, height=0, alpha=0.5, size=0.5, stroke=0.5) +
    geom_line(data=summary.df, aes(x=!!xvar, y=var.mean, group=animal) , color='#444444') +
    facet_grid(. ~ channelLocation + exp, scales = 'free') +
    xlab('Stimulation') +
    scale_x_discrete(labels = xlabels) +
    scale_color_manual(values=expPalette) +
    gtheme +
    theme(legend.position = 'none')
}

plot.samples.with.prediction = function(df, var, model) {
  var = enquo(var)
  summary.df = create.animal.summary(filter(df, stage_desc != 'after_stim'), !!var)
  xlabels = c('off', 'on')
  cbind(df, pred = predict(model)) %>%
    ggplot(aes(x=laserOn, y=!!var)) +
    #geom_jitter(aes(color=animal), shape=1, width=0.25, height=0) +
    geom_jitter(aes(color=laserOn), shape=1, width=0.25, height=0, alpha=0.7, size=0.5) +
    #geom_line(aes(x=laserOn, y=pred, group=animal), color='#444444') +
    geom_line(data=summary.df,
              aes(x=laserOn, y=var.mean, group=animal), 
              #linetype='dashed', 
              color='#444444') +
    facet_grid(. ~ channelLocation + exp, scales = 'free') +
    xlab('Stimulation') +
    scale_x_discrete(labels = xlabels) +
    scale_color_manual(values=c(nolaser.col, laser.col)) +
    #scale_color_manual(values=cbPalette) +
    gtheme
}

calc.laser.fixed.effects = function(tested.df, 
                                    var, 
                                    fixed.effects.str = 'laserOn * exp',
                                    random.effects.str = '(1 + laserOn | exp_animal)',
                                    print.mean.stats=TRUE) {
  var = enquo(var)
  var.name = quo_name(var)
  form = glue::glue('{varname} ~ {fixed.effects} + {random.effects}', 
                    varname=var.name,
                    fixed.effects=fixed.effects.str,
                    random.effects=random.effects.str)
  m.full = lmerTest::lmer(form, 
                          data=tested.df,
                          REML = TRUE) 
  print(summary(m.full))
  
  if (print.mean.stats) {
    print('Mean stats')
    group_by(tested.df, channelLocation, exp, laserOn) %>% 
      dplyr::summarise(mean(!!var), sem(!!var), n()) %>% print
  }
  
  #print(coef(m.full)$animal)
  print(anova(m.full, refit=FALSE, ddf='Satterthwaite'))
  g = plot.model.diagnostics(m.full, tested.df$animal, tested.df$laserOn)
  print(g)
  return(m.full)
}

pairwise.post.hoc = function(m.full, factor.interaction=c('laserOn:exp'), p.adjust.method='holm') {
  x = lsmeansLT(m.full, pairwise=TRUE, which=factor.interaction)
  x$`Pr(>|t|).adj` = p.adjust(x$`Pr(>|t|)`, method=p.adjust.method)
  x
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