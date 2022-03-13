R"library('ggplot2')

fit_km = survfit(Surv(tempo/365, cens) ~ 1, data = $dados)

plot_dados = data.frame(fit_km$time, fit_km$surv, fit_km$n.event)
colnames(plot_dados) = c('Tempo', 'Sobr', 'Evento')
"

nc_a = modelo_nc["a"]
nc_b = modelo_nc["b"]

nm_a = modelo_nm["a"] 
nm_b = modelo_nm["b"] 
nm_p = modelo_nm["prob"] 
cm_a = modelo_cm["a"] 
cm_b = modelo_cm["b"] 
cm_p = modelo_cm["prob"] 

vals_dens = [0:0.01:23.7;]

surv_nc = sgompertz(vals_dens, median(nc_a), median(nc_b))
surv_nm = median(nm_p).^(pgompertz(vals_dens, median(nm_a), median(nm_b)))
surv_cm = median(cm_p) .+ (1-median(cm_p)).*(sgompertz(vals_dens, median(cm_a), median(cm_b)))

R"
d1 = data.frame(x = $vals_dens, y = $surv_nc)
d2 = data.frame(x = $vals_dens, y = $surv_nm)
d3 = data.frame(x = $vals_dens, y = $surv_cm)

plot_gg = ggplot() +
	geom_line(aes(x = Tempo, y = Sobr, colour = 'a'), data = plot_dados) +
	theme(plot.title = element_text(hjust = 0.5)) + 
	labs(x = 'Tempo', y = 'S(t)') +
	ylim(0, 1) + 
	theme_minimal()  +
	geom_line(aes(x = x, y = y, colour = 'b'), data = d1) +
	geom_line(aes(x = x, y = y, colour = 'c'), data = d2) +
	geom_line(aes(x = x, y = y, colour = 'd'), data = d3) +
    scale_color_manual(name = 'Ajuste',
				 values = c(
				   'royalblue',
				   'springgreen3',
				   'gold',
				   'brown2'),
				 labels = c('Kaplan-Meier','Sem cura', 'Com Mistura', 'Sem Mistura')) +
  theme(legend.position = 'bottom')
  
plot_gg

ggsave('figuras/ajustes.pdf', plot_gg, units = 'in', width = 7, height = 5)
"