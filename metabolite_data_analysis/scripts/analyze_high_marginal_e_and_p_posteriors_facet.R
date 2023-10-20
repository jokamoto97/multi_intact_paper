#make facet plot to compare multi-INTACT-only and multi-INTACT + KBA sign concordance results

suppressPackageStartupMessages(require(data.table));
suppressPackageStartupMessages(require(dplyr));
suppressPackageStartupMessages(require(ggplot2));
suppressPackageStartupMessages(require(tidyr));
suppressPackageStartupMessages(require(tidyverse));
suppressPackageStartupMessages(require(ggpattern));



infile1 <- "results/ukbb_gtex_multi_intact/marginal_ep_posterior_triplets_0_5_sign_concordance.txt"

infile2 <- "results/ukbb_gtex_multi_intact/marginal_ep_posterior_triplets_0_5_sign_concordance_kba.txt"

outfile <- "results/ukbb_gtex_multi_intact/marginal_ep_posterior_triplets_0_5_sign_concordance_facet_order_prop_conc.pdf"



m_intact_rst <- fread(infile1,h=T) %>%
	mutate(analysis = "m_intact_only")

m_intact_kba_rst <- fread(infile2,h=T) %>%
	mutate(analysis = "m_intact_kba")


plotdf <- rbind.data.frame(m_intact_rst,m_intact_kba_rst)


#Tissue order based on total pairs
tissue_order_total <- m_intact_rst %>%
	group_by(Tissue) %>%
	summarise(n = sum(N_pairs)) %>%
	arrange(desc(n)) %>%
	dplyr::select(Tissue) %>% unlist() %>% as.vector()


#Tissue order based on proportion signs agree
tissue_order_prop_sign <- m_intact_rst %>%
	pivot_wider(names_from = Type,values_from=N_pairs) %>%
	mutate(prop_conc = Concordant_marginal_z_signs/(Discordant_marginal_z_signs + Concordant_marginal_z_signs)) %>%
	arrange(desc(prop_conc)) %>%
	dplyr::select(Tissue) %>% unlist %>% as.vector()




pdf(outfile)
plotdf  %>%
        pivot_wider(names_from = Type,values_from=N_pairs) %>%
        mutate(prop_conc = format(round(Concordant_marginal_z_signs/(Discordant_marginal_z_signs + Concordant_marginal_z_signs),digits = 2),nsmall=2)) %>%
	mutate(n_tot = Concordant_marginal_z_signs +Discordant_marginal_z_signs) %>%
        pivot_longer(cols = Concordant_marginal_z_signs:Discordant_marginal_z_signs,
                        names_to = "Type",
                        values_to = "N_pairs") %>%
	mutate(prop_conc = as.character(prop_conc)) %>%
	mutate(prop_conc = case_when(Type == "Concordant_marginal_z_signs" ~ prop_conc,
                                     TRUE ~ " ")) %>%
	mutate(Type = case_when(Type == "Concordant_marginal_z_signs" ~ "TWAS and PWAS signs agree",
                                Type == "Discordant_marginal_z_signs" ~ "TWAS and PWAS signs disagree")) %>%
        mutate(Type = factor(Type)) %>%
        mutate(Type = fct_relevel(Type, c("TWAS and PWAS signs agree","TWAS and PWAS signs disagree"))) %>%
        mutate(Tissue = factor(Tissue)) %>%
        mutate(Tissue = fct_relevel(Tissue, rev(tissue_order_prop_sign))) %>%
        mutate(analysis = case_when(analysis == "m_intact_only" ~ "Multi-INTACT Only",
                                           TRUE ~ "Multi-INTACT and KBA")) %>%
        mutate(analysis = factor(analysis)) %>%
        mutate(analysis = fct_relevel(analysis,rev(c("Multi-INTACT and KBA","Multi-INTACT Only")))) %>%
        ggplot(aes(x=Tissue,y=N_pairs,fill = Type,pattern=Type)) +
	geom_bar_pattern(stat="identity",
			 color = "black",
                         pattern_fill = "black",
                         pattern_angle = 45,
                         pattern_density = 0.1,
			 pattern_spacing = 0.025,
			 pattern_key_scale_factor = 0.6)+
	scale_pattern_manual(values = c("TWAS and PWAS signs disagree" = "stripe", "TWAS and PWAS signs agree" = "none")) +
#        geom_bar(stat = 'identity') +
        geom_text(aes(label = prop_conc,y=n_tot),size = 2.25,hjust = -0.1) +
	facet_wrap(~analysis) +
	ylim(0,500) +
        coord_flip() +
        scale_fill_manual(values=c("#00BFC4","#F8766D")) +
#        geom_text(aes(label=n_pairs),hjust = -0.15,size = 2.75) +
        ylab("Number of implicated gene-metabolite pairs") +
#        scale_y_continuous(limits = c(0,450)) +
        theme_bw() +
        theme(text = element_text(size = 10, face = 'bold'),
        legend.title=element_blank(),legend.position = 'bottom')
dev.off()

