






library(ecoCopula)
# fit marginal model
spider_pa <- stackedsdm(comm_rarfy, ~ 1,
                        data = data.frame(metadata.sub$sample_env),
                        family="negative.binomial")
# fit copula ordination 
spid_lv = cord(spider_pa)
# biplot
plot(spid_lv, biplot = TRUE)

envcol <- ifelse(metadata.sub$sample_env == "desert",
                 "#E69F00",
                 "#666666") 
colvec <- c("#E69F00", "#666666")

env <- data.frame(sample_env = metadata.sub$sample_env)

plot(spid_lv, site.col = envcol,
    xlim=c(-1, 3),
    ylim=c(-2, 2))
ordiellipse(spid_lv, env$sample_env, label=TRUE)
with(env,
     legend("topright", 
            legend = levels(as.factor(sample_env)),
            #bty = "n",
            col = colvec,
            pch = 21,
            pt.bg = colvec))

# Check residuals
plot(spider_pa)

# fit marginal model
spider_nb <- stackedsdm(comm_rarfy, ~ 1,
                        data = data.frame(metadata.sub$sample_env),
                        family = "negative.binomial") 
# fit copula ordination 
spid_gr = cgr(spider_nb, seed = 3)
# biplot
plot(spid_gr, pad = 1)


# Comment interpréter l'ASV? Est-ce que je regarde les ordres? espèces?


# NMDS
comm_hel <- decostand(comm_rarfy,method='hellinger')

# NMDS ordination based on Bray-Curtis distances
comm_NMDS <- metaMDS(comm_hel, distance="bray", trace=FALSE)
ordiplot(comm_NMDS, cex = 0.5, type = 'text', display='sites')
ordiellipse(comm_NMDS, metadata.sub$sample_env, label=TRUE)
plot(ef)